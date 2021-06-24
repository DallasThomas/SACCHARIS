#!/usr/bin/perl -w
#
# Script Name: cazy_extract.pl
# Created:     September 22, 2015
# Modified     February 15, 2017
#	       July 13, 2018
#
# Usage:       cazy_extract.pl [options]
#
# Description: This program - using a cazy family name and a qualifier
#               queries www.cazy.org for either genbank ID's or Protein data
#               ID's (determined by qualifier/group).  ID's are extracted
#               and used to extract sequence results from either NCBI or PDB
#
# Mods:        An issue was found where some of the first listed PDB ID's 
#               did not link to a Fasta entry, therefore to deal with this
#               need to hand structure extraction differently - will remove
#               the structure portion from the regular extraction method and
#               do everything in one step for each pdf entry from Cazy
#
# CPAN Lib:    Need to Install LWP::Simple and URI::Fetch
#
#              You also need to modify line 236 of TagParser.pm
#		--> Old Line:	bless $self, $package;
#		--> New Line: 	bless $self, ref($package) || $package;
############################################################################
use strict;
use Data::Dumper;
use LWP::Simple;
use HTML::TagParser;
use Carp qw( croak );

use File::Basename;

# NCBI - User defined variables - These varibles must be defined
my ($ncbi_api, $ncbi_name, $ncbi_email);

# User Must register for their own api id, etc to fill in these fields
$ncbi_api = '';
$ncbi_name = '';
$ncbi_email = '';

ncbi_error() if !defined $ncbi_api || !defined $ncbi_name || !defined $ncbi_email; 

sub ncbi_error {
   croak "NCBI values for api_key, tool name and email must be defined.\n";
}

# Set-up for Arguments
use Getopt::Long;
use Pod::Usage;
Getopt::Long::Configure ("no_ignore_case");

my ($opt_help, $opt_man, $opt_fam, $opt_grp, $opt_frag );

# Parse command line options
my $PROG_NAME = basename($0);

GetOptions(
    'help!'         => \$opt_help,
    'man!'          => \$opt_man,
    'family=s'	    => \$opt_fam,
    'group=s'       => \$opt_grp,
    'Fragment=s'    => \$opt_frag,
) or pod2usage(-verbose => 1) && exit;

# Exit program and output Help if...
pod2usage(-verbose => 1) && exit if defined $opt_help;
pod2usage(-verbose => 2) && exit if defined $opt_man;
pod2usage(-verbose => 1) && exit if !defined $opt_fam && !defined $opt_grp;

# Multiple families can be passed in - to deal with this we will split the
#  string into an array and then loop through the array
my @family = split(/,/,$opt_fam);

# Need to make sure the Group name is lowercase
$opt_grp = lc $opt_grp;

# Need to set fragment default to on - so default is remove fragments
my $remove_fragments = $opt_frag ? $opt_frag : 'true';
$remove_fragments = lc $remove_fragments;

# Need to make sure $remove_fragments is either "true" or "false"
if (($remove_fragments ne 'true') && ($remove_fragments ne 'false')) {
   print "Incorrect Format for Fragments\n\n";
   pod2usage(-verbose => 1) && exit;
}

# Set default variables
my $webbase = 'http://www.cazy.org/';
my @cazy;
my %testhash;
my $bigseq = 0;
my $SEQCOUNT = 0;
my ($FRAG, $DUP) = (0, 0);

# Time to start looping and extracting what we want
foreach (@family) {

   ###############################################################################################
   # Webpage URL & Sequence Count
   ###############################################################################################
   # Create the webpage URL for Cazy
   my $f = uc $_;
   my $webpage = $webbase . $f . "_" . $opt_grp . ".html";

   # If there are more than 1000 sequences, there is more than one page
   #  we need to determine the number of sequences and then proceed as such
   my @pages = &query_pages($webpage);

   # Determine the total number of Sequences
   &get_seqcount($webpage);


   ###############################################################################################
   # Sequence Extraction
   ###############################################################################################

   # Variable Declaration
   my (@results, @acc_list, $iterator);
   my ($arr_pos, $list_size) = (0, 0);

   # Extraction from Query - Method will differ if user selected group is structure as extract from
   #  PDB and not NCBI 
   if ( $opt_grp eq 'structure' ) {

      # Call Query Structure for main webpage
      my $pdbstuff = &query_structure($webpage);
      push (@results, $pdbstuff);

      # If there are more than pages of structure call Query structure for each of these
      if ($bigseq == 1) {
	 foreach (@pages) { 
	   my $pdb_res = &query_structure($_); 
	   push (@results, $pdb_res);
	 }
      }

   } else { 

      # Extract the ID list from Cazy Using the Webpage
      &query_cazy($webpage);

      # Extract ID list from array of other Cazy pages
      if ($bigseq == 1) {
         foreach (@pages) {
           &query_cazy($_);
         }
      }

      # Seem to really only be able to extract 500 records at a time - so for the larger data samples
      #  we need to split them into subsets

      $iterator = 1;

      foreach (@cazy) {
        ########################################################################################################
        # Ran into a case where the Description in Cazy had a Number following description kind of like an
        #  accession number that messed up the script - have tried two fixes, this commented one is the
        #  second fix, where I screen these results out after they have been added to the array, the first
        #  solution is found in &query_cazy, and prevents them from being added, using that method right now
        #  and hoping it works for all cases
        ########################################################################################################
        #      # Need to make sure I actually have just accession numbers and nothing else
        #      if ($_ !~ m/^([A-Z]{2,5}\d{2,7})[.]/ || $_ !~ m/^([A-Z]{2,5}\_\d{2,9})[.]/) {
        #	       next;
        #      }
        ########################################################################################################
        if ($iterator == 351) {
 	        $iterator = 1;
 	        $arr_pos++;
        }

        if ($iterator == 1) {
 	        $acc_list[$arr_pos] = $_ . ",";
        } elsif ( ($iterator == 350) || ($list_size == $#cazy) ) {
 	        $acc_list[$arr_pos] = $acc_list[$arr_pos] . $_;
        } else {
 	        $acc_list[$arr_pos] = $acc_list[$arr_pos] . $_ . ",";
        }
        $iterator++;
        $list_size++;
      }

      # Take this array of strings - that is the lists of 500 accension numbers and loop through array
      #  querying ncbi and writing to outfile
      $iterator = 0;

#print Dumper(@acc_list); exit;

      foreach (@acc_list) {
         my $ncbistuff = &query_ncbi($_);
         push (@results, $ncbistuff);
         $iterator++;
      }
    }

#exit;

   ###############################################################################################
   # Print Results
   ###############################################################################################

   # The last step is to print the results to file, and modify the header of the sequences in the file
   my $outfile = $_ . '_cazy.fasta';
   open (OUT, ">$outfile") or die "Can't open $outfile:$!\n";
   binmode OUT, ":encoding(UTF-8)";  

   # Need to split the ncbi data in each string array to another array form processing
   my @final_results;

   foreach (@results) {
      my @tmp_results = split(/\n/, $_);
      push (@final_results, @tmp_results);
   }

   # Print to file
   if ( $opt_grp eq 'structure' ) {

     foreach (@final_results) {
       my (@l, @v, $newline);
       if ( $_ =~ m/.*?\>/ ) {
         @v = split(/\>/, $_);
         @l = split(/:/, $v[1]);
         $newline = ">" . $l[0] . "_" . $f . " " . $v[1];
         print OUT $newline, "\n";
       } else {
         print OUT $_, "\n";
       }
     }

   } else {

#print Dumper(@final_results); exit;

     foreach (@final_results) {
        my (@l, @v, $newline);
        if ( $_ =~ m/.*?\>/ ) {
          if ( $_ =~ m/\>[bdeglprs].*?\|/ ) {
            @l = split(/\|/, $_);
  	    @v = split(/\>/, $_);
            unless($l[3]) { print Dumper($_), "\n"; exit; }
            $newline = ">" . $l[3] . "_" . $f . " " . $v[1];
  	  } else {
            @v = split(/\>/, $_);
   	    @l = split(/ /, $v[1]);
  	    $newline = ">" . $l[0] . "_" . $f . " " . $v[1];
  	  }
  	  print OUT $newline, "\n";
        } else {
	  print OUT $_, "\n";
        }
     }

   }
   close OUT;
}

exit;

################################################################################
# This subroutine - Determines from the original page the total number of
#  sequences in the group - sets this value to global count, subroutine also
#  parses the html to determine how many pages are included with each family
#  group and sets page array to url for each subsequent page with accession
#  numbers
################################################################################
sub query_pages {

    # Passed in Variables
    my $web = shift;

    # Variable Declaration and Retrieval from Website
    my @array;
    my ($newpage, $num);
    my $count = 0;
    my $html = new HTML::TagParser->new( $web );
    my @clist = $html->getElementsByTagName( "span" );
    my @plist = $html->getElementsByTagName( "a" );
    my $debut_val;

    # Determine how many pages there are from the source
    foreach my $elem ( @plist ) {
       my $attr = $elem->attributes;
       my $text = $elem->innerText;
       my (@text, $value);

       foreach my $key ( sort keys %$attr ) {

	  # Need to know how many sequences per page
	  if ($key eq "href") {
	     $value = $elem->getAttribute( $key );
	     if ( $value =~ m/[a-zA-Z_].*?\?debut_[A-Z].*?=(.*?)#pagination_[A-Z].*?/i) { 
		$num = $1;
		$debut_val = $1 if $value =~ m/[a-zA-Z_].*?\?debut_([A-Z]+)/i;		
	     }
	  }

	  # Need to know how many pages
    	  if ($key eq "rel") {
      	     @text = split(/ /, $text);
	     if ($text[0] =~ m/\d+/) {
	        $count = $text[0] > $count ? $text[0] : $count;
	     }
	  }
       }
    }

    # Use the count to create the new pages and push onto array
    for (my $i = 1; $i < $count; $i++) {
        $newpage = $web . "?debut_" . $debut_val . "=" . $i*$num . "#pagination_" . $debut_val;
	push(@array, $newpage);
	$bigseq = 1;
    }
    return @array;
}

################################################################################
# This subroutine - Determines the total number of sequences in the group and
#  outputs this number to a file
################################################################################
sub get_seqcount {

    # Passed in Variables
    my $web = shift;

    # Variable Declaration and Retrieval from Website
    my @array;
    my ($newpage, $num);
    my $count = 0;
    my $html = new HTML::TagParser->new( $web );
    my @clist = $html->getElementsByTagName( "span" );
    my @plist = $html->getElementsByTagName( "a" );

    # Determine the total Number of sequences in the family group
    foreach my $elem ( @clist ) {
       my $attr = $elem->attributes;
       my $text = $elem->innerText;
       my @text;

       foreach my $key ( sort keys %$attr ) {
          if ($key eq "class") {
             @text = split(/ /, $text);

 	     # Match with the group passed
             if ( $text[0] =~ m/^$opt_grp$/i ) {
                if ( $text[1] =~ m/;(.*?)&/ ) {
                   $SEQCOUNT = $1;
                } elsif ( $text[1] =~ m/;(.*?)$/ ) {
                   $SEQCOUNT = $1;
                } else {
                   $SEQCOUNT = 0;
                }
             }
          }
       }
    }

    my $modified = $SEQCOUNT - ($DUP + $FRAG);

    # Create the oufile and print this seqcount out
    open (SEQOUT, ">>sequence.count") or die "Cannnot open sequence.count for writing: $!";
    print SEQOUT "$SEQCOUNT, $DUP, $FRAG, $modified";
    close (SEQOUT);
}

################################################################################
# This subroutine - uses the supplied webpage to get a list of Genbank
#  Accession ID's - have to guarantee there are no duplicates
################################################################################
sub query_cazy {

    # Passed in Variables
    my $web = shift;

#print Dumper($web);

    # Variable Declaration and Retrieval from Website
    my $html = new HTML::TagParser->new( $web );
    my @list = $html->getElementsByTagName( "td" );
    my $fragment = 0;

#print Dumper(@list);

    # Parse HTML Data
    foreach my $elem ( @list ) {
        my $attr = $elem->attributes;
        my $text = $elem->innerText;
        my ( @text, $t1 );

        foreach my $key ( sort keys %$attr ) {

          # First identify if the following sequence is a fragment from Protein Name
          if ( ( $key eq 'id' ) && ( $text =~ /fragment/ ) && ( $remove_fragments eq 'true' ) ) {
             $fragment = 1;
             $FRAG++;
          }

          # If the sequence is a fragment skip over or else add to list
          ########################################################################################################
          # To deal with cases, where accession like numbers are found in the organism name, I need to make sure I
          #  ignore these or the program crashes as I am adding descriptions to my array.  A second fix I
          #  implemented for this has been commented out above.  Here I have replaced the conditional commented
          #  out below, with a similar conditional, the difference is it is looking for text that starts with
          #  the Accession Number - this will remain good as long as there are no other organism names or protein
          #  name that start with a simiar pattern.  If so we might need to look at implementing the solution
          #  above
          ########################################################################################################
          #
          #  if ( ( $key eq 'id' ) && ( $text =~ m/.*?([A-Z]{2,5}\d{2,7})[.]/ || $text =~ m/.*?([A-Z]{2,5}\_\d{2,9})[.]/ ) && !( $text =~ m/&nbsp/ ) )
          #
          ########################################################################################################

          if ( ( $key eq 'id' ) && ( $text =~ m/^([A-Z]{2,5}\d{2,7})[.]/ || $text =~ m/^([A-Z]{2,5}\_\d{2,9})[.]/ ) && !( $text =~ m/&nbsp/ ) ) {
             if ( ( $text =~ m/.*?[.]\d{1,2}[A-Z0-9]/ ) && ( $fragment == 0 ) ) {
                @text = split(/([.]\d{1})/, $text);
		if ( $text[1] =~ m/\.$/ ) {
		   $t1 = $text[0];
		} else {
		   $t1 = $text[0] . $text[1];
		}

                # Test for duplicate accession numbers
                $DUP++ if exists $testhash{$t1};
                next if exists $testhash{$t1};
                # If does not exist - set accession number to array and set value to hash
                push(@cazy,$t1);
                $testhash{$t1} = 1;
             } elsif ( $fragment == 0 ) {
                if ( $text =~ m/\.$/ ) {
                   @text = split(/\./, $text);
                   $t1 = $text[0];
                } else {
                   $t1 = $text;
                }

                # Test for duplicate accession numbers
                $DUP++ if exists $testhash{$t1};
                next if exists $testhash{$t1};
                # If does not exist - set accession number to array and set value to hash
                push(@cazy,$t1);
                $testhash{$t1} = 1;
             } else {
                $fragment = 0;
             }
          }
        }
    }
    return;
}

################################################################################
# This subroutine - uses the Genbank Accession Id List from cazy to extract
#  the fasta information from NCBI
################################################################################
sub query_ncbi{

    # Passed in Variables
    my $id_list = shift;

    # Variable Declaration and Set-up
    my ($ebase, $esearch, $esearch_result, $efetch, $efetch_result);
    my ($web, $key, $count);

    # Set-up the Query URL
    my $utils = 'http://www.ncbi.nlm.nih.gov/entrez/eutils';

    # Set-up a search out to the eSearch program:
    #    - db is protein and search term is left blank for now <term>
    #$esearch = $utils . '/esearch.fcgi?db=protein&term=';
    $ebase = $utils . '/esearch.fcgi?db=protein&email=' . $ncbi_email . '&tool='
		. $ncbi_name . '&api_key=' . $ncbi_api;
    $esearch = $ebase . '&term='; 

    # Submit the search to retrieve a count of total number of sequences
    $esearch_result = get( $esearch . $id_list );

    # Extract the count and submit search again to retrieve XML based results
    # - set the number of results we want to count <retmax>
    $esearch_result =~ m|.*<Count>(.*)</Count>.*|s;
    $count = $1;

    #$esearch = $utils . '/esearch.fcgi?db=protein&retmax=' . $count . '&term=';
    $esearch = $ebase . '&retmax=' . $count . '&term=';
    $esearch_result = get( $esearch . $id_list . '&usehistory=y' );

    # Extract the WebEnv and QueryKey
    $esearch_result =~ m|.*<QueryKey>(.*)</QueryKey>.*|s;
    $key = $1;

    $esearch_result =~ m|.*<WebEnv>(.*)</WebEnv>.*|s;
    $web = $1;

    # Fetch the Fasta data from NCBI using the esearch results
    #$efetch = $utils . '/efetch.fcgi?db=protein&query_key=' . $key . '&WebEnv='
    #             . $web . '&rettype=fasta&retmode=text';
    $ebase = $utils . '/efetch.fcgi?db=protein&email=' . $ncbi_email . '&tool='
                . $ncbi_name . '&api_key=' . $ncbi_api;
    $efetch = $ebase . '&query_key=' . $key . '&WebEnv=' . $web . '&rettype=fasta&retmode=text';

    $efetch_result = get( $efetch );
    
    # Remove Spaces between each of the sequences
    $efetch_result =~ s/\n+/\n/g;

    return $efetch_result;
}

################################################################################
# This subroutine - uses the supplied webpage to get the PDB entry and use that
#  PDB Id from cazy to extract the fasta information from PDB
################################################################################
sub query_structure{

    # Passed in Variables
    my $web = shift;

    # Variable Declaration and Retrieval from Website
    my $html = new HTML::TagParser->new( $web );
    my @list = $html->getElementsByTagName( "td" );
    my $fragment = 0;
    my %testhash;
    my ($pdb_fasta, @seqs);

    # Parse HTML Data
    foreach my $elem ( @list ) {
        my $attr = $elem->attributes;
        my $text = $elem->innerText;
        my @text;

        foreach my $key ( sort keys %$attr ) {

          # First identify if the following sequence is a fragment from Protein Name
          if ( ( $key eq 'id' ) && ( $text =~ /fragment/ ) && ( $remove_fragments eq 'true' ) ) {
             $fragment = 1;
             $FRAG++;
          }

          # Extracting PDB ID's then use to extract fasta
	  # If the sequence is a fragment skip over else, add to list - PDB IDs are 4 Character ID with first character a digit
          #  followed by three alphanumeric characters
          if ( ( $key eq 'id' ) && ( $text =~ m/^[1-9][0-9a-zA-z]{3}\[.*/ ) && !( $text =~ m/&nbsp/ ) ) {
            if ( ( $text =~ m/^[1-9][0-9a-zA-z]{3}\[.*\]$/ ) && ( $fragment == 0 ) ) {
              @text = split(/(\[.*?\])/, $text);

              # Test for duplicate accession numbers
              $DUP++ if exists $testhash{$text[0]};
              next if exists $testhash{$text[0]};

              # If does not exist - Query PDB with ID - If fasta extracted add to string else next
	      $pdb_fasta = &query_pdb($text[0]);
              if ($pdb_fasta eq '1') {
	         next;
	      } else {
	         push(@seqs,$pdb_fasta);
                 $testhash{$text[0]} = 1;
              }

            } elsif ( ( $text =~ m/.*?\[.*?\][1-9]/ ) && ( $fragment == 0 ) ) {
              @text = split(/(\[.*?\])/, $text);

              # Test for duplicate accession numbers
              $DUP++ if exists $testhash{$text[0]};
              next if exists $testhash{$text[0]};

              # If does not exist - Query PDB with ID - If fasta extracted add to string else next
              $pdb_fasta = &query_pdb($text[0]);
              if ($pdb_fasta eq '1') {
                 next;
              } else {
                 push(@seqs,$pdb_fasta);
                 $testhash{$text[0]} = 1;
              }

            } elsif ( $fragment == 0 ) {
              # Extract ID or IDs from this selection
              @text = split(/(\s+)/, $text);
              my $tc = 0;

              foreach (@text) {
                if ( $_ =~ m/^[1-9][0-9a-zA-z]{3}\[.*\]/ ) {
                  my @tt = split(/(\[.*?\])/, $_);

                  # Test for duplicate accession numbers
                  $DUP++ if exists $testhash{$tt[0]};
                  next if exists $testhash{$tt[0]};

                  # If does not exist - Query PDB with ID - If fasta extracted add to string else next
		  $pdb_fasta = &query_pdb($tt[0]);
                  if ( $tc == 0 ) {
	            if ($pdb_fasta eq '1') {
                      next;
                    } else {
                      push(@seqs,$pdb_fasta); 
                      $testhash{$tt[0]} = 1;
		      $tc++;
                    }
                  } else {
                    $testhash{$tt[0]} = 1;
                  }
                }
              }
            } else {
              $fragment = 0;
            }
	  }
	}
    }

    # Need to split array into a string separated with Newlines
    my $struct_seqs = join ( '\n', @seqs); 

    return $struct_seqs;
}


################################################################################
# This subroutine - uses the past in PDB ID to extract fasta information from
#  the Protein Data Bank
################################################################################
sub query_pdb{

    # Passed in Variables
    my $id = shift;

    # Variable Declaration and Set-up
    my ($pdb, @pdb, $pdb_result);

    # Set-up the Protein Data Bank Query URL
    my $utils = 'http://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList=';
    my $tail = "&compressionType=uncompressed";

    # Grab fasta information with URL and id_list
    $pdb = get( $utils . $id . $tail );

    # Verify that a fasta sequence was output
    if ( $pdb =~ m/^\>[1-9A-Za-z]*.?/ ) {
       @pdb = split( />/, $pdb );

       # We just need the top sequence
       $pdb_result = '>' . $pdb[1];
    } else {
       $pdb_result = '1'; 
    } 

    return $pdb_result;
}

__END__


=pod

=head1 NAME

 cazy_extract.pl

=head1 SYNOPSIS

 cazy_extract.pl [-Options [--] [Arguments...]]

=head1 DESCRIPTION

 Extracts Protein Sequence data from NCBI or PDB according to Cazy Family and
 group.

 Requires a Cazy Family and Group (subset summary).  This information
 is used to create the URL of the specific Cazy search, whereby Genbank ID's or
 Protein Databank ID's are extracted from Cazy.  These ID's are used to extract
 the Protein Sequence Data from either NCBI or PDB and output to files.

 A group selection of 'structure' uses PDB.

 Switches can be done in long or short form and are case sensitive
 eg:
   cazy_extract.pl --family	{ family or families }
   cazy_extract.pl -f
   cazy_extract.pl -F           { fragments - to include or not }

 Families - you can select more than one family - to select more than one
 you must use this format:

   cazy_extract.pl -f 'GH1,GH2,...'

 Do not forget the single quotes.  Note quotes are not needed for a single
 family.

 Group Names : There are only 7 group names allowed

   all, archaea, bacteria, eukaryota, unclassified, structure, characterized

=head1 ARGUMENTS

 --family    <family_list>   : Family Name, single or list separated by comma in quotes
 --group     <group_name>    : Filename of fasta or fastq file
 --Fragments <boolean>       : true or false - eg. false - keep fragments
 --help                      : print Options and Arguments instead of fetching data
 --man                       : print complete man page instead of fetching data

=head1 Options

 Allowable Group Names:
  all
  archaea
  bacteria
  eukaryota
  unclassified
  structure
  characterized

 Fragment Options:
  true (default)             : Fragments are removed
  false

=head1 AUTHOR

 Dallas Thomas

=head1 TESTED

 Perl    5.10.1
 Debian  5
 Debian  6

=head1 BUGS

 None that I know of

=head1 TODO

 Nothing yet...

=head1 UPDATES

 Nothing yet...

=cut

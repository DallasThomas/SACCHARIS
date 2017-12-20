#!/usr/bin/perl -w 
#
# Script Name: cazy_parse.pl
# Created:     February 16, 2016
#
# Usage:       cazy_parse.pl [options]
#
# Description: This script modifies the header for each fasta sequence
#               extracted from cazy - to remove the chance of duplicate 
#               errors after sequences run through mothur and piped into
#               raxml
#
############################################################################
use strict;
use Data::Dumper;
use Bio::Seq;
use Bio::SeqIO;

# Set-up for Arguments
use Getopt::Long;
use Pod::Usage;
Getopt::Long::Configure ("no_ignore_case");

my ($opt_help, $opt_man);
my ($opt_infile, $opt_outfile);

GetOptions(
    'help!'         => \$opt_help,
    'man!'          => \$opt_man,
    'infile=s'      => \$opt_infile,
    'outfile=s'     => \$opt_outfile,
) or pod2usage(-verbose => 1) && exit;

# Exit program and output Help if...
pod2usage(-verbose => 1) && exit if defined $opt_help;
pod2usage(-verbose => 2) && exit if defined $opt_man;
pod2usage(-verbose => 1) && exit if (!defined $opt_infile || !defined $opt_outfile);

# Create Seq Object(s)
my $in_seq_obj = Bio::SeqIO->new( -file => $opt_infile, -format => "fasta", );
my $out_seq_obj = Bio::SeqIO->new( -file => ">$opt_outfile", -format => "fasta", );

# Variable
my %id_list;
my @id_parts;
my $id_counter = 0;

# Go through each sequence taking the id and parsing it
while ( 1 ) {

  my $idcnt = sprintf("%09d", $id_counter);
  my $seq_obj = $in_seq_obj->next_seq || last;
  my $id = $seq_obj->primary_id;

  # We just need the accession number - no tails
  @id_parts = split(/_/, $id);
  my $new_id;
 
  # Dealing with special case caused by actual underscore in the accession number
  if (exists $id_parts[2]) {
    $new_id = $id_parts[0] . "_" . $id_parts[1];
  } else {
    $new_id = $id_parts[0];
  }

  # Set the id and use counter as key for the hash table 
  if ( $id_list{$idcnt} ) {
    print "This really should not happen - what have you done - contact Dallas\n";
    exit(1);
  } else {
    $id_list{$idcnt} = $new_id;
  }

  # As we are changing the display ID - I still want the old id printed, will append to the description
  my $desc = $seq_obj->display_id . " -- " . $seq_obj->desc;
  
  # Create new object and print out
  my $new_seq_obj = Bio::Seq->new ( -display_id => $idcnt, -primary_id => $seq_obj->primary_id, -desc => $desc, 
		      		    -seq => $seq_obj->seq, -alphabet => $seq_obj->alphabet);	
  &printout($new_seq_obj);
  $id_counter++;
}

# The last phase of the program is to output our temp file of key and id pairs
# Read Identified Sequences from Config File into Hash 
open (TEMP, ">key_id_pairs.txt") or die "Can't open key_id_pairs.txt:$!\n";
foreach (keys %id_list )
{
   print TEMP $_ . "\t" . $id_list{$_} . "\n";
}
close TEMP;
exit;

## Subroutines ##

sub printout {

   # Passed in Variables
   my $seq = shift;

   # Print Out Seqs
   $out_seq_obj->write_seq($seq);
}

__END__


=pod

=head1 NAME

 cazy_parse.pl

=head1 SYNOPSIS

 cazy_parse.pl [-Options [--] [Arguments...]]

=head1 DESCRIPTION

 Parses Headers from a Cazy/NCBI Fasta Files using the ascession numbers. 

 Requires a Cazy/NCBI extracted fasta file as input and the output file name.
 Using the Primary ID from the fasta file created by cazy_extract.pl, this script
 creates another primary ID with a 8-10 character tag appened to the front of the
 regular primary id.  This tag needs to be unique for downstream processing of the 
 cazy pipeline. 

 NOTE: THIS SCRIPT WAS WRITTEN TO WORK EXCLUSIVELY WITH FASTA FILES CREATED BY THE
       CAZY_EXTRACT.PL SCRIPT

 Switches can be done in long or short form 
 eg:
   cazy_parse.pl --infile  { fasta filename }
   cazy_parse.pl -i        { fasta filename}

=head1 ARGUMENTS

 --infile   <fasta_file>    : Filename of Cazy/NCBI Extracted sequences (cazy_extract.pl)
 --outfile  <fasta_file>    : Filename of modified output
 --help                     : print Options and Arguments instead of fetching data
 --man                      : print complete man page instead of fetching data

=head1 Options

 Works with only fasta files extracted using cazy_extract.pl

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

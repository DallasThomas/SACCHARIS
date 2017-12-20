#!/usr/bin/perl -w 
#
# Script Name: prune_seqs.pl
# Created:     May 20, 2015
#
# Usage:       prune_seqs.pl [options]
#
# Description: This script reads in Seq_Ids and start and stop values from 
#               a Debris_list file, or any file of sequence ID's and coord 
#               values whereby the sequence ID's match fasta or fastq sequence
#               ID's.  Using the start and stop locations and length of the
#               sequence, the sequence is pruned at the 5' and 3; ends. The 
#               final pruned sequence is output.
##############
# Update:      July 28th, 2016 - Need to account for mutlimodules of same
#               family.  Some cases a sequence has 2 modules, need to 
#               keep both and modify the identifier name
############################################################################
use strict;
use Data::Dumper;

use Bio::Seq;
use Bio::SeqIO;

use File::Basename;

# Set-up for Arguments
use Getopt::Long;
use Pod::Usage;
Getopt::Long::Configure ("no_ignore_case");

my ($opt_help, $opt_man, $opt_cfg, $opt_fin, $opt_qin, $opt_fmt);

# Parse command line options
my $PROG_NAME = basename($0);

GetOptions(
    'help!'         => \$opt_help,
    'man!'          => \$opt_man,
    'config=s'	    => \$opt_cfg,
    'fasta=s'       => \$opt_fin,
) or pod2usage(-verbose => 1) && exit;

# Exit program and output Help if...
pod2usage(-verbose => 1) && exit if defined $opt_help;
pod2usage(-verbose => 2) && exit if defined $opt_man;
pod2usage(-verbose => 1) && exit if !defined $opt_cfg && !defined $opt_fin;

# Create input objects and Output filenames
my ($fbase, $qbase, $in_seq_obj, $in_qual_obj);
$in_seq_obj = Bio::SeqIO->new( -file => $opt_fin, -format => "fasta", );
$fbase = basename($opt_fin);
$fbase =~ s/\.[^.]+$//;
 
# Create output objects for a seq file
my $out_seq_obj;
$out_seq_obj = Bio::SeqIO->new( -file => ">$fbase.pruned.fasta", -format => "fasta", );

# Read Identified Sequences and Coordinate values from Config File into Hash 
open (BEG, "<$opt_cfg") or die "Can't open $opt_cfg:$!\n";
my %BEGIN;
while ( <BEG> )
{
   chomp;
   my @list = split(' ', $_);
   my $keyB = $list[0];
   my $coord = $list[1] . ":" . $list[2];

   # Need to account for potential duplicate hash keys
   push @{$BEGIN{$keyB}}, $coord;
}
close BEG;

# Now go through Fasta file given, extract matching sequence, write to temp file,
#  run through cutadapt and then output to final output file
while ( 1 ) {

    my $seq_obj;
    $seq_obj  = $in_seq_obj->next_seq || last;

    # If id's match hash - printout sequence
    if ($BEGIN{$seq_obj->id}) {

	# Need the sequence - to prune 
	my $sequence = $seq_obj->seq;
	my $sequence_length = $seq_obj->length;

	# Prune the sequence and set to seq_obj->seq - then printout
	my @cd = @{$BEGIN{$seq_obj->id}};

        # There might be more than one set of coordinates so we need to deal
        #  with this possibility
	my $arrsize = @cd;
	my $iter = 1;

	if ($arrsize > 1) {
	   foreach (@cd) {
	     my $newid = &parseid($seq_obj->display_id, $iter);
	     $iter++;
	     &process($seq_obj, $sequence, $sequence_length, $_, $newid);
	   }
	} else {
	   &process($seq_obj, $sequence, $sequence_length, $cd[0], $seq_obj->display_id);
        }
    } 
}

exit;

## Subroutines ##

sub parseid {

  # Passed in Variables
  my ($id, $val) = @_;

  # Need to equate the val with a letter of the alphabet
  my $letter = chr(96+$val);

  return $id . $letter;
}

sub process {

   # Passed in Variables
   my ($s_obj, $seq, $seq_len, $coordinate, $id) = @_;

   # Process this stuff
   my @coord = split(':', $coordinate);
   my @seqarray = split('', $seq);
   my $start = $coord[0] - 1;
   my $stop = $coord[1] - 1;
   my $new_seq;
        
   if ($seq_len < $coord[1]) {
      print "There was an error - in coordinates for $s_obj->id, moving to next seq\n";
      next;
   }

   foreach my $i ($start..$stop) {
      $new_seq .= $seqarray[$i];
   }

   my $new_seq_obj = Bio::Seq->new ( -display_id => $id, -primary_id => $id, -desc => $s_obj->desc,
                                     -seq => $new_seq, -alphabet => $s_obj->alphabet);  

   &printout($new_seq_obj);
}

sub printout {

   # Passed in Variables
   my $seq = shift;

   # Print Out Seqs
   $out_seq_obj->write_seq($seq);

}

__END__


=pod

=head1 NAME

 prune_seqs.pl

=head1 SYNOPSIS

 prune_seqs.pl [-Options [--] [Arguments...]]

=head1 DESCRIPTION

 Prunes Sequence Data from Fasta Files Based on Start and Stop positions. 

 Requires a file of Sequence ID's and coordinates for which to extract and prune.  
 Uses this list of Sequence ID's to extract the Sequences from a Fast{a|q} file
 whereby it uses the coordinate information to chop bases from the 5' and 3' ends
 outputting the pruned sequence to the output file.  Needs proper format to know 
 which files are required and which are not.

 Switches can be done in long or short form and are case sensitive
 eg:
   prune_seqs.pl --config
   prune_seqs.pl -c
   prune_seqs.pl -f     { fasta filename}

=head1 ARGUMENTS

 --config   <config_file>   : Filename of Sequence ID list
 --fasta    <seq_file>      : Filename of fasta file
 --help                     : print Options and Arguments instead of fetching data
 --man                      : print complete man page instead of fetching data

=head1 Options

 Works with only fasta files at this time

 Config File Must Be formatted as such:

	 	<Sequence_ID>	<Start>	<Stop>

	Eg.

		gi|10336504|dbj|BAB13699.2|	468	595
		gi|10434341|dbj|BAB14227.1|	431	554
		gi|10435819|dbj|BAB14676.1|	197	335
		...

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

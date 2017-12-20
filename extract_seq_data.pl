#!/usr/bin/perl -w 
#
# Script Name: extract_seq_data.pl
# Created:     February 21, 2012
#
# Usage:       extract_seq_data.pl [options]
#
# Description: This script reads in Seq_Ids from a Debris_list file, or any
#               file of sequence ID's that match fasta or fastq sequence
#               ID's and outputs the matching sequence data in the format
#               specified.
############################################################################
# Modified:    May 30th, 2017 - Added the ability to have a user defined
#                output filename
############################################################################
use strict;
use Data::Dumper;

use Bio::SeqIO;
use Bio::Seq;

use File::Basename;

# Set-up for Arguments
use Getopt::Long;
use Pod::Usage;
Getopt::Long::Configure ("no_ignore_case");

my ($opt_help, $opt_man, $opt_cfg, $opt_fin, $opt_qin, $opt_out, $opt_fmt);

# Parse command line options
my $PROG_NAME = basename($0);

GetOptions(
    'help!'         => \$opt_help,
    'man!'          => \$opt_man,
    'config=s'	    => \$opt_cfg,
    'fastaq=s'      => \$opt_fin,
    'outfile=s'     => \$opt_out,
    'qual=s'        => \$opt_qin,
    'Format=s'      => \$opt_fmt,
) or pod2usage(-verbose => 1) && exit;

# Exit program and output Help if...
pod2usage(-verbose => 1) && exit if defined $opt_help;
pod2usage(-verbose => 2) && exit if defined $opt_man;
pod2usage(-verbose => 1) && exit if !defined $opt_cfg && !defined $opt_fin;

# Set-up variables based on Arguments Passed In
my $format =  $opt_fmt ? $opt_fmt : "fastq";

# Create input objects and Output filenames
my ($fbase, $qbase, $in_seq_obj, $in_qual_obj);
$in_seq_obj = Bio::SeqIO->new( -file => $opt_fin, -format => $format, );
$fbase = basename($opt_fin);
$fbase =~ s/\.[^.]+$//;

if ($format eq "fasta" && defined $opt_qin) {
   $in_qual_obj = Bio::SeqIO->new( -file => $opt_qin, -format => 'qual', );
   $qbase = basename($opt_qin);
   $qbase =~ s/\.[^.]+$//;
} 

# Create output objects for both a seq and qual file
my $outfilename = $opt_out ? $opt_out : "$fbase.mod";
my $outqual;

if (($format eq "fastq") && (!defined $opt_out)) {
   $outfilename = "$outfilename.fastq";
} elsif (!defined $opt_out) {
   $outfilename = "$outfilename.fasta";
   $outqual = "$outfilename.qual";
}

my ($out_seq_obj, $out_qual_obj);
if (!defined $qbase) {
   $out_seq_obj = Bio::SeqIO->new( -file => ">$outfilename", -format => $format, );
} else {
   $out_seq_obj = Bio::SeqIO->new( -file => ">$outfilename", -format => $format, );
   $out_qual_obj = Bio::SeqIO->new( -file => ">$outqual", -format => 'qual', -width  => 22, );
}


# Read Identified Sequences from Config File into Hash 
open (BEG, "<$opt_cfg") or die "Can't open $opt_cfg:$!\n";
my %BEGIN;
while ( <BEG> )
{
   chomp;
   my $keyB = $_;
   $BEGIN{$keyB} = ( 1 );
}
close BEG;

# Now go through Fast{a|q} file and qual file if given and extract matching sequences
my ($seq_obj, $qual_obj);

while ( 1 ) {

    my ($seq_obj, $qual_obj);
    $seq_obj  = $in_seq_obj->next_seq || last;

    if ($format eq "fasta" && defined $opt_qin) {
       $qual_obj = $in_qual_obj->next_seq;
       die "Id's don't match!\n" unless $seq_obj->id eq $qual_obj->id;
    }

    # If id's match hash - printout sequence
    if ($BEGIN{$seq_obj->id}) {
	&printout($seq_obj, $qual_obj);
    } 
}

exit;

## Subroutines ##

sub printout {

   # Passed in Variables
   my ($seq, $qual) = @_;

   # Print Out Seqs
   $out_seq_obj->write_seq($seq);

   if (defined $qual) {
      $out_qual_obj->write_seq($qual)
   }
}

__END__


=pod

=head1 NAME

 extract_seq_data.pl

=head1 SYNOPSIS

 extract_seq_data.pl [-Options [--] [Arguments...]]

=head1 DESCRIPTION

 Extracts Sequence and Quality Data from Fast{a|q} and Quality Files. 

 Requires a file of Sequence ID's for which to extract.  Uses this list
 of Sequence ID's to extract Sequences from a Fast{a|q} file to the output
 file.  Needs proper format to know which files are required and which are
 not. See below for proper formats.

 Switches can be done in long or short form and are case sensitive
 eg:
   extract_seq_data.pl --config
   extract_seq_data.pl -c
   extract_seq_data.pl -f     { fastaq filename}
   extract_seq_data.pl -F     { Format }

=head1 ARGUMENTS

 --config   <config_file>   : Filename of Sequence ID list
 --fastaq   <seq_file>      : Filename of fasta or fastq file
 --qual     <qual_file>     : Filename of quality file {optional}
 --outfile  <out_file>      : Filename for the output file - User selection
 --Format   <format>        : Format
 --help                     : print Options and Arguments instead of fetching data
 --man                      : print complete man page instead of fetching data

=head1 Options

 Allowable Formats:
  fasta
  fastq                   : Sanger variant {fastq-sanger} (Default)
  fastq-solexa            : Illumina 1.0
  fastq-illumina          : Illumina 1.3

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

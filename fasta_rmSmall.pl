#!/usr/bin/perl -w 
#
# Script Name: fasta_rmSmall.pl
# Created:     February 13, 2018
#
# Usage:       fasta_rmSmall.pl [options]
#
# Description: This script takes in minsize value for sequence and removes 
#               sequences that do not meet this minsize from a fasta file 
############################################################################
use strict;
use Data::Dumper;

use File::Basename;

# Set-up for Arguments
use Getopt::Long;
use Pod::Usage;
Getopt::Long::Configure ("no_ignore_case");

my ($opt_help, $opt_man, $opt_size, $opt_fasta);

# Parse command line options
my $PROG_NAME = basename($0);

GetOptions(
    'help!'         => \$opt_help,
    'man!'          => \$opt_man,
    'size=s'	    => \$opt_size,
    'fasta=s'       => \$opt_fasta,
) or pod2usage(-verbose => 1) && exit;

# Exit program and output Help if...
pod2usage(-verbose => 1) && exit if defined $opt_help;
pod2usage(-verbose => 2) && exit if defined $opt_man;
pod2usage(-verbose => 1) && exit if !defined $opt_size && !defined $opt_fasta;

# Create input objects and Output filenames
my $fbase;
$fbase = basename($opt_fasta);
$fbase =~ s/\.[^.]+$//;
 
# Create output filename
my $output_file = $fbase . ".pruned.fasta";

# Read in file and screen for size - output to file
open (INFILE, "<$opt_fasta") or die "Can't open $opt_fasta:$!\n";
open (OUTFILE, ">>$output_file") or die "Can't open $output_file:$!\n";

# Changing scope as using local to change vars
{
   # Want to stop reading at > so we can basically read an entire sequence into buffer at a time
   #  Note: Using local to change $/ can be dangerous as this resets $/ from \n, will need to change back later
   local $/ = ">";

   # Loop through fasta file and select sequences that meet our criteria
   while (<INFILE>) {
      chomp;
      next unless /\w/;         # Looking for white space or not
      s/>$//gs;                 # This is an in-case substitution condition - matches from > to end of string found in $_, globally including \n 
      my @sequence = split /\n/;
      my $header = shift @sequence;
      my $seq_length = length join "", @sequence;

      if ( $seq_length >= $opt_size ) { 
	print OUTFILE ">$_";    # Because we changed $/ above $_ stores the entire sequence
      } 
   }
   # Change $/ back
   local $/ = "\n";
}


__END__


=pod

=head1 NAME

 fasta_rmSmall.pl

=head1 SYNOPSIS

 fasta_rmSmall.pl [-Options [--] [Arguments...]]

=head1 DESCRIPTION

 Screens Fasta files and removes Sequence Data depending on Minimum Length Specified. 

 Requires for input a fasta file and the minimum sequence length.  Reads in sequences  
 into buffer one at a time and screens based on length of sequence, writing output to
 an output file.  Needs both a fasta file and a length cutoff.

 Switches can be done in long or short form and are case sensitive
 eg:
   fasta_rmSmall.pl --size
   fasta_rmSmall.pl -s
   fasta_rmSmall.pl -f     { fasta filename}

=head1 ARGUMENTS

 --size     <min_length>    : Minimum Sequence Length for Prunning
 --fasta    <seq_file>      : Filename of fasta file
 --help                     : print Options and Arguments instead of fetching data
 --man                      : print complete man page instead of fetching data

=head1 Options

 Works with only fasta files at this time

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

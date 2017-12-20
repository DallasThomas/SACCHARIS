#!/usr/bin/perl -w 
#
# Script Name: muscle_rename.pl
# Created:     November 7, 2016
#
# Usage:       muscle_rename.pl [options]
#
# Description: This script switches the modified header names in the 
#               phylip output file to proper accession numbers that were 
#               parsed during the cazy_parse phase.  Script makes output
#               more readable
#
############################################################################
use strict;
use Data::Dumper;
use File::Basename;

# Set-up for Arguments
use Getopt::Long;
use Pod::Usage;
Getopt::Long::Configure ("no_ignore_case");

my ($opt_help, $opt_man);
my ($opt_infile, $opt_keyfile);

GetOptions(
    'help!'         => \$opt_help,
    'man!'          => \$opt_man,
    'infile=s'      => \$opt_infile,
    'keyfile=s'     => \$opt_keyfile,
) or pod2usage(-verbose => 1) && exit;

# Exit program and output Help if...
pod2usage(-verbose => 1) && exit if defined $opt_help;
pod2usage(-verbose => 2) && exit if defined $opt_man;
pod2usage(-verbose => 1) && exit if (!defined $opt_infile || !defined $opt_keyfile);

# Read in the Key and ID pairs
open (BEG, "<$opt_keyfile") or die "Can't open $opt_keyfile:$!\n";
my %BEGIN;
while ( <BEG> )
{
   chomp;
   my $keyB = $_;
   my @keys = split(/\t/, $keyB);
   $BEGIN{$keys[0]} = ( $keys[1] );
}
close BEG;

# Create the new oufile name
my $fbase = basename($opt_infile);
$fbase =~ s/\.[^.]+$//;
$fbase = $fbase . "_mod.phyi";

# Read in from the Infile and use the key data to swap id information
open (IN, "<$opt_infile") or die "Can't open $opt_infile:$!\n";
open (OUT, ">>$fbase") or die "Can't open $fbase:$!\n";
my $iter = 0;

while (<IN>) {

  # Need to skip the first line
  if ($iter == 0) {
    $iter++;
    print OUT $_;
    next;
  }

  # Look for lines that match the keys, if they match substitute value for key and print out
  #  else just print line as is
  if ($_ =~ m/^\d+/) {
    # First need to get the digit to match
    my @list = split(/ /, $_);

    # Need to account for multi-modularity - ids have an (a,b,...) postfix
    if ($list[0] =~ m/^(\d+)([a-z])/) {
       my $val = $BEGIN{ $1 };
       $list[0] = $val . $2;
    } else { 
       my $val = $BEGIN{ $list[0] };
       $list[0] = $val;
    }
    my $line = join(' ', @list);
    print OUT $line;
  } else {
    print OUT $_;
  }
}
close IN;
close OUT;
exit;

__END__


=pod

=head1 NAME

 muscle_rename.pl

=head1 SYNOPSIS

 muscle_rename.pl [-Options [--] [Arguments...]]

=head1 DESCRIPTION

 Uses the Key-Value Pairs from cazy_parse to replace unique headers in the muscle
 phylip output file with headers that are more readable. 

 Requires a muscle output file in phylip format and a key-value-id list generated 
 by cazy_parse.pl.  This script reads through the phylip file and replacess the 
 unique header/id or key with the value to create a more legible file prior to 
 submission to the tree building phase.

 NOTE: THIS SCRIPT WAS WRITTEN TO WORK EXCLUSIVELY WITH PHYLIP FILES AND OUTPUT
       FROM THE CAZY_PARSE.PL SCRIPT

 Switches can be done in long or short form 
 eg:
   muscle_rename.pl --infile  { fasta filename }
   muscle_rename.pl -i        { fasta filename}

=head1 ARGUMENTS

 --infile   <fasta_file>    : Filename of Cazy/NCBI Extracted sequences (cazy_extract.pl)
 --keyfile  <key-val_file>  : Filename of the key value pairs list from cazy_parse.pl
 --help                     : print Options and Arguments instead of fetching data
 --man                      : print complete man page instead of fetching data

=head1 Options

 Works with Phylip files and a key-value id list file similar to cazy_parse.pl output

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

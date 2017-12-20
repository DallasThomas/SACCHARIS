#!/usr/bin/perl

my $pgm = $0;                   # name of program
$pgm =~ s#.*/##;                # remove part up to last slash
my $seed = 1;			# random number seed
my $line_size = 50;		# size of lines to print
my $off = 1;			# first position of sequence to print
my $len = -1;			# maximum length of sequence to print
my $rest = "";			# file to receive remainder of seqs
 
$usage = <<USAGE;               # usage message
  USAGE:
        $pgm <fasta> <n> [-seed <seed>] [-rest <rest>] 
		[-off <off>] [-len <len>]

		<fasta>		name of FASTA sequence file
		<n>		number of sequences to output
		[-seed <seed>]	random number seed; default: $seed
		[-rest <rest>]	name of file to receive the FASTA
				sequences not being output; default: none
		[-off <off>] 	print starting at position <off> in each
				sequence; default: $off
 		[-len <len>]	print up to <len> characters for each 
				sequence; default: print entire sequence

	Output a random subsample of size <n> of the sequences in 
	a FASTA sequence file.  The seed of the random generator can be 
	changed using -seed, otherwise the same subset of sequences will 
	always be output.  If requested, the remaining sequences will 
	be output to a file named <rest>, which is useful for 
	cross-validation.

	You can also choose to only output portions of each sequence
	using the -off and -len switches.  If the sequences have 
	UCSC BED file headers (e.g., ">chr1:0-99999"), the headers will
	be adjusted to reflect -off and -len.

	Writes to standard output.
USAGE

if ($#ARGV+1 < 2) {             # wrong number of arguments
  die $usage;
}

# get input arguments
my $fasta = shift;		# name of fasta file
my $n = shift;			# size of subsample
while ($#ARGV >= 0) {
  $_ = shift;
  if ($_ eq "-seed") {		# random seed
    $seed = shift;
  } elsif ($_ eq "-rest") {	# rest file name
    $rest = shift;
  } elsif ($_ eq "-off") {	# offset
    $off = shift;
  } elsif ($_ eq "-len") {	# maximum length
    $len = shift;
  } elsif ($_ eq "-rest") {	# file for remainder
    $rest = shift;
  } else {
    die $usage;
  }
}

# read in FASTA file and make index
open(FASTA, "<$fasta") || die "Couldn't open file `$fasta'.\n";
my $byte = 0;
my %index;			# ID-to-start index
my $id;				# sequence ID
my @rest;			# dummy
my @id_list;			# list of all IDs
while (<FASTA>) {
  if (/^>/) {
    ($id, @rest) = split;
    $index{$id} = $byte; # start of sequence record
    push @id_list, $id;
  } 
  $byte += length;
} # read FASTA file

# check that there are enough IDs
my $nseqs = @id_list;
die ("Not enough sequences ($nseqs); $n requested.\n") if ($nseqs < $n);

# shuffle the list of IDs
srand($seed);
shuffle(\@id_list);
#print join " ", @id_list, "\n";

# output the requested number of FASTA sequences to STDOUT
foreach $id (@id_list[0..$n-1]) { 
  print_fasta_seq_portion(*FASTA, *STDOUT, $id, $off, $len, \%index);
} # id

# output the remainder of the sequences if requested
# to the "rest" file
if ($rest) {
  open(REST, ">$rest") || die("Can't open file `$rest'.\n");
  foreach $id (@id_list[$n..$nseqs-1]) { 
    print_fasta_seq_portion(*FASTA, *REST, $id, $off, $len, \%index);
  }
}

################################################################################
# print (a portion of) a FASTA sequence
# Assumes FASTA file is open and the index contains
# the file byte offset for a given ID.
sub print_fasta_seq_portion {
  my ($fasta, $output, $id, $off, $len, $index) = @_;

  my $addr = $index{$id};		# address of sequence 
  die "Can't find target $id.\n" if ($addr eq undef);
  seek($fasta, $addr, 0);		# move to start of target

  # save ID line
  $id_line = <$fasta>;

  my $seq = "";
  # read in sequence lines for this sequence
  while (<$fasta>) {			# read sequence lines
    if (/^>/) {last}			# start of next sequence
    chop;
    $seq .= $_;
  }

  # get length of sequence
  $length = length($seq);

  # print ID of FASTA sequence
  $_ = $id_line;
  if (/^(>chr[\dXY]+):(\d+)-(\d+)/i) {	# handle BED format
    $chr = $1;
    $start = $2;
    $end = $3;
    $comment = $'; 
    $start = $start + ($off - 1);	# new start; 0-based
    # BED end is really "end+1"
    $end = ($len != -1) ? $start + $len : $start + $length;
    # print ID for the sequence in BED format, adjusting for offset and length
    printf($output "%s:%s-%s%s", $chr, $start, $end, $comment);
  } else {				# handle other formats
    printf($output $_);			# print ID for this sequence
  }

  # print sequence in lines of length $line_size
  # get portion of sequence to print if -off and/or -len given
  if ($off != 1 || $len != -1) {
    if ($len == -1) {
      $seq = substr($seq, $off-1);
    } else {
      $seq = substr($seq, $off-1, $len);
    }
  }
  for ($i=0; $i<length($seq); $i+=$line_size) {
    print $output substr($seq, $i, $line_size), "\n";
  }
} # print_fasta_seq_portion

# shuffle a list in place
sub shuffle (\@) { 
    my $r=pop; 
    $a = $_ + rand @{$r} - $_ 
      and @$r[$_, $a] = @$r[$a, $_] 
        for (0..$#{$r}); 
}

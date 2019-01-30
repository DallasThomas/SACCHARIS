#!/usr/bin/perl -w 
#
# Script Name: Saccharis.pl
# Created:     January 28, 2016
#
# Usage:       Saccharis.pl [options]
#
# Description: This script is a wrapper script for a Cazy pipeline.. 
#
#
############################################################################
use strict;
use Data::Dumper;
use File::Basename;
use File::Copy;
use File::chdir;
use threads;
use Benchmark;
use POSIX ;
use Switch;
#use Cwd;
use List::Util qw( reduce );
use Bio::SeqIO;

# Set Date
use Date::Calc qw/Date_to_Text Date_to_Text_Long Today/;
my @date_today = Today();

# Set-up for Arguments
use Getopt::Long;
use Pod::Usage;
Getopt::Long::Configure ("no_ignore_case");

my ($opt_help, $opt_man);
my ($opt_homedir, $opt_threads, $opt_grp, $opt_fam, $opt_frag, $opt_rax, $opt_userseq);
my ($opt_RAX, $opt_FAST);
my $FRAGS;

# Parse command line options
my $PROG_NAME = basename($0);
my $ARGS = join(" ", @ARGV);

GetOptions(
    'help!'         => \$opt_help,
    'man!'          => \$opt_man,
    'directory=s'   => \$opt_homedir,
    'threads=s'     => \$opt_threads,
    'group=s'       => \$opt_grp,
    'family=s'      => \$opt_fam,
    'Fragments!'    => \$opt_frag,
    'RAXML!'        => \$opt_RAX,
    'raxml=s'       => \$opt_rax,
    'seqfile=s'     => \$opt_userseq,
) or pod2usage(-verbose => 1) && exit;

# Exit program and output Help if...
pod2usage(-verbose => 1) && exit if defined $opt_help;
pod2usage(-verbose => 2) && exit if defined $opt_man;
pod2usage(-verbose => 1) && exit if (!defined $opt_grp || !defined $opt_fam);

# Need to validate if UserSeqs are included they are in proper fasta format
if (defined $opt_userseq) {
  my $seqio = Bio::SeqIO->new(-file => $opt_userseq, -format => 'fasta');
  while (my $nextseq = $seqio->next_seq) { 
    next;
  }
}
my $userfile = $opt_userseq ? $opt_userseq : undef;

# Set the Tree Building Program - Default is FastTree -> 1
my $TreeBuild = 1;

if (defined $opt_RAX) {
  $TreeBuild = 0;
}

# Set Default Variables
my $home = $opt_homedir ? $opt_homedir : getcwd();
my $threads = $opt_threads ? $opt_threads : 2;
my $rax = $opt_rax ? $opt_rax : 'raxmlHPC-PTHREADS-SSE3';

# Need to know the number of sequences in the fasta file
my $NumberSeqs;

# Global Variable
my $prot_muscle;

# Fragments - if Flag is given - implies that fragments are to be included
if ( defined $opt_frag ) { $FRAGS = 'false'; } else { $FRAGS = 'true'; }

# Allowing for multiple groups per family
my @group = split(/,/,$opt_grp);

# Screen output - reads different based on input flags
print "==============================================================================\n\n";
print "Date: ", Date_to_Text(@date_today), "\n";
print "\nStarting Cazy Pipeline using the following options...\n\n";
print "\t - Family: \t" . $opt_fam . "\n";
print "\t - Group(s): \t" , $opt_grp . "\n";
print "\t - Threads: \t" . $threads . "\n";
print "\t - HomeDir: \t" . $home . "\n";
my $uf = $userfile ? $userfile : "...";
print "\t - UserFile: \t" . $uf . "\n";
if ( defined $opt_RAX ) {
   print "\t - RaxML: \t" . $rax . "\n";
} else {
   print "\t - Tree Generation: \tFastTree\n";
}
if ( defined $opt_frag ) { 
   print "\t - Fragments are included...\n";
} else {
   print "\t - Fragments are not included...\n"; 
}
print "\n";
print "Program Command : $PROG_NAME $ARGS \n\n";
print "==============================================================================\n\n";

###########################################################################################
# Starting of the pipeline - this of course is looped based on number of groups
###########################################################################################
my $fam_dir = $home . "/" . $opt_fam;
unless ( -d $fam_dir ) { mkdir($fam_dir, 0755); }
$CWD = $fam_dir;

foreach my $i (@group) {

  # Create a Directory for each specific group
  my $grp_dir = $fam_dir . "/" . $i;
  unless( -d $grp_dir ) { mkdir($grp_dir, 0755); }

  # Set current working directory
  local $CWD = $grp_dir;

  # Proceed with the Pipeline
  print "Begin pipeline analysis for group: $i of family $opt_fam\n";

  #######################################
  # Step One - Cazy Extract
  #######################################
  
  # FileName Delcarations
  my ($cazy_ext_file, $cazy_file);

  # Start of Benchmarks
  my $t0 = Benchmark->new;

  # Call Cazy Extract Subroutine
  print "Cazy Extract is proceeding for $i of family $opt_fam\n";
  $cazy_ext_file = &cazy_extract($opt_fam, $i, $FRAGS, $grp_dir);
  print "Completed Cazy Extract\n\n";
  print "==============================================================================\n\n";

  # Set the Sequence Total Variable - read this from the sequence.count file in Cazy Folder
  my ($TSEQCOUNT, $FRAGMENTS, $DUPLICATES, $MSEQCOUNT);
  my @CNTDATA;
  open (SQ, "<$CWD/cazy/sequence.count") or die "Cannnot open $CWD/cazy/sequence.count for reading: $!";
  while (<SQ>) {
    chomp $_;
    @CNTDATA = split(/,/, $_);
    $TSEQCOUNT = $CNTDATA[0];
    $FRAGMENTS = $CNTDATA[2];
    $DUPLICATES = $CNTDATA[1];
    $MSEQCOUNT = $CNTDATA[3]; 
  }
  close(SQ);
  print "\n";
  print "Notice\n\n";
  print "\tThe following Duplicate and Modified Count may not be completely accurate, it is recommended you,\n";
  print "\t confirm the Modified count with the total number of sequences in the extracted fasta file.\n\n";
  print "\t\t***Total Sequence Count     -> $TSEQCOUNT\n";
  print "\t\t***Total Duplicate Count    -> $DUPLICATES\n";
  print "\t\t***Total Fragmment Count    -> $FRAGMENTS\n";
  print "\t\t***Modified Sequence Count  -> $MSEQCOUNT\n";  

  # Validate Cazy Sequences match Total Count
  &validate_seqfile($cazy_ext_file, $TSEQCOUNT);

  # Call the Cazy Parser to setup fasta header lines to prevent identical hits
  print "Cazy parsing is proceeding - need to prevent duplicate headers\n";
  $cazy_file = &cazy_parse($cazy_ext_file, $grp_dir);
  print "Completed Cazy Parse\n\n";
  print "==============================================================================\n\n";

  # Validate Parsed Sequences match Total Count
  &validate_seqfile($cazy_file, $TSEQCOUNT);

  # Run First Benchmark
  my $t_extract = Benchmark->new;
  my $td = timediff($t_extract, $t0);
  print "*********************************************\n";
  print "* Cazy Extraction Takes\n";
  print "*      --> ", timestr($td), " to run\n";
  print "*********************************************\n";

  #######################################
  # Step Two - Combine User and Cazy
  #######################################

  # Need a way to prevent user files from being added if Saccharis is restarted
  my $user_validation = "user.data.check";

  if ( -f $user_validation ) {
    print "User Data has already been added to Cazy data, continuing with script\n";
  } else {
    if (defined $opt_userseq) {
      print "User has defined a Fasta file, which will be appeneded to the Cazy Extract file\n";
      print "Need to first determine if the User Fasta file has a proper header\n";
    
      # Need to determine the validity of the user file
      my $user_file = "grep ^\\> $userfile \| awk '{if (\$1 !~ \/>U[0-9]{1,9}\$\/) exit(1)}'";
      system($user_file);

      if ($? == 0) {
         print "User File is valid format - proceeding...\n\n";
         print "Appending $userfile to $cazy_file\n";
         open ( CAZY, ">>", $cazy_file ) or die "Could not open file $cazy_file: $!";
         open ( USER, "<", $userfile ) or die "Could not open file $userfile: $!";

         while ( my $line = <USER> ) { print CAZY $line; }
         close(CAZY);
         close(USER);
         print "Files have been merged\n\n";
         print "==============================================================================\n\n";
      } else {
         print "User File Format is Invalid - Exiting - Please see Manual for correct formatting\n\n";
         exit(1);
      }
    }
    # Create the User Validation File and enter whatever
    open ( USERV, ">", $user_validation ) or die "Could not open file $user_validation: $!";
    print USERV "User data added - Validated";
    close(USERV);
  } 

  # Set NumberSeq count for all sequences
  $NumberSeqs = `grep -c ^\\> $cazy_file`;
  print "Final Sequence Count prior to analysis --> \t$NumberSeqs\n\n";

  # Run Second Benchmark
  my $t_merge = Benchmark->new;
  $td = timediff($t_merge, $t_extract);
  print "*********************************************\n";
  print "* Fasta file Concatenation Takes\n";
  print "*      --> ", timestr($td), " to run\n";
  print "*********************************************\n";

  #######################################
  # Step Three - dbCAN, extract & prune
  #######################################

  print "dbCAN processing of $cazy_file is underway\n";
  my $dbcan_file = &dbcan($cazy_file, $opt_fam, $grp_dir);
  print "Completed dbCAN Processing\n\n";
  print "==============================================================================\n\n";

  # Set NumberSeq count to pruned sequences
  $NumberSeqs = `grep -c ^\\> $dbcan_file`;
  print "Pruned Sequence Count --> \t$NumberSeqs\n\n";

  # Run Third Benchmark
  my $t_dbcan = Benchmark->new;
  $td = timediff($t_dbcan, $t_merge);
  print "*********************************************\n";
  print "* dbCAN Processing Takes\n";
  print "*      --> ", timestr($td), " to run\n";
  print "*********************************************\n";

  #######################################
  # Step Four - Muscle
  #######################################

  # FastTree variant file
  my $MFast = $grp_dir . "/muscle/" . $opt_fam . ".muscle_aln_mod_fast.phyi";

  print "Muscle alignment of $dbcan_file is underway\n";
  my $muscle_file = &muscle($dbcan_file, $opt_fam, $grp_dir, $MFast, 0);
  print "Completed Muscle Alignment\n\n";
  print "==============================================================================\n\n";

  # Run Fourth Benchmark
  my $t_muscle = Benchmark->new;
  $td = timediff($t_muscle, $t_dbcan);
  print "*********************************************\n";
  print "* Muscle Alignment Takes\n";
  print "*      --> ", timestr($td), " to run\n";
  print "*********************************************\n";

  #######################################
  # Step Five - Prottest3
  #######################################

  # Prottest does not like phylip files with headers greater than 10
  print "Prottest3 tree modeling of $prot_muscle is underway\n";
  my $TreeModel = &prottest($prot_muscle, $cazy_file, $dbcan_file, $opt_fam, $grp_dir, $threads, $MFast);
  print "Best model found via Prottest\n\n";
  print "==============================================================================\n\n";

  # Run Five Benchmark
  my $t_prot = Benchmark->new;
  $td = timediff($t_prot, $t_muscle);
  print "*********************************************\n";
  print "* Prottest3 Tree Modeling Takes\n";
  print "*      --> ", timestr($td), " to run\n";
  print "*********************************************\n";
  
  #######################################
  # Step Six - Build Tree
  #######################################

  # Adding FastTree to the Mix of Programs here - this will change...
  my $best_tree;

  if ( $TreeBuild == 0 ) {
    print "RaxML - Tree building of $muscle_file is underway\n";
    $best_tree = &raxml($muscle_file, $TreeModel, $opt_fam, $grp_dir, $rax, $threads);
  } else {
    print "FastTree - Tree building of $MFast is underway\n";
    $best_tree = &fasttree($MFast, $TreeModel, $opt_fam, $grp_dir);
  }
  print "Completed Building of Tree\n\n";
  print "==============================================================================\n\n";

  # Run Last Benchmark
  my $t_rax = Benchmark->new;
  $td = timediff($t_rax, $t_prot);
  print "*********************************************\n";
  print "* Tree Building Takes\n";
  print "*      --> ", timestr($td), " to run\n";
  print "*********************************************\n";

  # Copy the Best tree to the current working directory
  my $final_tree_file = $opt_fam . "_" . $i . ".tree";
  copy($best_tree, $final_tree_file);

  # Final Benchmark tests
  my $t1 = Benchmark->new;
  $td = timediff($t1, $t0);
  print "*********************************************\n";
  print "* Cazy Pipeline Took in Total\n";
  print "*      --> ", timestr($td), " to finish\n";
  print "*********************************************\n\n";

  print "Finished Cazy pipeline analysis for group: $i of family $opt_fam\n";
  print "==============================================================================\n";
  print "==============================================================================\n\n";
}

print "Cazy Pipeline Finished\n\n";

exit;

###############################################################################
# Cazy Extract - Using the Family and Group Names Extracts the Fasta sequence 
#                data the family in that group from the Cazy Website
# Notes:  Fragment removal is default in the perl code so will only pass the  
#          Fragment flag if they are to be kept in
###############################################################################
sub cazy_extract {

  # Passed in Variables
  my ($fam, $grp, $frag, $dir) = @_;

  # Varible Declaration
  my ($cmd1, $cmd2) = ("", "");
  my $cazybase = $fam . "_cazy.fasta";

  # Create Directory for Cazy and change to this directory
  my $local = $dir . "/cazy";
  unless ( -d $local ) { mkdir($local, 0755); }
  local $CWD = $local;

  # Run Cazy Extract if file does not already exist
  print "Running Cazy Extract Script\n";
  if ( -f $cazybase ) {
    print "Cazy Extract has already been run, continuing with script\n";
  } else {
    if ($frag eq 'false') {
      $cmd1 = "cazy_extract.pl -f $fam -g $grp -F $frag; ";
    } else {   
      $cmd1 = "cazy_extract.pl -f $fam -g $grp; ";
    }
    &run_cmd($cmd1, $cmd2);
    print "Cazy Sequences have been extracted\n\n";
  }

  # Create the name of the cazy output
  my $fullcazy = $local . "/" . $cazybase;

  return $fullcazy;
}

###############################################################################
# Cazy Parse - Cazy sometimes uses an 'a,b' variant in some of the accession 
#               numbers in the fasta headers - which kind of creates duplicates
#               when muscle cuts headers to 10 characters.  Need to fix this
###############################################################################
sub cazy_parse {

  # Passed in Variables
  my ($file, $dir) = @_;

  # Varible Declaration
  my ($cmd1, $cmd2) = ("", "");

  # Create Directory for Cazy and change to this directory
  my $local = $dir . "/cazy";
  local $CWD = $local;

  # Need to create the output name
  my $basename = basename($file);
  $basename =~ s/\.[^.]+$//;
  my $cazyparse = $basename . ".mod.fasta";

  # Run Cazy Parse
  if ( -f $cazyparse ) {
    print "Cazy Parse has already been run, continuing with script\n";
  } else {
    print "Running Cazy Parse \n";
    $cmd1 = "cazy_parse.pl -i $file -o $cazyparse; ";
    &run_cmd($cmd1, $cmd2);
    print "Cazy sequence headers have been modified\n\n";
  }

  return $local . "/" . $cazyparse;
}

###############################################################################
# Validate Seqfile - Before proceeding we need to determine whether the cazy 
#                     extract and parse subroutines had output sequence files
#                     with a sequence count withing 10% of the total sequence
#                     count obtained from Cazy - best to figure this out now
###############################################################################
sub validate_seqfile {

  # Passed in Variables
  my ($file, $MSCNT) = @_;

  # Variables
  my ($seqnum, $tenpercent);
  my ($cmd1, $cmd2) = ("", "");  

  print "Validating the total number or sequences within: $file\n\n";

  # Count Number of Sequences in file
  $seqnum = `grep -c ^\\> $file`;
  chomp $seqnum;

  # Determine the Ten Percent Number of Total Seqs
  $tenpercent = ceil($MSCNT * .10);
  $tenpercent = $MSCNT - $tenpercent;

  # Print out the Information
  print "\n---------------------------------------------------------------------\n";
  print "Sequence Count in $file\n\n";
  print " Count -> $seqnum\n";
  print " Total -> $MSCNT\n";
  print " 10%    -> $tenpercent\n";

  # Set the Warnings and the errors
  if ($seqnum < $tenpercent) {
    print "\nWARNING - There is a Discrepancy in total number of Fasta Sequences and Sequence count, verify and run again!!\n\n";
    #exit(1);
  } elsif ( $seqnum < $MSCNT ) {
    print "\nWARNING - Counts do not match - this is caused by duplicates and fragments.  As duplicates are hard to keep\n";
    print "\ttrack of especially when the group is 'structure' it is best to work with total count for validation\n\n"; 
  } else {
    print "\nNumbers look good - Proceed as normal\n\n";
  }
  print "---------------------------------------------------------------------\n\n"; 
}

###############################################################################
# dbCAN - Uses dbcan to identify specific pul/family regions in the sequences
#          which is a bit mute since this is already identified in the
#          sequences extracted from Cazy - but is usefull if user adds their 
#          own sequences
#
# Notes:  To run properly need to create sym links in working directory to 
#          a few dbcan files
############################################################################### 
sub dbcan {

  # Get passed in Variables
  my ($cazy, $fam, $dir) = @_;

  # Variable Declaration
  my ($cmd1, $cmd2) = ("", "");
  my $prune = $fam . ".extracted.pruned.fasta";

  # Create Directory for dbcan and change to this directory
  my $local = $dir . "/dbcan";
  unless ( -d $local ) { mkdir($local, 0755); }
  local $CWD = $local;

  # Create basname of $cazy
  my $fbase = basename($cazy);
  $fbase =~ s/\.[^.]+$//;

  # Do not need to do if already done before
  if ( -f $prune ) {
    print "dbCAN has been run before, continue with script\n";
  } else {

    # Need to create symbolic links 
    print "Creating Symbolic links to dbcan files for running on cazy sequences\n";
    $cmd1 = "ln -s /usr/local/dbcan/all.hmm.ps.len .; ";
    $cmd2 = "ln -s /usr/local/dbcan/hmmscan-parser.sh .; ";
    &run_cmd($cmd1, $cmd2);
    $cmd1 = "";
    $cmd2 = "";
    print "Links created, proceed with dbcan\n\n";

    # Declare variables
    my $hmmfile = $fam . ".hmm.out";
    my $hmmdm = $hmmfile . ".dm";
    my $dbcanfile = $fam . ".dbcan.ps";

    # Run dbCan
    print "Running hmmscan\n";
    $cmd1 = "hmmscan --domtblout $hmmdm /usr/local/dbcan/dbCAN-fam-HMMs.txt $cazy > $hmmfile; ";
    &run_cmd($cmd1, $cmd2);
    $cmd1 = "";
    print "hmmscan complete\n\n";
    print "Running hmmscan-parser\n";
    $cmd1 = "sh hmmscan-parser.sh $hmmdm > $dbcanfile; ";
    &run_cmd($cmd1, $cmd2);
    $cmd1 = "";
    print "hmmscan-parser complete\n\n";

    # Proceed with Extracting Sequences and Pruning Sequences Before returning filename

    ## Extract Sequences
    print "Extract Sequences that hit family: $fam\n";
    my $tag = $fam;
    $cmd1 = "grep $tag'[_0-9]*\.hmm' $dbcanfile > $fam.dbcan.final.txt; ";
    &run_cmd($cmd1, $cmd2);
    $cmd1 = "";
  
    $cmd1 = "awk -F\" \" '{print \$3}' $fam.dbcan.final.txt > $fam.dbcan.final.seqlist; ";
    &run_cmd($cmd1, $cmd2);
    $cmd1 = "";

    $cmd1 = "extract_seq_data.pl --config $fam.dbcan.final.seqlist --fastaq $cazy --Format fasta; ";
    &run_cmd($cmd1, $cmd2);
    $cmd1 = "";

    my $extract_out = $fbase . ".mod.fasta";
    $cmd1 = "mv $extract_out $fam.extracted.fasta; ";
    &run_cmd($cmd1, $cmd2);
    $cmd1 = "";
    print "Extraction of sequences is completed\n\n";

    ## Prune Sequences
    print "Pruning extracted sequences based on dbCAN hit locations\n";
    $cmd1 = "awk -F\" \" '{print \$3 \"\\t\" \$8 \"\\t\" \$9}' $fam.dbcan.final.txt > seqid_coord.txt; ";
    &run_cmd($cmd1, $cmd2);
    $cmd1 = "";

    $cmd1 = "prune_seqs.pl -c seqid_coord.txt -f $fam.extracted.fasta; ";
    &run_cmd($cmd1, $cmd2);
    $cmd1 = "";
    print "Prunning completed\n\n";
  }

  # Create filename to return
  my $prunedfile = $local . "/" . $prune; 
   
  return $prunedfile;
}

###############################################################################
# Muscle - Run the pruned sequence file through a muscle alignment and output 
#           results in phylip interleaved format
############################################################################### 
sub muscle {

  # Passed in Variables
  my ($dbcan, $fam, $dir, $muscle_fast, $R) = @_;

  # Variable Declaration
  my ($cmd1, $cmd2) = ("", "");

  # Create Directory for muscle and change to this directory
  my $local = $dir . "/muscle";
  unless ( -d $local ) { mkdir($local, 0755); }
  local $CWD = $local;

  # Muscle output filename
  my $muscle_file = $fam . ".muscle_aln.phyi";

  # Run Muscle Alignment
  if ( -f $muscle_file ) {
    print "Muscle has already been run, continuing script\n";
  } else {
    print "Running the muscle alignment on the pruned fasta data\n";
    $cmd1 = "muscle -in $dbcan -phyiout $muscle_file; ";
    &run_cmd($cmd1, $cmd2);
    print "Muscle Alignment completed\n\n";
  }

  # Set prot_muscle
  $prot_muscle = $local . "/" . $muscle_file;

  # Muscle rename output filename
  my $muscle_ren_file = $fam . ".muscle_aln_mod.phyi";
  
  # $dir will be different depending on if this is first run or subsample
  my $id_file;
  if ($R == 0) {
    $id_file = $dir . "/cazy/key_id_pairs.txt";
  } else {
    $id_file = "../../cazy/key_id_pairs.txt"; 
  }

  # Run the Muscle Renamer script
  if ( -f $muscle_ren_file ) {
    print "Muscle Rename has already been run, continuing script\n";
  } else {
    print "Running the muscle rename script on the muscle output\n";
    $cmd1 = "muscle_rename.pl -i $muscle_file -k $id_file; ";
    &run_cmd($cmd1, $cmd2);
    print "Muscle Rename Completed\n\n";
  }

  # Need to account for if RAxML is chosen or if FastTree is Used
 
  # Create the FastTree Variant
  if ( -f $muscle_fast ) {
    print "Muscle - FastTree Variant has been created, continuing script\n";
  } else {
    print "Creating a FastTree variant file of the muscle output\n";
    $cmd1 = "awk '{if (\$1 \~ \/[0-9]\/) print \$0 \; else print \"\",\$0}' $muscle_ren_file | sed 's/^ \$//' > $muscle_fast; ";
    &run_cmd($cmd1, $cmd2);
    print "FastTree Variant Created\n\n";
  }

  # The plain muscle file is easiest for prottest so will us it with subsample run
  if ( $R == 0) {
    return $local . "/" . $muscle_ren_file;
  } else {
    return $local . "/" . $muscle_file;
  }
}

###############################################################################
# Prottest 3 - This program is used to determine the best model for RaxML to  
#               build the tree with - uses built in parser
############################################################################### 
sub prottest {

  # Get passed in Variables
  my ($muscle, $cazy, $dbcan, $fam, $dir, $thr, $MF) = @_;

  # Variable Declaration
  my ($cmd1, $cmd2) = ("", "");
  my $prot_file = $fam . "_prot_model.txt";
  my $RaxML_Tree;

  # Create Directory for prottest and change to this directory
  my $local = $dir . "/prottest";
  unless ( -d $local ) { mkdir($local, 0755); }
  local $CWD = $local;

  # Need to confirm file has not been made already and if so get model from file
  if ( -f $prot_file ) {
    print "Prottest has be run already, continue with script\n";

    # need to set variable from what is in file
    my @ar;
    open (IN, "<$prot_file") or die "Cannot open $prot_file $!\n";
    while (<IN>) {
      push(@ar, $_);
    }
    close(IN);
    $RaxML_Tree = $ar[0];
  } else {

    # Before Actually running prottest we need to make sure there are not more 
    #  than 4000 taxa in the fasta file or prottest will not run, if there are 
    #  more than 4000 taxa, we need to make a subsample dataset to run with
    if ( -d "subsample" ) {
      print "Subsample step has been run prior, Continuing\n\n";
      $muscle = $local . "/subsample/muscle/" . $fam . ".muscle_aln.phyi";
    } else {
      if ($NumberSeqs >= 4000) {
      	# Must use the dbcan pruned file as if large user dataset could grab seqs that will get pruned
        $muscle = &subsample($dbcan, $fam, $local, $thr, $MF);
      }
    }

    # Copy the properties file for prottest to the cwd
    my $prottest_properties = $local . "/prottest.properties";
    if ( -f $prottest_properties ) {
      print "Properties File has already been copied, proceeding with prottest run\n\n";
    } else {
      $cmd1 = "cp /usr/local/prottest3/prottest.properties .; ";
      &run_cmd($cmd1, $cmd2);
      $cmd1 = "";
    }

    # Run prottest on the muscle results
    print "Running prottest - search for best model\n";
    my $outfile = $fam . ".prottest.out";
    $cmd1 = "java -jar /usr/local/prottest3/prottest3.jar -i $muscle -o $outfile -S 0 -all-distributions -AIC -AICC -BIC -DT -threads $thr; ";
    system $cmd1;
    print "Prottest finished - proceeding with parsing prottest results\n";

    # Parse the prottest results to obatain the model for raxML 
    my (@model, %raxmodel, $matrix);
    my ($i, $g) = (0,0);
    my $file = $local . "/" . $outfile;

    open (IN, "<$file") or die "Cannot open $file $!\n";
    while (<IN>) {
      chop;
      # Only want lines that give the best model
      next if (!/^Best/);
      # Push only the name of the model into the array, discard preceding text
      push(@model,( split /:/ )[ -1 ]);
    }
    close(IN);

    # Use the models parsed from the file to create the raxml modelname and push to hash incrementing identical values
    foreach (@model) {
  
      # Deterine if matrix model uses I and G
      $i = 1 if /\+I/;
      $g = 1 if /\+G/;
  
      # Set the Matrix name
      if ( $_ =~ /^\s+(.*?)\+/ || $_ =~ /^\s+(.*?)$/ ) { $matrix = $1 }

      # Set Tree ModelName based on RAxML or FastTree
      my $rxm;
      if ( $TreeBuild == 0 ) {
        # Create the RaxML ModelName
        if ($g > 0) {
          if ($i > 0) {
            $rxm = 'PROTGAMMAI' . $matrix;
          } else {
            $rxm = 'PROTGAMMA' . $matrix;
          }
        } else {
          if ($i > 0) {
            $rxm = 'PROTCATI' . $matrix;
          } else {
            $rxm = 'PROTCAT' . $matrix;
          }
        }
      } else {
	# Create the FastTree ModelName
	my $mod;

	if ($g > 0) {
	  $mod = 'gamma';
	} else {
	  $mod = 'cat';
	}

	if ($matrix eq 'WAG') {
	  $rxm = $mod . "-wag";
	} elsif ($matrix eq 'LG') {
	  $rxm = $mod . "-lg";
	} else {
	  $rxm = $mod . "-jtt";
	}
      }  
      $i = 0; $g = 0;

      # Set the Tree ModelName to a hash
      if (exists $raxmodel{ $rxm }) {
        $raxmodel{ $rxm } = $raxmodel{ $rxm } + 1;
      } else { 
        $raxmodel{ $rxm } = 1;
      }
    }

    # Set - final matrix name to the hash key with the largest count
    $RaxML_Tree = reduce { $raxmodel{$a} > $raxmodel{$b} ? $a : $b } keys %raxmodel;
    print "Parsing of prottest files completed and have best model: $RaxML_Tree for Raxml Run\n\n";
   
    # Print out Result to file
    open (OUT, ">$prot_file") or die "Can't open $prot_file:$!\n";
    print OUT $RaxML_Tree;
    close OUT;
  }

  return $RaxML_Tree;
}

###############################################################################
# RaxML - Now it is time to use muscle alignment and the model chosen by 
#          prottest to build the best tree 
############################################################################### 
sub raxml {

  # Passed in Variables
  my ($muscle, $tree, $fam, $dir, $raxml, $threads) = @_;

  # Variable Declaration
  my ($cmd1, $cmd2) = ("", "");
  my $rax_tree = "RAxML_bipartitions.A1";
  my $bootstrap = 100;

  # Calculate an optimal number of threads to use
  my $muscle_depth = `head -1 $muscle`;
  my @md = split(/ /, $muscle_depth);
  my $opt_thr = ceil( $md[1] / 300 );

  # Now set threads used to the lower value of $opt_thr and $threads
  if ( $threads < $opt_thr ) {
    $opt_thr = $threads;
  } elsif ($threads >= 16) {
    $opt_thr = 16;
  } 

  # Create Directory for muscle and change to this directory
  my $local = $dir . "/raxml";
  unless ( -d $local ) { mkdir($local, 0755); }
  local $CWD = $local;

  # Lets get raxml running
  if ( -f $rax_tree ) {
    print "\n\nEntire script has been run before, why are you running this again???\n\n";
    exit;
  } else {
    print "Building best tree - using RaxML\n";
    if ($opt_thr == 1) {
      $cmd1 = "$raxml -f a -p 0111 -s $muscle -n A1 -m $tree -x 0123 -# $bootstrap; ";
    } else {
      $cmd1 = "$raxml -f a -T $opt_thr -p 0111 -s $muscle -n A1 -m $tree -x 0123 -# $bootstrap; ";
    }
    &run_cmd($cmd1, $cmd2);
    print "RaxML has finished\n\n";
  }

  return $local . "/" . $rax_tree;
}

###############################################################################
# FastTree - Use muscle alignment and the model chosen by prottest to build 
#             tree - faster than RAxML
############################################################################### 
sub fasttree {

  # Passed in Variables
  my ($muscle, $tree, $fam, $dir) = @_;

  # Variable Declaration
  my ($cmd1, $cmd2) = ("", "");
  my $fast_tree = "FastTree_bootstrap.tree";

  # Create Directory for muscle and change to this directory
  my $local = $dir . "/fasttree";
  unless ( -d $local ) { mkdir($local, 0755); }
  local $CWD = $local;

  # Lets parse the Tree Model to set proper flags in run
  my @model = split(/-/, $tree);

  # Lets get fasttree running
  if ( -f $fast_tree ) {
    print "\n\nEntire script has been run before, why are you running this again???\n\n";
    exit;
  } else {
    print "Building best tree - using FastTree\n";
    
    if ( $model[0] eq 'cat' ) {

      if ( $model[1] eq 'wag' ) {
	$cmd1 = "fasttree -wag -out $fast_tree $muscle; ";
      } elsif ( $model[1] eq 'lg' ) {
	$cmd1 = "fasttree -lg -out $fast_tree $muscle; ";
      } else {
	$cmd1 = "fasttree -out $fast_tree $muscle; ";
      }
    } else {

      if ( $model[1] eq 'wag' ) {
	$cmd1 = "fasttree -wag -gamma -out $fast_tree $muscle; ";
      } elsif ( $model[1] eq 'lg' ) {
	$cmd1 = "fasttree -lg -gamma -out $fast_tree $muscle; ";
      } else {
	$cmd1 = "fasttree -gamma -out $fast_tree $muscle; ";
      }
    }
    &run_cmd($cmd1, $cmd2);
    print "FastTree has finished\n\n";
  }

  return $local . "/" . $fast_tree;
}

###############################################################################
# Subsample - If there are >4000 taxa in the Cazy fasta file, we need to create 
#              a smaller subset to then test with prottest, as prottest will
#              crash on alignments with >4000 taxa
############################################################################### 
sub subsample {

  # Passed in Variables
  my ($pfile, $fam, $dir, $threads, $mf) = @_;

  # Variable Declaration
  my ($cmd1, $cmd2) = ("", ""); 

  # Create Directory for muscle and change to this directory
  my $local = $dir . "/subsample";
  unless ( -d $local ) { mkdir($local, 0755); }
  local $CWD = $local;

  print "#######################################################################################\n";
  print "#######################################################################################\n";
  print "#  Creating a Subsample of the Cazy Fasta Sequences to run through prottest as the    #\n";
  print "#   original fasta file has over 4000 taxa, which is beyound the bounds of prottest   #\n";
  print "#######################################################################################\n";
  print "#######################################################################################\n\n";

  # Extract the Subsample of Sequences - uses a script from MIT
  print "Extracting subsample...\n\n";
  my $cazy_file = $fam . "_subsample.fasta";
  $cmd1 = "fasta_subsample.pl $pfile 1500 -seed 7 > $cazy_file; ";
  &run_cmd($cmd1, $cmd2);

  $cazy_file = $local . "/" . $cazy_file;

  # Call dbcan
  print "Running dbcan...\n\n";
  my $dbcan_file = &dbcan($cazy_file, $fam, $local);

  # Call muscle
  print "Running muscle...\n\n";
  my $muscle_file = &muscle($dbcan_file, $fam, $local, $mf, 1);

  print "#######################################################################################\n";
  print "#######################################################################################\n";
  print "#  Subsample process has been created and returning new subsample muscle filename     #\n";
  print "#   for prottest testing                                                              #\n";
  print "#######################################################################################\n";
  print "#######################################################################################\n\n";
  
  # Return muscle file
  return $muscle_file;
}

###############################################################################
# Run Command - This runs the system commands from each of the 4 steps in qc
#                utilizing threads (multi-core) 
############################################################################### 
sub run_cmd {

   # Get passed in Variables
   my ($cmd1, $cmd2) = @_;

   # Print Out Job Process
   print "Commands To Run - \n\t$cmd1\n";
   print "\t$cmd2\n" if defined $cmd2;
   print "\n";

   # Running Command
   if (defined $cmd2) {
      my $thr1 = threads->new(\&thrsub, $cmd1);
      my $thr2 = threads->new(\&thrsub, $cmd2);
      $thr1->join();
      $thr2->join();
   } else {
      my $thr1 = threads->new(\&thrsub, $cmd1);
      $thr1->join();
   } 
}

###############################################################################
# Threads - This threads the commands (multi-core) 
############################################################################### 
sub thrsub {

   # Passed in Variables
   my ($cmd) = @_;

   # Run the command
   print "Threading: $cmd\n";
   system $cmd;
   print "************************** Threading Ended\n";
}

__END__


=pod

=head1 NAME

 Saccharis.pl

=head1 SYNOPSIS

 Saccharis.pl [-Options [--] [Arguments...]]

=head1 DESCRIPTION

 Wrapper Script for the Cazy Pipeline

 Uses a number of different scripts and programs to perform Cazy family 
 identification and linkage through a phylogenetic tree. 
 Extracts sequences from the Cazy website, based on family and group name
 and will combine with User sequences. Process then proceeds through
 pruning, alignment, tree model determination and building of tree.

 Uses User defined flags as shown below. A logfile can be created, via
 redirecting STDOUT to a file.

 Switches can be done in long or short form and are case sensitive
 eg:
   Saccharis.pl --raxml
   Saccharis.pl -r
   Saccharis.pl -F     
   Saccharis.pl -f     { family }

 Please use the manual option to get more information on each of the flags.

=head1 ARGUMENTS

 --directory   <path>
 --threads     <threads>              
 --group       <group_name>          
 --family      <family_name>
 --RAXML                   : Flag - to choose RAXML over FastTree
 --raxml       <version_raxml>        
 --seqfile     </path/filename>
 --Fragments               : Boolean Flag - if added means keep Fragments
 --help                    : print Options and Arguments instead of fetching data
 --man                     : print complete man page instead of fetching data

=head1 Options
 
 --directory | -d  <path>

     : You can set a predefined output directory with this flag, though you
       must give it's path.  Default is Current Working Directory (CWD)

 --threads | -t    <threads>
 
     : This is basically straightforward.  Some scripts allow the use of
       multi-core processing.  Set a number in here from 1 to > <max_cores>
       The default is set at 2.

 --group | -g      <group_name>

     : This can be a single group or a list of groups to run.  Lists are comma
       delimited in single quotes. They are case-sensitive - lowercase please...

       -> eg. characterized
       -> eg. 'characterized, all, bacteria'

     : Allowable GroupNames

       -> all, archaea, bacteria, eukaryota, unclassified, characterized,
          subfamilies

 --family | -f     <family_name>
 
     : This is a single family name. 

       -> eg. GH43
       
       NOTE: You can only run a single family per run...

 --raxml | -r      <version_raxml>

     : There are 3 different raxml programs that you can run, which one you 
       choose will be based off of the age and type of your CPU.
 
       -> Older/Slower CPU           -- raxml
       -> Somewhat Older/Faster CPU  -- raxmlHPC-PTHREADS-SSE3 (default)
       -> Newest Procesors           -- raxmlHPC-PTHREADS-AVX

 --seqfile | -s    </path/filename>

     : If you would like to add your own sequences to this run - this is your 
       chance.  Sequences MUST be in FASTA FORMAT - if they are not the script
       will fail.  Make sure to include path with filename.
     : Sequence Header lines must start with a uniqure 9 character string with
       the following format:

       -> >U[0-9]{1,8} - That is a >U followed by 8 digits

	ex. U00000001

 --Fragments | -F

     : This is a boolean value flag that by default is set to TRUE, which means
       fragments are left out by default.  If you would like to include fragment
       sequences from CAZY, include this flag in our call

=head1 AUTHOR

 Dallas Thomas

=head1 CREDITS

 Visual Genomics Lab at the University of Calgary

=head1 TESTED

 Perl    5.10.1
 Debian  5
 Debian  6

=head1 BUGS

 None that I know of

=head1 TODO

 The Duplicate Removal Phase is rather poor and I do not know how effective it
  really is.  This phase needs a complete overhaul and better algorithms.
 I am working on this issue.

=head1 UPDATES

 Nothing yet...

=cut

# SACCHARIS
Sequence Analysis and Clustering of CarboHydrate Active enzymes for Rapid Informed prediction of Specificity (SACCHARIS), is a PERL based pipeline designed to improve functional predictions of uncharacterized sequences for any CAZyme or CBM family currently maintained on the CAZy website or within user-defined datasets.
# Citation
When using SACCHARIS please site the following paper:
> Jones DR, Thomas DK, Alger N, Ghavidel A, Inglis GD, Abbott DW. SACCHARIS: An automated pipeline to inform discover of new carbohydrate active enzyme activities within polyspecific families and de novo sequence datasets. Biotechnology for Biofuels, (out for review, BBIO-D-17-00286).
# License
This software is distributed under the terms of the GPL, version 2 or later, excepting that:
- The third party programs and scripts used by SACCHARIS are covered by the terms of their respective licenses
# Additional Scripts
With this package I have included a copy of:
  - fasta_rmSmall.pl
    - This script screens a fasta file and will remove sequence data where length of sequence is smaller than a user definined minimum length
# Requirements
- Perl Libraries
  - Bio::Seq, Bio::SeqIO
  - Date::Calc
  - File::chdir
  - GetOpt::Long
  - HTML::TagParser  (*see note below*)
  - List::Util
  - LWP::Simple
  - Threads
- Third Party Software
  - dbCAN
  - HMMER 3.1
  - MUSCLE
  - ProtTest 3
  - RAxML
  - FastTree (Requires version 2.1.10 or greater)
  - Fasta_subsample.pl
- Notices
  - If you experience an error pertaining to an uninitialized $esearch value, confirm you have the following packages installed
    -  libwww-perl (linux)
    -  LWP::Protocol::https (OSX)
# Installation
1. Install all Requirements
2. Clone Repository `git clone`
3. Copy Scripts to a location in the *PATH*
## HMMER Installation
1. Download [HMMER](http://hmmer.org/download.html)
2. Extract archive
3. Copy or Move folder to `/usr/local/hmmer`
4. Add binaries directory to your *Path*
## dbCAN Installation
1. `mkdir /usr/local/dbcan`
2. Download [dbCAN](http://csbl.bmb.uga.edu/dbCAN/download.php)
   - dbCAN-fam-HMMs.txt
   - hmmscan-parser.sh
3. Format HMM db 
   - `hmmpress dbcan-fam-HMMs.txt`
## ProtTest Installation
1. Clone Repository `git clone`
2. Move directory to `/usr/local/prottest3`
## MUSCLE Installation
1. Download [MUSCLE](http://www.drive5.com/muscle/downloads.htm)
2. Copy `muscle` binary to location in the path like `/usr/local/bin`
## FastTree Installation
1. Install as per directions [Here](http://www.microbesonline.org/fasttree/#Install)
## RAxML Installation
1. Clone Repository `git clone`
2. Follow Directions in `README` to create Executables
3. Move executables to a location in the *Path*
## Fasta_subsample.pl Installation
1. Script is included with SACCHARIS
2. Script was written by Timothy L. Bailey and William Noble
3. Script is part of the [MEME Suite](http://web.mit.edu/meme_v4.11.4/share/doc/overview.html)
# Modifications - HTML::TagParser
- Tagparser.pm throws an error on line 236 - **Fix - Alter Line to**
  - `bless $self, ref($package) || $package;`
# NCBI Registration - cazy_extract.pl
- NCBI E-Utilities Registration is Required for running of the cazy_extract.pl script
  - Steps Involved
    - Send email to eutilities@ncbi.nlm.nih.gov that includes the desired values for email address and tool name
        - eg. tool = SaccharisTool, email = your.name@domain
    - Create and account on NCBI (https://www.ncbi.nlm.nih.gov/account/) in the Settings page create and API Key
    - Uncomment and Modify Lines 31-33 of cazy_extract.pl - use the information from steps 1 and 2 to modify the script
- **cazy_extract.pl will not run without this information**
# Running SACCHARIS
- In terminal follow Usage as given by
  - `Saccharis.pl` ; or
  - `Saccharis.pl -m`

# SACCHARIS
Sequence Analysis and Clustering of CarboHydrate Active enzymes for Rapid Informed prediction of Specificity (SACCHARIS), is a PERL based pipeline designed to improve functional predictions of uncharacterized sequences for any CAZyme or CBM family currently maintained on the CAZy website or within user-defined datasets.
# Citation
When using SACCHARIS please site the following paper:
> Jones DR, Thomas DK, Alger N, Ghavidel A, Inglis GD, Abbott DW. SACCHARIS: An automated pipeline to inform discover of new carbohydrate active enzyme activities within polyspecific families and de novo sequence datasets. Biotechnology for Biofuels, (out for review, BBIO-D-17-00286).
# License
This software is distributed under the terms of the GPL, version 2 or later, excepting that:
- The third party programs and scripts used by SACCHARIS are covered by the terms of their respective licenses
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
  - FastTree
  - Fasta_subsample.pl
# Installation
1. Install all Requirements
2. Clone Repository `git clone`
3. Copy Scripts to a location in the *PATH*
## HMMER Installation
1. Download [HMMER](http://hmmer.org/download.html)
2. Extract archive
3. Copy or Move folder to `/usr/local/hmmer`
4. Add binaries directory to your *Path*

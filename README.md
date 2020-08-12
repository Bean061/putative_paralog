# Putative Paralogs Detection (PPD)

This method is based on shared heterozygous information to detect putative paralogs for Hyb-Seq data. Highly recommended for Hyb-Seq downstream analysis.

## Prerequisites:
Run the main [HybPiper](https://github.com/mossmatters/HybPiper) script and check the results in [individual]/[gene]/[individual]/sequences/intron/ and [individual]/[gene]/[individual]/sequences/FNA/ directories.

All the intron and exon sequences are considered as input for Step1.

## Software/dependencies

1. [Picard](http://broadinstitute.github.io/picard/) for Step1.

2. [GATK](https://software.broadinstitute.org/gatk/download/) for Step1.

3. [python3](https://www.python.org/downloads/) and its dependencies for Step2.
* [Biopython](https://biopython.org/wiki/Packages). Easy installation from [conda](https://biopython.org/wiki/Packages) 
* [numpy](https://numpy.org/doc/stable/user/whatisnumpy.html). Easy installation from [conda](https://anaconda.org/anaconda/numpy)

## Environment
Examples can be run on Mac or Linux.

## Steps
It includes three steps as follows:

1. Concatenate all the supercontigs into a single file. All supercontigs are stored in two new files named supercontig and exon.
  This script is modified from [Mossmatters github "Alleles from HybSeq Data"](https://github.com/mossmatters/phyloscripts/tree/master/alleles_workflow).
  
  ```
  Please type in the following arguments after Step1.sh in order: 
  * The direcotory HybPiper.
  * full name of namelist.
  example ./Step1.sh ./HybPiper-master/[Hybpiper_result_file]/ /[full_path]/namelist.txt
  ```
2.

3.


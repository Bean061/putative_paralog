# Putative Paralogs Detection (PPD)

This method is based on shared heterozygous information to detect putative paralogs for Hyb-Seq data. Highly recommended for Hyb-Seq downstream analysis.

## Prerequisites:
Run the main [HybPiper](https://github.com/mossmatters/HybPiper) script and check the results in [individual]/[gene]/[individual]/sequences/intron/ and [individual]/[gene]/[individual]/sequences/FNA/ directories.

All the intron and exon sequences are considered as input for Step1.

## Software/dependencies

1. [Picard](http://broadinstitute.github.io/picard/)

2. [GATK](https://software.broadinstitute.org/gatk/download/)

3. python3 dependencies
* [Biopython](https://biopython.org/wiki/Packages). Easy installation from [conda](https://biopython.org/wiki/Packages) 
* [numpy](https://numpy.org/doc/stable/user/whatisnumpy.html). Easy installation from [conda](https://anaconda.org/anaconda/numpy)

## Environment
examples can be run on Mac or Linux.

It includes three steps as follows:



1. Concatenate all the supercontigs into a single file. Given a directory gerated by HybPiper called prefix:

2.

3.


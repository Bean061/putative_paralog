#!/bin/bash

#### plese cite 
# Zhou et al., unpublished
# Kates, H.R., Johnson, M.G., Gardner, E.M., Zerega, N.J. and Wickett, N.J., 2018. Allele phasing has minimal impact on phylogenetic reconstruction from targeted nuclear gene sequences in a case study of Artocarpus. American journal of botany, 105(3), pp.404-416.
# Please check the oringinal website https://github.com/mossmatters/phyloscripts/tree/master/alleles_workflow
# one major change is adding -k 100 in BWA.

#########CHANGE THESE PATHS AS NEEDED###########
## If you install the picard and gatk, you don't need to set these PATHs.

#gatkpath=/opt/Software/GenomeAnalysisTK.jar
#picardpath=/opt/Software/picard/build/libs/picard.jar


################Prepared files #################

###### make output directory for iupac consensus sequences, which will be used for putative paralogs detection.
### mkdir iupac_supercontig
### mkdir iupac_exon

### 1. namelist, it is as same as the one used for hybpiper
### 2. paired-end raw reads fastq.gz files after quality control
### 3. cataneate all supercontigs into one $prefix.supercontigs.fasta for every species.
### you can use the following lines from Hybpiper website. Change the directory if needed

# while read prefix
# do

# cat ${prefix}/*/${prefix}/sequences/intron/*_supercontig.fasta > ${prefix}.supercontigs.fasta
# cat ${prefix}/*/${prefix}/intronerate.gff > ${prefix}.intronerate.fasta

# done < namelist.txt


#############COMMAND LINE ARGUMENTS############
 


###change dir to current working directory
echo "This script is used for generating the degenerated seuqences.
"

echo "Please type in the following arguments in order: 
1. BWA -k matching base length (if raw reads are 150 bp, we recommand use 100 here). 
2. ploidy number
3. a direcotory containing consensus sequences from Step 1.
4. a directory containing all raw reads file (keep raw reads in fastq.gz format)
5. output directory
6. full name of namelist.
example ./Step2.sh 100 2 ../
"

echo "CHECKING................................"
echo "BWA -k is $1"
echo "ploidy is $2"
echo "input direcotory for consensus sequences from Step 1: $3."
echo "input directory for raw reads (fastq.gz) is $4/"
echo "output directory is $5/"
echo "The namelist is $6"

while read prefix
do

### change the following path to the raw read fastq files.
read1fq=$4/${prefix}/${prefix}_QCP_R1.fastq.gz
read2fq=$4/${prefix}/${prefix}_QC_R2.fastq.gz

mkdir $prefix
cd $prefix

### copy the supercontigs under the new coresponding species folder.
supercontig=./${prefix}.supercontigs.fasta
cp $3/${prefix}.supercontigs.fasta ./


#####STEP ZERO: Make Reference Databases

picard CreateSequenceDictionary \
R=$supercontig 
bwa index $supercontig
samtools faidx $supercontig

#####STEP ONE: Map reads

echo "Mapping Reads"

### change the -k if you have longer reads.
bwa mem -k $1 $supercontig $read1fq $read2fq | samtools view -bS - | samtools sort - -o $supercontig.sorted.bam


### add 
picard FastqToSam  \
F1=$read1fq \
F2=$read2fq \
O=$supercontig.unmapped.bam \
SM=$supercontig

picard MergeBamAlignment \
ALIGNED=$supercontig.sorted.bam \
UNMAPPED=$supercontig.unmapped.bam \
O=$supercontig.merged.bam \
R=$supercontig

#####STEP TWO: Mark duplicates

echo "Marking Duplicates"
picard MarkDuplicates \
I=$supercontig.merged.bam \
O=$supercontig.marked.bam \
M=$supercontig.metrics.txt

#######STEP THREE: Identify variants, select only SNPs

echo "Identifying variants"

samtools index $supercontig.marked.bam
#samtools mpileup -B -f $supercontig $supercontig.marked.bam -v -u > $supercontig.vcf

echo "$2"
gatk HaplotypeCaller \
-R $supercontig \
-I $supercontig.marked.bam \
-ploidy $2 \
-O $supercontig.vcf



time gatk SelectVariants \
-R $supercontig \
-V $supercontig.vcf \
-select-type SNP \
-O $supercontig.snps.vcf 


######STEP FOUR: Output new supercontig FASTA with ambiguity codes

echo "Generating IUPAC FASTA file"

gatk FastaAlternateReferenceMaker \
-R $supercontig \
-O $prefix.degenerated.fasta \
-V $supercontig.snps.vcf \
--use-iupac-sample $supercontig


#### change the directory name
#### supercontig result
cp -r $prefix.degenerated.fasta $5

#### exon result
# cp -r $prefix.degenerated.fasta /Users/xianglab/tools/HybPiper-master/iupac2/iupac_exon/

cd ..


done < $6




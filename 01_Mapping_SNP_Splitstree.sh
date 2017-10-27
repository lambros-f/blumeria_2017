#!/bin/bash

# This is the bash script used to map all genomic reads to the DH14 genome
# and subsequently extract converage per SPs or SNPs for the phylogeny analysis



### Software used and other global variables that have to be set
gatk='java -jar ~/utils/gatk/GenomeAnalysisTK.jar'
bwa=~/utils/bin/bwa
star=~/utils/bin/STAR
samtools=~/utils/samtools-1.3.1/samtools
picard='java -jar ~/utils/bin/picard-2.8.2.jar'
snpeff='java -Xmx4g -jar ~/utils/snpEff/snpEff.jar eff'

JAVA_TOOL_OPTIONS='-Xmx90000m'
###

### This is a list with all the datasets downloaded from NCBI
### including the SRA id, isolate name and formae speciales they belong to
runs=/hpcwork/rwth0146/published_data/01_blumeria_other_genomic_reads/runslist.txt
y="$(wc -l runslist.txt | cut -f1)"

### Here you basically start a loop to process each genome individually,
### you can turn this to multiple scripts and submit them seperately in you HPC

x=1
while [ $x -le $y ]
do
genome=~/genomes/bgh_dh14_v4.fa
run_name="$(sed "$x!d" $runs | cut -f3 -d ' ')"
fsp="$(sed "$x!d" $runs | cut -f2 -d ' ')"
isolate="$(sed "$x!d" $runs | cut -f1 -d ' ')"
readx="$(sed "$x!d" $runs | cut -f3 -d ' ')"
read1=$readx\_1.fastq.gz
read2=$readx\_2.fastq.gz


cd $x-$fsp-$isolate #go to where the reads are stored
mkdir trimmed

### Trimming starts here, some dataset are single end reads
### so the script has to be adjusted accordingly

java -jar ~/utils/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
-threads 16 \
$read1 \
$read2 \
trimmed-$read1.fq \
unpaired-trimmed-$read1.fq \
trimmed-$read2.fq \
unpaired-trimmed-$read2.fq \
ILLUMINACLIP:~/utils/Trimmomatic-0.36/adapters/TruSeq2-PE.fa:5:30:10 \
SLIDINGWINDOW:3:18 \
TRAILING:18 \
MINLEN:20 

### With the paired-end datasets, I didn't keep the broken pairs
read1=$x-$fsp-$isolate/trimmed-$read1.fq 
read2=$x-$fsp-$isolate/trimmed-$read2.fq 

### Here the mapping with BWA starts, assuming you 
### have already generated an index file for your genome
### using bwa index $genome

$bwa mem -t 12 ../dh14_genome_index/bgh_dh14_v4.fa $read1 $read2 > $isolate-$fsp.sam

### Here Picard sorts and generates bam file

$picard \
AddOrReplaceReadGroups \
I=$isolate-$fsp.sam \
O=rg_added_sorted.bam \
SO=coordinate \
RGID=$isolate \
RGLB=$isolate-$fsp \
RGPL=ILLUMINA \
RGPU= \
RGSM=$isolate

$picard \
MarkDuplicates \
I=rg_added_sorted.bam \
O=$isolate-$fsp-dedupped.bam \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
M=output.metrics

### As described in the TSL Pathogenomics Workshop 2017 (https://github.com/SUSTAIN-COST-Action-Norwich)
### Keeping all sites, SNPs and non-SNP variants (indels)

$gatk \ " >> vcf_$x.sh
-T HaplotypeCaller \ " >> vcf_$x.sh
-R $genome \ " >> vcf_$x.sh
-ploidy 1 \ " >> vcf_$x.sh
--input_file $isolate-$fsp-dedupped.bam \ " >> vcf_$x.sh
--emitRefConfidence BP_RESOLUTION \ " >> vcf_$x.sh
--variant_index_type LINEAR \ " >> vcf_$x.sh
--variant_index_parameter 128000 \ " >> vcf_$x.sh
-o $isolate-$fsp\_sorted.g.vcf " >> vcf_$x.sh

done

### Now we join all VCF and filter variants

$gatk \
-T GenotypeGVCFs \
-R $genome \
--variant $isolate-$fsp #here include all VCFs by reusing the --variant 
--includeNonVariantSites \
-o variants.vcf

java -jar ~/gatk/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $genome \
-selectType SNP \
--variant variants.vcf \
-o variants-snps.vcf

java -jar /home/lf216591/utils/gatk/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R $genome \
--variant variants-snps.vcf \
--filterExpression "QD < '+5.0+'" --filterName "QDFilter" \
--filterExpression "QUAL < 5000.0" --filterName "QualFilter" \
--filterExpression "MQ < 20.0" --filterName "MQ" \
--filterExpression "ReadPosRankSum < -2.0+" \
--filterName "ReadPosRankSum" \
--filterExpression "ReadPosRankSum > 2.0" \
--filterName "ReadPosRankSum" \
--filterExpression "MQRankSum < -2.0" \
--filterName "MQRankSum" \
--filterExpression "MQRankSum > 2.0" \
--filterName "MQRankSum" \
--filterExpression "BaseQRankSum < -2.0" \
--filterName "BaseQRankSum" \
--filterExpression "BaseQRankSum > 2.0" \
--filterName "BaseQRankSum" \
-o variants-snps-filtered.vcf

### Now we use vcftools to generate select SNPs
### We keep only common SNP position, and allow no missing data

~/vcftools \
--vcf variants-snps-filtered.vcf \
--max-missing 1 \
--remove-filtered-all \
--recode \
--out variants-snps-filtered.vcf

### You can use the python script vcf_to_fasta_bp_resolution.py
### to exctract sequences per isolate from the vcf file
### https://github.com/SUSTAIN-COST-Action-Norwich

vcf_to_fasta_bp_resolution.py variants-snps-filtered.vcf

### This script extracts sequences per scaffold
### which is pretty handy when you have a chromosome level assembly
### but in our case you want one merged file so you can import everything
### in Splitstree, this should do the work

for isolate in `grep '>' variants-snps-filtered.vcf_scaffold_1.fasta | tr -d '>' | head -n -1 `; do
echo "Working on $isolate"
for file in `ls *.fasta`; do
#echo "Working on $file"
sed -n -e "/$isolate\$/,/>/ p" $file | head -n -1  >> $isolate-fasta.fa
done
grep -v '>' $isolate-fasta.fa | tr -d '\n' > 1.fa
ex -sc "1i|>$isolate" -cx 1.fa
awk '{ print length($0); }' 1.fa
cat 1.fa >> merged_snps.fa
rm 1.fa $isolate-fasta.fa
done

isolate=dh14
echo 'Working on DH14'
for file in `ls *.fasta`; do
sed -n -e "/BR32/,/>/ p" $file >> $isolate-fasta.fa
done
grep -v '>' $isolate-fasta.fa | tr -d '\n' > 1
ex -sc "1i|>DH14" -cx 1
cat 1 >> merged_snps.fa
rm 1 $isolate-fasta.fa

### The merged_snps.fa file can be imported to Splitstree
### and you can generate a network tree or a UPGMA tree





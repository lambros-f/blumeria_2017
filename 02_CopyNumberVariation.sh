#--- Chunk included in the R documentation ---

### First you need the per base coverage of the reference genome (DH14 in this case)
### based on the mapping done with BWA (see 01 script)

$bedtools genomecov -ibam $isolate.sorted.bam -bga > $isolate.sorted.bam.bedgraph

### Then you intersect the coverage file with the positions of the SP-coding genes
### for this you can generate a bed file based on the gff file of the annotation with
### the ranges of the CDS annotation (not the mRNA or the gene, since the include UTRs
### so your calculation will be off).

$bedtools \
intersect \
-wb -b secreted.bed -a $isolate.sorted.bam.bedgraph > 903-hordei-K1_ONT_genomic.bedgraph_overlap.all




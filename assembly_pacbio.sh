#!/bin/sh
# GENOME ASSEMBLY USING PACBIO SEQUENCING DATA ----------------
# This is the pipeline for the genome assembly using long-reads from PacBio

# Outline ------------------------
## 0. Stats on FASTQ using fastq-stats
## 1. Mapping 
##      1a. minimap2
##      1b. ngmlr
## 2. Mapping stats using samtools
## 3. Overall coverage evaluation using bedtools genomecov
## 4. Evaluation of mapping at fixed coverage (200X)
## 5. Assembly with Canu  
## 6. Assembly stats using FastaSeqStats
## 7. Extraction longest assembly
## 8. Annotation reconstructed genome using GeSeq


root="/data/99_genomicclass/griecomc/project/pacbio"
REF="/data/99_genomicclass/griecomc/project/reference/prunus_avium.fasta"
REF2="/data/99_genomicclass/griecomc/project/reference/prunus_apetala.fasta"
FASTQ="${root}/SRR4280451_1.fastq.gz"

BAM="${root}/PB_SRR4280451_minimap_pavium.bam"
BAM_mapped="${root}/PB_SRR4280451_minimap_pavium_mapped.bam"
BAM_ngmlr="${root}/PB_SRR4280451_ngmlr_pavium.bam"
BAM_ngmlr_mapped="${root}/PB_SRR4280451_ngmlr_pavium_mapped.bam"


BAM200x="${root}/PB_SRR4280451_minimap_pavium_200x.bam"
BAM200x_mapped="${root}/PB_SRR4280451_minimap_pavium_200x_mapped.bam"
BAM200x_ngmlr="${root}/PB_SRR4280451_ngmlr_pavium_200x.bam"
BAM200x_ngmlr_mapped="${root}/PB_SRR4280451_ngmlr_pavium_200x_mapped.bam"


BAM2="${root}/PB_SRR4280451_minimap_papetala.bam"
BAM2_ngmlr="${root}/PB_SRR4280451_ngmlr_papetala.bam"
BAM2_mapped="${root}/PB_SRR4280451_minimap_papetala_mapped.bam"
BAM2_ngmlr_mapped="${root}/PB_SRR4280451_ngmlr_papetala_mapped.bam"

FASTQ_minimap_mapped="${root}/PB_SRR4280451_minimap_pavium_mapped.fastq"
FASTQ_ngmlr_mapped="${root}/PB_SRR4280451_ngmlr_pavium_mapped.fastq"


# 0. Stats on FASTQ *******************
fastq-stats $FASTQ > "${root}/SRR4280451.fastq.stats.txt"


# 1. Mapping *******************

## a. minimap2 ------------------
echo "Running minimap2 using prunus avium reference..."
minimap2 -ax map-pb $REF $FASTQ | samtools view -hSb > $BAM

echo "Running minimap2 using prunus apetala reference..."
minimap2 -ax map-pb $REF2 $FASTQ | samtools view -hSb > $BAM2

## b. ngmlr ---------------------
echo "Running ngmlr using prunus avium reference..."
ngmlr -t 4 -r $REF -q $FASTQ | samtools view -hSb > $BAM_ngmlr

echo "Running ngmlr using prunus apetala reference..."
ngmlr -t 4 -r $REF2 -q $FASTQ | samtools view -hSb > $BAM2_ngmlr


# 2. Mapping stats using samtools ---------
echo "Running samtools on prunus avium mapped with minimap2..."
samtools sort $BAM -o $BAM
samtools index -b $BAM
samtools view -Sb -F 4 -o $BAM_mapped $BAM
samtools sort $BAM_mapped -o $BAM_mapped

samtools stats $BAM -r $REF > "${root}/PB_SRR4280451_minimap_pavium_stats.txt"
samtools stats $BAM_mapped -r $REF > "${root}/PB_SRR4280451_minimap_pavium_mapped_stats.txt"

echo "Running samtools on prunus avium mapped with ngmlr..."
samtools sort $BAM_ngmlr -o $BAM_ngmlr
samtools index -b $BAM_ngmlr
samtools view -Sb -F 4 -o $BAM_ngmlr_mapped $BAM_ngmlr
samtools sort $BAM_ngmlr_mapped -o $BAM_ngmlr_mapped

samtools stats $BAM_ngmlr -r $REF > "${root}/PB_SRR4280451_ngmlr_pavium_stats.txt"
samtools stats $BAM_ngmlr_mapped -r $REF > "${root}/PB_SRR4280451_ngmlr_pavium_mapped_stats.txt"


echo "Running samtools on prunus apetala with minimap2..."
samtools sort $BAM2 -o $BAM2
samtools index -b $BAM2
samtools view -Sb -F 4 -o $BAM2_mapped $BAM2
samtools sort $BAM2_mapped -o $BAM2_mapped

samtools stats $BAM2 -r $REF2 > "${root}/PB_SRR4280451_minimap_papetala_stats.txt"
samtools stats $BAM2_mapped -r $REF2 > "${root}/PB_SRR4280451_minimap_papetala_mapped_stats.txt"

cho "Running samtools on prunus apetala mapped with ngmlr..."
samtools sort $BAM2_ngmlr -o $BAM2_ngmlr
samtools index -b $BAM2_ngmlr
samtools view -Sb -F 4 -o $BAM2_ngmlr_mapped $BAM2_ngmlr
samtools sort $BAM2_ngmlr_mapped -o $BAM2_ngmlr_mapped

samtools stats $BAM2_ngmlr -r $REF2 > "${root}/PB_SRR4280451_ngmlr_papetala_stats.txt"
samtools stats $BAM2_ngmlr_mapped -r $REF2 > "${root}/PB_SRR4280451_ngmlr_papetala_mapped_stats.txt"


echo "Running bedtools bamtofastq"
## I convert the BAM (mapped to p.avium) into fastq in order to perform the subsampling using seqtk sample
bedtools bamtofastq -i $BAM_mapped -fq $FASTQ_minimap_mapped
bedtools bamtofastq -i $BAM_ngmlr_mapped -fq $FASTQ_ngmlr_mapped


# 3. Bedtools genome coverage *******************
## Mapping coverage of the whole BAM mapped
samtools sort $BAM_mapped -o $BAM_mapped
samtools sort $BAM_ngmlr_mapped -o $BAM_ngmlr_mapped
#samtools sort $BAM2_mapped -o $BAM2_mapped
#samtools sort $BAM2_ngmlr -o $BAM2_ngmlr

bedtools genomecov -d -ibam $BAM_mapped -g "$REF.fai" > "${root}/PB_SRR4280451_minimap_pavium_mapped_ref.covbed.txt"
bedtools genomecov -d -ibam $BAM_ngmlr_mapped -g "$REF.fai" > "${root}/PB_SRR4280451_ngmlr_pavium_mapped_ref.covbed.txt"


# 4. Evaluation of the coverage of the BAM made up with subset *******************
#       for a common depth coverage among the different technologies

#Extract reads using seqtk sample
# --------- FORMULA ----------------------------
# Coverage=Length_mean_read * N_reads / Genome size
# N_reads = Genome_size * Coverage / Length_mean_read

# genome_size=158k
# Coverage=200k
# Length_mean_read=7,697
# -----------------------------------------------

seqtk sample -s1000 $FASTQ_minimap_mapped 4105 > "${root}/PB_SRR4280451_minimap_pavium_mapped_200x.fastq"
seqtk sample -s1000 $FASTQ_ngmlr_mapped 4105 > "${root}/PB_SRR4280451_ngmlr_pavium_mapped_200x.fastq"

# To get the BAM i have to map again to the reference
# using minimap2
minimap2 -ax map-ont $REF "${root}/PB_SRR4280451_minimap_pavium_mapped_200x.fastq" | samtools view -hSb > $BAM200x
samtools view -Sb -F 4 -o $BAM200x_mapped $BAM200x

# using ngmlr
ngmlr -t 4 -r $REF -q "${root}/PB_SRR4280451_ngmlr_pavium_mapped_200x.fastq" -x ont | samtools view -hSb > $BAM200x_ngmlr
samtools view -Sb -F 4 -o $BAM200x_ngmlr_mapped $BAM200x_ngmlr


# If the input is in BAM (-ibam) format, the BAM file must be sorted by position. 
# Using samtools sort aln.bam aln.sorted will suffice.

samtools sort $BAM200x_mapped -o $BAM200x_mapped
bedtools genomecov -d -ibam $BAM200x_mapped -g "$REF.fai" > "${root}/PB_SRR4280451_minimap_pavium_mapped_ref_200x.covbed.txt"

samtools sort $BAM200x_ngmlr_mapped -o $BAM200x_ngmlr_mapped
bedtools genomecov -d -ibam $BAM200x_ngmlr_mapped -g "$REF.fai" > "${root}/PB_SRR4280451_ngmlr_pavium_mapped_ref_200x.covbed.txt"



# 5. Assembly with canu *****************************
mkdir ${root}/canu

echo "Running seqtk sample"
#Then I perform the subsampling of the fastq 
# 1k,2k,3k,5k,10k,20k,25k,30k,40k,50k

for l in 1 2 3 5 10 20 25 30 40 50; do  
    seqtk sample -s1000 $FASTQ_minimap_mapped ${l}000 > "${root}/PB_SRR4280451_minimap_pavium_mapped_sub${l}k.fastq"
    seqtk sample -s1000 $FASTQ_ngmlr_mapped ${l}000 > "${root}/PB_SRR4280451_ngmlr_pavium_mapped_sub${l}k.fastq"

done


# Running canu 
echo "Running canu"

# canu [-haplotype|-correct|-trim] \
#    [-s <assembly-specifications-file>] \
#    -p <assembly-prefix> \
#    -d <assembly-directory> \
#    genomeSize=<number>[g|m|k] \
#    [other-options] \
#    [-trimmed|-untrimmed|-raw|-corrected] \
#    [-pacbio|-nanopore|-pacbio-hifi] *fastq
# The -p option, to set the file name prefix of intermediate and output files, is mandatory. 
# If -d is not supplied, canu will run in the current directory, 
# otherwise, Canu will create the assembly-directory and run in that directory. 
# It is _not_ possible to run two different assemblies in the same directory.

# on the whole sample 
canu -d "${root}/canu/canu_minimap_whole" -p PB_SRR4280451_minimap_whole genomeSize=158k -pacbio $FASTQ_minimap_mapped
canu -d "${root}/canu/canu_ngmlr_whole" -p PB_SRR4280451_ngmlr_whole genomeSize=158k -pacbio $FASTQ_ngmlr_mapped

# on the subsamples
for l in 1 2 3 5 10 20 25 30 40 50; do
    canu -d "${root}/canu/canu_minimap_${l}k" -p PB_SRR4280451_minimap_sub${l}k genomeSize=158k -pacbio "${root}/PB_SRR4280451_minimap_pavium_mapped_sub${l}k.fastq"
    canu -d "${root}/canu/canu_ngmlr_${l}k" -p PB_SRR4280451_ngmlr_sub${l}k genomeSize=158k -pacbio "${root}/PB_SRR4280451_ngmlr_pavium_mapped_sub${l}k.fastq"

done


# 6. Assembly stats using FastaSeqStats *******************
mkdir ${root}/canu_stats

FastaSeqStats -i "${root}/canu/canu_minimap_whole/PB_SRR4280451_minimap_whole.contigs.fasta" > "${root}/canu_stats/PB_SRR4280451_minimap_whole.contigs.stats.txt"
FastaSeqStats -i "${root}/canu/canu_ngmlr_whole/PB_SRR4280451_ngmlr_whole.contigs.fasta" > "${root}/canu_stats/PB_SRR4280451_ngmlr_whole.contigs.stats.txt"

for l in 1 2 3 5 10 20 25 30 40 50; do
    FastaSeqStats -i "${root}/canu/canu_minimap_${l}k/PB_SRR4280451_minimap_sub${l}k.contigs.fasta" > "${root}/canu_stats/PB_SRR4280451_minimap_sub${l}k.contigs.stats.txt"
    FastaSeqStats -i "${root}/canu/canu_ngmlr_${l}k/PB_SRR4280451_ngmlr_sub${l}k.contigs.fasta" > "${root}/canu_stats/PB_SRR4280451_ngmlr_sub${l}k.contigs.stats.txt"
done


touch "${root}/canu_pacbio_stats.txt"

for f in $(ls ${root}/canu_stats | grep PB);
do
j=$(basename $f);
echo $j >> "${root}/canu_pacbio_stats.txt";
cat ${root}/canu_stats/$f >> "${root}/canu_pacbio_stats.txt";
echo "          ****************        " >> "${root}/canu_pacbio_stats.txt";
done


# 7. Extraction longest assembly using FastaExtract *******************

# 8. Annotation reconstructed genome using GeSeq *******************
# https://chlorobox.mpimp-golm.mpg.de/geseq.html








#!/bin/sh
# GENOME ASSEMBLY USING ILLUMINA SEQUENCING DATA ----------------
# This is the pipeline for the genome assembly using short-reads from ILLUMINA 


# Outline ------------------------
## 0. Removing evantual adapters
## 1. Stats on FASTQ using fastq-stats
## 2. Mapping using bowtie
## Note: there is the code to run BWA but here it was not run
## 3. Mapping stats using samtools
## 4. Overall coverage evaluation using bedtools genomecov
## 5. Evaluation of mapping at fixed coverage (200X)
## 6. Assembly with AbySS
##  6a. using SE reads
##  6b. using PE (separated withsamtools view) 
##  6c. using PE reads (separated with samtools view) 
## 7. Assembly stats using FastaSeqStats
## 8. Extraction longest assembly
## 9. Annotation reconstructed genome using GeSeq



root="/data/99_genomicclass/griecomc/project/illumina"
REF="/data/99_genomicclass/griecomc/project/reference/prunus_avium.fasta"
REF2="/data/99_genomicclass/griecomc/project/reference/prunus_apetala.fasta"
FASTQ1="${root}/SRR10362958_1.fastq.gz"
FASTQ2="${root}/SRR10362958_2.fastq.gz"
FASTQ1_trim="${root}/SRR10362958_1.trim.fastq.gz"
FASTQ2_trim="${root}/SRR10362958_2.trim.fastq.gz"
BAM="${root}/ILLU_SRR10362958.bowtie.ref.bam"
BAM2="${root}/ILLU_SRR10362958.bowtie.ref_apetala.bam"

BAM_mapped="${root}/ILLU_SRR10362958.map.bowtie.ref.bam"
BAM2_mapped="${root}/ILLU_SRR10362958.map.bowtie.ref_apetala.bam"

BAM200x="${root}/ILLU_SRR10362958.bowtie.ref.200x.bam"
BAM200x_mapped="${root}/ILLU_SRR10362958.map.bowtie.ref.200x.bam"

FASTQ_mapped="${root}/ILLU_SRR10362958.map.bowtie.ref.fastq"


# BAM_bwa="${root}/ILLU_SRR10362958.bwa.ref.bam"
# BAM_bwa_mapped="${root}/ILLU_SRR10362958.bwa.map.ref.bam"


# 0. Remove eventual adapters *******************
echo "Running fastq-mcf"
fastq-mcf /data/99_genomicclass/00_shared_data/Genomics/IlluminaAdapters_V2.fasta SRR10362958_1.fastq.gz -o SRR10362958_1.trim.fastq.gz
fastq-mcf /data/99_genomicclass/00_shared_data/Genomics/IlluminaAdapters_V2.fasta SRR10362958_2.fastq.gz -o SRR10362958_2.trim.fastq.gz
#no adapters found

# 1. Stats on FASTQ using fastq-stats *******************
echo "Running fastq-stats"
fastq-stats $FASTQ1 > "${root}/SRR10362958_1.fastq.stats.txt"
fastq-stats $FASTQ2 > "${root}/SRR10362958_2.fastq.stats.txt"
fastq-stats $FASTQ1_trim > "${root}/SRR10362958_1.trim.fastq.stats.txt"
fastq-stats $FASTQ2_trim > "${root}/SRR10362958_2.trim.fastq.stats.txt"


# 2. Mapping using bowtie *******************

# Create an index using Bowtie2-build
echo "Creating an index of the reference"
bowtie2-build $REF prunus_avium_ref
bowtie2-build $REF2 prunus_apetala_ref

# I can not use the trim.fastq since the two files have a different number of reads and so bowtie get error:
# "fewer reads in file specified with -2"
# I used the raw reads since it were cleaned and no had adapters

bowtie2 -p 8 -x ${root}/prunus_avium_ref -1 $FASTQ1 -2 $FASTQ2 | samtools view -Sb > $BAM
bowtie2 -p 8 -x ${root}/prunus_apetala_ref -1 $FASTQ1 -2 $FASTQ2 | samtools view -Sb > $BAM2

# 3. Mapping stats using samtools *******************

samtools view -Sb -F 4 -o $BAM_mapped $BAM
samtools sort $BAM_mapped -o $BAM_mapped

samtools sort $BAM -o $BAM
samtools stats $BAM -r $REF > "${root}/ILLU_SRR10362958_bowtie_pavium_stats.txt"
samtools stats $BAM2 -r $REF2 > "${root}/ILLU_SRR10362958_bowtie_papetala_stats.txt"

samtools view -Sb -F 4 -o $BAM2_mapped $BAM2
samtools sort $BAM2_mapped -o $BAM2_mapped
samtools stats $BAM2 -r $REF2 > "${root}/ILLU_SRR10362958_bowtie_papetala_stats.txt"

#to convert the BAM into SAM
#samtools view -h -o $BAM_mapped.sam $BAM


# #Mapping using BWA ----- NOT DONE!!!!
# bwa index $REF 
# bwa mem $REF $FASTQ1 $FASTQ2 | samtools view -Sb > $BAM_bwa

# samtools view -Sb -F 4 -o $BAM_bwa_mapped $BAM_bwa
# samtools sort $BAM_bwa_mapped -o $BAM_bwa_mapped

# samtools sort $BAM_bwa -o $BAM_bwa
# samtools stats $BAM_bwa -r $REF > "${root}/ILLU_SRR10362958_bwa_pavium_stats.txt"



# 4. Overall coverage evaluation using bedtools genomecov *******************
# Mapping coverage of the whole BAM mapped
samtools sort $BAM_mapped -o $BAM_mapped

bedtools genomecov -d -ibam $BAM_mapped -g "$REF.fai" > "${root}/ILLU_SRR10362958_bowtie_pavium_mapped_ref.covbed.txt"


# 5. Evaluation of the coverage of the BAM made up with subset *******************
#       for a common depth coverage among the different technologies

# Extract reads using seqtk sample
# --------- FORMULA ----------------------------
# Coverage=Length_mean_read * N_reads / Genome size
# N_reads = Genome_size * Coverage / Length_mean_read

# genome_size=158k
# Coverage=200k
# Length_mean_read=140
# N_reads= 158k*200k / 140 = 225000 reads (to be splitted in R1 and R2)
# -----------------------------------------------

# subsamples generated using seqtk sample
seqtk sample -s1000 $FASTQ_mapped.R1.fastq 112500 > "${root}/ILLU_SRR10362958_bowtie_pavium_mapped_200x_R1.fastq"
seqtk sample -s1000 $FASTQ_mapped.R2.fastq 112500 > "${root}/ILLU_SRR10362958_bowtie_pavium_mapped_200x_R2.fastq"

# To get the BAM i have to map again to the reference
bowtie2 -p 8 -x ${root}/prunus_avium_ref -1 "${root}/ILLU_SRR10362958_bowtie_pavium_mapped_200x_R1.fastq" -2 "${root}/ILLU_SRR10362958_bowtie_pavium_mapped_200x_R2.fastq" | samtools view -Sb > $BAM200x
samtools view -Sb -F 4 -o $BAM200x_mapped $BAM200x


# If the input is in BAM (-ibam) format, the BAM file must be sorted by position. 
# Using samtools sort aln.bam aln.sorted will suffice.

samtools sort $BAM200x_mapped -o $BAM200x_mapped
bedtools genomecov -d -ibam $BAM200x_mapped -g "$REF.fai" > "${root}/ILLU_SRR10362958_bowtie_pavium_mapped_ref_200x.covbed.txt"


# 6. Assembly **************************

## 6a. ABYSS USING READS AS SINGLE END 

## converting the BAM_mapped in the fastq
bedtools bamtofastq -i $BAM_mapped -fq $FASTQ_mapped

### on the whole BAM mapped with bowtie 
for k in `seq 17 8 96`; do
    mkdir ${root}/abyss_assembly/abyss_whole/k$k
    abyss-pe -C ${root}/abyss_assembly/abyss_whole/k$k name=pavium k=$k in=$BAM_mapped
done

### on the subsamples 
#### creating subsets
for l in 10 20 25 50 100 200 300 400 500 900; do  
    seqtk sample -s1000 $FASTQ_mapped ${l}000 > "${root}/ILLU_SRR10362958_pavium_bowtie_mapped_sub${l}k.fastq"
done

#### assembly
for s in 10 20 25 50 100 200 300 400 500 900; do
    mkdir ${root}/abyss_assembly/abyss_sub${s}k

    for k in `seq 17 8 96`; do
    mkdir ${root}/abyss_assembly/abyss_sub${s}k/sub${s}k_k$k
    abyss-pe -C ${root}/abyss_assembly/abyss_sub${s}k/sub${s}k_k$k name=pavium_sub${s}k k=$k in="${root}/ILLU_SRR10362958_pavium_bowtie_mapped_sub${s}k.fastq"
    done
done


### 6b. ABYSS using PE reads (meaning two fastq files) ---------------

#### a. Reads separated with samtools view
#### try to separate R1 reads from the R2 ones, in order to give two input files to abyss and perform an assembly PE
samtools view -hbf 64 $BAM_mapped > $BAM_mapped.R1.bam
samtools view -hbf 128 $BAM_mapped > $BAM_mapped.R2.bam
samtools stats $BAM_mapped.R1.bam -r $REF > "${root}/ILLU_SRR10362958_mapped_R1_bowtie_pavium_stats.txt"
samtools stats $BAM_mapped.R2.bam -r $REF > "${root}/ILLU_SRR10362958_mapped_R2_bowtie_pavium_stats.txt"

#### converting the BAM_mapped in the fastq using bedtools bamtofastq
bedtools bamtofastq -i $BAM_mapped.R1.bam -fq $FASTQ_mapped.R1.fastq
bedtools bamtofastq -i $BAM_mapped.R2.bam -fq $FASTQ_mapped.R2.fastq


#### assembly with whole reads
for k in `seq 17 8 96`; do
    mkdir ${root}/abyss_assemblyPE/abyss_whole/k$k
    abyss-pe -C ${root}/abyss_assemblyPE/abyss_whole/k$k name=pavium k=$k in="${FASTQ_mapped}.R1.fastq ${FASTQ_mapped}.R2.fastq"
done

#### on the subsamples
#### subsampling
for l in 5 10 12 25 50 100 150 200 250 450; do  
    seqtk sample -s1000 $FASTQ_mapped.R1.fastq ${l}000 > "${root}/ILLU_SRR10362958_pavium_bowtie_mapped_sub${l}k_R1.fastq"
    seqtk sample -s1000 $FASTQ_mapped.R2.fastq ${l}000 > "${root}/ILLU_SRR10362958_pavium_bowtie_mapped_sub${l}k_R2.fastq"

done

#### assembly
for s in 5 10 12 25 50 100 150 200 250 450; do
    mkdir ${root}/abyss_assemblyPE/abyss_sub${s}k

    for k in `seq 17 8 96`; do
    mkdir ${root}/abyss_assemblyPE/abyss_sub${s}k/sub${s}k_k$k
    abyss-pe -C ${root}/abyss_assemblyPE/abyss_sub${s}k/sub${s}k_k$k name=pavium_sub${s}k k=$k in="${root}/ILLU_SRR10362958_pavium_bowtie_mapped_sub${s}k_R1.fastq ${root}/ILLU_SRR10362958_pavium_bowtie_mapped_sub${s}k_R2.fastq"
    done
done



### 6b. Reads separated in R1 and R2 using bedtools bamtofastq (another way to separate reads R1 and R2)
#### separating reads into fastq1 and fastq2
samtools sort -n $BAM_mapped -o $BAM_mapped

#### converting the BAM_mapped in the fastq
bedtools bamtofastq -i $BAM_mapped \
                      -fq $FASTQ_mapped.end1.fq \
                      -fq2 $FASTQ_mapped.end2.fq

#### assembly on the whole
for k in `seq 17 8 96`; do
    mkdir ${root}/abyss_assemblyPEbis/abyss_whole/k$k
    abyss-pe -C ${root}/abyss_assemblyPEbis/abyss_whole/k$k name=pavium k=$k in="$FASTQ_mapped.end1.fq $FASTQ_mapped.end2.fq"
done

#### subsampling
for l in 5 10 12 25 50 100 150 200 250 450; do  
    seqtk sample -s1000 $FASTQ_mapped.end1.fq ${l}000 > "${root}/ILLU_SRR10362958_pavium_bowtie_mapped_sub${l}k_R1.fastq.end1.fq"
    seqtk sample -s1000 $FASTQ_mapped.end2.fq ${l}000 > "${root}/ILLU_SRR10362958_pavium_bowtie_mapped_sub${l}k_R2.fastq.end2.fq"
done

#### assembly on subsamples
for s in 5 10 12 25 50 100 150 200 250 450; do
    mkdir ${root}/abyss_assemblyPEbis/abyss_sub${s}k

    for k in `seq 17 8 96`; do
    mkdir ${root}/abyss_assemblyPEbis/abyss_sub${s}k/sub${s}k_k$k
    abyss-pe -C ${root}/abyss_assemblyPEbis/abyss_sub${s}k/sub${s}k_k$k name=pavium_sub${s}k k=$k in="${root}/ILLU_SRR10362958_pavium_bowtie_mapped_sub${s}k_R1.fastq.end1.fq ${root}/ILLU_SRR10362958_pavium_bowtie_mapped_sub${s}k_R2.fastq.end2.fq"
    done
done






# 7. Assembly stats using FastaSeqStats

## 7a. Stats on assembly using SE reads
### on the whole
for k in `seq 17 8 96`; do
    FastaSeqStats -i "${root}/abyss_assembly/abyss_whole/k$k/pavium-scaffolds.fa" > "${root}/abyss_stats/pavium-whole_k${k}-scaffold.stats.txt"
done

### on the subsamples
for s in 10 20 25 50 100 200 300 400 500 900; do
    
    for k in `seq 17 8 96`; do
    FastaSeqStats -i "${root}/abyss_assembly/abyss_sub${s}k/sub${s}k_k$k/pavium_sub${s}k-scaffolds.fa" > "${root}/abyss_stats/pavium-sub${s}k_k${k}-scaffold.stats.txt"
    done
done

touch "${root}/abyss_illuminaSE_stats.txt"

for f in $(ls ${root}/abyss_stats | grep pavium);
do
j=$(basename $f);
echo $j >> "${root}/abyss_illuminaSE_stats.txt";
cat ${root}/abyss_stats/$f >> "${root}/abyss_illuminaSE_stats.txt";
echo "          ****************        " >> "${root}/abyss_illuminaSE_stats.txt";
done



## 7b. Stats on assembly using PE reads
### on the whole
for k in `seq 17 8 96`; do
    FastaSeqStats -i "${root}/abyss_assemblyPE/abyss_whole/k$k/pavium-scaffolds.fa" > "${root}/abyss_statsPE/pavium-whole_k${k}-scaffold.stats.txt"
done

### on the subsamples
for s in 5 10 12 25 50 100 150 200 250 450; do
    
    for k in `seq 17 8 96`; do
    FastaSeqStats -i "${root}/abyss_assemblyPE/abyss_sub${s}k/sub${s}k_k$k/pavium_sub${s}k-scaffolds.fa" > "${root}/abyss_statsPE/pavium-sub${s}k_k${k}-scaffold.stats.txt"
    done
done

touch "${root}/abyss_illuminaPE_stats.txt"

for f in $(ls ${root}/abyss_statsPE | grep pavium);
do
j=$(basename $f);
echo $j >> "${root}/abyss_illuminaPE_stats.txt";
cat ${root}/abyss_statsPE/$f >> "${root}/abyss_illuminaPE_stats.txt";
echo "          ****************        " >> "${root}/abyss_illuminaPE_stats.txt";
done



## 7c. Stats on assembly using PE reads bis
### on the whole
for k in `seq 17 8 96`; do
    FastaSeqStats -i "${root}/abyss_assemblyPEbis/abyss_whole/k$k/pavium-scaffolds.fa" > "${root}/abyss_statsPEbis/pavium-whole_k${k}-scaffold.stats.txt"
done

### on the subsamples
for s in 5 10 12 25 50 100 150 200 250 450; do
    
    for k in `seq 17 8 96`; do
    FastaSeqStats -i "${root}/abyss_assemblyPEbis/abyss_sub${s}k/sub${s}k_k$k/pavium_sub${s}k-scaffolds.fa" > "${root}/abyss_statsPEbis/pavium-sub${s}k_k${k}-scaffold.stats.txt"
    done
done

touch "${root}/abyss_illuminaPEbis_stats.txt"

for f in $(ls ${root}/abyss_statsPEbis | grep pavium);
do
j=$(basename $f);
echo $j >> "${root}/abyss_illuminaPEbis_stats.txt";
cat ${root}/abyss_statsPEbis/$f >> "${root}/abyss_illuminaPEbis_stats.txt";
echo "          ****************        " >> "${root}/abyss_illuminaPEbis_stats.txt";
done


# 8. Extraction longest assembly
# Select the longest contig assembled with FastaExtract
FastaExtract -f "${root}/abyss_assemblyPE/abyss_sub100k/sub100k_k73/pavium_sub100k-scaffolds.fa" -l 100000 -o "${root}/abyssPE_sub100k_k73.fa"
FastaExtract -f "${root}/abyss_assemblyPE/abyss_sub100k/sub100k_k81/pavium_sub100k-scaffolds.fa" -l 100000 -o "${root}/abyssPE_sub100k_k81.fa"
FastaExtract -f "${root}/abyss_assemblyPE/abyss_sub100k/sub100k_k89/pavium_sub100k-scaffolds.fa" -l 100000 -o "${root}/abyssPE_sub100k_k89.fa"
FastaExtract -f "${root}/abyss_assemblyPE/abyss_sub150k/sub150k_k57/pavium_sub150k-scaffolds.fa" -l 100000 -o "${root}/abyssPE_sub150k_k57.fa"
FastaExtract -f "${root}/abyss_assemblyPE/abyss_sub150k/sub150k_k57/pavium_sub150k-scaffolds.fa" -l 130000 -o "${root}/abyssPE_sub150k_k57_longest.fa"


# 9. Annotation reconstructed genome using GeSeq
# https://chlorobox.mpimp-golm.mpg.de/geseq.html




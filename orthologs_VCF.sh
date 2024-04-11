#!/bin/bash
source /local/anaconda3/bin/activate

# run using nohup bash orthologs_VCF.sh > orthologs_VCF.out &!


MAPPED_PTH="/home/cody/ETHIOPIAN_CALOSOMA/data_2023/2mapped/HIC_genome"
ORTHO_CONS="/home/cody/ETHIOPIAN_CALOSOMA/data_2023/ortho_3consensus"
GENOME="HIC_genome.fasta" # can also be changed to be first bash input GENOME=$1
# THREADS=$2 # gets added to the native-pair-hmm-threads flag for gatk

cd /home/cody/Calosoma_phylo/ETHIOPIAN_CALOSOMA/data_2023

mkdir -p ortho_4snps/${GENOME}

conda activate picard_gatk4

# create list of samples as a variable
SAMPLES=""
for DATA in $(ls ${ORTHO_CONS}/*_q20_sorted.bam); do 
    SAMPLES=${SAMPLES}" -I "$DATA; 
done

gatk CreateSequenceDictionary -R ${GENOME}

# index all bam files

for SAMPLE in $(ls ${ORTHO_CONS}/*_q20_sorted.bam); do 
    samtools index -@ 4 ${SAMPLE}
done

# *_sorted_keep_rg_pcrdup.bam
gatk HaplotypeCaller \
    --native-pair-hmm-threads 4 \
    -R ${GENOME} \
    -O ortho_4snps/${GENOME}_snpcalling.vcf \
    ${SAMPLES}
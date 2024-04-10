#!/bin/bash
source /local/anaconda3/bin/activate
# Modified from Ines Carrrasquer "ortho_consenus_lastpart.sh" script by Cody Raul Cardenas 2024.04
# file is run like so
# modified_ortho_consensus_after.sh genome.fa mincov3 30 100 3 100

# where:
REFERENCE=$(echo $1) # Full path to the genome/reference fasta. The fasta header should be only one word, e.g. CM067376.1
sam2bed_mincov=$(echo $2) # Minimum read coverage to produce bed file per sample. Recommended value: 3
sharedbed_dist=$(echo $3) # Maximum distance between regions while producing the shared-by-all bed file. It prevents the bed file from being massive. Recommended value: 30
sharedbed_minlen=$(echo $4) # Minimum lenght of the region in the shared-by-all bed file to be kept. Recommended value: 100
consensus_mincov=$(echo $5) # Minimum read coverage to produce consensus sequence. Recommended value: 3
consensus_minlen=$(echo $6) # Minimum lenght to produce consensus sequence. Recommended value: 100

# requires mapped reads (see data in 2mapped/HIC_genome dir) generated using the "mapping.sh" script
MAPPED_PTH="/home/cody/ETHIOPIAN_CALOSOMA/data_2023/2mapped/HIC_genome"
ORTHO_CONS="/home/cody/ETHIOPIAN_CALOSOMA/data_2023/ortho_3consensus"

###########################################################
# since we didn't clean up by mapping quality, do that now
# explicitly indicate WHERE to do the work
###########################################################
cd /home/cody/Calosoma_phylo/ETHIOPIAN_CALOSOMA/data_2023
conda activate sam-bam-bedtools

for SAMPLE in $(ls 0clean_demux/ | cut -d "_" -f 1 | sort -u); do
# with a mapped sam/bam file use:
printf "mapping filter for ${SAMPLE}\n"
    samtools view -bF 4 -q 20 \
        ${MAPPED_PTH}/${SAMPLE}_sorted_keep_rg.bam \
        > ${ORTHO_CONS}/${SAMPLE}_q20.bam
printf "sorting ${SAMPLE}\n"
    samtools sort ${ORTHO_CONS}/${SAMPLE}_q20.bam \
        > ${ORTHO_CONS}/${SAMPLE}_q20_sorted.bam

    rm ${ORTHO_CONS}/${SAMPLE}_q20.bam;

# REFGROUP MAY NEED READDED IN THIS PIPELINE

done


##########
# bam2bed
##########
# get bed file with min coverage

for SAMPLE in $(ls 0clean_demux/ | cut -d "_" -f 1 | sort -u); do

printf "getting genome coverage for ${SAMPLE}\n"
    bedtools genomecov -bga -ibam \
        ${ORTHO_CONS}/${SAMPLE}_q20_sorted.bam | \
        awk '$4 > '"${sam2bed_mincov}"'' | \
        sort -V -k1,1 -k2,2n | \
        bedtools merge > ${ORTHO_CONS}/${SAMPLE}_q20_sorted.bed;
done

# create a shared bedfile
printf "creating shared bed file \n"
cat ${ORTHO_CONS}/*q20_sorted.bed | \
    sort -k1,1 -k2,2n | \
    bedtools merge -d ${sharedbed_dist} | \
    awk ' $3-$2 > '"${sharedbed_minlen}"' {print $0"\t"$1":"$2"-"$3 }' \
        > ${ORTHO_CONS}/SHARED_BY_ALL.bed

conda deactivate
#####################
# generate consensus
#####################

for SAMPLE in $(ls 0clean_demux/ | cut -d "_" -f 1 | sort -u); do
    # Get consensus sequence with same lenght than the complete chromosomes 
    # bcftools & vcfutils.pl are native on this environment
printf "generating consensus file ${SAMPLE}\n"
    samtools mpileup -A -uf ${REFERENCE} \
        ${ORTHO_CONS}/${SAMPLE}_q20_sorted.bam | \
        bcftools call -c | \
        bcftools filter -e "DP<${consensus_mincov}" | \
        vcfutils.pl vcf2fq \
            > ${ORTHO_CONS}/tmp.fastq;

        #From fastq to fasta
        conda activate seqtk
printf "converting tmp.fastq to tmp2.fasta file ${SAMPLE}\n"
        seqtk seq -A -n N ${ORTHO_CONS}/tmp.fastq \
            > ${ORTHO_CONS}/tmp2.fasta;
        conda deactivate

    #Extract the regions for this sample and give them the shared-by-all region code
    conda activate sam-bam-bedtools
printf "generating intersect file ${SAMPLE}\n"
    bedtools intersect \
        -a ${ORTHO_CONS}/SHARED_BY_ALL.bed \
        -b ${ORTHO_CONS}/${SAMPLE}_q20_sorted.bed | \
        awk '{print $4 }' | \
        sort -u \
            > ${ORTHO_CONS}/tmp_${SAMPLE}_ALL.faidx;
    conda deactivate 
printf "indexing ${SAMPLE} tmp2.fasta\n"
    samtools faidx ${ORTHO_CONS}/tmp2.fasta;

    conda activate sam-bam-bedtools
printf "cleaning uip indexing ${SAMPLE} tmp2.fasta\n"
    samtools faidx \
        -r ${ORTHO_CONS}/tmp_${SAMPLE}_ALL.faidx \
        ${ORTHO_CONS}/tmp2.fasta \
            > ${ORTHO_CONS}/tmp3_clean.fasta;
    conda deactivate
    #Keep only those region longer than consensus_minlen 
    conda activate seqtk
printf "keeping regions longer than ${consensus_minlen} in ${SAMPLE} tmp2.fasta\n"
        seqtk comp ${ORTHO_CONS}/tmp3_clean.fasta | \
        awk '$2-$9 > '"${consensus_minlen}"' { print $1}' \
            > ${ORTHO_CONS}/tmp_seq;
    conda deactivate

printf "generate orthologous fasta file for ${SAMPLE} \n"
    perl fastaselect.pl ${ORTHO_CONS}/tmp3_clean.fasta \
        ${ORTHO_CONS}/tmp_seq \
            > ${ORTHO_CONS}/orthologous_${SAMPLE}.fasta;

done

rm ${ORTHO_CONS}/tmp*
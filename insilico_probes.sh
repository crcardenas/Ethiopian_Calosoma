#! /bin/bash
#CARAP3_L1.consens_uniq.fa are probes generated using ipyrad...
# this script maps the probes to a genome, and extracts a fasta file of unique mappings
WORKGINDIRECT=${1}
BEDSCORE=${2} # 30
cd ${WORKGINDIRECT}
# process probe references

# only do once
source /local/anaconda3/bin/activate /home/jeremy/local/envgatk
java -jar -mx128G /home/jeremy/local/envpicard/share/picard-2.20.2-0/picard.jar \
        CreateSequenceDictionary \
        R=../CARAP3_L1.consens_uniq.fa \
        O=../CARAP3_L1.consens_uniq.dict
conda deactivate

#These need to be indexed and a dictionary created for them.
for SAMPLE in $(ls *.fasta | cut -d "." -f1); do
break # if already done kill it.... worry about if then logic later
# make index for mapping
        bwa index ${SAMPLE}.fasta;
# create picard index of probe references & faidx file
# R= input fasta file; O= output file, appended out put with dict
        picard CreateSequenceDictionary \
                -R ${SAMPLE}.fasta \
                -O ${SAMPLE}.dict;
        samtools faidx ${SAMPLE}.fasta;
       
done

for SAMPLE in $(ls *.fasta); do
        ID=$(echo ${SAMPLE} | cut -d "." -f 1);
         # if already done kill it.... worry about if then logic later
        # map
        bwa mem ${SAMPLE} ../CARAP3_L1.consens_uniq.fa \
                > probes_on_${ID}.sam

        # sort and keep non-dupe maps
        samtools sort probes_on_${ID}.sam -o temp_1_sorted.bam
        samtools view -bF 4 temp_1_sorted.bam > \
        #        probes_on_${ID}_sorted_keep.bam
       
        # index
        samtools index probes_on_${ID}_sorted_keep.bam

        source /local/anaconda3/bin/activate /home/cody/.conda/envs/sam-bam-bedtools
        #use bedtools to slice out our sequences!
        bedtools bamtobed -i probes_on_${ID}_sorted_keep.bam \
                > ${ID}_mapped.bed
        # have to pick something, 30 seemed reasonable here.
        awk -v BEDSCORE=${BEDSCORE} '$5 >= BEDSCORE {print $0}' ${ID}_mapped.bed \
         > ${ID}_score-filter.bed
        # find overlaps in bed file
        bedtools merge -i ${ID}_score-filter.bed \
                -c 1 -o count \
                > input_overlap_count.tmp
        # filter rows that do not overlap
        awk '/\t1$/{print}' input_overlap_count.tmp \
                > input_no_overlap.tmp
        # intersect original input and keep rows that were found after filtering
        bedtools intersect -a ${ID}_mapped.bed \
                -b input_no_overlap.tmp -wa \
                > ${ID}_no_overlap.bed

        # brute force first, refine later, there might be more quantitative filters we can use.
        bedtools getfasta -fi ${SAMPLE} -bed ${ID}_no_overlap.bed - \
        > ${ID}_mapped.fasta
done
# once pipeline is working... 

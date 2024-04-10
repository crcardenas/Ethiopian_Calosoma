#!/bin/bash
source /local/anaconda3/bin/activate

cd /data/work/Calosoma_phylo/ETHIOPIAN_CALOSOMA/data_2023/1probes

for REFERENCE in HIC_genome; do #for REFERENCE in CarInt CARAP3_L1.consens_uniq GCA022063505_gen; do
break
# break if reference tools already run
# make index for mapping
	bwa index ${REFERENCE}.fa;

conda activate picard_gatk;
# create picard index of probe references
# R= input fasta file
# O= output file, appended with dict
	picard CreateSequenceDictionary \
		-R ${REFERENCE}.fa \
		-O ${REFERENCE}.dict;

conda deactivate;

	samtools faidx ${REFERENCE}.fa;
	
done

cd /home/cody/Calosoma_phylo/ETHIOPIAN_CALOSOMA/data_2023

#start a for loop to map and sort reads
for REFERENCE in HIC_genome; do #for REFERENCE in CarInt CARAP3_L1.consens_uniq GCA022063505_gen; do 
	for SAMPLE in CBX0213 CBX0214 CBX0215 CBX0216 CBX0217 CBX0218 CBX0219 CBX0220 CBX0221 CBX0222 CBX0223 CBX0225 CBX0226 CBX0227 CBX0228 CBX0230 CBX0234 CBX0235 CBX0236 CBX0237 CBX0238 CBX0239 CBX0240 CBX0229 CBX0231 CBX0224 CBX0232 CBX0233 DRR295720 CBX0058 CBX0066 CBX0068 CBX0035 CBX0022 CBX0061 CBX0038 CBX0020 CBX0028 CBX0040 CBX0062 CBX0039 CBX0071 CBX0032 CBX0034 CBX0033 SRR2083640 CBX0075 CBX0309 CBX0583 CBX0308 CBX0879 CBX0310 CBX0367 CBX0269 CBX0603; do # $(ls 0clean_demux/ | cut -d "_" -f 1 | sort -u); do

        echo -e "\nstarting ${SAMPLE}\n";
        echo -e "\nMapping ${SAMPLE}\n";
        # use bwa meme to map sequence data to hyradx probes
        bwa mem -t 3 1probes/${REFERENCE}.fa \
            0clean_demux/${SAMPLE}_R1.fastq.paired.fq \
            0clean_demux/${SAMPLE}_R2.fastq.paired.fq \
            > 2mapped/${REFERENCE}/${SAMPLE}_on_ref.sam;

        echo -e "\nSorting ${SAMPLE}\n";
        # use samtools to sort samples
        # first by sorting
        samtools sort 2mapped/${REFERENCE}/${SAMPLE}_on_ref.sam \
            -o 2mapped/${REFERENCE}/tmp_${SAMPLE}_sorted.bam;

        echo -e "\nPulling ${SAMPLE} Mapped Sequences\n"
        # then by selecting only those that map
        samtools view -bF 4 \
            2mapped/${REFERENCE}/tmp_${SAMPLE}_sorted.bam \
            > 2mapped/${REFERENCE}/${SAMPLE}_sorted_keep.bam;

	conda activate picard_gatk;

        echo -e "\nCreating read groups from ${SAMPLE} bam file\n"
        # create readgroups from sorted bam files
        picard AddOrReplaceReadGroups \
            -I 2mapped/${REFERENCE}/${SAMPLE}_sorted_keep.bam \
            -O 2mapped/${REFERENCE}/${SAMPLE}_sorted_keep_rg.bam \
            -RGID ${SAMPLE} \
            -RGLB ID -RGPL PL -RGPU PU \
            -RGSM ${SAMPLE}

        echo -e "\nCreating index of ${SAMPLE}_sorted_keep_rg.bam\n"
        # create index for downstream workflow
		 samtools index 2mapped/${REFERENCE}/${SAMPLE}_sorted_keep_rg.bam

# these throw an OOM error with large genomes? yikes. all notes on the broad institute suggest that this is an uncessary step these days? IDK about for capture methods like hyradx though... we will see!
        # echo -e "\nFinding indels in ${SAMPLE}\n"
        # find indels
      	#GenomeAnalysisTK -T RealignerTargetCreator \
       #     -I 2mapped/${REFERENCE}/${SAMPLE}_sorted_keep_rg.bam \
       #     -R 1probes/${REFERENCE}.fa \
       #     -o 2mapped/${REFERENCE}/${SAMPLE}.intervals

       # echo -e "\nAligining indels in ${SAMPLE}\n"
        # indel realignment
       # GenomeAnalysisTK -T IndelRealigner \
       #     -I 2mapped/${REFERENCE}/${SAMPLE}_sorted_keep_rg.bam \
       #     -R 1probes/${REFERENCE}.fa \
       #     -targetIntervals 2mapped/${REFERENCE}/${SAMPLE}.intervals \
       #     -o 2mapped/${REFERENCE}/${SAMPLE}_sorted_keep_rg_realign.bam

        # realign indels
        # java -jar -Xmx128G /home/jeremy/local/envgatk/opt/gatk-3.8/GenomeAnalysisTK.jar \
         #   -T IndelRealigner \
         #   -I 2mapped/${REFERENCE}/tmp_${SAMPLE}_sorted_keep_rg_realign.bam \
         #   -R 1probes/${REFERENCE}.fa \
         #   -targetIntervals 2mapped/${REFERENCE}/tmp_${SAMPLE}.intervals \
         #   -o 2mapped/${REFERENCE}/${SAMPLE}_sorted_keep_rg_realign.bam

        echo -e "\nPCR duplicate filtering of ${SAMPLE}\n"

        # remove PCR duplicates
		# Jeremy advised against removing PCR dupes here
        picard MarkDuplicates \
            -I 2mapped/${REFERENCE}/${SAMPLE}_sorted_keep_rg.bam \
            -O 2mapped/${REFERENCE}/${SAMPLE}_sorted_keep_rg_pcrdup.bam \
            -REMOVE_DUPLICATES true \
            -M 2mapped/${REFERENCE}/${SAMPLE}_PCR-dupe_metrics.txt

        echo -e "\nIndexing PCR dupe cleaned ${SAMPLE}\n"
        # index *_sorted_keep_rg_realign_pcrdup.bam file
        samtools index 2mapped/${REFERENCE}/${SAMPLE}_sorted_keep_rg_pcrdup.bam


        # have to remove damaged/ancient DNA in a seperate step
        # R is broken
        # see: http://ginolhac.github.io/mapDamage/
        # remove temp files
 #       echo -e "\nRemoving ${SAMPLE} temporary files\n"
 #       rm 2mapped/${REFERENCE}/tmp_*
 #       echo -e "\n${SAMPLE} finished\n"
	conda deactivate

#rm 2mapped/${REFERENCE}/tmp_*
	done;
done


## get mapping statistics
#for REFERENCE in CARAP3_L1.consens_uniq GCA022063505_gen; do
#for i in $(ls 2mapped/${REFERENCE}/*_sorted_keep.bam); do
#TOTAL=$(samtools flagstat $i | grep "QC-passed" | cut -d"+" -f1);
#SAMPLE=$(echo ${i} | cut -d"/" -f3)
#printf "${REFERENCE}\t${SAMPLE}\t${TOTAL}\n"; 
#done;
#done > mapstats.tsv
#! /bin/bash
# classic phylohyrad pipeline
#CARAP3_L1.consens_uniq.fa

# process probe references
for REFERENCE in CARAP3_L1.consens_uniq; do
        # make index for mapping
	bwa index ${REFERENCE}.fa;

         # if you see this line, you will need this package in your own conda environment
        source /local/anaconda3/bin/activate /home/jeremy/local/envpicard/;
        # create picard index of probe references
	/home/jeremy/local/envpicard/share/picard-2.20.2-0/picard.jar \
        CreateSequenceDictionary \
		-R ${REFERENCE}.fa \
		-O ${REFERENCE}.dict;
        conda deactivate;
	samtools faidx ${REFERENCE}.fa;
done

# Jeremy's chunk; if used see and cite https://github.com/JeremyLGauthier/PHyRAD

for i in `ls *_R1.fastq.paired.fq`; do
        sample=`echo $i | sed -e 's/_R1.fastq.paired.fq//g' -e 's/clean_demux-//g'`
        # map
        bwa mem CARAP3_L1.consens_uniq.fa ${sample}_R1.fastq.paired.fq ${sample}_R2.fastq.paired.fq > ${sample}_on_ref.sam
        # sort 
        samtools sort ${sample}_on_ref.sam -o temp_1_sorted.bam
        samtools view -bF 4 temp_1_sorted.bam > temp_1_sorted_keep.bam
        
        source /local/anaconda3/bin/activate /home/jeremy/local/envpicard/
        java -jar -mx128G /home/jeremy/local/envpicard/share/picard-2.20.2-0/picard.jar CreateSequenceDictionary R=CARAP3_L1.consens_uniq.fa O=CARAP3_L1.consens_uniq.dict
        java -jar -mx128G /home/jeremy/local/envpicard/share/picard-2.20.2-0/picard.jar AddOrReplaceReadGroups I=temp_1_sorted_keep.bam O=temp_1_sorted_keep_rg.bam ID=[${sample}] RGLB=[id] PL=[pl] PU=[pu] SM=[${sample}]
        conda deactivate

        samtools index temp_1_sorted_keep_rg.bam

        source /local/anaconda3/bin/activate /home/jeremy/local/envgatk
        GenomeAnalysisTK -T RealignerTargetCreator -I temp_1_sorted_keep_rg.bam -R CARAP3_L1.consens_uniq.fa -o temp.intervals
        java -jar -Xmx128G /home/jeremy/local/envgatk/opt/gatk-3.8/GenomeAnalysisTK.jar -T IndelRealigner -I temp_1_sorted_keep_rg.bam -R CARAP3_L1.consens_uniq.fa -targetIntervals temp.intervals -o ${sample}_sorted_keep_rg_realign.bam
        conda deactivate

        rm temp*
        rm CARAP3_L1.consens_uniq.dict
done


# to do later
# find freebayes path and format this for your pipeline!
# freebayes -f "$reference"_locus_uniq.fasta *.rescaled.bam > snpset_freebayes.vcf

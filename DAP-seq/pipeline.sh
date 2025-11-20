#!/bin/bash
#PBS -q regular
#PBS -N MYBS_rep2_mapping
#PBS -l nodes=1:ppn=10
#PBS -l walltime=96:00:00
#PBS -l mem=40gb

cd /scratch/hbgao/MYBS_DAP-seq_rep2/

#############################################
# 1. Quality trimming with fastp
#############################################
fastp \
  -i /project3/gen230/gen230-2/gen230/gen230_R1.fq.gz \
  -I /project3/gen230/gen230-2/gen230/gen230_R2.fq.gz \
  -o MYBS_all_rep2.trim.fastq.gz \
  -O MYBS_all_rep2.R.trim.fastq.gz

#############################################
# 2. Mapping reads with Bowtie2
#############################################
bowtie2 \
  -x /home/zefulu/genomes/Ta-ensemble/Ta \
  -1 MYBS_all_rep2.trim.fastq.gz \
  -2 MYBS_all_rep2.R.trim.fastq.gz \
  -S MYBS_all_rep2.sam \
  -p 10 \
  -X 1000

#############################################
# 3. Convert SAM to BAM, filter, sort, mark duplicates
#############################################
samtools view -q 10 -b MYBS_all_rep2.sam -o MYBS_all_rep2.bam
samtools sort MYBS_all_rep2.bam -o MYBS_all_rep2.sort.bam

java -jar /opt/picard-2.21.9/picard.jar MarkDuplicates \
  INPUT=MYBS_all_rep2.sort.bam \
  OUTPUT=MYBS_all_rep2.clean.bam \
  METRICS_FILE=MYBS_all_rep2.dup_metrics.txt \
  REMOVE_DUPLICATES=true \
  VALIDATION_STRINGENCY=LENIENT

samtools index -m 14 MYBS_all_rep2.clean.bam

#############################################
# 4. Generate coverage tracks
#############################################
bamCoverage --bam MYBS_all_rep2.clean.bam -o MYBS_all_rep2.clean.bw
bedtools bamtobed -i MYBS_all_rep2.clean.bam > MYBS_all_rep2.clean.bed
bedtools genomecov \
  -i MYBS_all_rep2.clean.bed \
  -split \
  -bg \
  -g /scratch/hbgao/wheat-ref/Ta.genome.sizes \
  > MYBS_all_rep2.clean.bg

#############################################
# 5. Normalize coverage by total intersected reads
#############################################
for i in MYBS_all_rep2.clean
do
    x=$(bedtools intersect -a ${i}.bed -b At_narrowPeak.bed -wa | bedtools sort -i - | wc -l)
    awk -v total=$x '{print $1 "\t" $2 "\t" $3 "\t" $4*1000000/total}' ${i}.bg > ${i}.norm.bg
    sort -k1,1 -k2,2n ${i}.norm.bg > ${i}.sorted.bg
    bedGraphToBigWig ${i}.sorted.bg /scratch/hbgao/wheat-ref/Ta.genome.sizes ${i}.norm.bw
done

#############################################
# 6. Call peaks with MACS2
#############################################

# Joint peak calling: replicate1 + replicate2
macs2 callpeak \
  -t /scratch/hbgao/MYBS_DAP-seq_rep2/MYBS_all_rep2.clean.bed \
     /scratch/hbgao/MYBS_DAP-seq/all_clean/MYBS_all.clean.bed \
  -g 1.6e10 \
  -q 0.05 \
  -n MYBS_all_rep1_rep2 \
  --outdir together_call_peak

# Individual replicate peak calling
macs2 callpeak \
  -t /scratch/hbgao/MYBS_DAP-seq_rep2/MYBS_all_rep2.clean.bed \
  -g 1.6e10 -q 0.05 \
  -n MYBS_all_rep2 \
  --outdir rep2_call_peak

macs2 callpeak \
  -t /scratch/hbgao/MYBS_DAP-seq/all_clean/MYBS_all.clean.bed \
  -g 1.6e10 -q 0.05 \
  -n MYBS_all_rep1 \
  --outdir rep1_call_peak

#############################################
# 7. Filter peaks by score (>6)
#############################################
for i in MYBS_all_rep1_peaks MYBS_all_rep1_rep2_peaks MYBS_all_rep2_peaks
do
    awk '{if($7 > 6) print $0}' ${i}.narrowPeak > ${i}.Score6.narrowPeak
done

#############################################
# 8. Further filter by column thresholds
#############################################
for i in MYBS_all_rep1_peaks.Score6 MYBS_all_rep1_rep2_peaks.Score6 MYBS_all_rep2_peaks.Score6
do
    for j in {8..14}
    do
        awk -v col=$j '{if($8 > col) print $0}' ${i}.narrowPeak > ${i}.p.${j}.narrowPeak
    done
done

echo "DAP-seq processing completed."

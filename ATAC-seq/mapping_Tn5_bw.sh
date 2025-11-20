#######【1. Mapping & BAM processing】######
#!/bin/bash
#PBS -q regular
#PBS -N ATAC2_all_W235_mapping
#PBS -l nodes=1:ppn=4
#PBS -l walltime=96:00:00
#PBS -l mem=40gb

cd /scratch/hbgao/ATAC2-seq/jiace

# ===========================
# Mapping + BAM Processing
# ===========================

REF=/scratch/hbgao/wheat_all_ref/bowtie-index/Ta_all
SAMPLES="W2_all"

for i in $SAMPLES
do
    echo "Processing $i ..."

    # Bowtie2 mapping
    bowtie2 -x $REF \
        -1 ${i}.trim.fastq.gz \
        -2 ${i}.R.trim.fastq.gz \
        -S ${i}.all.sam \
        -p 6 \
        -X 1000

    # SAM → BAM → sort
    samtools view -q 10 -b ${i}.all.sam -o ${i}.all.bam
    samtools sort ${i}.all.bam -o ${i}.allsorted.bam

    # Remove duplicates
    java -jar /opt/picard-2.21.9/picard.jar MarkDuplicates \
        INPUT=${i}.allsorted.bam \
        OUTPUT=${i}.allclean.bam \
        METRICS_FILE=${i}.dedup.metrics.txt \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT

    # Index BAM
    samtools index -m 14 ${i}.allclean.bam

    # bigwig
    bamCoverage --bam ${i}.allclean.bam -o ${i}.allclean.bw

    echo "Done: $i"
done

#######【2. remove chloroplast/mitochondria reads】#######
# Filtering mito/chloro and keep proper paired reads
samtools view -h sample.srtdup.bam \
    | grep -v chloroplast \
    | grep -v mitochondrinon \
    | samtools view -q 10 -F 1804 -f 0x2 - \
    > sample.filt.bam

#######【3. Peak calling with MACS2】#######
#!/bin/bash
#PBS -q regular
#PBS -N W_callpeak
#PBS -l nodes=1:ppn=4
#PBS -l walltime=96:00:00
#PBS -l mem=40gb

cd /scratch/hbgao/ATAC_all/bow
OUTDIR=/scratch/hbgao/ATAC_all/peak

mkdir -p $OUTDIR

SAMPLES="N0-2_W1_all N1-1_W2_all N1-2_W3_all N3-2_W5_all D2-1_W4_all"

for i in $SAMPLES
do
    echo "Calling peaks for $i ..."
    macs2 callpeak \
        -t ${i}.bed \
        -g 1.6e10 \
        --nomodel \
        --shift -100 \
        --extsize 200 \
        --keep-dup all \
        -n $i \
        --outdir $OUTDIR
done

#######【4. Peak length >200 statistics and sort based on p-value】#######

cd /scratch/hbgao/ATAC_all/peak

SAMPLES="N0-1 N0-2_W1_all N1-1_W2_all N1-2_W3_all N2-1 N2-2 N3-1 N3-2_W5_all D1-1 D1-2 D2-1_W4_all D2-2 D3-1 D3-2"

for i in $SAMPLES
do
    echo "Filtering and sorting: $i"

    awk 'BEGIN{OFS="\t"} {
        len = $3 - $2;
        if(len > 200) {
            print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10
        }
    }' $i.narrowPeak \
    | sort -k8,8n \
    > ${i}.L200_Psort.narrowPeak
done

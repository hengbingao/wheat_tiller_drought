
#####   1. fastp trim adaptor   ####
#!/bin/bash
#PBS -q regular
#PBS -N fastp
#PBS -l nodes=1:ppn=4
#PBS -l walltime=96:00:00
#PBS -l mem=40gb

cd /project3/GHB_RNA_seq/221215_A01047_0158_BHF5W5DSX5

RAW_OUT=/scratch/hbgao/MYBS_RNAseq/raw
mkdir -p $RAW_OUT

for i in $(ls *fastq.gz | sed 's/\.R[12].fastq.gz//g' | uniq)
do
    fastp \
        -i ${i}.R1.fastq.gz \
        -I ${i}.R2.fastq.gz \
        -o ${RAW_OUT}/${i}.trim.fastq.gz \
        -O ${RAW_OUT}/${i}.R.trim.fastq.gz \
        -w 20
done

####2. mapping + dedup + bigwig + count   ############


#!/bin/bash
#PBS -q regular
#PBS -N mapping_counts
#PBS -l nodes=1:ppn=4
#PBS -l walltime=10000:00:00
#PBS -l mem=40gb

cd /scratch/hbgao/MYBS_RNAseq/raw

GENOME=/project4/genomes/Ta-ensemble/Hisat2_v1.0/iwgsc1_tran
GTF=/home/zefulu/genomes/Ta-ensemble/Ta.gtf

sample_list="
G-74-D1_L4_393X93
G-74-D2_L4_378X78
G-74-D3_L4_379X79
G-74-N1_L4_380X80
G-74-N2_L4_381X81
G-74-N3_L4_394X94
OX4-D1_L4_382X82
OX4-D2_L4_395X95
OX4-D3_L4_383X83
OX4-N1_L4_384X84
OX4-N2_L4_385X85
OX4-N3_L4_386X86
WT-D1_L4_388X88
WT-D2_L4_389X89
WT-D3_L4_390X90
WT-N1_L4_391X91
WT-N2_L4_396X96
WT-N3_L4_392X92
"

for i in $sample_list
do
    echo "Processing $i ..."

    ## 1) Mapping
    hisat2 -p 20 -x $GENOME \
        -1 ${i}.trim.fastq.gz \
        -2 ${i}.R.trim.fastq.gz \
        -S ${i}.sam

    ## 2) Convert + Filter
    samtools view -q 10 -b ${i}.sam -o ${i}.bam
    samtools sort ${i}.bam -o ${i}.sort.bam

    ## 3) Remove duplicates
    java -jar /opt/picard-2.21.9/picard.jar MarkDuplicates \
        INPUT=${i}.sort.bam \
        OUTPUT=${i}.clean.bam \
        METRICS_FILE=${i}.metrics.txt \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT

    samtools index -m 14 ${i}.clean.bam

    ## 4) bigwig
    bamCoverage --bam ${i}.clean.bam -o ${i}.clean.bw

    ## 5) Count
    featureCounts -T 20 -a $GTF \
        -o ${i}.count \
        -p -B -C -g gene_id \
        ${i}.clean.bam

    ## optional cleanup
    rm ${i}.sam
    rm ${i}.bam
    rm ${i}.sort.bam

done

###### 3.extract read count of each samples  ####################
for i in *.count
do
    out=$(basename $i .count)
    grep -v '#' $i | cut -f7 > ${out}.readcounts
done

# 合并成一个矩阵
paste names *.readcounts > all.readc

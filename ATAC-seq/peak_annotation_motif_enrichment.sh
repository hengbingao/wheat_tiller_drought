#!/bin/bash
#PBS -q regular
#PBS -N peak_annotation_and_fimo
#PBS -l nodes=1:ppn=4
#PBS -l walltime=96:00:00
#PBS -l mem=40gb

########################################
# Working directory
########################################
cd /scratch/peaks_file


########################################
# 1. Peak annotation with HOMER
########################################
for i in $(ls *.narrowPeak | sed 's/.narrowPeak//')
do
    annotatePeaks.pl ${i}.narrowPeak \
        /scratch/hbgao/wheat-ref/Ta.fa \
        -gtf /scratch/hbgao/wheat-ref/Ta.gtf \
        -annStats ${i}.stats.txt \
        > ${i}.anno.xls
done


########################################
# 2. Generate shuffled control regions (within merged ACR regions)
########################################
for i in $(ls *.bed | sed 's/.bed//')
do
    bedtools shuffle \
        -i ${i}.bed \
        -g /scratch/hbgao/wheat-ref/Ta.chrom.size \
        -incl /scratch/hbgao/ATAC_all/Tn5/DEseq2/diff_peak/sizefactor/final_peak/merged_peak_delMYBS.bed \
        > ${i}.control.inACR.bed
done


########################################
# 3. Run FIMO motif scanning
########################################
for i in $(ls *.bed | sed 's/.bed//')
do
    bedtools getfasta \
        -fi /scratch/hbgao/wheat-ref/Ta.fa \
        -bed ${i}.bed \
        > ${i}.fa

    fimo --oc ${i} \
        /scratch/hbgao/wheat-ref/jasper_meme/root_cluster.meme \
        ${i}.fa
done


########################################
# 4. Summarize motif enrichment (motif hit counts)
########################################
FIMO_OUT=/scratch/hbgao/ATAC_all/Tn5/DEseq2/diff_peak/sizefactor/final_peak/Kmeans/fimo

for i in $(ls *.bed | sed 's/.bed//')
do
    cd ${FIMO_OUT}/${i}

    cut -f1 fimo.tsv \
        | sort \
        | uniq -c \
        | grep 'cluster' \
        | grep -v 'fimo' \
        > ${FIMO_OUT}/${i}.txt
done

cd ${FIMO_OUT}
paste *.txt > all_fimo.txt
ls *.txt | sed 's/.txt//' | tr '\n' '\t' > header.txt


########################################
# 5. FIMO on regions excluding ACRs
########################################
FIMO_NOACR=/scratch/hbgao/ATAC_all/Tn5/DEseq2/diff_peak/sizefactor/final_peak/Kmeans/fimo_noACR

for i in $(ls *.bed | sed 's/.bed//')
do
    bedtools shuffle \
        -i ${i}.bed \
        -g /scratch/hbgao/wheat-ref/Ta.chrom.size \
        -excl /scratch/hbgao/ATAC_all/Tn5/DEseq2/diff_peak/sizefactor/final_peak/merged_peak.bed \
        > ${i}.control.noACR.bed
done

for i in $(ls *.bed | sed 's/.bed//')
do
    bedtools getfasta \
        -fi /scratch/hbgao/wheat-ref/Ta.fa \
        -bed ${i}.bed \
        > ${i}.fa

    fimo --oc ${i} \
        /scratch/hbgao/wheat-ref/jasper_meme/root_cluster.meme \
        ${i}.fa
done

for i in $(ls *.bed | sed 's/.bed//')
do
    cd ${FIMO_NOACR}/${i}

    cut -f1 fimo.tsv \
        | sort \
        | uniq -c \
        | grep 'cluster' \
        | grep -v 'fimo' \
        > ${FIMO_NOACR}/${i}.txt
done

cd ${FIMO_NOACR}
paste *.txt > all_fimo.txt
ls *.txt | sed 's/.txt//' | tr '\n' '\t' > header.txt


########################################
# 6. Subgenome distribution of differential peaks
########################################
for i in $(ls *.bed | sed 's/.bed//')
do
    cut -f1 ${i}.bed \
        | sed -E 's/[1-7]//g' \
        | sort \
        | uniq -c \
        > ${i}.ABD
done


########################################
# 7. Prepare BED file with peak lengths
########################################
for i in $(ls *.bed | sed 's/.bed//')
do
    awk '{ print $1 "\t" $2 "\t" $3 "\t" ($3 - $2) }' ${i}.bed \
        > ${i}.peak_length.bed
done

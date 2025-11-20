#####peak annotation and TF root motif enrichment#####

#!/bin/bash 
#PBS -q regular 
#PBS -N peak_anno
#PBS -l nodes=1:ppn=4 
#PBS -l walltime=96:00:00 
#PBS -l mem=40gb
cd /scratch/peaks_file
for i in (ls *.narrowPeak|sed 's/.narrowPeak//g')
do
annotatePeaks.pl $i.narrowPeak /scratch/hbgao/wheat-ref/Ta.fa -gtf /scratch/hbgao/wheat-ref/Ta.gtf -annStats $i.txt > $i.anno.xls
done

#########################motif fimo #################
###################shuffle 一些control序列
for i in $(ls *.bed|sed "s/.bed/\t/g")
do
bedtools shuffle -i $i.bed -g /scratch/hbgao/wheat-ref/Ta.chrom.size -incl /scratch/hbgao/ATAC_all/Tn5/DEseq2/diff_peak/sizefactor/final_peak/merged_peak_delMYBS.bed > $i.nochange.control.bed
done

###################用fimo来寻找相关的位点#################
for i in $(ls *.bed|sed "s/.bed/\t/g")
do
bedtools getfasta -fi /scratch/hbgao/wheat-ref/Ta.fa -bed $i.bed > $i.fa
fimo --oc ${i} /scratch/hbgao/wheat-ref/jasper_meme/root_cluster.meme ${i}.fa
done


####################整理相关的富集次数#############
for i in $(ls *.bed|sed "s/.bed/\t/g")
do
cd /scratch/hbgao/ATAC_all/Tn5/DEseq2/diff_peak/sizefactor/final_peak/Kmeans/fimo/$i
cut -f1 fimo.tsv|sort|uniq -c|grep -e 'cluster'|grep -v 'fimo' > /scratch/hbgao/ATAC_all/Tn5/DEseq2/diff_peak/sizefactor/final_peak/Kmeans/fimo/$i.txt
done

paste *.txt > all_fimo.txt

ls *.txt|sed 's/.txt//g'|sed 's/\n/\t/g' > header

19508	12312	61871	44183	139959	139494	178940	223905	37479	27168	62932	56081	37302	37680	41453	42546

cat fimo_file.txt|awk '{print $1 "\t" $1/($1-19508)}'






#########################fimo noACR 
for i in $(ls *.bed|sed "s/.bed/\t/g")
do
bedtools shuffle -i $i.bed -g /scratch/hbgao/wheat-ref/Ta.chrom.size -excl /scratch/hbgao/ATAC_all/Tn5/DEseq2/diff_peak/sizefactor/final_peak/merged_peak_delMYBS.bed > $i.nochange.control.bed
done


for i in $(ls *.bed|sed "s/.bed/\t/g")
do
bedtools getfasta -fi /scratch/hbgao/wheat-ref/Ta.fa -bed $i.bed > $i.fa
fimo --oc ${i} /scratch/hbgao/wheat-ref/jasper_meme/root_cluster.meme ${i}.fa
done



for i in $(ls *.bed|sed "s/.bed/\t/g")
do
cd /scratch/hbgao/ATAC_all/Tn5/DEseq2/diff_peak/sizefactor/final_peak/Kmeans/fimo_noACR/$i
cut -f1 fimo.tsv|sort|uniq -c|grep -e 'cluster'|grep -v 'fimo' > /scratch/hbgao/ATAC_all/Tn5/DEseq2/diff_peak/sizefactor/final_peak/Kmeans/fimo_noACR/$i.txt
done

paste *.txt > all_fimo.txt

ls *.txt|sed 's/.txt//g'|sed 's/\n/\t/g' > header



################diff peak subgenome distribution
for i in $(ls *.bed|sed "s/.bed/\t/g")
do
cut -f1 $i.bed|sed -re 's/[1-7]//p'|sort|uniq -c > $i.ABD
done



for i in $(ls *.bed|sed "s/.bed/\t/g")
do
awk '{print $1 "\t" $2 "\t" $3 "\t" $3-$2}' $i.bed > $i.peak.bed
done

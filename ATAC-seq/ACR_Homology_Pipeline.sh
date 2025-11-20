#!/bin/bash
set -e

############################################################
### Step 1. Extract wheat homologous genes
############################################################

cat Ta.collinearity | grep -v '#' | cut -f2,3 > Ta.homolog
sed 's/\./\t/g' Ta.homolog | cut -f1,3 > Ta.homolog_gene


############################################################
### Step 2. Process ACR subgenome coordinate files (.gene)
############################################################

for i in $(ls *.txt | sed 's/.txt//'); do
    sed 's/gene://g' $i.txt \
    | awk '{print $1 ":" $2 "-" $3 "\t" $4 }' > $i.gene
done

for i in $(ls *.gene | sed 's/.gene//'); do
    awk '{print $1 "\t" $2 "\t" "'"$i"'" }' $i.gene > $i.gene2
done


############################################################
### Step 3. Generate ACR BED files
############################################################

cat ACR_hom_fasta.txt \
    | sed 's/:/\t/g' | sed 's/-/\t/g' \
    | bedtools sort -i - | uniq \
    > ACR_hom_fasta_loci.bed


############################################################
### Step 4. Extract sequences for all BED files
############################################################

for i in $(ls *.bed | sed 's/.bed//'); do
    bedtools getfasta -fi /scratch/hbgao/wheat-ref/Ta.fa \
                      -bed $i.bed > $i.fa
done


############################################################
### Step 5. Build BLAST databases (subA, subB, subD)
############################################################

makeblastdb -in ACR_hom_fasta_subA.fa -dbtype nucl -parse_seqids -out ACR_subA
makeblastdb -in ACR_hom_fasta_subB.fa -dbtype nucl -parse_seqids -out ACR_subB
makeblastdb -in ACR_hom_fasta_subD.fa -dbtype nucl -parse_seqids -out ACR_subD

blastn -query ACR_hom_fasta_subA.fa -out subA_B.blast -db ACR_subB -outfmt 6 -evalue 1e-5
blastn -query ACR_hom_fasta_subA.fa -out subA_D.blast -db ACR_subD -outfmt 6 -evalue 1e-5
blastn -query ACR_hom_fasta_subB.fa -out subB_D.blast -db ACR_subD -outfmt 6 -evalue 1e-5


############################################################
### Step 6. Extract p-values & bitscores
############################################################

for i in $(ls *.blast | sed 's/.blast//'); do
    cut -f1,2,11 $i.blast > $i.pvalue
    cut -f1,2,12 $i.blast > $i.score
done

cat subA_B.score subA_D.score subB_D.score > ACR_homo_score


############################################################
### Step 7. Compute conserved peak fragments
############################################################

for i in $(ls *.blast | sed 's/.blast//'); do
    cut -f1,7,8 $i.blast \
    | sed 's/:/\t/g' | sed 's/-/\t/g' \
    | awk '{print $1 "\t" $2 + $4 "\t" $2 + $5 "\t" $5 - $4}' \
    > $i.ABD_conserved.bed
done


############################################################
### Step 8. Identify pACR and dACR peaks
############################################################

grep -e 'pC' ACR_HC_hom.txt \
    | sed 's/:/\t/g' | sed 's/-/\t/g' | cut -f1,2,3 \
    | bedtools sort -i - | uniq > pACR_peak.bed

grep -e 'dC' ACR_HC_hom.txt \
    | sed 's/:/\t/g' | sed 's/-/\t/g' | cut -f1,2,3 \
    | bedtools sort -i - | uniq > dACR_peak.bed

# Split by subgenome
for t in pACR dACR; do
    grep -e 'A' ${t}_peak.bed > ${t}_subA_peak.bed
    grep -e 'B' ${t}_peak.bed > ${t}_subB_peak.bed
    grep -e 'D' ${t}_peak.bed > ${t}_subD_peak.bed
done


############################################################
### Step 9. Extract sequences for pACR/dACR peaks
############################################################

for i in $(ls *.bed | sed 's/.bed//'); do
    bedtools getfasta -fi /scratch/hbgao/wheat-ref/Ta.fa \
                      -bed $i.bed > $i.fa
done


############################################################
### Step 10. pACR BLAST database building + BLAST
############################################################

makeblastdb -in pACR_subA_peak.fa -dbtype nucl -parse_seqids -out pACR_subA
makeblastdb -in pACR_subB_peak.fa -dbtype nucl -parse_seqids -out pACR_subB
makeblastdb -in pACR_subD_peak.fa -dbtype nucl -parse_seqids -out pACR_subD

blastn -query pACR_subA_peak.fa -out subA_B.blast -db pACR_subB -outfmt 6 -evalue 1e-5
blastn -query pACR_subA_peak.fa -out subA_D.blast -db pACR_subD -outfmt 6 -evalue 1e-5
blastn -query pACR_subB_peak.fa -out subB_D.blast -db pACR_subD -outfmt 6 -evalue 1e-5
blastn -query pACR_subD_peak.fa -out subD_B.blast -db pACR_subB -outfmt 6 -evalue 1e-5

for i in $(ls *.blast | sed 's/.blast//'); do
    cut -f1,2,11,12 $i.blast > $i.score
done

cat subA_B.score subA_D.score subB_D.score subD_B.score > pACR_homo_score


############################################################
### Step 11. dACR BLAST database building + BLAST
############################################################

makeblastdb -in dACR_subA_peak.fa -dbtype nucl -parse_seqids -out dACR_subA
makeblastdb -in dACR_subB_peak.fa -dbtype nucl -parse_seqids -out dACR_subB
makeblastdb -in dACR_subD_peak.fa -dbtype nucl -parse_seqids -out dACR_subD

blastn -query dACR_subA_peak.fa -out subA_B.blast -db dACR_subB -outfmt 6 -evalue 1e-5
blastn -query dACR_subA_peak.fa -out subA_D.blast -db dACR_subD -outfmt 6 -evalue 1e-5
blastn -query dACR_subB_peak.fa -out subB_D.blast -db dACR_subD -outfmt 6 -evalue 1e-5
blastn -query dACR_subD_peak.fa -out subD_B.blast -db dACR_subB -outfmt 6 -evalue 1e-5

for i in $(ls *.blast | sed 's/.blast//'); do
    cut -f1,2,12 $i.blast > $i.score
done

cat subA_B.score subA_D.score subB_D.score > dACR_homo_score


############################################################
### Step 12. Generate conserved peak fragments for pACR/dACR
############################################################

for i in $(ls *.blast | sed 's/.blast//'); do
    cut -f1,7,8 $i.blast \
    | sed 's/:/\t/g' | sed 's/-/\t/g' \
    | awk '{print $1 "\t" $2 + $4 "\t" $2 + $5 "\t" $5 - $4}' \
    | awk '{print $1 ":" $2 "-" $3 "\t" $4}' \
    > $i.fragment
done

for i in $(ls *.blast | sed 's/.blast//'); do
    paste $i.score $i.fragment > $i.ACR_homo
done


############################################################
### Step 13. Filter high-confidence (HC) ACRs
############################################################

awk '{if($7 == $12) print $0}' dACR_homo_score > dACR_homo_HC_score

for i in $(ls *.csv | sed 's/.csv//'); do
    sed 's/"//g' $i.csv | sed 's/,/\t/g' > $i
done

awk '{if($23 == $24) print $0}' pACR_ABD_gene_anno \
| awk '{if($24 == $25) print $0}' \
> pACR_HC_ABD_bitscore_anno

awk '{if($22 == $23) print $0}' dACR_ABD_gene_anno \
| awk '{if($23 == $24) print $0}' \
> dACR_HC_ABD_bitscore_anno


############################################################
### Step 14. Merge conserved peak loci
############################################################

for i in $(ls *.txt | sed 's/.txt//'); do
    sed 's/:/\t/g' $i.txt \
    | sed 's/-/\t/g' \
    | bedtools sort -i - \
    | bedtools merge -i - \
    | awk '{print $1 "\t" $2 "\t" $3 "\t" $3-$2}' \
    > $i.bed
done

for type in pACR_conserved_peak dACR_conserved_peak; do
    grep -e 'A' $type.bed > $type.A.bed
    grep -e 'B' $type.bed > $type.B.bed
    grep -e 'D' $type.bed > $type.D.bed
done


############################################################
### Step 15. Annotate conserved peaks
############################################################

for i in $(ls *.bed | sed 's/.bed//'); do
    annotatePeaks.pl $i.bed /scratch/hbgao/wheat-ref/Ta.fa \
        -gtf /scratch/hbgao/wheat-ref/Ta.gtf \
        -annStats $i.txt \
        > $i.anno.xls
done


############################################################
### Step 16. Extract intergenic distances
############################################################

for i in $(ls *.xls | sed 's/.anno.xls//'); do
    grep -e 'Intergenic' $i.anno.xls \
    | cut -f10 \
    | sed 's/-//g' \
    | awk '{print $1 "\t" "'"$i"'" }' \
    | sed 's/\./_/g' \
    > $i.intergenic.distance
done


############################################################
### Step 17. All steps completed
############################################################

echo "All steps completed successfully."

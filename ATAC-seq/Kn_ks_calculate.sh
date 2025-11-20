#!/bin/bash

###############################################
# 1. Compute Kn/Ks for AXT files
###############################################

# Example single-file run
KnKs -i Triticum_urartu.Triticum_aestivum.sort.axt \
     -j 0.006 \
     -o noncoding.axt.knks

# Batch run for all *.axt
for i in $(ls *.axt | sed 's/.axt//' | sort)
do
    KnKs -i ${i}.axt -j 0.05 -o ${i}.knks
done


###############################################
# 2. Extract useful columns from Kn/Ks results
#    and generate subgenome-specific files
###############################################

# AB comparison
grep -v 'Kn' Ta_self.chrAB.knks \
    | awk '{print $2 "\t" $3 "\t" $4 "\t" $10}' \
    > Ta.chrAB_chrA.knks

# BD comparison
grep -v 'Kn' Ta_self.chrBD.knks \
    | awk '{print $2 "\t" $3 "\t" $4 "\t" $10}' \
    > Ta.chrBD_chrB.knks

# AD comparison
grep -v 'Kn' Ta_self.chrAD.knks \
    | awk '{print $5 "\t" $6 "\t" $7 "\t" $10}' \
    > Ta.chrAD_chrD.knks


###############################################
# 3. Merge all Kn/Ks data and convert to bedGraph
###############################################

cat Ta.chrAB_chrA.knks Ta.chrBD_chrB.knks Ta.chrAD_chrD.knks \
    | sort -k1,1 -k2,2n \
    > Kn.merge.bg

# Optional: generate raw genome coverage (if needed)
# bedtools genomecov -i Kn.bed -split -bg -g /scratch/hbgao/wheat-ref/Ta.genome.sizes > Kn.bg


###############################################
# 4. Map Kn values to specific loci
###############################################

bedtools map \
    -a Kn.loci.bed \
    -b Kn.merge.bg \
    -c 4 -o sum \
    > Kn.uniq.bg


###############################################
# 5. Convert to bigWig
###############################################

bedGraphToBigWig \
    Kn.uniq.bg \
    /scratch/hbgao/wheat-ref/Ta.genome.sizes \
    Kn.bw


###############################################
# 6. deepTools: signal alignment on ACR categories
###############################################

# Differential ACR vs. unchanged ACR
computeMatrix reference-point \
    -S Kn.bw \
    -R pACR_diffpeak.bed pACR.nochange.bed dACR_diffpeak.bed dACR.nochange.bed \
    --referencePoint center \
    -a 2000 -b 2000 \
    -o diffpeak_Kn.gz \
    -p 30 \
    -bs 20

plotProfile \
    -m diffpeak_Kn.gz \
    -out diffpeak_Kn.pdf \
    --perGroup \
    --numPlotsPerRow 1


# All ACR vs shuffled background
computeMatrix reference-point \
    -S Kn.bw \
    -R all_dACR.bed all_dACR.shuffle.bed all_pACR.bed all_pACR.shuffle.bed \
    --referencePoint center \
    -a 2000 -b 2000 \
    -o all_ACR_Kn.gz \
    -p 30 \
    -bs 20

plotProfile \
    -m all_ACR_Kn.gz \
    -out allACR_Kn.pdf \
    --perGroup \
    --numPlotsPerRow 1


###############################################
# 7. Extract Kn values overlapping peak categories
###############################################

bedtools intersect -a pACR_diffpeak.bed -b Kn.merge.bg -wb \
   | awk '{print $7 "\t" "pACR_diffpeak"}' \
   > pACR_diff_Kn

bedtools intersect -a pACR.nochange.bed -b Kn.merge.bg -wb \
   | awk '{print $7 "\t" "pACR_nochange"}' \
   > pACR_nochange_Kn

bedtools intersect -a dACR_diffpeak.bed -b Kn.merge.bg -wb \
   | awk '{print $7 "\t" "dACR_diffpeak"}' \
   > dACR_diffpeak_Kn

bedtools intersect -a dACR.nochange.bed -b Kn.merge.bg -wb \
   | awk '{print $7 "\t" "dACR_nochange"}' \
   > dACR_nochange_Kn

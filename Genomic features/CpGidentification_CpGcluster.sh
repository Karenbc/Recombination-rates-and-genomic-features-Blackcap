#!/bin/bash
#


#####CpGisland identification from blackcap reference genome using the program CpGcluster#####



#Tools:
module load perl/5.24.1
samtools
CpGcluster

#Folders

#Reference genome path
IN=/home/bascon.c/Reference_genome/NEW_Reference_genome

#The output should be a directory  containing sequence files in FASTA format for each chromosmes
OUT=/home/bascon.c/CpG_GC/CpGcluster/NEW_Ref/chr_ss_bSylAtr1_1_fasta




## Extract chromosomes in FASTA format from the reference genome 

for a in  1 2 3 Z 4 5 6 7 8 9 10 11 12 W 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 w_unkown 33 

do
samtools faidx $IN/bSylAtr1.1.fasta chr_${a} > $OUT/chr${a}.fa

echo "chr_${a} done"

done


##Run the program with recommended parameters (d=50, p-value=1e-5)

#perl CpGcluster.pl <assembly>  <d>  <P-value>

    #d is the threshold distance on basis of a given percentile. It is recommendable to use 50 (median distance)
    #P-value is the maximal P-value under which a CpG cluster is considered as a CpG island. The recommended limit is 1E-5


perl ~/Bin/CpGcluster/CpGcluster.pl $OUT 50 1E-5

echo "CpGcluster chromosomes done"





##MERGE CHROMOSOMES##

for a in  1 2 3 Z 4 5 6 7 8 9 10 11 12 W 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 w_unkown 33 

do

#Remove first line (header) of CpG cluster output
tail -n +2 $OUT/chr${a}_p50_CpGcluster.txt > $OUT/chr${a}_p50_CpGcluster_nohead.txt

#Add chromosome name

awk -v OFS='\t' {'print "chr_'"${a}"'",$2,$3,$4,$5,$6,$7,$8,$9,$10,$11'} $OUT/chr${a}_p50_CpGcluster_nohead.txt > $OUT/chr${a}_p50_CpGcluster.chr

echo "chr_${a} name done"
done

#Merge all chromosomes 

awk 'NF' $OUT/chr*_p50_CpGcluster.chr > All_chr_p50_CpGcluster_chrname_merged_head





##FILTER CpGi FROM CPGcluster OUTPUT##


#Filter out CpGs shorter or equal to 50bp
awk '$5 >= 50 {print $0}' All_chr_p50_CpGcluster_chrname_merged_head > All_chr_p50_CpGcluster_chrname_merged_head_Filt_50Lng

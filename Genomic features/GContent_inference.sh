#!/bin/bash
#
# 






#Tools:
module load perl/5.24.1
#geecee
#GCProfile (http://tubic.tju.edu.cn/GC-Profile/download/readme.txt)


#my folders



# Directory containing sequence files in FASTA format for each chromosmes
IN=/home/bascon.c/CpG_GC/GContent/GCProfile/chr_ss_bSylAtr1_1_fasta
#Output path
OUT=/home/bascon.c/CpG_GC/GContent/Geecee_Out
#Reference genome
In_Ref=/home/bascon.c/Reference_genome/NEW_Reference_genome
#
OUT_GCprofile=/home/bascon.c/CpG_GC/GContent/GCProfile/Ref_WholeGenome


##CALCULATE GC content for each chromosome##

for a in  1 2 3 Z 4 5 6 7 8 9 10 11 12 W 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 w_unkown 33 

do
geecee -sequence $IN/chr${a}.fasta    -outfile $OUT/${a}
echo "chromosome_${a} GC content done"
done



##CALCULATE GC content for the whole genome with GCprofile##

#The reference genome should be wthin the same folder as the output
cp $In_Ref/bSylAtr1.1.fasta $OUT_GCprofile
cd $OUT_GCprofile
GCProfile bSylAtr1.1.fasta -t 300 -m
#t : threshold for segmentation
#m : multiplot mode,  plots are placed on the same page

##The output are values for the whole genome , not specifing the chromosomes. 
#To create a data table with the inferences specifiying each chromosome:

for a in  1 2 3 Z 4 5 6 7 8 9 10 11 12 W 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33
do
#create column with 2 positions

paste -sd "\t\n" $IN/chr${a}.SegGC > $IN/GCProfile_output_2positions/chr${a}.SegGC_2pos

#add chromosome names

awk -v OFS='\t' {'print "chr_'"${a}"'",$1,$3,$2'} $IN/GCProfile_output_2positions/chr${a}.SegGC_2pos > $IN/GCProfile_output_2positions/chr_${a}.SegGC_2pos_chr

done

#merge all chr files in one

awk 'NF' $IN/GCProfile_output_2positions/chr*.SegGC_2pos_chr > $IN/GCProfile_output_2positions/chr.SegGC_2pos_chr.WholeGenome

#remove missing values "

grep -Eiv '(")' < $IN/GCProfile_output_2positions/chr.SegGC_2pos_chr.WholeGenome > $IN/GCProfile_output_2positions/chr.SegGC_2pos_chr.WholeGenome_nowWZ.nomissing.bed


##Calculate GC% in bins --> weighted average GC content in different winow sizes


GCbin=/home/bascon.c/CpG_GC/GContent/GCProfile/chr_ss_bSylAtr1_1_fasta/GCProfile_output_2positions/GContent_bins/Weght_Avrg_GC
WinGenome=/home/bascon.c/Genome_Annotation/New_Ref/windows_bedtools


# weighted mean (average taking into account physical distance between region where GC was calculated)
script=/home/bascon.c/scripts/Recombination_Analisis/Rec_Rates

$script/meanweighted_recombination_inwindows.py  $IN/GCProfile_output_2positions/chr.SegGC_2pos_chr.WholeGenome_nowWZ.nomissing.bed 200000 > $GCbin/GC_AvrgWgth_200kb_WholeGenome.bed

$script/meanweighted_recombination_inwindows.py  $IN/GCProfile_output_2positions/chr.SegGC_2pos_chr.WholeGenome_nowWZ.nomissing.bed 500000 > $GCbin/GC_AvrgWgth_500kb_WholeGenome.bed

$script/meanweighted_recombination_inwindows.py  $IN/GCProfile_output_2positions/chr.SegGC_2pos_chr.WholeGenome_nowWZ.nomissing.bed 1000000 > $GCbin/GC_AvrgWgth_1Mb_WholeGenome.bed






#This (following) was not used int he end
#In 1Mb500kbOvW
#bedtools intersect -a $WinGenome/genome.1Mb_windows.500kbovW.nowWZ.bed  -b $IN/GCProfile_output_2positions/chr.SegGC_2pos_chr.WholeGenome_nowWZ.nomissing.bed -wa -wb > $GCbin/GContent_1Mb500kbOvW.wawb.bed

#the meanGCcontent per window

#sort -k1,1 -k2,2n $IN/GCProfile_output_2positions/chr.SegGC_2pos_chr.WholeGenome_nowWZ.nomissing.bed > $IN/GCProfile_output_2positions/chr.SegGC_2pos_chr.WholeGenome_nowWZ.nomissing.sort.bed

#bedtools map -a $WinGenome/genome.1Mb_windows.500kbovW.nowWZ.sort.bed -b $IN/GCProfile_output_2positions/chr.SegGC_2pos_chr.WholeGenome_nowWZ.nomissing.sort.bed -c 4 -o mean > $GCbin/GContent_1Mb500kbOvW.mean.bed

#!/bin/bash
#
# 
###ESTIMATION OF GC CONTENT IN THE BLACKAP GENOME######


#Tools:
module load perl/5.24.1
#geecee (https://www.bioinformatics.nl/cgi-bin/emboss/help/geecee)
#GCProfile (http://tubic.tju.edu.cn/GC-Profile/download/readme.txt)


#my folders

# Directory containing sequence files in FASTA format for each chromosmes
IN=/home/bascon.c/CpG_GC/chr_ss_bSylAtr1_1_fasta

#Output path for Geecee
OUT=/home/bascon.c/CpG_GC/GContent/Geecee_Out

#Reference genome
In_Ref=/home/bascon.c/Reference_genome/NEW_Reference_genome

#Output path for GCProfile
OUT_GCprofile=/home/bascon.c/CpG_GC/GContent/GCProfile/Ref_WholeGenome



##CALCULATE GC CONTENT PER CHROMOSOME##

for a in  1 2 3 Z 4 5 6 7 8 9 10 11 12 W 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 

do
geecee -sequence $IN/chr${a}.fasta    -outfile $OUT/${a}
echo "chromosome_${a} GC content done"
done



##CALCULATE GC CONTENT FOR THE WHOLE GENOME WITH GCprofile##

#The reference genome should be wIthin the same folder as the output
cp $In_Ref/bSylAtr1.1.fasta $OUT_GCprofile
cd $OUT_GCprofile
GCProfile bSylAtr1.1.fasta -t 300 -m
   #t : threshold for segmentation
   #m : multiplot mode,  plots are placed on the same page

##The output are values for the whole genome , not specifing the chromosomes. 
#To create a whole genome data table specifiying chromosomes:

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

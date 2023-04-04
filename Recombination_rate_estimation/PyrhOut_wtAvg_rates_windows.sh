#Calculate from pyrho output the weighted average recombination rate in non-overlapping windows using a python script

#Input file (rec rates converted to cM/Mb)
In_rec=/home/Cont_Res0
#Path for output file
Rec_bin=/home/bascon.c/Recombination_Analysis/New_Ref/Rec_Bins_AvrgWeight
#Script to calculate rec rates in windows
script=/home/scripts/Recombination_Analisis/Rec_Rates

#Modules:
module load python/3.6.0 


#Code:

#50kb windows

$script/Rec_windows.py  $In_rec/WholeGenome_chr_ContRes0_n38_Pen20_W50 50000 > $Rec_bin/OutPyrho_r_ContRes0_pen20W50_AvrgWgth_50kb


#200kb windows

$script/Rec_windows.py  $In_rec/WholeGenome_chr_ContRes0_n38_Pen20_W50 200000 > $Rec_bin/OutPyrho_r_ContRes0_pen20W50_AvrgWgth_200kb


#1Mb windows

$script/Rec_windows.py  $In_rec/WholeGenome_chr_ContRes0_n38_Pen20_W50 1000000 > $Rec_bin/OutPyrho_r_ContRes0_pen20W50_AvrgWgth_1Mb


#With pyrho output calculate weighted average recombination rate in non-overlapping windows using python script

In_rec=/home/bascon.c/PYRHO/Optimize_RecRates/NEW_Ref/Continent_Residents/Cont_Res0
Rec_bin=/home/bascon.c/Recombination_Analysis/New_Ref/Rec_Bins_AvrgWeight
script=/home/bascon.c/scripts/Recombination_Analisis/Rec_Rates


module load python/3.6.0 

#200kb windows

$script/Rec_windows.py  $In_rec/WholeGenome_chr_ContRes0_n38_Pen20_W50 200000 > $Rec_bin/OutPyrho_r_ContRes0_pen20W50_AvrgWgth_200kb


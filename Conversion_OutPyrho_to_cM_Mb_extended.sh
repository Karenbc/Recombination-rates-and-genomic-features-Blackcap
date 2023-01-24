#!/bin/bash


#my folders

IN=Optimize_RecRates_Output
OUT=ContRes0_n38_Pen20_W50_convert_2Pos
Script=/home/bascon.c/scripts/Pyrho/Coversions




#The format of pyrho output is (chr_3 example): 

head chr_3_ContRes0_n38_Pen20_W50.rmap

74      1817    4.770579636540578e-09
1817    4057    4.770579636540578e-09
4057    9555    4.770579636540578e-09
9555    9573    4.770579636540578e-09
9573    9575    4.770579636540578e-09
9575    9585    4.770579636540578e-09
9585    9590    4.770579636540578e-09
9590    9593    4.770579636540578e-09
9593    9595    4.770579636540578e-09
9595    9597    4.770579636540578e-09

first position: start position or interval
second column: end position or interval
third column: r (per-site per-generation recombination rate)



#The conversion of pyrho output (r) to cM and cM/Mb is implemented in the script called convert_Output_pyrho_cM_Mb.awk. 
#The details of the conversion:

        W=$2-$1 #interval or window size
	r=$3 #r
	cM=100*W*r # due to cM
	cMb=r*1e+8 #cM/Mb is calculated as 100*W*r/1e-6*w , simplfied is equal to r*1e+8
	Map+=100*W*r #cumulative cM 
	




## Command to convert Pyrho output to cM/Mb and calculate genetic cumulative map (cM)  (using script : $Script/convert_Output_pyrho_cM_Mb.awk) for each chromosome.

for a in  1 2 3 Z 4 5 6 7 8 9 10 11 12 W 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 
do


awk -f $Script/convert_Output_pyrho_cM_Mb.awk $IN/chr_${a}_ContRes0_n38_Pen20_W50.rmap > $OUT/Continent_Residents/Cont_Res0/cM_GenMap/Conv_startEndpos/chr_${a}_ContRes0_n38_Pen20_W50_convert_2Pos_W_Final.txt


echo "chr${a} from r to cM/Mb and Genetic Map (cM)"

#The output contains first column: start position, second column: end position, third column: cM/Mb, fourth column: genetic Map, fifth column: interval where recombination rate was calculated  





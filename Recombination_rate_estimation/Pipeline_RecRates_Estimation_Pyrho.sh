#!/bin/bash
#
#SBATCH --job-name=ContRes_Recombination_Pyrho
#SBATCH --ntasks=4
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=120:05:00
#  maximum requested memory
#SBATCH --mem=70G
#  write std out and std error to these files
#SBATCH --error=/home/bascon.c/stdout/LookUptable_Hyperparam_Pyrho_ContRes_n38_N48.%J.err
#SBATCH --output=/home/bascon.c/stdout/LookUptable_Hyperparam_Pyrho_ContRes_n38_N48.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bascon.c@evolbio.mpg.de
#SBATCH --partition=standard

#  add your code here:


#mytools

 module load python/3.5.3


#my folders
demogTable=/home/bascon.c/demography
OUTPUT_table=/home/bascon.c/Pyrho/LookupTable
Out_Hyperparam=/home/bascon.c/Rec_LocalPCA_Jun/Islands/Rec_Pyrho/Azores/Hyperparam


#Load source

source /home/bascon.c/tutorial-env/bin/activate

#Step 1# Pre-compute a lookup table for our 19 diploid individuals (38 haploids) 

pyrho make_table -n 38 -N 48 --mu 4.6e-9 --outfile $OUTPUT_table/Continental_Residents_n38_N48_lookuptable.hdf --approx --msmc_file $demogTable/residentcontinent_sample1.final.txt --decimate_rel_tol 0.1


echo "LookupTable n38 N48 done"


#Step 2# Hyperparam: find parameter settings that are optimal for the demography data. This command simulates data and performs optimization under a number of different hyperparameter settings and outputs the accuracy in terms of a number of different metrics.

pyrho hyperparam -n 38 --mu 4.6e-9  --ploidy 2 --blockpenalty 10,20,50,100 --windowsize 20,25,50,100 --logfile . --tablefile $OUTPUT_table/Continental_Residents_n38_N48_lookuptable.hdf --num_sims 3 --msmc_file $demogTable/residentcontinent_sample1.final.txt  --outfile $OUTPUT_Hyperparam/Continental_Residents_n38_hyperparam_results_1.txt 

echo "Hyperparam n38 N48 done"

#Optimize

#Recombination rates inference run with perl script: ContRes0_Pyrho_Optimize.pl



#converted the Output of pyrho to cM/Mb and genetic map (cM)




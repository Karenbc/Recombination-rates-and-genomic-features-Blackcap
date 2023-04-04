This repository holds the scripts used to estimate recombination rates and annotate the genome. 
It contains three main folders:

i)   Calling Variants: It contains the scripts used to call SNPs from re-sequencing data.

ii)  Recombination rates: The rates of recombination were inferred using pyrho (Kamm et al., 2016; Spence & Song, 2019). 
The script called Pipeline_RecRates_Estimation_Pyrho.sh contains the computation of the lookup-table and the simulation step, named hyperparam, to find optimal parameters (window size and block penalty) for the estimation. 
The inference of recombination rates from real data is included in the script named RecRates_inference_optimize.pl. This script generates several parallel jobs corresponding to each chromosome from where recombination rate is inferred. 

The output of pyrho (r,recombination rate sper site per generation) for each chromosome is converted to cM/Mb with the script called Pipeline_OutPyrho_cM_Mb.sh using a command named convert_Output_pyrho_cM_Mb.awk.

For the calculation of recombination rates in non-overlapping windows with different sizes (50 kb ,100 kb,200 kb and 1Mb), a python script under the name Rec_windows.py was used. This script calculates the average recombination rates weighted by the physical distance betwen each pair of sites where recombination was estimated. The script is implemented in PyrhOut_wtAvg_rates_windows.sh.

iii) Genome annotation: 

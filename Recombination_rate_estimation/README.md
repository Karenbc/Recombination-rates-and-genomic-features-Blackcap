**Recombination rate estimation**

The recombination rates were inferred using pyrho (Kamm et al., 2016; Spence & Song, 2019).

The script called RecRates_inference_I_Pyrho.sh contains the computation of the lookup-table and the simulation step, named hyperparam, to find optimal parameters (window size and block penalty) for the estimation. The inference of recombination rates from real data (VCF files) is included in the script named RecRates_inference_II_pyrho.pl. This script generates several parallel jobs corresponding to each chromosome where the recombination rate is inferred.

The demographic data used for the recombination rate estimation, the computed lookup-table, and hyperparam output for the blackcap are located in the folder called "Dataset". 


**Recombination rates conversion**

The output of pyrho (r, recombination rates per site per generation) for each chromosome is converted to cM/Mb with the script called Conversion_OutPyrho_to_cM_Mb_extended.sh using a script named convert_Output_pyrho_cM_Mb.awk.

**Recombination rates calculated in windows**

For the calculation of recombination rates in non-overlapping windows with different sizes (50 kb ,100 kb, 200 kb, and 1 Mb), a python script under the name Rec_windows.py was used. This script calculates the average recombination rates weighted by the physical distance between each pair of sites where recombination was estimated. The script is implemented in PyrhOut_wtAvg_rates_windows.sh.

The recombination rates calculated in windows can be found in the folder "Dataset".

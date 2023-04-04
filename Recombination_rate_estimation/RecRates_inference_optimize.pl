#!/usr/bin/perl -w
use strict;

# author: Karen Bascon  Cardozo
# date: January 2020
# usage: Estimate recombination rates from VCFs for each chromosome with Optimize tool in Pyrho

my $qsub_template = do { local $/; <DATA> };
my @files = qw/ 1 2 3 Z 4 5 6 7 8 9 10 11 12 W 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 w_unkown 33 /;

for my $file ( @files ) {

    (my $sbatch_cmd = $qsub_template) =~ s/!file!/$file/g;
    my $sbatch_out = "$file.Pyrho_Optimize_NewRef_chr_pen20w50.sh";
    open my $OUT, ">", $sbatch_out or die "Could not open file $sbatch_out: $!";

    print $OUT $sbatch_cmd;

    close $OUT;

    system "sbatch < $sbatch_out";
    system "mv *Pyrho_Optimize_NewRef_chr_pen20w50.sh /home/bascon.c/scripts/Pyrho/Perl_Optimize/NEW_Ref/SH_Optimize_chr_SS";
} 


__DATA__
#!/bin/bash
#
#SBATCH --job-name=chr!file!CR0_Pyrho_Optimize_RecRate
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=72:05:00
#SBATCH --mem=60G
#SBATCH --error=/home/bascon.c/stdout/PyRho/Optimize_RecRates/NEW_Reference/Continental_Residents/Cont_Res0/chr_!file!_ContRes0_n38_pen100_w100.%J.err
#SBATCH --output=/home/bascon.c/stdout/PyRho/Optimize_RecRates/NEW_Reference/Continental_Residents/Cont_Res0/chr_!file!_ContRes0_n38_pen100_w100.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bascon.c@evolbio.mpg.de
#SBATCH --partition=standard






#mytools

 module load python/3.5.3

#my folders

OUTPUT_table=/home/bascon.c/PYRHO/Lookup_Tables
VCF=/home/bascon.c/VCF/NEWref_VCF/Population_Subsets/Continental_Residents/Cont_Res0/BioFiltered_VCF
Out_RecRate=/home/bascon.c/PYRHO/Optimize_RecRates/NEW_Ref/Continent_Residents/Cont_Res0

#Code:

source /home/bascon.c/tutorial-env/bin/activate

#export NUMEXPR_MAX_THREADS=10 

##3rd step in pyrho: Inference of recombination-rates with optimize



bgzip $VCF/Chr_scaffolds/chr_!file!_ContRes0_HFilt_19samples_GQ20_DPmin5_DPmax60_gatkrmfilt_maxmiss07_Maf3_biallelic_1.vcf

pyrho optimize --vcffile $VCF/Chr_scaffolds/chr_!file!_ContRes0_HFilt_19samples_GQ20_DPmin5_DPmax60_gatkrmfilt_maxmiss07_Maf3_biallelic_1.vcf.gz --windowsize 50 --blockpenalty 20 --tablefile $OUTPUT_table/Continental_Residents_n38_N48_lookuptable.hdf --ploidy 2 --outfile $Out_RecRate/chr_!file!_ContRes0_n38_Pen20_W50.rmap --logfile .

echo "chr!file! Continental_Resident done"


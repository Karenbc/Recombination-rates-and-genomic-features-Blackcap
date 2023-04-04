#!/usr/bin/perl -w
use strict;

# usage: Estimate recombination rates from VCF files for each chromosome with optimize tool in Pyrho.

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
############################

#Tools

 module load python/3.5.3

#Folders

#Lookup table 
OUTPUT_table=/home/Lookup_Tables
#SNPs in VCF file
VCF=/home/NEWref_VCF/Population_Subsets/Continental_Residents/Cont_Res0/BioFiltered_VCF

#path for output
Out_RecRate=/home/bascon.c/PYRHO/Optimize_RecRates/NEW_Ref/Continent_Residents/Cont_Res0

#Code:

source /home/bascon.c/tutorial-env/bin/activate



##3rd step in pyrho: Inference of recombination-rates with optimize



bgzip $VCF/Chr_scaffolds/chr_!file!_ContRes0_HFilt_19samples_GQ20_DPmin5_DPmax60_gatkrmfilt_maxmiss07_Maf3_biallelic_1.vcf

pyrho optimize --vcffile $VCF/Chr_scaffolds/chr_!file!_ContRes0_HFilt_19samples_GQ20_DPmin5_DPmax60_gatkrmfilt_maxmiss07_Maf3_biallelic_1.vcf.gz --windowsize 50 --blockpenalty 20 --tablefile $OUTPUT_table/Continental_Residents_n38_N48_lookuptable.hdf --ploidy 2 --outfile $Out_RecRate/chr_!file!_ContRes0_n38_Pen20_W50.rmap --logfile .

echo "chr!file!  done"


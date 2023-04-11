#!/usr/bin/perl -w
use strict;

# author: Andrea Bours 
# date, last modification: Feb 2021
# usage: perl multiple.pl
# submits sbatch file for each item in list
# adjust project name, paths and directories as needed in the script
# ensure you have the sbatch folder created
# modified by Andrea based on GATK best practices 2019
# the first part of the script is to make multiple scripts in one go from a template, in which different variables can be provided dependent on desired output files. The bash script is written and placed into a scripts folder specified in line 31
# variables are specified in line 19, contains per line from left to right (separated by tabs): the parts of the genome (1st...10th), on which the script should operate.
# based on the variables one needs to make interval lists, which specifies the chromosome/scaffolds the variable represents.
# additionally a list of all samples should be made, see READ.me. 
# line 19 specifies by what the variables above are replaced
# above variables can be adjusted according to the needs of the data 
my $qsub_template = do { local $/; <DATA> };

my @files = qw/ 1st 2nd 3rd 4th 5th 6th 7th 8th 9th 10th /;

for my $file ( @files ) {
    (my $sbatch_cmd = $qsub_template) =~ s/!file!/$file/g;
    my $sbatch_out = "$file.genotypegvcfs_all_sites.sh";
    open my $OUT, ">", $sbatch_out or die "Could not open file $sbatch_out: $!";

    print $OUT $sbatch_cmd;

    close $OUT;

    system "sbatch < $sbatch_out";
    system "mv *genotypegvcfs_all_sites.sh <path_to_dir>/scripts";
} 

# Here the actual bash script starts as will be made by the above perl command
# Here the genome is split into 10 chunks (equal), as provided from the previous script: 3.combinegvcfs_script.pl 
# memory requirement is specified at 50G, be aware that for back ground processes the actual memory usage specified in the commands is lower. NOTE: REQUESTING MORE MEMORY WILL NOT SPEED UP THE PROGRAM!
# runtime requirement is 100 hours, this should provide enough time for 1/10th of the genome to be genotyped (buffer already added).
# checkpoints are there to help in restarting the script at controlled points. NOTE: occasionally before restarting the last file may need to be removed, to ensure a new file is created. 
# checkpoints control if a file is made and contains data, and often also if anything went wrong while running the program (found by exception in log file).
# adjust to your own email for slurm to report to and adjust the paths for standard out and error to create files in.
# To aid in computation the beegfs directory, ensure you make a folder in this directory! This should have been done already based on the previous script.
# this is computationally light, time is more limiting factor than memory.
# standard partition is sufficient for this part of the pipeline.


__DATA__
#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=genotype_all_sites_!file!
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, important if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime
#SBATCH --time=100:00:00
#  maximum requested memory
#SBATCH --mem=50G
#  write std out and std error to these files
#SBATCH --error=<path_to_dir>/stdout/genotype_all_sites_!file!.%J.err
#SBATCH --output=<path_to_dir>/stdout/genotype_all_sites_!file!.%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<email>
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

#add your code here:

# first some housekeeping to run the programs without any problems, additionally variables are set to clean up the code. 
# NOTE: always load the modules before hand to ensure that they are in the partition of the cluster that you're using. If you need to switch between modules up-date the script accordingly.

# load necessary modules
# module needed for: GATK v4.1.7.0; picard v2.21.9 
module load java/x64/8u121

#create variables to ease writing and reading code
# variables to easy use programs, we use the most recent version keep in mind that there may have been some important changes and up-date programs and script accordingly
picard=/data/biosoftware/Picard/v2.21.9/picard.jar
gatk4=/data/biosoftware/GATK/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar


# personal adjustment, this is not necessary
name=!file!

# variables to directories used in this pipeline, make these folders in advance. Here I go deeper into subsequent directories that will contain files of the next step. more directories can be added, but here I choose to compartmentalize based on major products and subsequent analyses on them.
dir=~/variant_calling
gvcf=${dir}/gvcf
vcf=${gvcf}/vcf
pipeline=${dir}/pipeline
beegfs=/mnt/beegfs/<path_to_dir>

# the log directory is made to be able to trouble shoot easier as std.err can be rerouted here and is then searchable by the check-points in this script which are written in stdout, note make a new log dir or empty old one before rerunning pipeline as similar named files will be over written in the current set up.
log=${dir}/log

# origin of the reference, note this is the new renamed_reordered blackcap reference in my personal folder, change the reference accordingly. It's inmportant to prepare the reference for use, with bwa index, Picard CreateSequenceDictionary and Samtools faidx (these files should be in same folder as reference)
ref=<path_to_dir>/reference.fasta

#######
# start of the actual process

# Step 1:
# genotypegvcfs is best run from the beegfs, this will speed up the process therefore the files were created in this directory, with the previous script. the output will be directly placed into the home directory. this run included non-variant sites, via: --include-non-variant-sites option

date
echo 'Step 1: genotypegvcfs all sites, running in the beegfs, output in personal home directory'

# first go to beegfs directory
cd ${beegfs}


java -Xms40G -Xmx40G -jar ${gatk4} GenotypeGVCFs -R ${ref} -V ${name}.g.vcf.gz -O ${vcf}/${name}_all_sites.vcf.gz --include-non-variant-sites 2> ${log}/${name}_genotype_all_sites.log

#change to the directory in which the output of above command is placed.
cd ${vcf}

# Checkpoint 1
date
echo 'Checkpoint 1'

if grep -i exception ${log}/${name}_genotype_all_sites.log ; then 
    echo "failed"
    exit 1
else echo "no problem occurred while running genotypegvcfs"
fi

if test -s ${name}_all_sites.vcf.gz ; then
	echo "combined all sites vcf exists and contains data"
else echo "failed"
	exit 1
fi

# Step 2:
# validate the compiled files
date
echo 'Step 2: validating the combined all sites vcf'

# this version of GATK has an issue in which the vcf contains sites that are considered erroneous (alternative allele not found -> will be fixed by filtering the files downstream), therefore validatevariants is run with --warn-on-errors

java -Xms40G -Xmx40G -jar ${gatk4} ValidateVariants -R ${ref} -V ${name}_all_sites.vcf.gz -L ${pipeline}/${name}_intervals.list --warn-on-errors 2> ${log}/${name}_genotype_all_sites_valivar.log

# Checkpoint: 2
date
echo 'Checkpoint 2'

if grep -i exception ${log}/${name}_genotype_all_sites_valivar.log ; then 
    echo "failed"
    exit 1
else echo "the variants are valid"
fi


echo 'all sites vcf has been made and validatevariants has run, check manually the warnings'

# In general these sites to which the warnings refer will be filtered out during subsequent filtering of the sites.

# REMINDER TO COPY THE GVCFS FROM THE BEEGFS TO HOME DIRECTORY!

# next step is to gather the vcfs into one, see script 5.gathervcf_script.sh or 5a.gathervcf_all_sites_script.sh

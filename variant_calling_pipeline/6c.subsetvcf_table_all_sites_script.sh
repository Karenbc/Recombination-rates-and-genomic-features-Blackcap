#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=subset_all_sites
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, important if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime
#SBATCH --time=72:00:00
#  maximum requested memory
#SBATCH --mem=100G
#  write std out and std error to these files
#SBATCH --error=<path_to_dir>/stdout/subset_all_sites.%J.err
#SBATCH --output=<path_to_dir>/stdout/subset_all_sites.%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<email>
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

# author: Andrea Bours 
# date, last modification: June 2021

# memory requirement is specified at 100G, be aware that for back ground processes the actual memory usage specified in the commands is lower. 
# runtime requirement is 72 hours, this should provide enough time for removing the INDELS and making the table
# checkpoints are there to help in restarting the script at controlled points. NOTE: occasionally before restarting the last file may need to be removed, to ensure a new file is created. 
# checkpoints control if a file is made and contains data, and often also if anything went wrong while running the program (found by exception in log file).
# adjust to your own email for slurm to report to and adjust the paths for standard out and error to create files in.
# standard partition is sufficient for this part of the pipeline.


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

# variables to directories used in this pipeline, make these folders in advance (I'll create these folders again on my rerun to reduce my current clutter in my first round folders)
dir=~/variant_calling
pipeline=${dir}/pipeline
gvcf=${dir}/gvcf
vcf=${gvcf}/vcf

# the log directory is made to be able to trouble shoot easier as std.err can be rerouted here and is then searchable by the check-points in this script, note to self make a new log dir or empty old one before rerunning pipeline 
log=${dir}/log

# origin of the reference, note this is the new blackcap reference in my personal folder, change the reference accordingly. It's important to prepare the reference for mapping, with bwa index, Picard CreateSequenceDictionary and Samtools faidx (these files should be in same folder as reference)
ref=<path_to_dir>/reference.fasta

#######
# Start of the actual process

# Step 1:
# subset the vcf made in the previous script by SNP

date
echo 'Step 1: subsetting final_all_sites_vcf by removing INDELs from vcf'

# change into the correct directory
cd ${vcf}

java -Xms64G -Xmx64G -jar ${gatk4} SelectVariants -R ${ref} -V final_all_sites.vcf.gz --select-type-to-exclude INDEL -O final_all_sites_no_indel.vcf.gz 2> ${log}/final_all_sites_subset_no_indel.log

# Checkpoint: 1
date
echo 'Checkpoint 1'

if grep -i exception ${log}/final_all_sites_subset_no_indel.log ; then 
    echo "failed"
    exit 1
else echo "no problem occurred when removing INDELs"
fi

if test -s final_all_sites_no_indel.vcf.gz ; then
	echo "the final all_sites_no_indel vcf exists and contains data"
else echo "failed"
	exit 1
fi

# Step 2:
# check the subsetted INDEL vcf with validatevariants
date
echo 'Step 2: validating the final all_sites_no_indel vcf'

# this version of GATK has an issue in which the vcf contains sites that are considered erroneous (alternative allele not found -> will be fixed by filtering the files downstream), therefore validatevariants is run with --warn-on-errors
java -Xms64G -Xmx64G -jar ${gatk4} ValidateVariants -R ${ref} -V final_all_sites_no_indel.vcf.gz -L ${pipeline}/whole_genome_intervals.list --warn-on-errors 2> ${log}/final_all_sites_no_indel_valivar.log

# Checkpoint: 2
date
echo 'Checkpoint 2'

if grep -i exception ${log}/final_all_sites_no_indel_valivar.log ; then 
    echo "failed"
    exit 1
else echo "the variants are valid"
fi

# Step 3:
# get table of the info fields for each indel from the indel vcf to hard filter it in R
date
echo "Step 3: making table to check how to filter in R"

java -Xms64G -Xmx64G -jar ${gatk4} VariantsToTable -R ${ref} -V final_all_sites_no_indel.vcf.gz -O final_all_sites.table -F ID -F CHROM -F POS -F FILTER -F QD -F QUAL -F DP -F AD -F FS -F MQ -F MQRankSum -F ReadPosRankSum -F SOR -F NO-CALL -F TRANSITION -F HET -F NCALLED -F MULTI-ALLELIC -raw

# Checkpoint: 3
date
echo 'Checkpoint 3'

if test -s final_all_sites.table ; then
	echo "the final_all_sites table exists and contains data"
else echo "failed"
	exit 1
fi


echo "subsetted the vcf by removing INDELs and created a table for the next step of checking hardfilters in R"

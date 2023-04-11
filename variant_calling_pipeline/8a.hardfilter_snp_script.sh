#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=hardfilter_subset_snp
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, important if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime
#SBATCH --time=100:00:00
#  maximum requested memory
#SBATCH --mem=50G
#  write std out and std error to these files
#SBATCH --error=<path_to_dir>/stdout/hardfilter_subset_snp.%J.err
#SBATCH --output=<path_to_dir>/stdout/hardfilter_subset_snp.%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<email>
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

# author: Andrea Bours 
# date, last modification: June 2021

# memory requirement is specified at 50G, be aware that for back ground processes the actual memory usage specified in the commands is lower. 
# runtime requirement is 100 hours, this should provide enough time for hardfiltering based on the parameters found in previous step.
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

# origin of the reference, note this is the old blackcap reference in my personal folder, change the reference accordingly. It's important to prepare the reference for mapping, with bwa index, Picard CreateSequenceDictionary and Samtools faidx (these files should be in same folder as reference)
ref=<path_to_dir>/reference.fasta

#######
# Start of the actual process

# Step 1:
# reset the filters according to the decisions made through the, snp specific, report made by 7.hard_filter.Rmd

date
echo 'Step 1: reset filters based on report'
# NOTE: gatk warns to provide the filters as a single attribute with it's own criteria and NOT as a string like done here

java -Xms42G -Xmx42G -jar ${gatk4} VariantFiltration -R ${ref} -V final_snp.vcf.gz -O hard_filter_snp.vcf.gz --filter-name "Low_QD" --filter-expression "QD < 2.5" --filter-name "high_FS" --filter-expression "FS > 45.0" --filter-name "High_SOR" --filter-expression "SOR > 3.0" --filter-name "Low_MQ" --filter-expression "MQ < 40.0" --filter-name "Low_MQRankSum" --filter-expression "MQRankSum < -12.5" --filter-name "Low_ReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0" 2> ${log}/final_snp_hard_filter.log

# Checkpoint 1
date
echo 'Checkpoint 1'

if grep -i exception ${log}/final_snp_hard_filter.log ; then 
    echo "failed"
    exit 1
else echo "no problem occurred when resetting the filters"
fi

if test -s hard_filter_snp.vcf.gz ; then
	echo "the hard filter snp vcf exists and contains data"
else echo "failed"
	exit 1
fi

# Step 2:
# check the reset vcf with validatevariants
date 
echo 'Step 2: validating the hard_filter_snp vcf'

java -Xms42G -Xmx42G -jar ${gatk4} ValidateVariants -R ${ref} -V hard_filter_snp.vcf.gz -L ${pipeline}/whole_genome_intervals.list --warn-on-errors 2> ${log}/hard_filter_snp_valivar.log

# Checkpoint: 2
date 
echo 'Checkpoint 2'

if grep -i exception ${log}/hard_filter_snp_valivar.log ; then 
    echo "failed"
    exit 1
else echo "the variants are valid"
fi

# Step 3:
# parse the vcf with the filters set into select variants, to remove all sites that do not meet the hard filters.
date 
echo 'Step 3: selecting the SNPs that passed the hardfilters'

java -Xms42G -Xmx42G -jar ${gatk4} SelectVariants -R ${ref} -V hard_filter_snp.vcf.gz --exclude-filtered -O snp.vcf.gz 2> ${log}/snp_remove_hard_filter.log

# Checkpoint 3
date 
echo 'Checkpoint 3'

if grep -i exception ${log}/snp_remove_hard_filter.log ; then 
    echo "failed"
    exit 1
else echo "no problem occured when removing variants"
fi

if test -s snp.vcf.gz ; then
    echo "the snp vcf exists and contains data"
else echo "failed"
    exit 1
fi

# Step 4:
# check the hardfiltered vcf with validatevariants
date
echo 'Step 4: validating the snp vcf'

java -Xms42G -Xmx42G -jar ${gatk4} ValidateVariants -R ${ref} -V snp.vcf.gz -L ${pipeline}/whole_genome_intervals.list --warn-on-errors 2> ${log}/snp_valivar.log

# Checkpoint 4
date
echo 'Checkpoint 4'

if grep -i exception ${log}/snp_valivar.log ; then 
    echo "failed"
    exit 1
else echo "the variants are valid"
fi

echo 'the SNP vcf is hard-filtered and ready for analyses'

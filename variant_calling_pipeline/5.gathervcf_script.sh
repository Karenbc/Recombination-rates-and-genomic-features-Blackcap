#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=gathervcfs
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, important if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=24:00:00
#  maximum requested memory
#SBATCH --mem=50G
#  write std out and std error to these files
#SBATCH --error=<path_to_dir>/stdout/gathervcfs.%J.err
#SBATCH --output=<path_to_dir>/stdout/gathervcfs.%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<email>
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

# author: Andrea Bours 
# date, last modification: Feb 2021

# memory requirement is specified at 50G, be aware that for back ground processes the actual memory usage specified in the commands is lower. 
# runtime requirement is 24 hours, this should provide enough time for gathering the vcfs
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


# variables to directories used in this pipeline, make these folders in advance. Here I go deeper into subsequent directories that will contain files of the next step. more directories can be added, but here I choose to compartmentalize based on major products and subsequent analyses on them.
dir=~/variant_calling
gvcf=${dir}/gvcf
vcf=${gvcf}/vcf
pipeline=${dir}/pipeline

# the log directory is made to be able to trouble shoot easier as std.err can be rerouted here and is then searchable by the check-points in this script which are written in stdout, note make a new log dir or empty old one before rerunning pipeline as similar named files will be over written in the current set up.
log=${dir}/log

# origin of the reference, note this is the new renamed_reordered blackcap reference in my personal folder, change the reference accordingly. It's important to prepare the reference for use, with bwa index, Picard CreateSequenceDictionary and Samtools faidx (these files should be in same folder as reference)
ref=<path_to_dir>/reference.fasta

#######
# Start of the actual process

# Step 1:
# gather the vcfs made previously into one big final vcf.

date
echo 'Step 1: gathering the vcfs into one final vcf'

# first change to the vcf directory
cd ${vcf}

java -Xms42G -Xmx42G -jar ${picard} GatherVcfs I=1st.vcf.gz I=2nd.vcf.gz I=3rd.vcf.gz I=4th.vcf.gz I=5th.vcf.gz I=6th.vcf.gz I=7th.vcf.gz I=8th.vcf.gz I=9th.vcf.gz I=10th.vcf.gz O=./final.vcf.gz R=${ref} 2> ${log}/final_vcf_gather.log

# Checkpoint: 1
date 
echo 'Checkpoint 1'

if grep -i exception ${log}/final_vcf_gather.log ; then 
    echo "failed"
    exit 1
else echo "no problem while running gathervcfs"
fi

if test -s final.vcf.gz ; then
	echo "the final vcf exists and contains data"
else echo "failed"
	exit 1
fi

# Step 2:
# as gathervcfs is a picard program we need to index the final vcf in order to be able to continue running gatk programs.
date
echo 'Step 2: indexing the final vcf'

java -Xms42G -Xmx42G -jar ${gatk4} IndexFeatureFile -I final.vcf.gz 2> ${log}/final_vcf_index.log

# Checkpoint: 2
date
echo 'Checkpoint 2'

if grep -i exception ${log}/final_vcf_index.log ; then 
    echo "failed"
    exit 1
else echo "no problem creating the index of the vcf"
fi

if test -s final.vcf.gz.tbi ; then
	echo "the final vcf index exists and contains data"
else echo "failed"
	exit 1
fi

# Step 3:
# validate the gathered file
date 
echo 'Step 3: validating the final vcf'

# this version of GATK has an issue in which the vcf contains sites that are considered erroneous (alternative allele not found -> will be fixed by filtering the files downstream), therefore validatevariants is run with --warn-on-errors
# NOTE MAKE A INTERVAL LIST OF THE COMPLETE GENOME 

java -Xms42G -Xmx42G -jar ${gatk4} ValidateVariants -R ${ref} -V final.vcf.gz -L ${pipeline}/whole_genome_intervals.list --warn-on-errors 2> ${log}/final_vcf_valivar.log

# Checkpoint: 3
date
echo 'Checkpoint 3'

if grep -i exception ${log}/final_vcf_valivar.log ; then 
    echo "failed"
    exit 1
else echo "the variants are valid"
fi

echo 'final vcf is complete, next steps is subsetting and filtering'

# Next steps are subsetting of the completed vcf into separate SNP and INDEL vcfs (and create tables for "hardfiltering" of variants), see script 6a.subsetvcf_table_snp_script.sh and 6b.subsetvcf_table_indel_script.sh

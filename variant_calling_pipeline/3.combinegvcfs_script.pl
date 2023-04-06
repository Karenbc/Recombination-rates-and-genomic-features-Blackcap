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
    my $sbatch_out = "$file.combinegvcfs.sh";
    open my $OUT, ">", $sbatch_out or die "Could not open file $sbatch_out: $!";

    print $OUT $sbatch_cmd;

    close $OUT;

    system "sbatch < $sbatch_out";
    system "mv *combinegvcfs.sh /home/bours/scripts";
} 

# Here the actual bash script starts as will be made by the above perl command
# Here the genome is split into 10 chunks (approx. equal) for the program to run on, this speeds up the process. However depending on the amount of samples or sizes of the g.vcfs this may not be needed.
# NOTE: always subset the reference genome, never the samples (this pipeline gets it's strength from using all samples available).
# memory requirement is specified at 420G (yes this much is needed at least), be aware that for back ground processes the actual memory usage specified in the commands is lower.
# runtime requirement is based on trial and error.
# checkpoints are there to help in restarting the script at controlled points. NOTE: occasionally before restarting the last file may need to be removed, to ensure a new file is created. 
# checkpoints control if a file is made and contains data, and often also if anything went wrong while running the program (found by exception in log file).
# adjust to your own email for slurm to report to and adjust the paths for standard out and error to create files in.
# To aid in computation the beegfs directory is used, therefore ensure that you make a folder in this directory! the beegfs is used as the program will be opening and closing a lot of files.
# This program is computationally the heaviest, therefore the highmem partition should be used!
# in doubt confirm with IT use of beegfs and highmem

__DATA__
#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=!file!_combinegvcfs
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, important if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime
#SBATCH --time=120:00:00
#  maximum requested memory
#SBATCH --mem=420G
#  write std out and std error to these files
#SBATCH --error=/home/bours/stdout/!file!_combinegvcfs.%J.err
#SBATCH --output=/home/bours/stdout/!file!_combinegvcfs.%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bours@evolbio.mpg.de
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=highmem

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

# variables to directories used in this script, make these folders in advance. Here I go deeper into subsequent directories that will contain files of the next step. more directories can be added, but here I choose to compartmentalize based on major products and subsequent analyses on them.
dir=~/variant_calling
gvcf=${dir}/gvcf
pipeline=${dir}/pipeline
beegfs=/mnt/beegfs/bours

# the log directory is made to be able to trouble shoot easier as std.err can be rerouted here and is then searchable by the check-points in this script which are written in stdout, note make a new log dir or empty old one before rerunning pipeline as similar named files will be over written in the current set up.
log=${dir}/log

# origin of the reference, note this is the new renamed_reordered blackcap reference in my personal folder, change the reference accordingly. It's important to prepare the reference for use, with bwa index, Picard CreateSequenceDictionary and Samtools faidx (these files should be in same folder as reference)
ref=/home/bours/reference/renamed_reorder_new_reference.fasta

#######
# start of the actual process

# Step 1:
# this step is best run from the beegfs folder! NOTE: BEEGFS IS NOT BACKED UP, THEREFORE COPY THE FILES AFTER COMPLETION TO YOUR HOME DIRECTORY, BUT ONLY AFTER THE NEXT SCRIPTS 4/4a HAVE BEEN RUN.
# combining the g.vcfs per partition of the reference genome.
date
echo 'Step 1: combinegvcfs, output will be in personal beegfs folder'

# first go to beegfs directory
cd ${beegfs}

java -Xms375G -Xmx375G -jar ${gatk4} CombineGVCFs -R ${ref} --variant ${pipeline}/sample.list -O ${name}.g.vcf.gz -L ${pipeline}/${name}_intervals.list 2> ${log}/${name}_combine.log

# checkpoint: 1
date
echo 'Checkpoint 1'

if grep -i exception ${log}/${name}_combine.log ; then 
    echo "failed"
    exit 1
else echo "no problem running combinegvcfs"
fi

if test -s ${name}.g.vcf.gz ; then
	echo "combined g.vcf exists and contains data"
else echo "failed"
	exit 1
fi

# Step 2:
# Validate the created g.vcf, important to provide the -gvcf option
date
echo 'Step 2: validating the combined g.vcf'

java -Xms375G -Xmx375G -jar ${gatk4} ValidateVariants -R ${ref} -V ${name}.g.vcf.gz -gvcf -L ${pipeline}/${name}_intervals.list 2> ${log}/${name}_combine_valivar.log

# checkpoint: 2
date
echo 'Checkpoint 2'

if grep -i exception ${log}/${name}_combine_valivar.log ; then 
    echo "failed"
    exit 1
else echo "the variants are valid"
fi

echo 'the g.vcfs have been combined, output in beegfs and variants have been validated'

# REMEMBER THIS IS NOT BACKED UP AFTER THE NEXT SCRIPT 4 COPY THE FILES TO YOUR OWN HOME DIRECTORY!

# next step is to genotype the gvcs, meaning making vcfs of either only callable sites or all sites, see script 4.genotypegvcfs_script.pl or 4a.genotypegvcfs_all_sites_script.pl
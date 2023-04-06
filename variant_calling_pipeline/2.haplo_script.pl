#!/usr/bin/perl -w
use strict;

# author: Andrea Bours 
# date, last modification: Feb 2021
# usage: perl multiple.pl
# submits sbatch file for each item in list
# adjust project name, paths and directories as needed, in both the script as well as the variables provided in the txt
# ensure you have the sbatch folder created
# modified by Andrea based on GATK best practices 2019
# the first part of the script is to make multiple scripts in one go from a template, in which different variables can be provided dependent on input files, which is accessed through haplo_variables.txt. The bash script is written and placed into a scripts folder specified in line 39
# haplo_variables.txt, opened in line 17, contains on per line from left to right (separated by tabs): sample name; runtime estimate in hours. for example line see READ.me.
# line 20 specifies by what the variables above are replaced
# above variables can be adjusted according to the needs of the data 
my $qsub_template = do { local $/; <DATA> };

open IN, "<haplo_variables.txt" or die "Could not open file file.txt, $!";

my @files = (<IN>);
my @variables = qw/ !name! !runtime! /;

for my $line (@files) {
    my @line_array = split ' ', $line;
    my $sbatch_cmd = $qsub_template;

    for (my $i=0; $i < scalar @variables; $i++) {
        my $rex = qr/$variables[$i]/;
        $sbatch_cmd =~ s/${rex}/$line_array[$i]/g;
    }

    my $sbatch_out = "$line_array[0].haplo.sh";
    open my $OUT, ">", $sbatch_out or die "Could not open file $sbatch_out: $!";

    print $OUT $sbatch_cmd;

    close $OUT;

    system "sbatch < $sbatch_out";
    system "mv *haplo.sh /home/bours/scripts";
} 

close(IN); 

# Here the actual bash script starts as will be made by the above perl command
# for the submission specifics, the haplo_variables.txt provides names to distinguish jobs, and an estimated runtime necessary for the script to complete. 
# memory requirement is specified at 40G, be aware that for back ground processes the actual memory usage specified in the commands is lower. NOTE: it's not useful to increase the memory requested as the programs don't require this and this way you don't penalize yourself in the cluster.
# runtime requirement are based on trial and error, however calculations for estimates will be given in READ.me, buffer times are added in these calculations. This step takes time!
# checkpoints are there to help in restarting the script at controlled points. NOTE: occasionally before restarting the last file may need to be removed, to ensure a new file is created. 
# checkpoints control if a file is made and contains data, and often also if anything went wrong while running the program (found by exception in log file).
# adjust to your own email for slurm to report to and adjust the paths for standard out and error to create files in.
# standard partition is sufficient for this part of the pipeline

__DATA__
#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=!name!_haplo
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime
#SBATCH --time=!runtime!:00:00
#  maximum requested memory
#SBATCH --mem=40G
#  write std out and std error to these files
#SBATCH --error=/home/bours/stdout/!name!_haplo.%J.err
#SBATCH --output=/home/bours/stdout/!name!_haplo.%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bours@evolbio.mpg.de
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
name=!name!

# variables to directories used in this script, make these folders in advance(!). Here I go deeper into subsequent directories that will contain files of the next step. more directories can be added, but here I choose to compartmentalize based on major products and subsequent analyses on them.
dir=~/variant_calling
ubam=${dir}/ubam
xtbam=${ubam}/xtbam
finalbam=${xtbam}/finalbam
gvcf=${dir}/gvcf

# with the below variable one let the script use the specific file only, this enables the script the execute, without a separate command to distinguish between samples that have been merged. However it's important to name the merged final accordingly
bam=!name!_final.bam

# the log directory is made to be able to trouble shoot easier as std.err can be rerouted here and is then searchable by the check-points in this script which are written in stdout, note make a new log dir or empty old one before rerunning pipeline as similar named files will be over written in the current set up.
log=${dir}/log

# origin of the reference, note this is the new renamed_reordered blackcap reference in my personal folder, change the reference accordingly. It's important to prepare the reference for use, with bwa index, Picard CreateSequenceDictionary and Samtools faidx (these files should be in same folder as reference)
ref=/home/bours/reference/renamed_reorder_new_reference.fasta

#######
# start of the actual script to call haplotypes from the input.

# Step 1:
# run haplotypecaller to create a g.vcf which contains all bases that were callable by the mapper. 
date
echo 'Step 1: Haplotypecaller, to call haplotypes'

# first change into the directory where the bams are
cd ${finalbam}

java -Xms32G -Xmx32G -jar ${gatk4} HaplotypeCaller -R ${ref} -I ${bam} -ERC GVCF -O ${gvcf}/${name}.g.vcf.gz 2> ${log}/${name}_haplo.log

# Checkpoint: 1
date
echo 'Checkpoint 1'

if grep -i exception ${log}/${name}_haplo.log ; then 
    echo "failed"
    exit 1
else echo "no problem making the g.vcf"
fi

if test -s ${gvcf}/${name}.g.vcf.gz ; then
    echo "g.vcf exists and contains data"
else echo "failed"
    exit 1
fi

# Step 2:
# validate the compiled files

# change into the directory containing the g.vcf 
cd ${gvcf}

# Validate the created g.vcf, important to provide the -gvcf option
date
echo 'Step 2: validating the g.vcf'

java -Xms32G -Xmx32G -jar ${gatk4} ValidateVariants -R ${ref} -V ${name}.g.vcf.gz -gvcf 2> ${log}/${name}_valivar.log

# Checkpoint: 2
date
echo 'Checkpoint 2'

if grep -i exception ${log}/${name}_valivar.log ; then 
    echo "failed"
    exit 1
else echo "the variants are valid"
fi

#note one still has to check the std.out

echo "gzipped g.vcf file created and variants have been validated"

# next step is to combine the g.vcfs into one g.vcf, see script 3.combinegvcfs_script.pl

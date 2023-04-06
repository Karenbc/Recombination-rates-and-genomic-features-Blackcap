#!/usr/bin/perl -w
use strict;

# author: Andrea Bours 
# date, last modification: Feb 2021
# usage: perl multiple.pl
# submits sbatch file for each item in list
# adjust project name, paths and directories as needed, in both the script as well as the variables provided in the txt
# check/change readgroup information according to the samples
# ensure you have the sbatch folder created
# modified by Andrea based on GATK best practices 2019
# the first part of the script is to make multiple scripts in one go from a template, in which different variables can be provided dependent on input files, which is accessed through merge_bam.txt. The bash script is written and placed into a scripts folder specified in line 40
# merge_bam.txt, opened in line 18, contains per line from left to right (seperated by tabs): the bams to be merged from 1..n (be aware that this script does not accommodate to difference in the amount of files to be merged); last entry is the name of the merged file. for example line see READ.me.
# line 21 specifies by what the variables above are replaced
# above variables can be adjusted according to the needs of the data 
my $qsub_template = do { local $/; <DATA> };

open IN, "<merge_bam.txt" or die "Could not open file file.txt, $!";

my @files = (<IN>);
my @variables = qw/ !bam_1! !bam_2! !merged! /;

for my $line (@files) {
    my @line_array = split ' ', $line;
    my $sbatch_cmd = $qsub_template;

    for (my $i=0; $i < scalar @variables; $i++) {
        my $rex = qr/$variables[$i]/;
        $sbatch_cmd =~ s/${rex}/$line_array[$i]/g;
    }

    my $sbatch_out = "$line_array[0].merge_bam.sh";
    open my $OUT, ">", $sbatch_out or die "Could not open file $sbatch_out: $!";

    print $OUT $sbatch_cmd;

    close $OUT;

    system "sbatch < $sbatch_out";
    system "mv *merge_bam.sh /home/bours/scripts";
} 

close(IN); 

# Here the actual bash script starts as will be made by the above perl command
# for the submission specifics, the merge_bam.txt provides names to distinguish jobs. 
# memory requirement is specified at 35G, be aware that for back ground processes the actual memory usage specified in the commands is lower. 30 hours should be enough time to merge the files and run quality control programs 
# checkpoints are there to help in restarting the script at controlled points. NOTE: occasionally before restarting, the last file may need to be removed ensuring a new file is created. 
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
#SBATCH --job-name=!merged!_merge_bam
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, important if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime
#SBATCH --time=30:00:00
#  maximum requested memory
#SBATCH --mem=35G
#  write std out and std error to these files
#SBATCH --error=/home/bours/stdout/!merged!_merge_bam.%J.err
#SBATCH --output=/home/bours/stdout/!merged!_merge_bam.%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bours@evolbio.mpg.de
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

#add your code here:

# first some housekeeping to run the programs, additionally variables are set to clean up the code. 
# NOTE: always load the modules before hand to ensure that they are in the partition of the cluster that you're using. If you need to switch between modules up-date the script accordingly.

# load necessary modules
# module needed for: multiqc v1.7
module load python/3.6.0

# module needed for: fastqc v0.11.8; picard v2.21.9; qualimap v2.2.1
module load java/x64/8u121

# module loads just in case
module load perl/5.24.1
module load R/3.5.3

#create variables to ease writing and reading code
# variables to easy use programs, we use the most recent version keep in mind that there may have been some important changes and up-date programs and script accordingly
picard=/data/biosoftware/Picard/v2.21.9/picard.jar

# personal adjustment, this is not necessary 
merged=!merged!

# names to distinguish the different runs (here 2)
bam_1=!bam_1!
bam_2=!bam_2!

# variables to directories used in this script, make these folders in advance(!). Here I go deeper into subsequent directories that will contain files of the next step. more directories can be added, but here I choose to compartmentalize based on major products and subsequent analyses on them.
dir=~/variant_calling
ubam=${dir}/ubam
xtbam=${ubam}/xtbam
finalbam=${xtbam}/finalbam
metrix=${xtbam}/metrix
qualimap=${finalbam}/qualimap
metrics=${finalbam}/metrics
fast_qc=${ubam}/fast_qc
multiqc_data=${finalbam}/multiqc_data

# to aid writing the script when running final quality checks a variable for the merged final bam file is set
gbam=${merged}_final.bam

# the log directory is made to be able to trouble shoot easier as std.err can be rerouted here and is then searchable by the check-points in this script which are written in stdout, note make a new log dir or empty old one before rerunning pipeline as similar named files will be over written in the current set up.
log=${dir}/log

# origin of the reference, note this is the new renamed_reordered blackcap reference in my personal folder, change the reference accordingly. It's important to prepare the reference for use, with bwa index, Picard CreateSequenceDictionary and Samtools faidx (these files should be in same folder as reference)
ref=/home/bours/reference/renamed_reorder_new_reference.fasta

#######
# start of the actual merging pipeline
# before starting the next part of the pipeline metrics and quality checks have to be evaluated.

# Step 1:
# Merge of the bam files produced by the 1.preproc_map_clean_script.pl that have been sequenced more than once.
date 
echo 'Step 1: merging the seperate bams into one bam'

# change into the last final bam directory.
cd ${finalbam}

# merge the two bam files into one
java -Xms20G -Xmx20G -jar ${picard} MergeSamFiles I=${bam_1}_final.bam I=${bam_2}_final.bam O=${gbam} SO=coordinate CREATE_INDEX=true 2> ${log}/${merged}_merge.log

# Checkpoint: 1
date
echo 'Checkpoint 1'

if grep -i exception ${log}/${merged}_merge.log ; then 
    echo "failed"
    exit 1
else echo "no problem occurred with merging the bams"
fi

if test -s ${gbam} ; then
    echo "merged final bam exists and contains data"
else echo "failed"
    exit 1
fi

echo 'the merged final bam file is created, next steps are quality control of the created file' 

#######
# after production of the merged final bam, quality control of the file needs to be performed. 
# most of these feature can be summarized with multiqc. Here a multiqc report is produced for each sample, but after all samples have run one can make a summary report in multiqc from all samples.

# Step 2:
# Validation of the bam file, here the WARNING for the NM is ignored as this will most likely come up, but isn't important to be corrected for downstream purposes
date 
echo 'Step 2: quality control'
echp 'validation of merged bam'
java -Xms20G -Xmx20G -jar ${picard} ValidateSamFile I=${gbam} O=${metrics}/${merged}_validatesamfile R=${ref} IGNORE=MISSING_TAG_NM

# collect metrics to evaluate the quality of the alignment 
# CollectMultipleMetrics default setting collects: CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, QualityScoreDistribution, MeanQualityByCycle, CollectBaseDistributionByCycle, here there is an addition of the following metric: CollectGcBiasMetrics
date
echo 'collection of metrics of alignment quality'
java -Xms20G -Xmx20G -jar ${picard} CollectMultipleMetrics I=${gbam} O=${metrics}/${merged}_multiple_metrics R=${ref} Program=CollectGcBiasMetrics

java -Xms20G -Xmx20G -jar ${picard} CollectRawWgsMetrics I=${gbam} O=${metrics}/${merged}_raw_wgs_metrics R=${ref} INCLUDE_BQ_HISTOGRAM=true

java -Xms20G -Xmx20G -jar ${picard} CollectWgsMetrics I=${gbam} O=${metrics}/${merged}_wgs_metrics R=${ref} INCLUDE_BQ_HISTOGRAM=true

# the program versions used for Picard and multiqc have an incompatibility, requires a small adjustment to the output file of CollectWgsMetrics, as multiqc will not recognise the output.

sed -i -e 's/picard.analysis.WgsMetrics/picard.analysis.CollectWgsMetrics$WgsMetrics/' ${metrics}/${merged}_wgs_metrics 

# qualimap bamqc, the quality control of the bam
date
echo 'qualimap quality control merged bam'

# change in to the folder in which qualimap results should be deposited
cd ${qualimap}

# running qualimap with bamqc function
qualimap --java-mem-size=20G bamqc -bam ${finalbam}/${gbam} -c -outdir ${merged}

#Step 8:
# merge everything produced under step 7 into on per sample multiqc report, this is for only the merged files. Any quality control files of the separate runs are redundant (no fastqc reports as well as the ignore setting)
date
echo 'Step 8: compile multiqc'

# change back to the last final bam directory.
cd ${finalbam}

multiqc ${metrics}/${merged}* ${qualimap}/${merged} --ignore ${metrics}/${bam_1}* --ignore ${metrics}/${bam_2}* --ignore ${qualimap}/${bam_1} --ignore ${qualimap}/${bam_2} -i ${merged} -f -o ${multiqc_data}

date
echo 'merged bam is produced and quality control has run, evaluate the metrics'


# NOTE: evaluate the metrics

# When all samples have run go to next step haplotypecaller, see script 2.haplo_script.pl

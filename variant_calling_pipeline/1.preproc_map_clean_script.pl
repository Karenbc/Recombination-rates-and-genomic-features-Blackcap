#!/usr/bin/perl -w
use strict;

# author: Andrea Bours 
# date, last modification: Feb 2021
# usage: perl multiple.pl
# submits sbatch file for each item in list
# adjust project name, paths and directories as needed, in both the script as well as the variables provided in the txt read by the perl part of the script
# check/change readgroup information according to the samples
# ensure you have the sbatch folder created, to which the scripts are written
# modified by Andrea based on GATK best practices 2019
# the first part of the script is to make multiple scripts in one go from a template, in which different variables can be provided dependent on input files, which is accessed through preproc_map_clean_variables.txt. The bash script is written and placed into a scripts folder specified in line 40
# preproc_map_clean_variables.txt, opened in line 18, contains per line from left to right (separated by tabs): the base of resequencing files; the name of the sample; the path to R1; the path to R2; library; readgroup; parameter for optical duplicate pixel distance; runtime estimate in hours. for example line see README.md.
# line 21 specifies by what the variables above are replaced
# above variables can be adjusted according to the needs of the data
my $qsub_template = do { local $/; <DATA> };

open IN, "<preproc_map_clean_variables.txt" or die "Could not open file file.txt, $!";

my @files = (<IN>);
my @variables = qw/ !reseq_name! !sample_name! !read_R1! !read_R2! !lib! !rgroup! !optic! !runtime! /;

for my $line (@files) {
    my @line_array = split ' ', $line;
    my $sbatch_cmd = $qsub_template;

    for (my $i=0; $i < scalar @variables; $i++) {
        my $rex = qr/$variables[$i]/;
        $sbatch_cmd =~ s/${rex}/$line_array[$i]/g;
    }

    my $sbatch_out = "$line_array[0].preproc_map_clean.sh";
    open my $OUT, ">", $sbatch_out or die "Could not open file $sbatch_out: $!";

    print $OUT $sbatch_cmd;

    close $OUT;

    system "sbatch < $sbatch_out";
    system "mv *preproc_map_clean.sh <path_to_dir>/scripts";
} 

close(IN); 

# Here the actual bash script starts as will be made by the above perl command
# for the submission specifics, the preproc_map_clean_variables.txt provides names to distinguish jobs, and an estimated runtime necessary for the script to complete. 
# the scripts specifies the use of 6 tasks, which are only used for the mapping step, thus when in resubmitting the script without mapping performed, adjust the amount of tasks to 1.
# memory requirement is specified at 35G, be aware that for background processes the actual memory usage specified in the commands is lower. NOTE: it's not useful to increase the memory requested as the programs don't require this, of course depends on cluster set-up.
# runtime requirement are based on trial and error, however calculations for estimates will be given in README.md, buffer times are already added in these calculations. 
# checkpoints are there to help in restarting the script at controlled points. NOTE: occasionally before restarting, the last file may need to be removed ensuring a new file is created by the script. 
# checkpoints control if a file is made and contains data, and often also if anything went wrong while running the program (found by exception in log file).
# adjust to your own email for Slurm to report to and adjust the paths for standard out and error to create files in.
# standard partition is sufficient for this part of the pipeline

__DATA__
#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=!reseq_name!_preproc_map_clean
#  how many cpus are requested
#SBATCH --ntasks=6
#  run on one node, important if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime
#SBATCH --time=!runtime!:00:00
#  maximum requested memory
#SBATCH --mem=35G
#  write std out and std error to these files
#SBATCH --error=<path_to_stdout_directory>/!reseq_name!_preproc_map_clean.%J.err
#SBATCH --output=<path_to_stdout_directory>/!reseq_name!_preproc_map_clean.%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<email>
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

#add your code here:

# first some housekeeping to run the programs, additionally variables are set to clean up the code. 
# NOTE: always load the modules beforehand to ensure that they are in the partition of the cluster that you're using. If you need to switch between modules up-date the script accordingly.

# load necessary modules
# module needed for: multiqc v1.7
module load python/3.6.0

# module needed for: fastqc v0.11.8; picard v2.21.9; qualimap v2.2.1
module load java/x64/8u121

# module loads just in case
module load perl/5.24.1
module load R/3.5.3

# create variables to ease writing and reading code:
# variables to ease use of programs, we use the most recent version keep in mind that there may have been some important changes and up-date programs and scripts accordingly
picard=/data/biosoftware/Picard/v2.21.9/picard.jar

# personal adjustment, this is not necessary
name=!reseq_name!
optic=!optic!

# this variable directly feeds the path of the forward and reverse reads
read_R1=!read_R1!
read_R2=!read_R2!

# variables to directories used in this script, make these folders in advance(!). Here I go deeper into subsequent directories that will contain files of the next step. More directories can be added, but here I choose to compartmentalize based on major products and subsequent analyses on them.
dir=~/variant_calling
ubam=${dir}/ubam
xtbam=${ubam}/xtbam
finalbam=${xtbam}/finalbam
metrix=${xtbam}/metrix
qualimap=${finalbam}/qualimap
metrics=${finalbam}/metrics
fast_qc=${ubam}/fast_qc
multiqc_data=${finalbam}/multiqc_data

# to aid writing the script, for running final quality checks a variable for the final bam file is set
gbam=${name}_final.bam

# the log directory is made to be able to trouble shoot easier as std.err can be rerouted here and is then searchable by the check-points in this script which are written in stdout, note make a new log dir or empty old one before rerunning pipeline as similar named files will be over written in the current set up.
log=${dir}/log

# information for readgroup, this is necessary for merging samples that have multiple runs/lanes (have been resequenced multiple times or on multiple lanes for the same run), lib and rgroup can be found in the blackcap key excel sheet (if needed rgroup are in the raw reads), adjustments may need to be made dependent on the sequencing platform used.
lib=!lib!
rgroup=!rgroup!

# this is a way to ease the use of the perl pipeline, it also eases the merge of multiple runs of a sample, later in the pipeline, as they need to have an identical "sample_name", otherwise the merge fails. It's no problem to have different readgroups in a merged file.
sample_name=!sample_name!

# origin of the reference, note this is the new renamed_reordered blackcap reference in my personal folder, change the reference accordingly. It's inmportant to prepare the reference for use, with bwa index, Picard CreateSequenceDictionary and Samtools faidx (these files should be in same folder as reference)
ref=<path_to_>/reference.fasta

#######
# start of the actual preprocessing, mapping and cleaning part of the pipeline.
# before starting the next part of the pipeline metrics and quality checks have to be evaluated.

# Step 1:
# quality control of the raw fastq reads, with fastqc. these output files can go into multiqc, however manual check is need for features not supported in multiqc. NOTE: the output name cannot be controlled in the command, thus make sure that the variables are pointing to the correct files.
date 
echo 'Step 1: quality control of the raw fastq reads, with fastqc'
echo 'read_R1'
fastqc ${read_R1} --outdir=${fast_qc}

echo 'read_R2'
fastqc ${read_R2} --outdir=${fast_qc}

# Checkpoint 1, checks whether the correct files are made and contain data
date
echo 'Checkpoint 1'

if test -s ${fast_qc}/${name}*html; then
    echo "fastqc of R1 done"
else echo "failed"
    exit 1
fi

if test -s ${fast_qc}/${name}*html; then
    echo "fastqc of R2 done"
else echo "failed"
    exit 1
fi

# Step 2:
# compiling the reads from the fastq files in unaligned bam files, in this step we set several variables inside the file, taken from the preproc_map_clean_variables.txt.
date 
echo 'Step 2: creating ubam'

#change directory to the directory in which the ubams are deposited
cd ${ubam}

java -Xms20G -Xmx20G -jar ${picard} FastqToSam F1=${read_R1} F2=${read_R2} O=${name}_unmapped.bam RG=${rgroup} SM=${sample_name} LB=${lib} PU=${rgroup} PL=illumina 2> ${log}/${name}_unmapped.log

# Checkpoint 2, check-points for 1) checking if "exception" appeared in the log, which means problems and 2) controlling if the file is made and contains data.
date
echo 'Checkpoint 2'

if grep -i exception ${log}/${name}_unmapped.log ; then 
    echo "failed"
    exit 1
else echo "no problem creating ubam"
fi

if test -s ${name}_unmapped.bam; then
    echo "ubam exists and contains data"
else echo "failed"
    exit 1
fi

# Step 3:
# mark the illumina adapters, with a xt flag. NOTE: the input bam must be query sorted. The following adapters are checked, the default: indexed, dual-indexed and paired.
date
echo 'Step 3: adapter cleaning'

# change directory to the directory in which the xtbams are deposited
cd ${xtbam}

java -Xms20G -Xmx20G -jar ${picard} MarkIlluminaAdapters I=${ubam}/${name}_unmapped.bam O=${name}_unmapped_xt.bam M=${metrix}/${name}_adapter_metrics.txt 2> ${log}/${name}_adapterxt.log

# Checkpoint 3
date
echo 'Checkpoint 3'

if grep -i exception ${log}/${name}_adapterxt.log ; then 
    echo "failed"
    exit 1
else echo "no problem creating the xtbam"
fi

if test -s ${name}_unmapped_xt.bam; then
    echo "adapter filtered bam exists and contains data"
else echo "failed"
    exit 1
fi

# NOTE: the adapter_metrics.txt are used to make histograms for each sample in R, there is a separate script for this, after all samples have completed running this can be merged into 1 graph.

# Step 4:
# Pipe that first makes the produced xtbams into a interleaved fastq (as bwa mem only takes this) -> mapping program, bwa mem -> merge info from the ubam to the mapped bam, with MergeBamAlignment.
# A pipe is used here to reduce the amount of files made. 
# Here the first program sets the bases flagged with xt to be part of adapters to have a quality of 2, thus they are not hard trimmed but soft trimmed, when they are mapped. 
# the resulting bam is query sorted to provide the correct input for markduplicates in order to have the desired output (see explanation below).
date 
echo 'Step 4: mapping and cleaning'

java -Xms20G -Xmx20G -jar ${picard} SamToFastq INPUT=${name}_unmapped_xt.bam FASTQ=/dev/stdout CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true 2> ${log}/${name}_pipe.log | 
bwa mem -M -t 6 -p ${ref} /dev/stdin 2>> ${log}/${name}_pipe.log |
java -Xms20G -Xmx20G -jar ${picard} MergeBamAlignment ALIGNED_BAM=/dev/stdin UNMAPPED_BAM=${ubam}/${name}_unmapped.bam OUTPUT=${name}_clean.bam R=${ref} ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true PRIMARY_ALIGNMENT_STRATEGY=MostDistant INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 ATTRIBUTES_TO_RETAIN=XS SO=queryname 2>> ${log}/${name}_pipe.log

# Checkpoint 4
date
echo 'Checkpoint 4'

if grep -i exception ${log}/${name}_pipe.log ; then 
    echo "failed"
    exit 1
else echo "no problem creating the clean bam"
fi

if test -s ${name}_clean.bam; then
    echo "clean bam exists and contains data"
else echo "failed"
    exit 1
fi

# Step 5:
# mark duplicates, in this pipeline it is also done on secondary and supplementary reads, hence the importance of a query sorted input.
# This step only marks the duplicate reads, it does not remove them from the bam file. Haplotypecaller will take account of the duplicate flags and ignores them. 
# this step might be excluded when having a PCR-free library.
# Based on the sequencing platform there are differences in OPTICAL_DUPLICATE_PIXEL_DISTANCE: nextseq is the default setting at 100; hiseq and novaseq both should have 2500. 
# settings are a result of an investigation in march 2019.
date
echo 'Step 5: mark duplicates'

java -Xms20G -Xmx20G -jar ${picard} MarkDuplicates I=${name}_clean.bam O=${name}_ddupl_clean.bam OPTICAL_DUPLICATE_PIXEL_DISTANCE=${optic} M=${metrix}/${name}_ddupl.txt 2> ${log}/${name}_markdupl.log

#Checkpoint 5
date
echo 'Checkpoint 5'

if grep -i exception ${log}/${name}_markdupl.log ; then 
    echo "failed"
    exit 1
else echo "no problem marking duplicates"
fi

if test -s ${name}_ddupl_clean.bam; then
    echo "ddupl bam exists and contains data"
else echo "failed"
    exit 1
fi

# Step 6:
# sorting the bams by coordinates, creating thereby the finished file.
date 
echo 'Step 6: sorting the final bam'

# change into the last final bam directory.
cd ${finalbam}

#SortSam to get coordinate sorted bam, this is needed as the next steps use a coordinate sorted bam
java -Xms20G -Xmx20G -jar ${picard} SortSam I=${xtbam}/${name}_ddupl_clean.bam O=${gbam} SO=coordinate CREATE_INDEX=true

# Checkpoint 6
date 
echo 'Checkpoint 6'

if test -s ${gbam} ; then
    echo "final bam exists and contains data"
else echo "failed"
    exit 1
fi


echo 'the final bam file is created, next steps are quality control of the created file' 

#######
# after production of the final bam, quality control of the created file needs to be performed. 
# most of these outputs can be summarized with multiqc. Here a multiqc report is produced for each sample, but after all samples have run one can make a summary report in multiqc from all samples.

# Step 7:
# Validation of the bam file, here the WARNING for the NM is ignored as this will most likely come up, but isn't important to be corrected for downstream purposes. 
# However if there is an error with this one can run Picard's SetNmMDAndUqTags, after the files has been sorted with SortSam, most likely the tag has been incorrectly computed because of specifying a query sorted output from MergBamAlignement
date 
echo 'Step 7: quality control'
echo 'validation of bam'
java -Xms20G -Xmx20G -jar ${picard} ValidateSamFile I=${gbam} O=${metrics}/${name}_validatesamfile R=${ref} IGNORE=MISSING_TAG_NM

# collect metrics to evaluate the quality of the alignment 
# CollectMultipleMetrics default setting collects: CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, QualityScoreDistribution, MeanQualityByCycle, CollectBaseDistributionByCycle, here there is an addition of the following metric: CollectGcBiasMetrics
date
echo 'collection of metrics of alignment quality'
java -Xms20G -Xmx20G -jar ${picard} CollectMultipleMetrics I=${gbam} O=${metrics}/${name}_multiple_metrics R=${ref} Program=CollectGcBiasMetrics

java -Xms20G -Xmx20G -jar ${picard} CollectRawWgsMetrics I=${gbam} O=${metrics}/${name}_raw_wgs_metrics R=${ref} INCLUDE_BQ_HISTOGRAM=true

java -Xms20G -Xmx20G -jar ${picard} CollectWgsMetrics I=${gbam} O=${metrics}/${name}_wgs_metrics R=${ref} INCLUDE_BQ_HISTOGRAM=true

# the program version used for Picard and multiqc have an incompatibility, requires a small adjustment to the output file of CollectWgsMetrics, as multiqc will not recognise the output.

sed -i -e 's/picard.analysis.WgsMetrics/picard.analysis.CollectWgsMetrics$WgsMetrics/' ${metrics}/${name}_wgs_metrics 

# qualimap bamqc, the quality control of the bam
date
echo 'qualimap quality control bam'

# change into the folder in which the qualimap results should be deposited
cd ${qualimap}

# running qualimap with bamqc function
qualimap --java-mem-size=20G bamqc -bam ${finalbam}/${gbam} -c -outdir ${name}

#Step 8:
# compile everything produced under step 7 into one multiqc report, this is per run! 
date
echo 'Step 8: compile multiqc'

# change back to the finalbam directory.
cd ${finalbam}

multiqc ${fast_qc}/${name}* ${metrix}/${name}* ${metrics}/${name}* ${qualimap}/${name}* -i ${name} -f -o ${multiqc_data}

date
echo 'bam is produced and quality control has run, evaluate the metrics'


# NOTE: if there have been multiple runs of a sample, these finalbams will have to be merged after this script is finished, continue with script 1a.merge_bam_script.pl, this has to be done before haplotypecaller is used. On the merged bam quality control will also be run.

# NOTE: evaluate the adapter metrics and the mark duplicate metrics

# When all samples have run go to next step haplotypecaller, see script 2.haplo_script.pl

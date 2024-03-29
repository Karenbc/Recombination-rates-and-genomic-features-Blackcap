README.md for variant calling pipeline (SNPs and INDELs): from preprocessing to quality control of variants
by Andrea Bours, 2021/06


This README informs in which order to run the scripts in the pipeline, and which additional files (within one example line) need to be made to run the scripts in the correct order. Any other information that is not provided in the script can be found here. 


1) preprocessing of reads, mapping to reference genome and quality control of mapped reads

	script: 1.preproc_map_clean_script.pl
	additional txt: preproc_map_clean_variables.txt

example line for preproc_map_clean_variables.txt (tab separated):

sample sample <path_to_dir>/sample_R1.fastq.gz <path_to_dir>/sample_R2.fastq.gz <library_name> <read_group> 100 72

rgroup information you get from the reads -> important if more lanes are used.
runtime estimates, are based on running the programs at 20G max, while requesting 35G; and based on mapping the same species:
					if R1 and R2 each are ~5 GB (total ~10GB) runtime with buffer = 54h 
					if R1 and R2 each are ~10 GB (total  ~20GB) runtime with buffer = 72h
					if R1 and R2 each are ~15 GB (total ~30 GB) runtime with buffer = 90h

calculations were based on 10GB per R1/R2 (total 20GB) runs taking about 24hx3 to provide buffer time = 72h,
increase or decrease of 5GB per R1/R2 (total in- or decrease of 10GB) is 18h more or less runtime

Evaluate the quality of the produced files before proceeding to step 2 (unless samples have been resequenced multiple times see step 1a first)

The only intermediate file that should go into the archive directory alongside the final.bams from this pipeline is the clean.bams. This as the other intermediate files can be made easily, but bwa mem has a random part in the mapping that could result in slightly different files upon rerun of the mapping step

NOTE when having a PCR-free library, the MarkDuplicates step can be excluded


1a) when samples have been resequenced multiple times or run in multiple lanes (at the same time) the resulting bams need to be merged
NOTE: in this pipeline initial bams need to be distinguishable from each other, the final bam will have only the sample name
this currently only merges 2 bams into one, can be expanded upon

	script: 1a.merge_bam_script.pl
	additional txt: merge_bam.txt

example line for merge_bam.txt (tab separated):

sample_1_1 sample_1_2 sample_1

evaluation of quality control on merged bams is important for downstream analyses


2) calling haplotypes from the mapped reads -> making g.vcf per sample (be aware that this g.vcf is specific to GATK)

	script: 2.haplo_script.pl
	additional txt: haplo_variable.txt

example line for haplo_variable.txt (tab separated):

A251015 170

runtime estimates, are based on running the programs at max 32G, while requesting 40G; and based on mapping on the same species:
					if bam is 15GB, runtime with buffer = 160h
					if bam is 20GB, runtime with buffer = 170h
					if bam is 25GB, runtime with buffer = 180h

calculations were based on 15GB bam taking about 110h to provide buffer time add 50 hours, 
increase or decrease of 5GB in bam size, amounts to a 10h more or less runtime
(NOTE the program doesn't speed up when adding more memory)


3) combine per sample g.vcfs into one g.vcfs, but splitting the reference into 10 parts
 NOTE the reference genome is as evenly divided into ten parts as can be done without losing the order of chromosomes/scaffolds and no splitting of chromosomes
 this is done to speed up the pipeline, however in this step and the next never separate based on samples always separated by reference.

	script: 3.combinegvcfs_script.pl
	additional txt: sample.list
				[1st-10th]_intervals.list

example line for sample.list:
<path_to_dir>/variant_calling/gvcf/1204.g.vcf.gz

example line for [1]_intervals.list
chr_1

before running this script a personal folder in the beegfs has to be made


4) genotype the g.vcfs, resulting in vcfs

	script: 4.genotypegvcfs_script.pl
	additional txt: [1st-10th]_intervals.list (same file as above)


4a) genotype the g.vcfs, but leaving in the non-variant sites

	script: 4a.genotypegvcfs_all_sites_script.pl
	additional txt: [1st-10th]_intervals.list

After these scripts have been run copy the files from beegfs to your home directory as beegfs is not backed up, ensure that your copy finsihed correctly by running md5sum on the files


5) gather the vcfs 1st-10th into 1 vcf 

	script: 5.gathervcf_script.sh
	additional txt: whole_genome_intervals.list

The whole_genome_intervals.list is all chromosomes and scaffolds


5a) gather the all_sites_vcfs 1st-10th into 1 all_sites_vcf

	script: 5a.gathervcf_all_sites_script.sh
	additional txt: whole_genome_intervals.list (same as file above)


6) subset the vcf into the different variant categories as needed for analyses and hard filtering
	
6a) subset the vcf into only SNP variants, in order to filter later on

	script: 6a.subsetvcf_table_snp_script.sh
	additional txt: whole_genome_intervals.list (same as file above)

script creates a SNP.table to load into R to base hardfilters upon


6b) subset the vcf into only (short)-INDEL variants, in order to filter later on
 NOTE In general INDEL vcf will be much smaller, therefore less memory and time requested

	script: 6b.subsetvcf_table_indel_script.sh
	additional txt: whole_genome_intervals.list (same as file above)

script creates a INDEL.table to load into R to base hardfilters upon
	

6c) prepare the table for filtering the all_sites.vcf, for this the INDELs will have to be removed 
 NOTE In general this vcf will be much bigger, therefore more memory and time requested

	script: 6c.subsetvcf_table_all_sites_script.sh
	additional txt: whole_genome_intervals.list (same as file above)

script creates a all_sites.table to load into R to check hardfilters upon


7) render html reports to evaluate possible hard filter criteria per vcf

	scripts: 7.hard_filter.Rmd and 7a.run_hard_filter_Rmd.R

The .Rmd file creates a report based on the specifics parameters set when rendering the markdown report via the .R file.

The gatk hard filters differ depending on the type of variants to be filtered, these parameters (params) can be provided as the ..._gatk. While evaluating different values for these filters can be evaluated by changing the ..._sugg parameters (params) in the .R file (or in the .Rmd).

					at the time of hard filtering gatk recommended the following filters:
					SNP: QD<2.0, FS>60.0, SOR>3.0, MQ<40.0, MQRankSum<-12.5, ReadPosRankSum<-8.0
					INDEL:QD<2.0, FS>200.0, ReadPosRankSum<-20.0 

					I adjusted the filters as follows:
					SNP: QD<2.5, FS>45, SOR>3.0, MQ<40.0, MQRankSum<-12.5, ReadPosRankSum<-8.0
					INDEL: QD<2.5, FS>200.0, ReadPosRankSum<-20.0, with additional filter MQ<40.0

					For the all_sites vcf I used the same filters as for the SNPs to ensure that the same SNPs are present between the two vcfs


This file can be used for the different variant tables created in step 6a/b/c.


8) Hard filter the different vcfs to ensure that only "high" quality variants are in the vcf
	
8a) hardfilter the SNP vcf

	script: 8a.hardfilter_snp_script.sh
	additional txt: whole_genome_intervals.list (same as file above)
	

8b) hardfilter the INDEL vcf

	script: 8b.hardfilter_indel_script.sh
	additional txt: whole_genome_intervals.list (same as file above)
	

8c) hardfilter the all_sites vcf
 NOTE here I filter with the SNP specific filter, as this is better suited for the specific analyses I will do on this data.

	script: 8c.hardfilter_all_sites_script.sh
	additional txt: whole_genome_intervals.list (same as file above)



The VCFs are ready for analyses, be aware that additional filtering may need to be performed depending on the analysis

#!/bin/bash
#
#  find homologs in swissprot database to found proteins
# 
#SBATCH --job-name=blast
#  how many cpus are requested
#SBATCH --nodes=1
#SBATCH --ntasks=12
#  maximum walltime, here 10min
#SBATCH --time=100:00:00
#  maximum requested memory
#SBATCH --mem=20G
#  write std out and std error to these files
#SBATCH --error=blastp.%J.err
#SBATCH --output=blastp.%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<your name>@evolbio.mpg.de
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard
 
blastp -query bSylAtr.filtered.proteins.fasta -db swissprot -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out output.blastp

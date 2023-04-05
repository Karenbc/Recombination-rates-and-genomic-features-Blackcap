#!/bin/bash
#
#  find interpro domains, pfam domains and Gene Ontologies of the found proteins
# 
#SBATCH --job-name=iprscan
#  how many cpus are requested
#SBATCH --nodes=1
#SBATCH --ntasks=12
#  maximum walltime, here 10min
#SBATCH --time=100:00:00
#  maximum requested memory
#SBATCH --mem=20G
#  write std out and std error to these files
#SBATCH --error=iprscan.%J.err
#SBATCH --output=iprscan.%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<your name>@evolbio.mpg.de
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

module load java/x64/11u1 

/data/biosoftware/InterProScan/interproscan-5.48-83.0/interproscan.sh -appl pfam -dp -f TSV -goterms -iprlookup -pa -t p -i bSylAtr.all.maker.proteins.fasta -o output.iprscan -cpu 12

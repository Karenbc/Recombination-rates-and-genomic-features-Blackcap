#!/bin/bash
#
#  run maker
# 
#SBATCH --job-name=maker
#  how many cpus are requested
#SBATCH --nodes=12
#  maximum walltime, here 10min
#SBATCH --time=200:00:00
#  maximum requested memory
#SBATCH --mem=20G
#  write std out and std error to these files
#SBATCH --error=maker.%J.err
#SBATCH --output=maker.%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<your name>@evolbio.mpg.de
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

mpiexec /data/biosoftware/maker/maker3mpi/maker/bin/maker -fix_nucleotides -base bSylAtr maker_opts.ctl maker_bopts.ctl maker_exe.ctl

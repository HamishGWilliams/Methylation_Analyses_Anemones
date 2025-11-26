#!/bin/bash
#SBATCH --mem 8G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=david.fisher@abdn.ac.uk
#SBATCH --time=1-00:00:00

# sbatch /uoa/home/s02df9/equina_meth_analysis/scripts/run_loop_gunzip_cpg.sh

cd /uoa/home/s02df9/equina_meth_analysis/'NEOF Methyl seq'/Trimmed

for d in ./*/ ; do
        (cd "$d" && gunzip -c bedgraph_output_cpg.gz.bismark.cov.gz  > meth_CpG_cov_reads );  done



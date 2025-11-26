#!/bin/bash
#SBATCH --mem 96G
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=david.fisher@abdn.ac.uk
#SBATCH --time=5-00:00:00

# sbatch /uoa/home/s02df9/equina_meth_analysis/scripts/run_bismark_extract_loop.sh

module load  bowtie2/2.4.2
module load  bismark/0.23.0


cd /uoa/home/s02df9/equina_meth_analysis/'NEOF Methyl seq'/Trimmed

for d in ./*/ ; do
        (cd "$d" && bismark_methylation_extractor    -p --no_overlap --comprehensive --report --multicore 8  *.bam ) ;
                    done

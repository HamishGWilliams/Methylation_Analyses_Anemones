#!/bin/bash
#SBATCH --mem 24G
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=david.fisher@abdn.ac.uk
#SBATCH --time=1-00:00:00

# sbatch /uoa/home/s02df9/equina_meth_analysis/scripts/run_loop_bismark2bedgraph_cphg.sh

module load  bowtie2/2.4.2
module load  bismark/0.23.0


cd /uoa/home/s02df9/equina_meth_analysis/'NEOF Methyl seq'/Trimmed

for d in ./*/ ; do
        (cd "$d" && bismark2bedGraph --ample_memory --CX --cutoff 5 --output bedgraph_output_cpg CpG_*_pe.txt );

                    done

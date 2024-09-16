#!/bin/bash
#SBATCH --mem 96G
#SBATCH --partition uoa-compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=h.williams.22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/slurm_outputs
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/Methylation_Analyses/Methylation_Analyses_Anemones/slurm_errors/%x_%j.err

# sbatch /uoa/home/r02hw22/Equina_Methylation_Analysis/Scripts/run_loop_bismark2bedgraph_chg.sh

module load  bowtie2/2.4.2
module load  bismark/0.23.0

cd /uoa/home/r02hw22/Equina_Methylation_Analysis/Data2/Trimmed

# Try running this with only 1 sample first
# Try running each of the configurations in an interactive sessions at a time and see if it works

for d in ./*/ ; do
	# Process CHG methylated cystosines
    (cd "$d" && bismark2bedGraph --ample_memory --CX --cutoff 5 --output bedgraph_output_chg CHG_*_pe.deduplicated.txt);

    done

# Try hashing out this second loop first to make sure bismark2bedgraph is working.

for d in ./*/ ; do
    # Unzip the CHG bedgraph file
    (cd "$d" && gunzip -c bedgraph_output_chg.gz.bismark.cov.gz > meth_CHG_cov_reads );

    done


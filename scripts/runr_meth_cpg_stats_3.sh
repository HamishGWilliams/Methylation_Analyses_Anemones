#!/usr/bin/env Rscript
#SBATCH --mem 24G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=david.fisher@abdn.ac.uk
#SBATCH --time=1-00:00:00

#Rscript /uoa/home/s02df9/equina_meth_analysis/scripts/runr_meth_cpg_stats_3.sh

library(dplyr)

setwd("/uoa/home/s02df9/equina_meth_analysis/outputs")

3_cpgs_raw = data.table::fread("../NEOF Methyl seq/Trimmed/Sample_3-3_D/meth_cov_cpg_reads", header = F, skip = 1)

colnames(3_cpgs_raw) = c("chrom", "start", "end", "pc_meth", "count_meth", "count_un")

3_cpgs = 3_cpgs_raw %>%
           mutate(coverage = count_meth + count_un) %>%
           filter(coverage %in% (5:29), count_meth>0) %>%  #only test CpGs with any methylated reads. 
      #So can't just take mean later, need to sum this and divide by total Cs from raw dataset
           rowwise() %>% #necessary to make binom test work on each row rather than whole columns
           mutate(pv_meth = binom.test(count_meth, coverage, 0.01, alternative = "greater")$p.value) %>%
           #what we're testing is whether the observed freq is sig HIGHER from the error rate of 0.01.
           ungroup() %>%
           mutate(pv_meth_c = p.adjust(pv_meth, method = "fdr"),  #then correct p-value for number of tests
                  is_meth = ifelse(pv_meth_c < 0.05, 1, 0))    #assign CpG as "methylated" if still clearly higher than error rate
                  #note that for Cs with low coverage (e.g. 5) need 40+ pc_meth to be assigned meth status
                  #while for Cs with higher coverage (e.g. 28) 7.14 % is enough
                  
mean_meth_rate =  sum(3_cpgs$is_meth)  / nrow(3_cpgs_raw)     #this appears to work

write.table(mean_meth_rate,
           "3_mean_cpg_meth.txt") #overall methylation rate of Cs
           
n_c_chrom =  3_cpgs_raw %>%
             group_by(chrom) %>%
             summarise(n = n())

mean_meth_rate_chrom = 3_cpgs %>%
                group_by(chrom) %>%
                summarise(total_meth = sum(is_meth)) %>%
                     left_join(n_c_chrom, by = "chrom") %>%
                     mutate(mean = total_meth / n)

total_meth_chrom = 3_cpgs %>%
                group_by(chrom) %>%
                summarise(total_meth = sum(is_meth))

mean_meth_rate_chrom =  3_cpgs_raw %>%
             mutate(coverage = count_meth + count_un) %>%
             filter(coverage %in% (5:29))  %>%
             group_by(chrom) %>%
             summarise(n = n()) %>%
             left_join(total_meth_chrom, by = "chrom") %>%
             mutate(mean = total_meth / n)

write.table(mean_meth_rate_chrom,
           "3_mean_cpg_meth_bychrom.txt")


#per 1000 base pairs

total_meth_1tbps = 3_cpgs %>%  
                 group_by(chrom, group = cut(start, breaks = seq(0, max(start), 1000))) %>%
                 summarise(total_meth = sum(is_meth))

mean_meth_rate_1tbps = 3_cpgs_raw %>%
           mutate(coverage = count_meth + count_un) %>%
           filter(coverage %in% (5:29)) %>%
           group_by(chrom, group = cut(start, breaks = seq(0, max(start), 1000))) %>%
           summarise(n = n() )   %>%
           left_join(total_meth_1tbps, by = c("group", "chrom")) %>%
           mutate(mean = total_meth / n)

write.table(mean_meth_rate_1tbps,
           "3_mean_cpg_meth_by1tbps.txt")

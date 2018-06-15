#!/bin/bash
snakemake --use-conda --cluster-config cluster_config.yaml --cluster "qsub -pe smp {cluster.num_cores} -l h_vmem={cluster.maxvmem}" -j 50 -p $@ 2> snakemake.log
mail -s "snakemake finished" heath.obrien@gmail.com < snakemake.log

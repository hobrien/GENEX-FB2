#!/bin/bash

snakemake --use-conda --cluster-config cluster_config.yaml --cluster "qsub -pe smp {cluster.num_cores} -l h_vmem={cluster.maxvmem}" -j 20 $@ 2> snakemake.log
mail -s "snakemake finished" heath.obrien@gmail.com < snakemake.log

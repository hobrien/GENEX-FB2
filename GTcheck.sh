#!/bin/bash

snakemake -p -s Snakefiles/GTcheck.smk --use-conda --cluster-config cluster_config.yaml --cluster "qsub -pe smp {cluster.num_cores} -l h_vmem={cluster.maxvmem}" -j 20 $@ 2> GTcheck.log
mail -s "GTcheck finished" $EMAIL < GTcheck.log

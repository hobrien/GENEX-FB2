#!/bin/bash

snakemake -s Snakefiles/eQTL.smk --use-conda --cluster-config cluster_config.yaml --cluster "qsub -pe smp {cluster.num_cores} -l h_vmem={cluster.maxvmem}" -j 20 $@ 2> eQTL.log
mail -s "eQTL analysis finished" $EMAIL < eQTL.log

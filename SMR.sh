#!/bin/bash

snakemake -p -s Snakefiles/SMR.smk --use-conda --cluster-config cluster_config.yaml --cluster "qsub -pe smp {cluster.num_cores} -l h_vmem={cluster.maxvmem}" -j 20 $@ 2> SMR.log
mail -s "SMR finished" $EMAIL < SMR.log

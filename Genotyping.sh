#!/bin/bash

snakemake -p -s Snakefiles/Genotyping.smk --use-conda --cluster-config cluster_config.yaml --cluster "qsub -pe smp {cluster.num_cores} -l h_vmem={cluster.maxvmem}" -j 20 $@ 2> genotyping.log
mail -s "genotyping QC finished" $EMAIL < genotyping.log

#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash

source activate py35
snakemake -p --cluster "qsub" -j 50

#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -l h_vmem=40G
#

Rscript R/MatrixEQTL.R  --genotypes Genotypes/Plink/recoded.traw \
  --counts Data/MalevsFemale.complete.txt --snps Genotypes/Plink/genotypes.map 
  --genes Data/geneloc.txt --cofactors Data/SampleInfo.txt \
  --cis MatrixEQTL/cis_eqtl.txt --trans MatrixEQTL/trans_eqtl.txt \
  --image MatrixEQTL/results.RData

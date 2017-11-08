#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -l h_vmem=30G
#

source activate peer2
Rscript /c8000xd3/rnaseq-heath/GENEX-FB2/R/PEER.R  -p diploid -n 25 \
  -e /c8000xd3/rnaseq-heath/GENEX-FB2/Data/counts_vst.txt \
  -o /c8000xd3/rnaseq-heath/GENEX-FB2/factors.txt \
  -r /c8000xd3/rnaseq-heath/GENEX-FB2/residuals.txt \
  -a /c8000xd3/rnaseq-heath/GENEX-FB2/alpha.txt \
  -c /c8000xd3/rnaseq-heath/GENEX-FB2/Data/SampleInfo.1.txt


# Workflow for eQTL analyses in the developing human brain

- This workflow consists of a series of [Snakemake](https://snakemake.readthedocs.io/en/stable) workflows
that can be executed using the supplied bash scripts:

- `Genotyping.sh`: run initial genotyping QC and prepare data for imputation
- `eQTL.sh`: QC and filtering of post-imputation genotypes and expression data; eQTL analysis
- `GTcheck.sh`: map reads to reference, call snps, and compare to imputed genotypes
    - This needs to be run after eQTL analysis because it requires processed vcf files from eTQL workflow
- `SMR.sh`: SMR analyses of eQTL data (requires config file that will be made available on publication)

- Expression data were created from the workflows detailed [here](https://github.com/hobrien/GENEX-FB1)

- Software used:
    - [bcftools](http://www.htslib.org/doc/bcftools.html) v1.6
    - [bedtools](http://bedtools.readthedocs.io/en/latest) v2.26.0
    - [bgzip](http://www.htslib.org/doc/bgzip.html) v1.3.1
    - [checkVCF.py](https://github.com/zhanxw/checkVCF) v1.4
    - [CrossMap](http://crossmap.sourceforge.net) v0.2.5
    - [fastQTL](http://fastqtl.sourceforge.net) v2.184_gtex
    - [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) v2.0.5
    - [pandas](https://pandas.pydata.org) v0.21.0
    - [PEER](https://github.com/PMBio/peer/wiki) v1.3
    - [PLINK](https://www.cog-genomics.org/plink2) v1.90b3.39
    - [samtools](http://www.htslib.org/doc/samtools.html) v1.3.1
    - [smr](http://cnsgenomics.com/software/smr) v0.69
    - [tabix](http://www.htslib.org/doc/tabix.html) v1.3.1
    - [tidyverse](https://www.tidyverse.org) v1.1.1
    - [vcf-sort](https://github.com/vcftools/vcftools/blob/master/src/perl/vcf-sort)


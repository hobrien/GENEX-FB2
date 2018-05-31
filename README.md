# Workflow for eQTL analyses in the developing human brain

- This workflow consists of a series of [Snakemake](https://snakemake.readthedocs.io/en/stable) workflows
that can be executed using the supplied bash scripts:

- `Genotyping.sh`: run initial genotyping QC and prepare data for imputation
- `eQTL.sh`: QC and filtering of post-imputation genotypes and expression data; eQTL analysis
- `GTcheck.sh`: map reads to reference, call snps, and compare to imputed genotypes
    - This needs to be run after eQTL analysis because it requires processed vcf files from eTQL workflow
- `SMR.sh`: SMR analyses of eQTL data (requires config file that will be made available on publication)

- Expression data were created from the workflows detailed [here](https://github.com/hobrien/GENEX-FB1)

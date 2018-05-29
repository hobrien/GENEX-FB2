
import warnings
import math

configfile: "config.yaml"

rule all:
    input:
        "Genotypes/Raw/PLINK_040517_1054/FB_new_excl_blank.bed"
        
rule make_bed:
    input:
        "Genotypes/Raw/PLINK_040517_1054/Nick_Bray_OmniExpress_230317.ped"
    output:
        "Genotypes/Raw/PLINK_040517_1054/FB_new.bed"
    params:
        input = "Genotypes/Raw/PLINK_040517_1054/Nick_Bray_OmniExpress_230317"
        output = "Genotypes/Raw/PLINK_040517_1054/FB_new"
    shell:
        "plink --file {params.input} --make-bed --out {parmas.output}"

rule exclude_blank:
    input: 
        rules.make_bed.output
    output:
        "Genotypes/Raw/PLINK_040517_1054/FB_new_excl_blank.bed"
    params:
        input = "Genotypes/Raw/PLINK_040517_1054/FB_new"
        output = "Genotypes/Raw/PLINK_040517_1054/FB_new_excl_blank"
        exclude = "<(echo 24; seq 45 48)"
    shell:
        "plink --bfile {params.input} --remove-fam {params.exclude} --make-bed --out {parmas.output}"

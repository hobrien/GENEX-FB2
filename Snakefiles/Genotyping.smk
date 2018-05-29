
import warnings
import math

configfile: "config.yaml"

rule all:
    input:
        "Genotypes/Raw/Merged/FB_all.het",
        "Genotypes/Raw/Merged/FB_all.imiss",
        "Genotypes/Raw/Merged/FB_all.sexcheck"
        
rule make_bed:
    input:
        "Genotypes/Raw/PLINK_040517_1054/Nick_Bray_OmniExpress_230317.ped"
    output:
        "Genotypes/Raw/PLINK_040517_1054/FB_new.bed"
    params:
        input = "Genotypes/Raw/PLINK_040517_1054/Nick_Bray_OmniExpress_230317",
        output = "Genotypes/Raw/PLINK_040517_1054/FB_new"
    shell:
        "plink --file {params.input} --make-bed --out {params.output}"

rule exclude_blank:
    input: 
        rules.make_bed.output
    output:
        "Genotypes/Raw/PLINK_040517_1054/FB_new_excl_blank.bed"
    params:
        input = "Genotypes/Raw/PLINK_040517_1054/FB_new",
        output = "Genotypes/Raw/PLINK_040517_1054/FB_new_excl_blank",
        exclude = "<(echo 24; seq 45 48)"
    shell:
        "plink --bfile {params.input} --remove-fam {params.exclude} --make-bed --out {params.output}"

rule merge_plink:
    input:
        set1 = "Genotypes/Raw/Set1/FB_Merged.bed",
        set2 = rules.exclude_blank.output
    output:
        "Genotypes/Raw/Merged/FB_all.bed"
    params:
        set1 = "Genotypes/Raw/Set1/FB_Merged",
        set2 = "Genotypes/Raw/PLINK_040517_1054/FB_new_excl_blank",
        output = "Genotypes/Raw/Merged/FB_all",
    shell:
        "plink --bfile {params.set2} --bmerge {params.set1} --make-bed --out {params.output}"

rule check_het:
    input:
        rules.merge_plink.output
    output:
        "Genotypes/Raw/Merged/FB_all.het"
    params:
        prefix = "Genotypes/Raw/Merged/FB_all"
    shell:
        "plink --bfile {params.prefix} --het --out {params.prefix}"
       
rule check_imiss:
    input:
        rules.merge_plink.output
    output:
        "Genotypes/Raw/Merged/FB_all.imiss"
    params:
        prefix = "Genotypes/Raw/Merged/FB_all"
    shell:
        "plink --bfile {params.prefix} --missing --out {params.prefix}"
        
rule sex_check:
    input:
        sample_info = config['sample_info'],
        genotypes = rules.merge_plink.output
    output:
        "Genotypes/Raw/Merged/FB_all.sexcheck"
    params:
        sex_info = "Genotypes/Raw/Merged/sex.txt",
        prefix = "Genotypes/Raw/Merged/FB_all",
    shell:
        "cat {input.sample_info} | cut -f 1,2,5 | perl -pe 's/Male/1/; s/Female/2/; s/NA/0/' > {params.sex_info}; "  
        "plink --bfile {params.prefix} --update-sex {params.sex_info} --make-bed --out {params.prefix}; "
        "plink --bfile {params.prefix} --check-sex --out {params.prefix}"

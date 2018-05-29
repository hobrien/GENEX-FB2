import yaml
import os
configfile: "config.yaml"
SMR=yaml.load(open('smr.yaml', 'r'))

rule all:
    input: 
       expand("SMR/mysmr_{level}_{gwas}_all.smr.gz", level = ['gene', 'transcript'], gwas=config['GWAS']),
       expand("SMR/plot/clozuk_{level}.{gene_id}.txt", level=['gene'], gene_id=SMR['gene']['clozuk']),
       expand("SMR/plot/clozuk_{level}.{gene_id}.txt", level=['transcript'], gene_id=SMR['transcript']['clozuk']),
       expand("SMR/plot/neuroticism_{level}.{gene_id}.txt", level=['gene'], gene_id=SMR['gene']['neuroticism']),
       expand("SMR/plot/neuroticism_{level}.{gene_id}.txt", level=['transcript'], gene_id=SMR['transcript']['neuroticism']),
       expand("SMR/plot/bip_{level}.{gene_id}.txt", level=['gene'], gene_id=SMR['gene']['bip']),
       expand("SMR/plot/mdd_{level}.{gene_id}.txt", level=['transcript'], gene_id=SMR['transcript']['mdd'])

rule prepare_smr:
    input:
        "FastQTL/all_eqtls_{level}.{chunk}_q05.txt.gz",
        "Genotypes/Combined/snp_positions.txt",
        "Data/{level}loc.txt"
    output:
        "SMR/myquery.{level}.{chunk}.txt.gz"
    shell:
        "Rscript R/PrepareSMR.R {input} {output}"

rule cat_eqtls_tr:
    input:
        expand("SMR/myquery.transcript.{chunk}.txt.gz", chunk=range(1,101))
    output:
        ["SMR/myquery.transcript.chr" + str(x+1) + ".txt" for x in range(23)]
    shell:
        "echo -e 'SNP\tChr\tBP\tA1\tA2\tFreq\tProbe\tProbe_Chr\tProbe_bp\tGene\tOrientation\tb\tse\tp' | tee {output}; "
        "zcat {input} | awk '{{print>>\"SMR/myquery.transcript.chr\"$2\".txt\"}}'"

rule cat_eqtls_gene:
    input:
        expand("SMR/myquery.gene.{chunk}.txt.gz", chunk=range(1,101))
    output:
        ["SMR/myquery.gene.chr" + str(x+1) + ".txt" for x in range(23)] 
    shell:
        "echo -e 'SNP\tChr\tBP\tA1\tA2\tFreq\tProbe\tProbe_Chr\tProbe_bp\tGene\tOrientation\tb\tse\tp' | tee {output}; "
        "zcat {input} | awk '{{print>>\"SMR/myquery.gene.chr\"$2\".txt\"}}'"

rule make_besd:
    input:
        "SMR/myquery.{level}.chr{chr}.txt"
    output:
        "SMR/mybesd.{level}.chr{chr}.besd"
    params:
        "SMR/mybesd.{level}.chr{chr}"
    shell:
        "smr --qfile {input} --make-besd --out {params}"

# rule prepare_gwas:
#     input:
#         config['gwas']
#     output:
#         "CLOZUK_recoded.txt"
#     shell:
#         "Rscript R/PrepareGWAS.R {input} {output}"

rule smr:
    input:
        gwas = lambda wildcards: config['GWAS'][wildcards.gwas],
        besd = rules.make_besd.output
    output:
        "SMR/mysmr_{level}_{gwas}.chr{chr}.smr"
    params:
        plink_prefix = "Genotypes/Plink/genotypes",
        besd_prefix = "SMR/mybesd.{level}.chr{chr}",
        out_prefix = "SMR/mysmr_{level}_{gwas}.chr{chr}"
    threads: 1
    shell:
        "smr --bfile {params.plink_prefix} --gwas-summary {input.gwas} "
        "--beqtl-summary {params.besd_prefix} --out {params.out_prefix} "
        "--thread-num {threads} --peqtl-smr 0.01" 

rule cat_smr:
    input:
        lambda wildcards: expand("SMR/mysmr_{level}_{gwas}.chr{chr}.smr", level = wildcards.level, gwas=wildcards.gwas, chr=range(1,23))
    output:
        "SMR/mysmr_{level}_{gwas}_all.smr.gz"
    shell:
        "tail -n+2 {input} | grep -v '==>' | gzip -c > {output}"

rule genloc_smr:
    input:
        "Data/{level}loc.txt"
    output:
        "Data/{level}loc_smr.txt"
    shell:
        "awk -v OFS=\"\t\" '{{sub(/\.[0-9]*/, \"\", $1); sub(/chr/, \"\", $2); print $2, $3, $4, $1, $5}}' {input} > {output}"

rule plot_smr:
    input:
        gwas = lambda wildcards: config['GWAS'][wildcards.gwas],
        genloc = rules.genloc_smr.output,
        besd = lambda wildcards: expand("SMR/mybesd.{level}.chr{chr_num}.besd", level=wildcards.level, chr_num=SMR[wildcards.level][wildcards.gwas][wildcards.gene_id])
    output:
        "SMR/plot/{gwas}_{level}.{gene_id}.txt"
    params:
        plink_prefix = "Genotypes/Plink/genotypes",
        besd_prefix = lambda wildcards: expand("SMR/mybesd.{level}.chr{chr_num}", level=wildcards.level, chr_num=SMR[wildcards.level][wildcards.gwas][wildcards.gene_id]),
        out_prefix = "SMR/{gwas}_{level}",
        genloc = "Data/{level}loc_smr.txt",
        gene_id = "{gene_id}"
    shell:
        "smr --bfile {params.plink_prefix} --gwas-summary {input.gwas} "
        "--beqtl-summary {params.besd_prefix} --out {params.out_prefix} --plot "
"--probe {params.gene_id} --probe-wind 500 --gene-list {params.genloc}"
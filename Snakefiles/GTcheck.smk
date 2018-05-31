from Helpers import get_sequences
configfile: "config.yaml"

files=get_sequences(config['seqfile'])
files = {k: v for k, v in files.items() if k not in ['17046','16385','17048','16024','16115','11449','16972']}

# to make DAG: snakemake -np --dag | dot -Tsvg > dag.svg
dag = 0  # set to 1 to produce DAG image without tons of duplicates 
if dag:
    k,v=files.popitem(); files = {k:v}
else:
   files = {k: v for k, v in files.items() if k not in ['17046','16385','17048','16024','16115','11449','16972']}

"""
rules:
- hisat: map reads to reference
- sort_bam: sort mapping file
- samtools_index: index mapping file
- call_snps: run samtools mpileup and bcftools call to call SNPs from mapping file
- index_vcf: index SNP data 
- gt_check: run bcftools GTcheck
- summarise_gt_check: plot Discordance for all samples
"""

rule all:
    input:
        "GTcheck/summary.png"

rule hisat:
    input:
        reads = lambda wildcards: files[wildcards.sample]
    output:
        temp("BAM/{sample}.bam")
    params:
        input1 = lambda wildcards: ','.join(files[wildcards.sample][::2]),
        input2 = lambda wildcards: ','.join(files[wildcards.sample][1::2]),        
        idx = config['reference']['index'],
        extra = '--known-splicesite-infile ' + config['reference']['splice_sites'],
        threads = 8
    benchmark:
        "Benchmarks/{sample}.hisat.benchmark.txt"
    log:
        "Logs/{sample}_hisat_map.txt"
    shell:
        "(hisat2 {params.extra} --threads {params.threads}"
        " -x {params.idx} -1 {params.input1} -2 {params.input2}"
        " | samtools view -Sbh -o {output} -)"
        " 2> {log}"

rule sort_bam:
    input:
        rules.hisat.output
    output:
        "BAM/{sample}.sort.bam"
    params:
        "-m 4G"
    threads: 8
    wrapper:
        "0.17.4/bio/samtools/sort"
        
rule samtools_index:
    input:
        rules.sort_bam.output
    output:
        "BAM/{sample}.sort.bam.bai"
    wrapper:
        "0.17.4/bio/samtools/index"

rule call_snps:
    input:
        bam=rules.sort_bam.output,
        index=rules.samtools_index.output
    output:
        "SNPcalls/{sample}.vcf.gz"
    params:
        fasta = config['reference']['genome']
    shell:
        "samtools mpileup -uf {params.fasta} {input.bam} | bcftools call -mv -Ov | "
        "awk '{{gsub(/##contig=<ID=chr/,\"##contig=<ID=\"); gsub(/^chr/,\"\"); print}}' | "
        "bcftools view -Oz - > {output}"

rule index_vcf:
    input:
        rules.call_snps.output
    output:
        "SNPcalls/{sample}.vcf.gz.csi"
    shell:
        "bcftools index {input}"

rule gt_check:
    input:
        vcf=rules.call_snps.output,
        index=rules.index_vcf.output
    output:
        "GTcheck/{sample}.tab"
    params:
        sample = "{sample}",
        prefix = "GTcheck/{sample}",
        genotyping="Genotypes/Combined/combined_filtered.vcf.gz"
    shell:
        "bcftools gtcheck -g {params.genotyping} -p {params.prefix} -S {params.sample} {input.vcf} "

rule summarise_gt_check:
    input:
        expand("GTcheck/{sample}.tab", sample=files.keys())
    output:
        plot="GTcheck/summary.png",
        table="GTcheck/summary.tsv"
    shell:
        "Rscript R/SummariseGTcheck.R -o {output.table} -p {output.plot} {input}"


#import yaml
#snakemake --use-conda --cluster-config cluster_config.yaml --cluster "qsub -pe smp {cluster.num_cores} -l h_vmem={cluster.maxvmem}" -j 20

configfile: "config.yaml"

qvals={'10': 0.1, '05': 0.05, '01': 0.01, '001': 0.001, '0001': 0.0001}

#SMR=yaml.load(open('smr.yaml', 'r'))
# to make DAG: snakemake -np --dag | dot -Tsvg > dag.svg
dag = 0
if dag:
    chr_num = 2
    num_permutations = 2
    config['gtex_samples'] = ['Brain_Cortex']
else:
   chr_num = 23
   num_permutations = 101

"""

New plan: 
- filter out non-overlapping SNPs from GTEx sig pairs
- find top SNP out of remaining eQTLs
- pull out p-vals from all_eqtls and from sig eQTLs
rules:
    rename_samples: replace genotyping well IDs with BrianBank IDs in imputed vcf files
    index_vcf: use bcftools to index vcf with renamed samples
    merge_vcf: merge vcf files from separate genotyping/imputation runs
    index_vcf2: use bcftools to index merged vcf files
    add_rsID: replace SNP IDs from imputation with rsIDs
    lift_over: convert from hg19 coordinates to GRCh38 coordinates 
    filter_chr: remove SNPs that do not map to chromosomes
    sort_vcf: coordinate sort vcf (sort order is changed during liftover)
    vcf_check: test for non-variable positions (?) and postions with >2 variants
    excluded_sites: list non variable sites (?)/ sites >2 variants for exclusion
    filter_dup: remove non variable sites (?)/ sites >2 variants
    combine_chromosomes: combine vcf files for individual chromosomes unto a single vcf
    sort_vcf2: coordinate sort combined vcf (no idea why it's unsorted)
    index_vcf3: use bcftools to index combined vcf files
    get_gene_positions: extract coordinates of genes from gtf file
    get_transcript_positions: extract coordinates of genes from gtf file
    filter_counts: filter low expression genes and remove excluded samples from count matrix
    filter_transcript_counts: filter low expression genes and remove excluded samples from count matrix
    bgzip_counts: compress BED file of  counts
    bgzip_transcript_counts: compress BED file of  counts
    index_counts: index BED file of counts
    select_samples: exclude sample from vcf and make sure order is the same as the counts file
    filter_tags: add info about MAF, HWE, etc. and filter SNPs
    snp_positions: extract genomic coordinates for each SNP (to be associated with eQTLs in q_values rule)
    plink_import: create (plink) bed file from vcf for PCA analysis
    plink_ld_prune: LD pruning analysis on plink-formatted genotype data
    plink_pca: PCA analysis on LD-pruned plink-formatted genotype data
    scz_ld: get list of SNPs tagged by schizophrenia index SNPs
    peer: PEER analysis on count data, using cofactors and PCA results
    format_cov:
    peer_nc: PEER analysis on count data, without cofactors or PCA results
    tabix_vcf: index final vcf file
    fast_qtl: run fast_qtl on each chunk of genome (calculate nominal p-values for all SNPs)
    cat_fast_qtl: concatinate fast_qtl nominal pass output
    rule dedup_fast_qtl:
    fast_qtl_permutations: run permutation analysis on each chunk of genome (corrected p-values for top SNPs)
    cat_permutations: concatinate fast_qtl permutation pass output
    q_values: calculate FDR for each eQTL
    rule gtex2bed:
    rule lift_over_bed:
    rule summarise_overlaps:

"""
rule all:
    input:
       "Peer/factors_nc.txt",
       "Genotypes/Plink/scz_ld.tags",
       "Results/snp_summary.txt",
       expand("FastQTL/sig_eqtls_{level}.{chunk}_q{fdr}.gz", level = ['gene', 'transcript'], chunk=range(1,num_permutations), fdr=['05', '01', '001', '0001']),
       expand("FastQTL/sig_snps_{level}_q05.gz", level=['gene', 'transcript']),
       expand("FastQTL/all_snps_{level}.all.txt.gz", level = ['gene', 'transcript']),
       expand("MatrixEQTL/{proximity}_eqtl_{level}.txt", proximity = ['cis', 'trans'], level = ['gene', 'transcript']),
       expand("FastQTL/all_eqtls_{level}.all_q{fdr}.gz", level = ['gene', 'transcript'], fdr=['10', '05']),
       expand("FastQTL/all_eqtls_{level}.all_gtex_{tissue}.txt.gz", level = ['gene'], tissue=['Brain_Frontal_Cortex_BA9']),
       expand("GTEx_Analysis_v7_eQTL/{tissue}_filtered.txt", tissue = config['gtex_samples']),

       
rule vcf_stats:
    input:
         "Genotypes/{run}/Renamed/chr{chr_num}.dose.renamed.vcf.gz"
    output:
        "Genotypes/{run}/Renamed/chr{chr_num}.dose.renamed.stats.txt"
    shell:
        "bcftools stats {input} > {output}"

rule summarise_stats:
    input:
        expand("Genotypes/{run}/Renamed/chr{chr_num}.dose.renamed.stats.txt", run=["Imputation3", "Imputation4"], chr_num=range(1,chr_num))
    output:
        "Results/snp_summary.txt"
    run:
        results = [] # list of 3-part tupples listing imputation run, chr_num and num_snps
        for filename in input:
            run = filename.split('/')[1]
            chr_num = filename.split('/')[3].split('.')[0]
            with open(filename, 'rt') as fh:
                for line in fh:
                    fields = line.split('\t')
                    if fields[0] == 'SN':
                      if fields[2] == 'number of SNPs:':
                         results.append((run, chr_num, fields[3].strip()))
                         break
                         
        with open(output[0], 'w') as out_fh:
            out_fh.write('\n'.join(['\t'.join(x) for x in results])+'\n')
 
rule rename_samples:
    """I need to run this before merging the two files because bcftools merge throws an error
    when the header does not include the following lines:
        ##FILTER=<ID=GENOTYPED,Description="Site was genotyped">
        ##FILTER=<ID=GENOTYPED_ONLY,Description="Site was genotyped only">"""
    input:
        sample_info=config["genotyping_info"],
        vcf="Genotypes/{run}/hg19/chr{chr_num}.dose.vcf.gz"
    output:
        "Genotypes/{run}/Renamed/chr{chr_num}.dose.renamed.vcf.gz"
    log:
        "Logs/Rename/chr{chr_num}{run}_rename.txt"
    run:
        from subprocess import Popen, PIPE
        import csv
        reader = csv.reader(open(input['sample_info'], 'r'), delimiter='\t')
        sample_info = dict((rows[2].encode('ascii'),rows[0].encode('ascii')) for rows in reader)
        p = Popen(['bcftools', 'view', '-h', input['vcf']], stdout=PIPE, stderr=PIPE)
        old_header, err = p.communicate()
        #old_header = old_header.replace(b"##contig=<ID=", b"##contig=<ID=chr")
        header_lines = old_header.splitlines()
        sample_names = header_lines[-1].split(b'\t')
        for i in range(9,len(sample_names)):
            old_name = b'_'.join(sample_names[i].split(b'_')[1:])
            sample_names[i] = sample_info[old_name]
        extra_lines = [b'##FILTER=<ID=GENOTYPED,Description="Site was genotyped">',
                       b'##FILTER=<ID=GENOTYPED_ONLY,Description="Site was genotyped only">']
        new_header = b'\n'.join(header_lines[:1] + extra_lines + header_lines[1:-1]+[b'\t'.join(sample_names)])
        p = Popen(['bcftools', 'reheader', '-h', '/dev/stdin', '-o', output[0], input['vcf']], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        (out, err) = p.communicate(input=new_header)
        with open(log[0], 'wb') as log:
            log.write(out)
            log.write(err)

rule index_vcf:
    input:
         "Genotypes/{run}/Renamed/chr{chr_num}.dose.renamed.vcf.gz"
    output:
        "Genotypes/{run}/Renamed/chr{chr_num}.dose.renamed.vcf.gz.csi"
    shell:
        "bcftools index {input}"

rule merge_vcf:
    input:
        file1="Genotypes/Imputation3/Renamed/chr{chr_num}.dose.renamed.vcf.gz", 
        file2="Genotypes/Imputation4/Renamed/chr{chr_num}.dose.renamed.vcf.gz",
        index1="Genotypes/Imputation3/Renamed/chr{chr_num}.dose.renamed.vcf.gz.csi", 
        index2="Genotypes/Imputation4/Renamed/chr{chr_num}.dose.renamed.vcf.gz.csi"
    output:
        "Genotypes/MergedImputations/chr{chr_num}.merged.bcf"
    log:
        "Logs/Merge/chr{chr_num}_merge.txt"
    shell:
       "(bcftools merge -Ov {input.file1} {input.file2} | bcftools filter -e 'GT =\".\"' -Ob -o {output} ) 2> {log}"

rule index_vcf2:
    input:
         "Genotypes/MergedImputations/chr{chr_num}.merged.bcf"
    output:
        "Genotypes/MergedImputations/chr{chr_num}.merged.bcf.csi"
    shell:
        "bcftools index {input}"

rule add_rsID:
    input:
        vcf="Genotypes/MergedImputations/chr{chr_num}.merged.bcf",
        index="Genotypes/MergedImputations/chr{chr_num}.merged.bcf.csi",
        rsID=config["reference"]["rsID"]
    output:
        "Genotypes/Annotated/chr{chr_num}.annotated.vcf.gz" 
    log:
        "Logs/Annotate/chr{chr_num}_annotate.txt"
    shell:
        "(bcftools annotate -c ID -Oz -a {input.rsID} -o {output} {input.vcf}) 2> {log}"

rule lift_over:
    input:
        vcf="Genotypes/Annotated/chr{chr_num}.annotated.vcf.gz", # CrossMap appears to require a vcf file, not bcf
        chain_file=config["reference"]["chain_file"],
        genome=config["reference"]["genome"]
    output:
        "Genotypes/GRCh38/chr{chr_num}.hg38.vcf"
    log:
        "Logs/Liftover/chr{chr_num}_liftover.txt"
    shell:
        "(CrossMap.py vcf {input.chain_file} {input.vcf} {input.genome} {output}) 2> {log}"

rule filter_chr:
    input:
        "Genotypes/GRCh38/chr{chr_num}.hg38.vcf"
    output:
        "Genotypes/FilterChr/chr{chr_num}.filter_chr.vcf.gz" 
    log:
        "Logs/FilterChr/chr{chr_num}_filter_chr.txt"
    params:
        chr =  "'^{chr_num}\\b'",
    shell:
        "(grep -e '^#' -e {params.chr} {input} | bcftools view -Oz -o {output}) 2> {log}"

rule sort_vcf:
    input:
        "Genotypes/FilterChr/chr{chr_num}.filter_chr.vcf.gz" # this is another tool that can't handle bcf
    output:
        "Genotypes/Sort/chr{chr_num}.sorted.vcf.gz"
    log:
        "Logs/Sort/chr{chr_num}_sort.txt"
    shell:
        "(vcf-sort {input} | bcftools view -Oz -o {output}) 2> {log}"
        
rule vcf_check:
    input:
        vcf="Genotypes/Sort/chr{chr_num}.sorted.vcf.gz", # this is another tool that can't handle bcf
        genome=config["reference"]["genome"]
    output:
        "Genotypes/Sort/chr{chr_num}.sorted.check.nonSnp",
        "Genotypes/Sort/chr{chr_num}.sorted.check.dup"
    params:
        prefix = "Genotypes/Sort/chr{chr_num}.sorted"    
    log:
        "Logs/CheckVCF/chr{chr_num}_vcf_check.txt"
    shell:
         "(Python/checkVCF.py -r {input.genome} -o {params.prefix} {input.vcf}) 2> {log}"

rule excluded_sites:
    input:
        nonSnp="Genotypes/Sort/chr{chr_num}.sorted.check.nonSnp",
        dup="Genotypes/Sort/chr{chr_num}.sorted.check.dup"
    output:
        "Genotypes/FilterDup/chr{chr_num}.excluded_sites.txt",
    shell:
        "cut -f 1,2 {input.nonSnp} > {output}; cut -f 2 {input.dup} | perl -pe 's/:/\t/' >> {output}"

rule filter_dup:
    input:
        vcf="Genotypes/FilterChr/chr{chr_num}.filter_chr.vcf.gz",
        excluded_sites="Genotypes/FilterDup/chr{chr_num}.excluded_sites.txt"
    output:
        "Genotypes/FilterDup/chr{chr_num}.filter_dup.bcf"
    log:
        "Logs/FilterDup/chr{chr_num}_filter_dup.txt"
    shell:
        "(if [ -s {input.excluded_sites} ] ; then bcftools filter -T ^{input.excluded_sites} -Ob -o {output} {input.vcf} ; else bcftools view -Ob -o {output} {input.vcf} ; fi) 2> {log}"

rule combine_chromosomes:
    input:
        expand("Genotypes/FilterDup/chr{chr_num}.filter_dup.bcf", chr_num=range(1,chr_num)) # change this to range(1,23) to run on all samples
    output:
        "Genotypes/Combined/combined.vcf.gz"
    log:
        "Logs/CombineChromosomes/combined.txt"
    shell:
        "(bcftools concat {input} -Oz -o {output}) 2> {log}"

rule sort_vcf2:
    input:
        rules.combine_chromosomes.output
    output:
        "Genotypes/Combined/combined_sort.vcf.gz"
    log:
        "Logs/CombineChromosomes/combined_sort.txt"
    shell:
        "(bcftools sort {input} -Oz -o {output}) 2> {log}"

rule index_vcf3:
    input:
         rules.sort_vcf2.output
    output:
        "Genotypes/Combined/combined_sort.vcf.gz.csi"
    shell:
        "bcftools index {input}"

rule get_gene_positions:
    input:
        gtf=config["reference"]["gtf"]
    output:
        "Data/geneloc.txt"
    params:
        level = "genloc"
    shell:
        "cat {input} | awk -v OFS=\"\t\" '{{if ($3 == \"{params.level}\") print $10, $1, $4, $5, $7}}' | sed 's/[\";]//g' > {output}"

rule get_transcript_positions:
    input:
        gtf=config["reference"]["gtf"]
    output:
        "Data/transcriptloc.txt"
    params:
        level = "transcript"
    shell:
        "cat {input} | awk -v OFS=\"\t\" '{{if ($3 == \"{params.level}\") print $12, $1, $4, $5, $7}}' | sed 's/[\";]//g' > {output}"

rule filter_counts:
    input:
        gene_counts = lambda wildcards: config["count_data"][wildcards.level],
        geneloc="Data/{level}loc.txt"
    output:
        "Data/expression_{level}.bed"
    params:
        min=5,
        num=10,
        excluded = "17046,16385,17048,16024,16115,11449,16972,13008"
    shell:
        "Rscript R/MakeBED.R --counts {input.gene_counts} --genes {input.geneloc} "
        "--min {params.min} --num {params.num} --out {output} --exclude {params.excluded}"

rule bgzip_counts:
    input:
        rules.filter_counts.output
    output:
        "Data/expression_{level}.bed.gz"
    shell:
        "bgzip {input}"

rule index_counts:
    input:
        rules.bgzip_counts.output
    output:
        "Data/expression_{level}.bed.gz.tbi"
    shell:
        "tabix -p bed {input}"

rule select_samples:
    input:
        expression = "Data/expression_gene.bed.gz",
        vcf = rules.sort_vcf2.output,
        index = rules.index_vcf3.output
    output:
        "Genotypes/Combined/combined_inc_samples.vcf.gz"
    shell:
        "bcftools view -s `zcat {input.expression} | head -1 | cut --complement -f 1-4 | perl -pe 's/\s+(?!$)/,/g'` "
        "{input.vcf} -Ou | bcftools sort -Oz -o {output} "

rule filter_tags:
    input:
        rules.select_samples.output
    output:
        "Genotypes/Combined/combined_filtered.vcf.gz"
    params:
        maf=.05,
        hwe=.0001,
        r2=.8
    shell:
        "bcftools +fill-tags {input} -Ou | bcftools view -e'MAF<{params.maf} || "
        "HWE<{params.hwe} || R2<{params.r2}' -Ou - | bcftools sort -Oz -o {output} - "

rule snp_positions:
    input:
        rules.filter_tags.output
    output:
        "Genotypes/Combined/snp_positions.txt"
    shell:
        "bcftools view -H {input} |cut -f 1,2,3,4,5 > {output}"

rule snp_positions_bed:
    input:
        rules.filter_tags.output
    output:
        "Genotypes/Combined/snp_positions.bed"
    shell:
        "bcftools view -H {input} | awk -v OFS=\"\t\" '{{print \"chr\" $1, $2-1, $2, $3, \".\", \"+\", $4, $5}}' > {output}"

rule plink_import:
    input:
        rules.filter_tags.output
    output:
        "Genotypes/Plink/genotypes.bed"
    params:
        prefix = "Genotypes/Plink/genotypes",         
    shell:
        "plink --vcf {input} --double-id --make-bed --out {params.prefix}"

rule plink_ld_prune:
    input:
        rules.plink_import.output
    output:
        "Genotypes/Plink/ld_prune.prune.in"
    params:
        input_prefix = "Genotypes/Plink/genotypes",
        output_prefix = "Genotypes/Plink/ld_prune"
    shell:
        "plink --bfile {params.input_prefix} --indep-pairwise 250 5 0.2 --out {params.output_prefix}"

rule plink_pca:
    input:
        bfile=rules.plink_import.output,
        included=rules.plink_ld_prune.output
    output:
        "Genotypes/Plink/pca.eigenvec",
    params:
        input_prefix = "Genotypes/Plink/genotypes",
        output_prefix = "Genotypes/Plink/pca",
        num_components = 3
    shell:
        "plink --bfile {params.input_prefix} --pca {params.num_components} --extract {input.included} --out {params.output_prefix}"

rule scz_ld:
    input:
        bfile=rules.plink_import.output,
        snps="Data/scz_snps.txt"
    output:
        "Genotypes/Plink/scz_ld.tags"
    params:
        input_prefix = "Genotypes/Plink/genotypes",
        output_prefix = "Genotypes/Plink/scz_ld",
    shell:
        "plink -bfile {params.input_prefix} --show-tags {input.snps} --list-all --out {params.output_prefix}"

rule tabix_vcf:
    input:
        rules.filter_tags.output
    output:
        "Genotypes/Combined/combined_filtered.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input}"

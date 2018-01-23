#snakemake --use-conda --cluster-config cluster_config.yaml --cluster "qsub -pe smp {cluster.num_cores} -l h_vmem={cluster.maxvmem}" -j 20

configfile: "config.yaml"


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
    rule prepare_smr:
    rule cat_smr:
    rule make_besd:
    rule prepare_gwas:
    rule smr:
    fast_qtl_permutations: run permutation analysis on each chunk of genome (corrected p-values for top SNPs)
    cat_permutations: concatinate fast_qtl permutation pass output
    q_values: calculate FDR for each eQTL
    rule gtex2bed:
    rule lift_over_bed:
    rule summarise_overlaps:

"""
rule all:
    input:
       "FastQTL/egenes.bed.gz",
       "FastQTL/FastQTL.all.txt.gz",
       "Peer/factors_nc.txt",
       "Genotypes/Plink/scz_ld.tags",
       "Results/GTExOverlaps.txt",
       expand("SMR/mysmr_{level}.smr", level = ['gene', 'transcript'])
       
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
        "Data/{level}loc.txt"
    params:
        level = "{level}"
    shell:
        "cat {input} | awk '{{if ($3 == \"{params.level}\") print $10, $1, $4, $5, $7}}' | sed 's/[\";]//g' > {output}"

rule filter_counts:
    input:
        gene_counts=config["count_data"]["{level}"],
        geneloc=rules.get_gene_positions.output
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

rule peer:
    input:
        pca = rules.plink_pca.output,
        counts = "Data/expression_gene.bed.gz"
    output:
        "Peer/factors.txt"
    params:
        sample_info=config["sample_info"],
        residuals = "Peer/residuals.txt",
        alpha = "Peer/alpha.txt",
        num_peer = 10
    log:
        "Logs/PEER/peer.txt"
    conda:
        "env/peer.yaml"
    shell:
        "(Rscript /c8000xd3/rnaseq-heath/GENEX-FB2/R/PEER.R -p {input.pca} "
        "-n {params.num_peer} -c {input.counts} -b {params.sample_info} "
        "-f {output} -r {params.residuals} -a {params.alpha}) > {log}"

rule format_cov:
    input:
        rules.peer.output
    output:
        "FastQTL/covariates.txt"
    run:
        import pandas as pd
        cov=pd.read_csv(input[0], sep='\t')
        cov['V1'] = 'sex' + cov['V1'].astype(str) # change sex to string (categorical)
        cov['V4'] = 'batch' + cov['V4'].astype(str) # change batch to string (categorical)
        cov.transpose().to_csv(output[0], sep='\t', header=False)
        
rule peer_nc:
    input:
        counts = "Data/expression_gene.bed.gz"
    output:
        "Peer/factors_nc.txt"
    params:
        residuals = "Peer/residuals_nc.txt",
        alpha = "Peer/alpha_nc.txt",
        num_peer = 10
    log:
        "Logs/PEER/peer_nc.txt"
    conda:
        "env/peer.yaml"
    shell:
        "(Rscript /c8000xd3/rnaseq-heath/GENEX-FB2/R/PEER.R  "
        "-n {params.num_peer} -c {input.counts} "
        "-r {params.residuals} -f {output}  -a {params.alpha}) > {log}"

rule tabix_vcf:
    input:
        rules.filter_tags.output
    output:
        "Genotypes/Combined/combined_filtered.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input}"

rule fast_qtl:
    input:
        counts = rules.bgzip_counts.output,
        count_index = rules.index_counts.output,
        genotypes = rules.filter_tags.output,
        genotype_index = rules.tabix_vcf.output,
        covariates = rules.format_cov.output
    output:
        "FastQTL/fastQTL.{level}.{chunk}.txt.gz"
    params:
        min = 1000,
        max = 10000,
        chunk = "{chunk}",
        num_chunks = 100
    log:
        "Logs/FastQTL/fastQTL_{level}_{chunk}.txt"
    shell:
        "fastQTL --vcf {input.genotypes} "
        "--bed {input.counts} --chunk {params.chunk} {params.num_chunks} "
        "--cov {input.covariates} --out {output} -- log {log} --normal"

# columns: gene_id, variant_id, tss_distance, ma_samples, ma_count, maf, pval_nominal, slope, slope_se
rule cat_fast_qtl:
    input:
        expand("FastQTL/fastQTL.{level}.{chunk}.txt.gz", chunk=range(1,num_permutations))
    output:
        "FastQTL/FastQTL_{level}.all.txt.gz"
    shell:
        "zcat {input} | gzip -c > {output}"

rule prepare_smr:
    input:
        rules.fast_qtl.output,
        rules.snp_positions.output,
        "Data/geneloc_{level}.txt"
    output:
        "SMR/myquery.{level}.{chunk}.txt.gz"
    shell:
        "Rscript R/PrepareSMR.R {input} {output}"

rule cat_smr:
    input:
        expand("SMR/myquery.{level}.{chunk}.txt.gz", chunk=range(1,num_permutations))
    output:
        "SMR/myquery_{level}.txt"
    run:
        from collections import defaultdict
        unique = defaultdict(set)
        import gzip
        with open(output[0], 'w') as smr_file:
            headers = ['SNP', 'Chr', 'BP', 'A1', 'A2', 'Freq', 'Probe', 'Probe_Chr', 'Probe_bp', 'Gene', 'Orientation', 'b', 'se', 'p']
            smr_file.write('\t'.join(headers) + '\n')
            for input_file in input:
              with gzip.open(input_file, 'rt') as f:
                for line in f.readlines():
                    fields = line.split('\t')
                    rsID = fields[0]
                    geneID = fields[6] 
                    if rsID in unique and geneID in unique[rsID]:
                        next
                    else:
                        smr_file.write(line)
                        unique[rsID].add(geneID)

rule make_besd:
    input:
        rules.cat_smr.output
    output:
        "SMR/mybesd_{level}.besd"
    params:
        "SMR/mybesd_{level}"
    shell:
        "smr --qfile {input} --make-besd --out {params}"

rule prepare_gwas:
    input:
        config['gwas']
    output:
        "SMR/CLOZUK_recoded.txt"
    shell:
        "Rscript R/PrepareGWAS.R {input} {output}"

rule smr:
    input:
        gwas = rules.prepare_gwas.output,
        besd = rules.make_besd.output
    output:
        "SMR/mysmr_{level}.smr"
    params:
        plink_prefix = "Genotypes/Plink/genotypes",
        besd_prefix = "SMR/mybesd_{level}",
        out_prefix = "SMR/mysmr_{level}"
    threads: 10    
    shell:
        "smr --bfile {params.plink_prefix} --gwas-summary {input.gwas} "
        "--beqtl-summary {params.besd_prefix} --out {params.out_prefix} --thread-num {threads}" 

rule fast_qtl_permutations:
    input:
        counts = rules.bgzip_counts.output,
        count_index = rules.index_counts.output,
        genotypes = rules.filter_tags.output,
        genotype_index = rules.tabix_vcf.output,
        covariates = rules.format_cov.output
    output:
        "FastQTL/permutations.{level}.{chunk}.txt.gz"
    params:
        min = 1000,
        max = 10000,
        chunk = "{chunk}",
        num_chunks = 100
    log:
        "Logs/FastQTL/fastQTL_{level}_{chunk}.txt"
    shell:
        "fastQTL --vcf {input.genotypes} --cov {input.covariates} --normal"
        " --bed {input.counts} --chunk {params.chunk} {params.num_chunks}"
        " --permute {params.min} {params.max} --out {output} -- log {log}"

#columns: gene_id, num_var, beta_shape1, beta_shape2, true_df, pval_true_df, variant_id, tss_distance, minor_allele_samples, minor_allele_count, maf, ref_factor, pval_nominal, slope, slope_se, pval_perm, pval_beta
rule cat_permutations:
    input:
        expand("FastQTL/permutations.{level}.{chunk}.txt.gz", chunk=range(1,num_permutations))
    output:
        "FastQTL/permutations_{level}.all.txt.gz"
    shell:
        "zcat {input} | gzip -c > {output}"

# columns: chr, start, end, gene_id, num_var, beta_shape1, beta_shape2, true_df, pval_true_df, variant_id, tss_distance, minor_allele_samples, minor_allele_count, maf, ref_factor, pval_nominal, slope, slope_se, pval_perm, pval_beta, qval, pval_nominal_threshold
rule q_values:
    input:
        eqtls=rules.cat_permutations.output,
        snp_pos=rules.snp_positions.output
    output:
        "FastQTL/egenes_{level}.bed.gz"
    params:
        fdr=.05
    log:
        "Logs/FastQTL/q_values_{level}.txt"
    shell:
        "(Rscript R/calulateNominalPvalueThresholds.R {input.eqtls} {input.snp_pos} {params.fdr} {output}) > {log}"

"""CrossMap is picky about the formatting of bed files, so I've had to discard some information
I was able to put the geneId in the score column and the nominal_p, slope, slope_se and q
into columns 7-10. Other columns can be recovered if needed by doing a join on the snp_id 
and gene_id columns
"""
rule gtex2bed:
    input:
        "GTEx_Analysis_v7_eQTL/{tissue}.v7.signif_variant_gene_pairs.txt.gz"
    output:
        "GTEx_Analysis_v7_eQTL/{tissue}_hg19.bed"
    run:
        import fileinput
        import warnings
        import gzip
        with open(output[0], 'w') as bed_file:
            with gzip.open(input[0], 'rt') as f:
                for line in f.readlines():
                    line = line.strip()
                    fields = line.split('\t')
                    try:
                        (chr, end, ref, alt, build) = fields[0].split('_')
                        bed_file.write('\t'.join(['chr' + chr, str(int(end)-1), end, fields[0], fields[1], '+']+fields[6:9]+fields[11:])+'\n')
                    except ValueError:
                        warnings.warn("skipping the following line:\n" + line)

rule lift_over_bed:
    input:
        bed=rules.gtex2bed.output,
        chain_file=config["reference"]["chain_file"],
    output:
        "GTEx_Analysis_v7_eQTL/{tissue}.bed"
    log:
        "Logs/LiftoverBED/{tissue}_liftover.txt"
    shell:
        "(CrossMap.py bed {input.chain_file} {input.bed} {output}) 2> {log}"

rule summarise_overlaps:
    input:
        sig_pairs = expand("GTEx_Analysis_v7_eQTL/{tissue}.bed", tissue = config['gtex_samples']),
        query = "FastQTL/egenes_gene.bed.gz"
    output:
        "Results/GTExOverlaps.txt"
    shell:
        "Rscript R/GetOverlaps.R --query {input.query} --outfile {output} {input.sig_pairs}"

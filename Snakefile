#snakemake --use-conda --cluster-config cluster_config.yaml --cluster "qsub -pe smp {cluster.num_cores} -l h_vmem={cluster.maxvmem}" -j 20

configfile: "config.yaml"

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
    plink_import: create (plink) bed file from vcf for PCA analysis
    plink_ld_prune: LD pruning analysis on plink-formatted genotype data
    plink_pca: PCA analysis on LD-pruned plink-formatted genotype data
    get_gene_positions: extract coordinates of genes from gtf file
    filter_counts: filter low expression genes and remove excluded samples from count matrix
    peer: PEER analysis on count data, using cofactors and PCA results
    make_bed: convert PEER residuals to BED file that can be used as input for FastQTL
    index_counts: index BED file of residualised counts
    select_samples: exclude sample from vcf and make sure order is the same as the counts file
    filter_tags: add info about MAF, HWE, etc. and filter SNPs
    rename_chromosomes: change chromosome names from 1,2,3,etc to chr1,chr2, etc. (this could be combined with rule above)
    fast_qtl: run permutation analysis on each chunk of genome (not sure if I also need to calculate nominal p-values)
    cat_eqtl: concatinate fast_qtl output
    q_values: calculate FDR for each eQTL

Todo:
- need to double-check what nonSNP files actually check for. I knew at one point, but I can't remember
- add genotyping run as cofactor?
- figure out why FastQTL keeps complaining that the tabix index files are older than the bed/vcf (I've just redone it manually for now)

Notes
- FastQTL runs on chunk1/1000 in just under 1 min (max vmem is 150MB)
- If things increase linearly, which seems reasonable, if I use 100 chunks, each will run in 90 min with 1.5 GB memory
- I'm currently limiting things to 20 simultaneous jobs (tho I could probably go higher)
- In principle, that means the job should take 5 * 1.5 = 7.5 hrs (technically 56*100*5/60/60=7.8)
"""
rule all:
    input:
       "FastQTL/results.txt"

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
        expand("Genotypes/FilterDup/chr{chr_num}.filter_dup.bcf", chr_num=range(1,23)) # change this to range(1,23) to run on all samples
    output:
        "Genotypes/Combined/combined.bcf"
    log:
        "Logs/CombineChromosomes/combined.txt"
    shell:
        "(bcftools concat -Ob -o {output} {input}) 2> {log}"

rule plink_import:
    input:
        rules.combine_chromosomes.output
    output:
        "Genotypes/Plink/genotypes.bed"
    params:
        prefix = "Genotypes/Plink/genotypes",         
    shell:
        "plink --bcf {input} --double-id --maf .05 --hwe .0001 --make-bed --out {params.prefix}"

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
        rules.plink_ld_prune.output
    output:
        "Genotypes/Plink/pca.eigenvec",
    params:
        input_prefix = "Genotypes/Plink/genotypes",
        output_prefix = "Genotypes/Plink/pca",
        num_components = 3
    shell:
        "plink --bfile {params.input_prefix} --pca {params.num_components} --extract {input} --out {params.output_prefix}"

rule get_gene_positions:
    input:
        gtf=config["reference"]["gtf"]
    output:
        "Data/geneloc.txt"
    shell:
        "cat {input} | awk '{{if ($3 == \"gene\") print $10, $1, $4, $5}}' | sed 's/[\";]//g' > {output}"

rule filter_counts:
    input:
        gene_counts=config["count_data"],
        geneloc="Data/geneloc.txt"
    output:
        "Data/expression.bed"
    params:
        min=5,
        num=10,
        excluded = "17046,16385,17048,16024,16115,11449"
    shell:
        "Rscript R/MakeBED.R --counts {input.gene_counts} --genes {input.geneloc} "
        "--min {params.min} --num {params.num} --out {output} --exclude {params.excluded}"

rule peer:
    input:
        pca=rules.plink_pca.output,
        counts = rules.filter_counts.output
    output:
        "Peer/residuals.txt"
    params:
        sample_info=config["sample_info"],
        factors = "Peer/factors.txt",
        alpha = "Peer/alpha.txt",
        num_peer = 10
    log:
        "Logs/PEER/peer.txt"
    conda:
        "env/peer.yaml"
    shell:
        "(Rscript /c8000xd3/rnaseq-heath/GENEX-FB2/R/PEER.R -p {input.pca} "
        "-n {params.num_peer} -c {input.counts} -b {params.sample_info} "
        "-r {output} -f {params.factors} -a {params.alpha}) > {log}"

rule make_bed:
    input:
        gene_counts=rules.peer.output,
        geneloc="Data/geneloc.txt"
    output:
        "Data/expression_residuals.bed"
    params:
        min=0,
        num=0
    shell:
        "Rscript R/MakeBED.R --counts {input.gene_counts} --genes {input.geneloc} "
        "--min {params.min} --num {params.num} --out {output}"

rule index_counts:
    input:
        rules.make_bed.output
    output:
        "Data/expression_residuals.bed.gz"
    shell:
        "bgzip {input} && tabix -p bed {output}" # this looks a bit awkward, but it gives the index for free

rule select_samples:
    input:
        expression=rules.make_bed.output,
        vcf=rules.combine_chromosomes.output
    output:
        "Genotypes/Combined/combined_inc_samples.bcf"
    shell:
        "bcftools view -s `head -1 {input.expression} | cut --complement -f 1-4 | perl -pe 's/\s+(?!$)/,/g'` "
        "{input.vcf} -Ob -o {output} "

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
        "&& tabix -p vcf {output}"

rule rename_chromosomes:
    input:
        rules.filter_tags.output
    output:
        "Genotypes/Combined/combined_filtered_renamed.vcf.gz"
    params:
        map = "Data/chr_map.txt"
    shell:
        "bcftools annotate --rename-chrs {params.map} -Oz -o {output} {input} && tabix -p vcf {output}"

rule fast_qtl:
    input:
        counts = rules.index_counts.output,
        genotypes = rules.rename_chromosomes.output
    output:
        "FastQTL/permutations.{chunk}.txt.gz"
    params:
        min = 1000,
        max = 10000,
        chunk = "{chunk}",
        num_chunks = 100
    log:
        "Logs/FastQTL/fastQTL_{chunk}.txt"
    shell:
        "~/src/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --vcf {input.genotypes} "
        " --bed {input.counts} --chunk {params.chunk} {params.num_chunks}"
        " --permute {params.min} {params.max} --out {output} -- log {log}"

rule cat_eqtls:
    input:
        expand("FastQTL/permutations.{chunk}.txt.gz", chunk=range(1,11))
    output:
        "FastQTL/permutations.all.txt.gz"
    shell:
        "zcat {input} | gzip -c > {output}"

rule q_values:
    input:
        rules.cat_eqtls.output
    output:
        "FastQTL/results.txt"
    params:
        fdr=.05
    shell:
        "Rscript ~/src/FastQTL-2.165.linux/scripts/calulateNominalPvalueThresholds.R "
        "{input} {params.fdr} {output}"

#snakemake --use-conda --cluster-config cluster_config.yaml --cluster "qsub -pe smp {cluster.num_cores} -l h_vmem={cluster.maxvmem}" -j 20

configfile: "config.yaml"

# this script currently adds sample names to vcf files, merges them, lifts them over to GRCh38, coordinate-sorts the output, and removes non-ref bases and duplicated sites (non-binary SNPs)
# previously I filtered SNPs with alt allele probabilities < 0.9
# there's also a bunch of stuff like removing SNPs with low MAF and testing for HWE that would probably be useful

rule all:
    input:
        "Peer/residuals.txt",
        "Genotypes/Combined/combined_inc_samples.bcf"

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

rule filter_prob:
    input:
        "Genotypes/FilterDup/chr{chr_num}.filter_dup.bcf"
    output:
        "Genotypes/FilterProb/chr{chr_num}.filter_prob.bcf"
    run:
        from pysam import VariantFile
        bcf_in = VariantFile(input[0])  # auto-detect input format
        bcf_out = VariantFile(output[0], 'wb', header=bcf_in.header)
        for site in bcf_in.fetch():
            keep_site = 0 # default option is to remove SNP
            for sample, rec in site.samples.items():
                if max(rec.get('GP')[1:]) > 0.9:
                    keep_site = 1 # do not remove SNP if either het or non-ref homo is greater than .9 for any sample
            if keep_site:
                bcf_out.write(site)

rule combine_chromosomes:
    input:
        expand("Genotypes/FilterProb/chr{chr_num}.filter_prob.bcf", chr_num=range(1,23)) # change this to range(1,23) to run on all samples
    output:
        "Genotypes/Combined/combined.bcf"
    log:
        "Logs/CombineChromosomes/combined.txt"
    shell:
        "(bcftools concat -Ob -o {output} {input}) 2> {log}"

rule fill_tags:
    input:
        rules.combine_chromosomes.output
    output:
        "Genotypes/Combined/combined_tags.bcf"
    shell:
        "bcftools +fill-tags {input} -Ob -o {output}"

rule filter_tags:
    input:
        rules.fill_tags.output
    output:
        "Genotypes/Combined/combined_filtered.bcf"
    params:
        maf=.05,
        hwe=.0001,
        r2=.6
    shell:
        "bcftools view -e'MAF<{params.maf} || HWE<{params.hwe} || R2<{params.r2}' {input} -Ob -o {output}"
        
rule plink_filter:
    input:
        "Genotypes/Combined/combined.bcf"
    output:
        "Genotypes/Plink/genotypes.bed"
    params:
        prefix = "Genotypes/Plink/genotypes",         
    shell:
        "plink --bcf {input} --double-id --maf .05 --hwe .0001 --make-bed --out {params.prefix}"

rule plink_ld_prune:
    input:
        rules.plink_filter.output
    output:
        "Genotypes/Plink/genotypes_ld_prune.prune.in"
    params:
        input_prefix = "Genotypes/Plink/genotypes",
        output_prefix = "Genotypes/Plink/genotypes_ld_prune"          
    shell:
        "plink --bfile {params.input_prefix} --indep-pairwise 250 5 0.2 --out {params.output_prefix}"

rule plink_pca:
    input:
        rules.plink_ld_prune.output
    output:
        "Genotypes/Plink/genotypes_pca.eigenvec",
    params:
        input_prefix = "Genotypes/Plink/genotypes",        
        output_prefix = "Genotypes/Plink/genotypes_pca",
        num_components = 3        
    shell:
        "plink --bfile {params.input_prefix} --pca {params.num_components} --extract {input} --out {params.output_prefix}"
 
rule plink_recode:
    input:
        "Genotypes/Combined/combined.bcf"
    output:
        "Genotypes/Plink/recoded.traw",
        "Genotypes/Plink/recoded.log"
    params:
        prefix = "Genotypes/Plink/recoded",
    shell:
        "plink --bcf {input} --double-id --maf .05 --hwe .0001 --recode A-transpose --out {params.prefix}"

rule get_gene_positions:
    input:
        gtf=config["reference"]["gtf"]
    output:
        "Data/geneloc.txt"
    shell:
        "cat {input} | awk '{{if ($3 == \"gene\") print $10, $1, $4, $5}}' | sed 's/[\";]//g' > {output}"

rule peer:
    input:
        rules.plink_pca.output
    output:
        "Peer/residuals.txt"
    params:
        sample_info=config["sample_info"],
        factors = "Peer/factors.txt",
        alpha = "Peer/alpha.txt",
        counts = "Data/counts_vst.txt",
        num_peer = 25,
        excluded = "17046,16385,17048,16024,16115"
    log:
        "Logs/PEER/peer.txt"
    conda:
        "env/peer.yaml"
    shell:
        "(Rscript /c8000xd3/rnaseq-heath/GENEX-FB2/R/PEER.R  -p {input} "
        "-n {params.num_peer} -c {params.counts} -b {params.sample_info} "
        "-e {params.excluded} -o {params.factors} -r {output} -a {params.alpha}) > {log}"

rule filter_counts:
    input:
        gene_counts=config["count_data"],
        geneloc="Data/geneloc.txt"
    output:
        "Data/expression.bed"
    params:
        min=6
    shell:
        "Rscript R/MakeBED.R --counts {input.gene_counts} --genes {input.geneloc} "
        "--min {params.min} --average --out {output}"

rule select_samples:
    input:
        expression=rules.filter_counts.output,
        vcf=rules.filter_tags.output
    output:
        "Genotypes/Combined/combined_inc_samples.bcf"
    params:
        min=6
    shell:
        "bcftools view -s `head -1 {input.expression} | cut --complement -f 1-4 | perl -pe 's/\t/,/g'` "
        "{input.vcf} -Ob -o {output} "
        
rule matrix_eqtl:
    input:
        genotypes="Genotypes/Plink/recoded.traw",
        snp_pos="Genotypes/Plink/genotypes.map",
        gene_counts=config["count_data"],
        gene_loc="Data/geneloc.txt",
        sample_info=config["sample_info"]
    output:
        image="MatrixEQTL/results.RData"
    params:
        cis="MatrixEQTL/cis_eqtl.txt",
        trans="MatrixEQTL/trans_eqtl.txt",
    log:
        "Logs/MatrixEQTL/matrix_eqtl.txt"
    shell:
        "(Rscript R/MatrixEQTL.R  --genotypes {input.genotypes} "
        "--counts {input.gene_counts} --snps {input.snp_pos} --genes {input.gene_loc} "
        "--cofactors {input.sample_info} --cis {params.cis} "
        "--trans {params.trans} --image {output.image}) 2> {log}"

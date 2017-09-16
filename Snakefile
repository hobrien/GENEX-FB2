configfile: "config.yaml"

# this script currently adds sample names to vcf files, merges them, lifts them over to GRCh38, coordinate-sorts the output, and removes non-ref bases and duplicated sites (non-binary SNPs)
# previously I filtered SNPs with alt allele probabilities < 0.9
# there's also a bunch of stuff like removing SNPs with low MAF and testing for HWE that would probably be useful

rule all:
    input:
        expand("Genotypes/FilterProb/chr{chr_num}.filter_prob.bcf", chr_num=range(21,23)) # change this to range(1,23) to run on all samples

rule rename_samples:
    """I need to run this before merging the two files because bcftools merge throws an error
    when the header does not include the following lines:
        ##FILTER=<ID=GENOTYPED,Description="Site was genotyped">
        ##FILTER=<ID=GENOTYPED_ONLY,Description="Site was genotyped only">"""
    input:
        sample_info=config["reference"]["sample_info"],
        vcf="Genotypes/{run}/hg19/chr{chr_num}.dose.vcf.gz"
    output:
        "Genotypes/{run}/Renamed/chr{chr_num}.dose.renamed.vcf.gz"
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
        p = Popen(['bcftools', 'reheader', '-Oz', '-h', '/dev/stdin', '-o', output[0], input['vcf']], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        p.communicate(input=new_header)

rule index_vcf:
    input:
         "Genotypes/{run}/Renamed/chr{chr_num}.dose.renamed.vcf.gz"
    output:
        "Genotypes/{run}/Renamed/chr{chr_num}.dose.renamed.vcf.gz.csi"
    shell:
        "bcftools index {input}"

rule merge_vcf:
   input:
       files=lambda wildcards: expand("Genotypes/{run}/Renamed/chr{chr_num}.dose.renamed.vcf.gz", chr_num=wildcards.chr_num, run=config["runs"]),
       index=lambda wildcards: expand("Genotypes/{run}/Renamed/chr{chr_num}.dose.renamed.vcf.gz.csi", chr_num=wildcards.chr_num, run=config["runs"])
   output:
       "Genotypes/MergedImputations/chr{chr_num}.merged.bcf"
   log:
       "Logs/Merge/chr{chr_num}_merge.txt"
   shell:
       "(bcftools merge -Ov {input.files} | bcftools filter -e 'GT =\".\"' -Ob -o {output} ) 2> {log}"
       
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
        chr =  "'^{chr_num}\\b'"
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

rule filter_dup:
    input:
        vcf="Genotypes/FilterChr/chr{chr_num}.filter_chr.vcf.gz",
        nonSnp="Genotypes/Sort/chr{chr_num}.sorted.check.nonSnp",
        dup="Genotypes/Sort/chr{chr_num}.sorted.check.dup"
    output:
        "Genotypes/FilterDup/chr{chr_num}.filter_dup.bcf"
    log:
        "Logs/FilterDup/chr{chr_num}_filter_dup.txt"
    run:
        from subprocess import call
        import pandas as pd
        try:
            nonSnp = pd.read_csv(input['nonSnp'], header=None, sep='\t')
        except pd.io.common.EmptyDataError:
            nonSnp = pd.DataFrame(pd.np.empty((0, 2)))
        try:
            dup_site = pd.read_csv(input['dup'], header=None, sep='\t')
            excluded_sites = '|'.join('POS='+snp for snp in set(dup_site[1].str.split(':', expand=True)[1].tolist() + nonSnp[1].tolist()))
        except pd.io.common.EmptyDataError:
            dup_site = pd.DataFrame(pd.np.empty((0, 2)))
            excluded_sites = '|'.join('POS=' + snp for snp in set(nonSnp[1].tolist()))
        with open(log[0], 'w') as log_file:
            if len(excluded_sites) > 0:
                log_file.write("filtering {} non-SNP sites and {} duplicate sites from {}\n".format(len(nonSnp), len(dup_site), input['vcf']))
                call(['bcftools', 'filter', '-Ob', '-e', excluded_sites, '-o', output[0], input['vcf']])
            else:
                log_file.write("No non SNP sites or duplicated sites in {}\n".format(input['vcf']))
                from shutil import copyfile
                call(['bcftools', 'view', '-Ob', '-o', output[0], input['vcf']])

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


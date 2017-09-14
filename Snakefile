configfile: "config.yaml"

rule all:
    input:
        expand("Genotypes/FilterDup/chr{chr_num}.dose.renamed.merged.annotated.hg38.sorted.filter_dup.vcf.gz", chr_num=range(21,23)) # change this to range(1,23) to run on all samples

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
        old_header = old_header.replace(b"##contig=<ID=", b"##contig=<ID=chr")
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
       "Genotypes/MergedImputations/chr{chr_num}.dose.renamed.merged.vcf.gz"
   log:
       "Logs/Merge/chr{chr_num}_merge.txt"
   shell:
       "(bcftools merge -Oz -o {output} {input.files}) 2> {log}"
       
rule index_vcf2:
    input:
         "Genotypes/MergedImputations/chr{chr_num}.dose.renamed.merged.vcf.gz"
    output:
        "Genotypes/MergedImputations/chr{chr_num}.dose.renamed.merged.vcf.gz.csi"
    shell:
        "bcftools index {input}"

rule add_rsID:
    input:
        vcf="Genotypes/MergedImputations/chr{chr_num}.dose.renamed.merged.vcf.gz",
        index="Genotypes/MergedImputations/chr{chr_num}.dose.renamed.merged.vcf.gz.csi",
        rsID=config["reference"]["rsID"]
    output:
        "Genotypes/Annotated/chr{chr_num}.dose.renamed.merged.annotated.vcf.gz"
    log:
        "Logs/Annotate/chr{chr_num}_annotate.txt"
    shell:
        "(bcftools annotate -c ID -Oz -a {input.rsID} -o {output} {input.vcf}) 2> {log}"

rule lift_over:
    input:
        vcf="Genotypes/Annotated/chr{chr_num}.dose.renamed.merged.annotated.vcf.gz",
        chain_file=config["reference"]["chain_file"],
        genome=config["reference"]["genome"]
    output:
        "Genotypes/GRCh38/chr{chr_num}.dose.renamed.merged.annotated.hg38.vcf"
    log:
        "Logs/Liftover/chr{chr_num}_liftover.txt"
    shell:
        "(CrossMap.py vcf {input.chain_file} {input.vcf} {input.genome} {output}) 2> {log}"

rule sort_vcf:
    input:
        "Genotypes/GRCh38/chr{chr_num}.dose.renamed.merged.annotated.hg38.vcf"
    output:
        "Genotypes/Sort/chr{chr_num}.dose.renamed.merged.annotated.hg38.sorted.vcf.gz"
    log:
        "Logs/Sort/chr{chr_num}_sort.txt"
    shell:
        "(vcf-sort {input} | bcftools view -Oz -o {output}) 2> {log}"
        
rule vcf_check:
    input:
        vcf="Genotypes/Sort/chr{chr_num}.dose.renamed.merged.annotated.hg38.sorted.vcf.gz",
        genome=config["reference"]["genome"]
    output:
        "Genotypes/Sort/chr{chr_num}.dose.renamed.merged.annotated.hg38.sorted.check.log"
    params:
        prefix = "Genotypes/Sort/chr{chr_num}.dose.renamed.merged.annotated.hg38.sorted"    
    log:
        "Logs/CheckVCF/chr{chr_num}_vcf_check.txt"
    shell:
         "(Python/checkVCF.py -r {input.genome} -o {params.prefix} {input.vcf}) 2> {log}"

rule filter_dup:
    input:
        vcf="Genotypes/Sort/chr{chr_num}.dose.renamed.merged.annotated.hg38.sorted.vcf.gz",
        nonSnp="Genotypes/Sort/chr{chr_num}.dose.renamed.merged.annotated.hg38.sorted.check.nonSnp",
        dup="Genotypes/Sort/chr{chr_num}.dose.renamed.merged.annotated.hg38.sorted.check.dup"
    output:
        "Genotypes/FilterDup/chr{chr_num}.dose.renamed.merged.annotated.hg38.sorted.filter_dup.vcf.gz"
    log:
        "Logs/FilterDup/chr{chr_num}_filter_dup.txt"
    run:
        from subprocess import Popen, PIPE
        import pandas as pd
        #input={'nonSnp': 'Genotypes/Imputation3/GRCh38_2/chr22.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.check.nonSnp', 'dup': 'Genotypes/Imputation3/GRCh38_2/chr22.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.check.dup'}
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
                log_file.write("filtering {} non-SNP sites and {} duplicate sites from {}".format(len(nonSnp), len(dup_site), input['vcf']))
                call(['bcftools', 'filter', '-Oz', '-e', excluded_sites, '-o', output[0], input['vcf']])
            else:
                log_file.write("No non SNP sites or duplicated sites in {}".format(input['vcf']))
                from shutil import copyfile
                copyfile(input['vcf'], output[0])


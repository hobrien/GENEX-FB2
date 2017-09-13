configfile: "config.yaml"

rule all:
    input:
        expand("Genotypes/GRCh38/chr{chr_num}.dose.renamed.merged.annotated.hg38.vcf", chr_num=range(21,23)) # change this to range(1,23) to run on all samples

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
    log:
        "Logs/Rename/chr{chr_num}_rename.txt"
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
        p = Popen(['bcftools', 'reheader', '-h', '/dev/stdin', '-o', output[0], input['vcf']], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        p.communicate(input=new_header)


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
         "Genotypes/MergedImputations/chr{chr_num}.dose.renamed.merged.vcf.gz",
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


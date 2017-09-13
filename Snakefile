configfile: "config.yaml"

rule all:
    input:
        expand("Genotypes/MergedImputations/chr{chr_num}.dose.renamed.merged.vcf.gz", chr_num=range(21,23)) # change this to range(1,23) to run on all samples

rule rename_samples:
    """I need to run this before merging the two files because bcftools merge throws an error
    when the header does not include the following lines:
        ##FILTER=<ID=GENOTYPED,Description="Site was genotyped">
        ##FILTER=<ID=GENOTYPED_ONLY,Description="Site was genotyped only">"""
    input:
        sample_info="Data/SampleInfo.txt",
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

rule index_vcf:
    input:
         "Genotypes/{run}/Renamed/chr{chr_num}.dose.renamed.vcf.gz"
    output:
        "Genotypes/{run}/Renamed/chr{chr_num}.dose.renamed.vcf.gz.csi"
    shell:
        "bcftools index {input}"

rule merge_vcf:
   """This is a bit fragile because it will break if there are more (or less) than 2 runs in the config files
      I'm also not 100% sure that the vcf files will always be the first and third argument
      I don't know of a better way to require the index files to be present tho"""
   input:
       files=lambda wildcards: expand("Genotypes/{run}/Renamed/chr{chr_num}.dose.renamed.vcf.gz", chr_num=wildcards.chr_num, run=config["runs"]),
       index=lambda wildcards: expand("Genotypes/{run}/Renamed/chr{chr_num}.dose.renamed.vcf.gz.csi", chr_num=wildcards.chr_num, run=config["runs"])
   output:
       "Genotypes/MergedImputations/chr{chr_num}.dose.renamed.merged.vcf.gz"
   log:
       "Logs/Merge/chr{chr_num}_merge.txt"
   shell:
       "(bcftools merge -o {output} {input.files}) 2> {log}"


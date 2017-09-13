configfile: "config.yaml"

rule all:
    input:
        expand("Genotypes/Imputation4/chr{chr_num}.merged.vcf", chr_num=range(1,22))
        
rule rename_samples:
    input:
        sample_info="Data/SampleInfo.txt"
        vcf="Genotypes/Imputation4/chr{chr_num}.dose.vcf.gz"
    output:
        "Genotypes/Imputation4/Renamed/chr{chr_num}.dose.renamed.vcf.gz"
    log:
        "logs/blast/{sample}.log"
    run:
        from subprocess import Popen, PIPE
        import csv
        reader = csv.reader(open(input['sample_info'], 'r'), delimiter='\t')
        sample_info = dict((rows[2].encode('ascii'),rows[0].encode('ascii')) for rows in reader)
        p = Popen(['bcftools', 'view', '-h', input['vcf']], stdout=PIPE, stderr=PIPE)
        old_header, err = p.communicate()
        header_lines = old_header.splitlines()
        sample_names = header_lines[-1].split(b'\t')
        for i in range(9,len(sample_names)):
            old_name = b'_'.join(sample_names[i].split(b'_')[1:])
            sample_names[i] = sample_info[old_name]
       new_header = b'\n'.join(header_lines[:-1]+[b'\t'.join(sample_names)])
       p = Popen(['bcftools', 'reheader', '-h', '/dev/stdin', '-o', output[0], input['vcf']], stdin=PIPE, stdout=PIPE, stderr=PIPE)
       p.communicate(input=new_header)


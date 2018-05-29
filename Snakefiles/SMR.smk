import yaml
import os
configfile: "config.yaml"
SMR=yaml.load(open('smr.yaml', 'r'))

rule all:
    input: 
       expand("SMR/mysmr_{level}_{gwas}_all.smr.gz", level = ['gene', 'transcript'], gwas=SMR['GWAS']),
       expand("SMR/plot/{gwas}_{level}.{gene_id}.txt", gwas=['clozuk'], level=['gene'], gene_id=SMR['gene']['clozuk']),

rule prepare_smr:
    input:
        "FastQTL/all_eqtls_{level}.{chunk}_q05.txt.gz",
        "Genotypes/Combined/snp_positions.txt",
        "Data/{level}loc.txt"
    output:
        "SMR/myquery.{level}.{chunk}.txt.gz"
    shell:
        "Rscript R/PrepareSMR.R {input} {output}"

rule cat_eqtls_tr:
    input:
        expand("SMR/myquery.transcript.{chunk}.txt.gz", chunk=range(1,101))
    output:
        ["SMR/myquery.transcript.chr" + str(x+1) + ".txt" for x in range(23)]
    shell:
        "echo -e 'SNP\tChr\tBP\tA1\tA2\tFreq\tProbe\tProbe_Chr\tProbe_bp\tGene\tOrientation\tb\tse\tp' | tee {output}; "
        "zcat {input} | awk '{{print>>\"SMR/myquery.transcript.chr\"$2\".txt\"}}'"

rule cat_eqtls_gene:
    input:
        expand("SMR/myquery.gene.{chunk}.txt.gz", chunk=range(1,101))
    output:
        ["SMR/myquery.gene.chr" + str(x+1) + ".txt" for x in range(23)] 
    shell:
        "echo -e 'SNP\tChr\tBP\tA1\tA2\tFreq\tProbe\tProbe_Chr\tProbe_bp\tGene\tOrientation\tb\tse\tp' | tee {output}; "
        "zcat {input} | awk '{{print>>\"SMR/myquery.gene.chr\"$2\".txt\"}}'"

rule make_besd:
    input:
        "SMR/myquery.{level}.chr{chr}.txt"
    output:
        "SMR/mybesd.{level}.chr{chr}.besd"
    params:
        "SMR/mybesd.{level}.chr{chr}"
    shell:
        "smr --qfile {input} --make-besd --out {params}"

rule plink_freqx:
    input:
        "LDSR/Reference/1000G_plinkfiles/1000G.mac5eur.{chr_num}.bed"
    params:
        prefix = "LDSR/Reference/1000G_plinkfiles/1000G.mac5eur.{chr_num}"
    output:
        "LDSR/Reference/1000G_plinkfiles/1000G.mac5eur.{chr_num}.frqx"
    shell:
        "plink -bfile {params.prefix} -freqx -out {params.prefix}"
        
rule cat_freq:
    input:
        expand("LDSR/Reference/1000G_plinkfiles/1000G.mac5eur.{chr_num}.frqx", chr_num=range(1,23))
    output:
        "LDSR/Reference/1000G_plinkfiles/1000G_AF.txt"
    shell:
        "tail -n+2 {input} | grep -v '==>' > {output}"
        
rule format_gwas:
    input:
        gwas = lambda wildcards: SMR['GWAS'][wildcards.gwas]['raw'],
        af = "LDSR/Reference/1000G_plinkfiles/1000G_AF.txt"
    output:
        "GWAS/{gwas}.smr"
    params:
        n = lambda wildcards: SMR['GWAS'][wildcards.gwas]['n'],
        cols = lambda wildcards: SMR['GWAS'][wildcards.gwas]['cols'],
        effect = lambda wildcards: SMR['GWAS'][wildcards.gwas]['effect']
    run:
        class AlleleFreq():
            """a simple container to store allele frequencies
            """
            def __init__(self, fields):
                (A1, A2, homo1_count, het_count, homo2_count) = fields
                self.AF ={}
                self.AF[A1] = (float(homo1_count) + float(het_count)*.5)/(float(homo1_count) + float(het_count) + float(homo2_count))
                self.AF[A2] = 1.0 - self.AF[A1]
            def af(self, allele):
                try:
                    return self.AF[allele]
                except KeyError:
                    raise KeyError("AlleleMismatch")
        
        assert(AlleleFreq(['A', 'T', 10, 10, 20]).af('A') == 0.375)
        assert(AlleleFreq(['A', 'T', 10, 10, 20]).af('T') == 0.625)
        try:
            AlleleFreq(['A', 'T', 10, 10, 20]).af('G')
        except KeyError as e:
            assert(str(e) == "'AlleleMismatch'")
        
        allele_frequencies = {}
        with open(input['af'], 'r') as af_fh:
            for line in af_fh:
                fields = line.split()
                try:
                    allele_frequencies[fields[1]] = AlleleFreq(fields[2:7])
                except ValueError:
                    warnings.warn("cannot get allele frequency for %s" % ' '.join(fields[2:7]))
                    
        with open(output[0], 'w') as output_fh:
                output_fh.write('\t'.join(('SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'n')) + '\n')
                header = 1
                with open(input['gwas'], 'r') as gwas_fh:
                    SNPs = set()
                    for line in gwas_fh:
                        line = line.strip()
                        if header:
                            header = 0
                        else:
                            fields = line.split()
                            SNP = fields[params[cols[SNP]]].split(':')[0] # remove info after rsID                            
                            if SNP[:3] != 'rs':
                                continue  # remove snps without rsIDs
                            if SNP in SNPs:
                                continue # remove duplicate SNPs
                            SNPs.add(SNP)
                            A1 = fields[params[cols[A1]]]
                            A2 = fields[params[cols[A2]]]
                            try:
                                freq = fields[params[cols[freq]]]
                            except IndexError:
                              try:
                                freq = str(allele_frequencies[fields[0]].af(fields[1]))
                              except KeyError as e:
                                if str(e) == "'AlleleMismatch'":
                                    warnings.warn("Allele mismatch for %s, skipping" % fields[0])
                                #else:
                                #    warnings.warn("Allele Freq missing for %s" % fields[0])
                                continue
                            if params['effect'] == 'or':
                                b = str(math.log(float(fields[params[cols[b]]])))
                            elif params['effect'] == '-beta':
                                b = str(-1*float(fields[params[cols[b]]]))
                            else:
                                b = fields[params[cols[b]]]                           
                            se = fields[params[cols[se]]] 
                            p = fields[params[cols[p]]] 
                            try:
                                n = fields[params[cols[n]]]
                            except IndexError:
                                try:
                                    n = str(params['n'])
                                except index_error:
                                    continue
                            output_fh.write('\t'.join(SNP, A1, A2, freq, b, se, p, n) + '\n')

rule smr:
    input:
        gwas = rules.format_gwas.output,
        besd = rules.make_besd.output
    output:
        "SMR/mysmr_{level}_{gwas}.chr{chr}.smr"
    params:
        plink_prefix = "Genotypes/Plink/genotypes",
        besd_prefix = "SMR/mybesd.{level}.chr{chr}",
        out_prefix = "SMR/mysmr_{level}_{gwas}.chr{chr}"
    threads: 1
    shell:
        "smr --bfile {params.plink_prefix} --gwas-summary {input.gwas} "
        "--beqtl-summary {params.besd_prefix} --out {params.out_prefix} "
        "--thread-num {threads} --peqtl-smr 0.01" 

rule cat_smr:
    input:
        lambda wildcards: expand("SMR/mysmr_{level}_{gwas}.chr{chr}.smr", level = wildcards.level, gwas=wildcards.gwas, chr=range(1,23))
    output:
        "SMR/mysmr_{level}_{gwas}_all.smr.gz"
    shell:
        "tail -n+2 {input} | grep -v '==>' | gzip -c > {output}"

rule genloc_smr:
    input:
        "Data/{level}loc.txt"
    output:
        "Data/{level}loc_smr.txt"
    shell:
        "awk -v OFS=\"\t\" '{{sub(/\.[0-9]*/, \"\", $1); sub(/chr/, \"\", $2); print $2, $3, $4, $1, $5}}' {input} > {output}"

rule plot_smr:
    input:
        gwas = rules.format_gwas.output,
        genloc = rules.genloc_smr.output,
        besd = lambda wildcards: expand("SMR/mybesd.{level}.chr{chr_num}.besd", level=wildcards.level, chr_num=SMR[wildcards.level][wildcards.gwas][wildcards.gene_id])
    output:
        "SMR/plot/{gwas}_{level}.{gene_id}.txt"
    params:
        plink_prefix = "Genotypes/Plink/genotypes",
        besd_prefix = lambda wildcards: expand("SMR/mybesd.{level}.chr{chr_num}", level=wildcards.level, chr_num=SMR[wildcards.level][wildcards.gwas][wildcards.gene_id]),
        out_prefix = "SMR/{gwas}_{level}",
        genloc = "Data/{level}loc_smr.txt",
        gene_id = "{gene_id}"
    shell:
        "smr --bfile {params.plink_prefix} --gwas-summary {input.gwas} "
        "--beqtl-summary {params.besd_prefix} --out {params.out_prefix} --plot "
        "--probe {params.gene_id} --probe-wind 500 --gene-list {params.genloc}"
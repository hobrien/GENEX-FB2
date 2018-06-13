import yaml
import os

SMR=yaml.load(open('smr.yaml', 'r'))

# to make DAG: snakemake -np --dag | dot -Tsvg > dag.svg
dag = 0  # set to 1 to produce DAG image without tons of duplicates 
if dag:
    gwas_list = ['clozuk']
    gene_ids = ['ENSG00000150967']
    chr_num = 2
    num_permutations = 2
    fdr_levels = ['05']
    expression_levels = ['gene']
else:
    gwas_list = SMR['GWAS']
    gene_ids = SMR['gene']['clozuk']
    chr_num = 23
    num_permutations = 101
    fdr_levels = ['100', '10', '05', '01', '001', '0001']
    expression_levels = ['gene', 'transcript']

"""
rules:
- prepare eQTL data for SMR:
    - prepare_smr: add allele info and gene position info to eQTL results
    - cat_eqtls_tr: combine formatted transcript-level eQTL data by chromosome
    - cat_eqtls_gene: combine formatted gene-level eQTL data by chromosome
    - make_besd: convert eQTL data to SMR input format
    
- prepare GWAS summary stats:
    - plink_freqx: extract allele frequency info from 1000 genomes project using plink
    - cat_freq: combine allele frequency from each chromosome
    - format_gwas: add allele freq to GWAS sumstats; transform betas
- SMR:
    - smr: run SMR on each chromosome for each GWAS
    - cat_smr: combine results from all chromosomes
    - genloc_smr: modify info about gene positions to format needed by SMR
    - plot_smr: create data to make effect size and locus plots for each SMR hit
        - info about SMR hits was added to smr.yaml after running previous steps
        - creating rules that depend on the results of the previous rules is challenging
        with Snakemake. My understanding is this is an area where Nextflow might be better
"""

rule all:
    input: 
       expand("SMR/mysmr_{level}_{gwas}_all.smr.gz", level = expression_levels, gwas=gwas_list),
       expand("SMR/plot/{gwas}_{level}.{gene_id}.txt", gwas=['clozuk'], level=['gene'], gene_id=gene_ids),


################################ prepare eQTL data for SMR ###############################
rule prepare_smr:
    input:
        "FastQTL/all_eqtls_{level}.{chunk}_q05.txt.gz",
        "Genotypes/Combined/snp_positions.txt",
        "Data/{level}loc.txt"
    output:
        "BESD/myquery.{level}.{chunk}.txt.gz"
    shell:
        "Rscript R/PrepareSMR.R {input} {output}"

rule cat_eqtls_tr:
    input:
        expand("BESD/myquery.transcript.{chunk}.txt.gz", chunk=range(1,num_permutations))
    output:
        ["BESD/myquery.transcript.chr" + str(x+1) + ".txt" for x in range(chr_num)]
    shell:
        "echo -e 'SNP\tChr\tBP\tA1\tA2\tFreq\tProbe\tProbe_Chr\tProbe_bp\tGene\tOrientation\tb\tse\tp' | tee {output}; "
        "zcat {input} | awk '{{print>>\"BESD/myquery.transcript.chr\"$2\".txt\"}}'"

rule cat_eqtls_gene:
    input:
        expand("BESD/myquery.gene.{chunk}.txt.gz", chunk=range(1,num_permutations))
    output:
        ["BESD/myquery.gene.chr" + str(x+1) + ".txt" for x in range(chr_num)] 
    shell:
        "echo -e 'SNP\tChr\tBP\tA1\tA2\tFreq\tProbe\tProbe_Chr\tProbe_bp\tGene\tOrientation\tb\tse\tp' | tee {output}; "
        "zcat {input} | awk '{{print>>\"BESD/myquery.gene.chr\"$2\".txt\"}}'"

rule make_besd:
    input:
        "BESD/myquery.{level}.chr{chr}.txt"
    output:
        "BESD/mybesd.{level}.chr{chr}.besd"
    params:
        "BESD/mybesd.{level}.chr{chr}"
    shell:
        "smr --qfile {input} --make-besd --out {params}"


############################### prepare GWAS summary stats ###############################
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
        expand("LDSR/Reference/1000G_plinkfiles/1000G.mac5eur.{chr_num}.frqx", chr_num=range(1,chr_num))
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
        import warnings
        import gzip
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
        
        def open_file(file_name):
            if file_name[-2:] == 'gz':
                return gzip.open(file_name, 'r')
            else:
                return open(file_name, 'r')
            
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
                try:
                    gwas_fh = open_file(input['gwas'][0])
                    SNPs = set()
                    for line in gwas_fh:
                        line = line.strip()
                        if header:
                            header = 0
                        else:
                            fields = line.split()
                            SNP = fields[params['cols']['SNP']].split(':')[0] # remove info after rsID                            
                            if SNP[:3] != 'rs':
                                continue  # remove snps without rsIDs
                            if SNP in SNPs:
                                continue # remove duplicate SNPs
                            SNPs.add(SNP)
                            A1 = fields[params['cols']['A1']]
                            A2 = fields[params['cols']['A2']]
                            try:
                                freq = fields[params['cols']['freq']]
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
                                b = str(math.log(float(fields[params['cols']['b']])))
                            elif params['effect'] == '-beta':
                                b = str(-1*float(fields[params['cols']['b']]))
                            else:
                                b = fields[params['cols']['b']]                           
                            se = fields[params['cols']['se']]
                            p = fields[params['cols']['p']]
                            try:
                                n = fields[params['cols']['n']]
                            except IndexError:
                                try:
                                    n = str(params['n'])
                                except index_error:
                                    continue
                            output_fh.write('\t'.join(SNP, A1, A2, freq, b, se, p, n) + '\n')
                finally:
                    close(gwas_fh)

########################################## SMR ##########################################
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
        lambda wildcards: expand("SMR/mysmr_{level}_{gwas}.chr{chr}.smr", level = wildcards.level, gwas=wildcards.gwas, chr=range(1,chr_num))
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
    output:
        "SMR/plot/{gwas}_{level}.{gene_id}.txt"
    params:
        plink_prefix = "Genotypes/Plink/genotypes",
        besd_prefix = lambda wildcards: expand("BESD/mybesd.{level}.chr{chr_num}", level=wildcards.level, chr_num=SMR[wildcards.level][wildcards.gwas][wildcards.gene_id]),
        out_prefix = "SMR/plot/{gwas}_{level}",
        genloc = "Data/{level}loc_smr.txt",
        gene_id = "{gene_id}"
    shell:
        "smr --bfile {params.plink_prefix} --gwas-summary {input.gwas} "
        "--beqtl-summary {params.besd_prefix} --out {params.out_prefix} --plot "
        "--probe {params.gene_id} --probe-wind 500 --gene-list {params.genloc}"
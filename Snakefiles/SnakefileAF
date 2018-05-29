"""
We need the allele frequencies from the 1000 genomes EUR population to run SMR on GWAS
where this is missing. As part of LDSC, I downloaded Plink file with names like '1000G.mac5eur.22.bed'
My guess, backed up somewhat by discussions like this (https://groups.google.com/forum/#!topic/ldsc_users/q367fQIOVwA)
is that these are plink files for 1000 genomes EUR samples with SNPs with minor allele counts
< 5 removed. Since there are 750 chromosomes in these files, this should work for any GWAS
that filtered on MAF of 1% or greater.

This command:

plink -bfile 1000G.mac5eur.22 -freq -out 1000G.mac5eur.22

gives MAF for each SNP, but SMR requires AF for the effect allele.
The plink manual says that A1 is *usually* the minor allele, but that's not a guarantee.

This command:

plink -bfile 1000G.mac5eur.22 -freqx -out 1000G.mac5eur.22 

gives counts for A1 and A2 homozygotes, which can be used to determine who is the minor allele.

This gives me two options:
    - I can use this command:
    
    awk '{if ($5 > $7) print $0}' 1000G.mac5eur.22.frqx
    
    to confirm that A1 is the minor allele, then run the command above 
    (assuming that the allele order will be the same) or
    - I can calculate the AF of A1 using ($5 + $6/2)/($5+$6+$7)
- I guess the second option is ultimately easier/safer.

- If I do this in R, I can then simply join on rsID, then use if else to determine if A1
is the effect allele (ifelse(A1_1000G==A1, AF, 1-AF)).

I also need to add N and do some filtering, possibly convert OR to beta
"""

import warnings
import math

configfile: "gwas_config.yaml"

rule all:
    input:
        expand("GWAS/{gwas}.smr", gwas=config['GWAS'].keys())
        
rule cat_freq:
    input:
        expand("LDSR/Reference/1000G_plinkfiles/1000G.mac5eur.{chr_num}.frqx", chr_num=range(1,23))
    output:
        "LDSR/Reference/1000G_plinkfiles/1000G_AF.txt"
    shell:
        "tail -n+2 {input} | grep -v '==>' > {output}"
        
rule plink_freqx:
    input:
        "LDSR/Reference/1000G_plinkfiles/1000G.mac5eur.{chr_num}.bed"
    params:
        prefix = "LDSR/Reference/1000G_plinkfiles/1000G.mac5eur.{chr_num}"
    output:
        "LDSR/Reference/1000G_plinkfiles/1000G.mac5eur.{chr_num}.frqx"
    shell:
        "plink -bfile {params.prefix} -freqx -out {params.prefix}"
        
rule add_af:
    input:
        gwas = lambda wildcards: config['GWAS'][wildcards.gwas]['raw'],
        af = "LDSR/Reference/1000G_plinkfiles/1000G_AF.txt"
    output:
        "GWAS/{gwas}.smr"
    params:
        n = lambda wildcards: config['GWAS'][wildcards.gwas]['n'],
        extra_cols = lambda wildcards: config['GWAS'][wildcards.gwas]['extra_cols'],
        effect = lambda wildcards: config['GWAS'][wildcards.gwas]['effect']
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
            for line in af_fh.readlines():
                fields = line.split()
                try:
                    allele_frequencies[fields[1]] = AlleleFreq(fields[2:7])
                except ValueError:
                    warnings.warn("cannot get allele frequency for %s" % ' '.join(fields[2:7]))
                    
        with open(output[0], 'w') as output_fh:
                output_fh.write('\t'.join(('SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'n')) + '\n')
                header = 1
                with open(input['gwas'], 'r') as gwas_fh:
                    for line in gwas_fh.readlines():
                        line=line.strip()
                        if header:
                            header = 0
                        else:
                            fields = line.split()
                            for i in sorted(params['extra_cols'], reverse=True):
                                 del fields[i]
                            try:
                                fields.insert(3, str(allele_frequencies[fields[0]].af(fields[1])))
                            except KeyError as e:
                                if str(e) == "'AlleleMismatch'":
                                    warnings.warn("Allele mismatch for %s, skipping" % fields[0])
                                #else:
                                #    warnings.warn("Allele Freq missing for %s" % fields[0])
                                continue
                            fields.append(str(params['n']))
                            if params['effect'] == 'or':
                                fields[4] = str(math.log(float(fields[4])))
                            elif params['effect'] == '-beta':
                                fields[4] = str(-1*float(fields[4]))                           
                            output_fh.write('\t'.join(fields) + '\n')

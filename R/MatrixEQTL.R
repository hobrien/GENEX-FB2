setwd("~/BTSync/FetalRNAseq/Github/GENEX-FB2/")

library(MatrixEQTL)
library(tidyverse)
library(stringr)
source(paste0('R/AnalyseDE.R'))

# Import normalised counts, filter, remove 'norm.', sort
counts <- read_delim("../Data/MalevsFemale.complete.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  filter(! is.na(padj)) %>%
  dplyr::select(Id, starts_with('norm'))
colnames(counts) <- str_replace(colnames(counts), 'norm\\.', '')
counts <- counts[ , order(names(counts))]
counts <- dplyr::select(counts, id=Id, everything())

gene_location_file_name = "../Data/geneloc.txt"
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos <- genepos %>% mutate(geneid= str_replace(geneid, '(ENSG\\d+)\\.\\d+', '\\1')) %>% rename(id=geneid)
genepos <- genepos %>% semi_join(counts, by=c("id"))

# import co-variates
# Seems that data are coerced into float/integer column-wiseduring import with 
# the LoadFile command even though the data are encoded row-wise. I think the 
# only solution to this is to round every thing to whole numbers.
target <- read_delim("../Data/SampleInfo.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(Sample=as.character(Sample))
target$Sex<-as.integer(factor(target$Sex))
target$Batch<-as.integer(factor(target$ReadLength))
target$RIN<- round(target$RIN, digits=0)
target <- arrange(target, Sample)
target <- as.data.frame(target)
rownames(target) <- target$Sample
target <- dplyr::select(target, -Sequencer, -ReadLength, -Sample)
target<- as.data.frame(t(target))
target$id=rownames(target)
target <- target[ , order(names(target))]
target <- dplyr::select(target, id, everything())

# convert vcf file to table of genotypes and list of SNP locations
#individual_vcf_files = "/Volumes/FetalRNAseq/Processed/chr*.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.gz"
#vcf_file_name = "/Volumes/FetalRNAseq/Processed/chr_all.dose.rename.filter_samples.filter_sites.rsID.recoded.GRCh38.sort.vcf.gz"
#python_script_name = "~/BTsync/FetalRNAseq/LabNotes/Python/VCF2numeric.py"
snps_location_file_name = "../MatrixEQTL/snp_pos.txt"
snps_file_name = "../MatrixEQTL/genotypes.txt"
alt_threshold = 10
# use file tests to avoid re-doing unnecessary shell commands, which are pretty time consuming
#if (! file_test("-f", vcf_file_name)) {
#  system(paste('bcftools concat', individual_vcf_files, '>', vcf_file_name))
#}

#if (! file_test("-f", snps_location_file_name) | ! file_test("-f", snps_file_name)) {
#  system(paste('bcftools view', vcf_file_name, '| python', python_script_name, snps_file_name, snps_location_file_name, alt_threshold))
#}

snp_tbl<- read_delim(snps_file_name, delim='\t')
#snp_tbl <- dplyr::rename(snp_tbl,  `18208` = `18121`)
snp_tbl <- snp_tbl[ , order(names(snp_tbl))]
snp_tbl <- rename(snp_tbl, id=IID)
snp_tbl <- dplyr::select(snp_tbl, id, everything())

# remove samples that are missing in one or more of the files
counts <- counts %>% dplyr::select(-one_of(setdiff(colnames(counts), colnames(target))))
counts <- counts %>% dplyr::select(-one_of(setdiff(colnames(counts), colnames(snp_tbl))))
target <- target %>% dplyr::select(-one_of(setdiff(colnames(target), colnames(counts))))
snp_tbl <- snp_tbl %>% dplyr::select(-one_of(setdiff(colnames(snp_tbl), colnames(counts))))

# save modified files and read in as matrix EQTL objects
counts_file_name = tempfile()
write_tsv(counts, counts_file_name)

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(counts_file_name);


covariates_file_name = tempfile()
write_tsv(target, covariates_file_name)

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

geneotype_file_name = "~/BTSync/FetalRNAseq/Github/GENEX-FB2/MatrixEQTL/genotypes_formatted.txt"
if (! file_test("-f", geneotype_file_name)) {
  write_tsv(snp_tbl, geneotype_file_name)
}

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile( geneotype_file_name );

#extract gene positions from GTF file
#if (! file_test("-f", gene_location_file_name)) {
#  system(paste('echo "geneid\tchr\ts1\ts2" >', gene_location_file_name))
#}
#system(paste("cat ~/BTSync/FetalRNAseq/Ref/genes.gtf | awk '{if ($3 == \"gene\") print $10, $1, $4, $5}' | sed 's/[\";]//g' >>", gene_location_file_name))
#gene_location_file_name = "~/BTSync/FetalRNAseq/Github/GENEX-FB2/MatrixEQTL/geneloc.txt"

snpspos = read_tsv(snps_location_file_name, col_names = c("snp",	"chr",	"pos"));



###############################################################
output_file_name_cis = "../MatrixEQTL/cis_eqtl.txt"
output_file_name_tra = "../MatrixEQTL/trans_eqtl.txt"

pvOutputThreshold_cis = 1e-5;
pvOutputThreshold_tra = 1e-10;

cisDist = 1e6;

useModel = modelLINEAR

errorCovariance = numeric()

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

write_tsv(me$cis$eqtls, "MatrixEQTL/cis_eqtl.txt")
write_tsv(me$trans$eqtls, "MatrixEQTL/trans_eqtl.txt")
write_tsv(filter(me$trans$eqtls, FDR <=.1e-10), "MatrixEQTL/trans_eqtl_e10.txt")

gene_info <- read_tsv("Data/genes.txt") %>%
  mutate(gene_id = sub("\\.[0-9]+", "", gene_id)) %>%
  dplyr::select(Id = gene_id, SYMBOL=gene_name, Chr=seqid)
genotypes <- read_delim("MatrixEQTL/genotypes_formatted.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
genotypes_filtered <- semi_join(genotypes, cis, by=c("id" = 'snps'))

cis <- me$cis$eqtls %>% 
   mutate(gene = sub("(ENSG[0-9]+)\\.[0-9]+", '\\1', gene),
       statistic = as.numeric(format(statistic, digits=2)), 
       pvalue = as.numeric(format(pvalue, digits=2)), 
       FDR = as.numeric(format(FDR, digits=2)),
       beta = as.numeric(format(beta, digits=2))) %>% 
  dplyr::rename(Id=gene, padj=FDR)

cis <- right_join(gene_info, cis)
cis <- data.frame(snps=genotypes_filtered$id, num_ref=rowSums(snps==0)) %>% inner_join(cis)

write_tsv(cis, "Shiny/GENEX-FB2/Data/cis_eqtl.txt")

trans <- me$trans$eqtls %>%
  mutate(gene = sub("(ENSG[0-9]+)\\.[0-9]+", '\\1', gene),
         statistic = as.numeric(format(statistic, digits=2)), 
         pvalue = as.numeric(format(pvalue, digits=2)), 
         FDR = as.numeric(format(FDR, digits=2)),
         beta = as.numeric(format(beta, digits=2))) %>%
  dplyr::rename(Id=gene, padj=FDR)

trans <- dplyr::select(cis, snps, cisId=Id, cisSYMBOL=SYMBOL) %>% left_join(trans)
trans <- right_join(gene_info, trans)
trans <- data.frame(snps=genotypes_filtered$id, num_ref=rowSums(snps==0)) %>% inner_join(trans)

write_tsv(trans, "Data/trans_eqtl.txt")

save.image(file="MatrixEQTL/results.RData")

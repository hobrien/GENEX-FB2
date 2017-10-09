#setwd("~/BTSync/FetalRNAseq/Github/GENEX-FB2/")

library(MatrixEQTL)
library(tidyverse)
library(stringr)
library("optparse")

option_list <- list(
  make_option(c("--genotypes"), type="character", default='../MatrixEQTL/genotypes.txt',
              help="Additive genotypes file from PLINK (recodeA)"),
  make_option(c("--counts"), type="character", default="../Data/MalevsFemale.complete.txt",
              help="Counts file from SARtools"),
  make_option(c("--snps"), type="character", default="../Genotypes/Plink/genotypes.map",
              help="SNP positions (plink map file)"),
  make_option(c("--genes"), type="character", default="../Data/geneloc.txt",
              help="Positions of genes (geneID, chr, start, end)"),
  make_option(c("--cofactors"), type="character", default="../Data/SampleInfo.txt",
              help="Information about each sample"),
  make_option(c("--cis"), type="character", default="../MatrixEQTL/cis_eqtl.txt",
              help="table of cis eQTLs"),
  make_option(c("--trans"), type="character", default="../MatrixEQTL/trans_eqtl.txt",
              help="table of transeQTLs"),
  make_option(c("--image"), type="character", default="../MatrixEQTL/results.RData",
              help="R image file with all outfiles"),
  make_option(c("--p_trans"), type="numeric", default=1e-4,
              help="R image file with all outfiles"),
  make_option(c("--p_cis"), type="numeric", default=1e-8,
              help="R image file with all outfiles")
 )
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


# Import normalised counts, filter, remove 'norm.', sort
print(paste("reading", opt$counts))
counts <- read_delim(opt$counts, "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  filter(! is.na(padj)) %>%
  dplyr::select(Id, starts_with('norm'))
colnames(counts) <- str_replace(colnames(counts), 'norm\\.', '')
counts <- counts[ , order(names(counts))]
counts <- dplyr::select(counts, id=Id, everything())

print(paste("reading", opt$genes))
genepos = read.table(opt$genes, header = TRUE, stringsAsFactors = FALSE);
genepos <- genepos %>% mutate(geneid= str_replace(geneid, '(ENSG\\d+)\\.\\d+', '\\1')) %>% rename(id=geneid)
genepos <- genepos %>% semi_join(counts, by=c("id"))

# import co-variates
# Seems that data are coerced into float/integer column-wiseduring import with 
# the LoadFile command even though the data are encoded row-wise. I think the 
# only solution to this is to round every thing to whole numbers.
print(paste("reading", opt$cofactors))
target <- read_delim(opt$cofactors, "\t", escape_double = FALSE, trim_ws = TRUE) %>%
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

print(paste("reading", opt$genotypes))
snp_tbl<- read_delim(opt$genotypes, delim='\t')
snp_tbl <- snp_tbl %>% dplyr::select(-CHR, -`(C)M`, -POS, -COUNTED, -ALT) %>%
  gather(key, value, -SNP) %>% 
  mutate(key=str_replace(key, '_.*', '')) %>% 
  spread(key, value)
#snp_tbl <- snp_tbl %>%  mutate(IID = str_replace(IID, '_[ACTG]', ''))
#snp_tbl <- dplyr::rename(snp_tbl,  `18208` = `18121`)
snp_tbl <- snp_tbl[ , order(names(snp_tbl))]
snp_tbl <- rename(snp_tbl, id=SNP)
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

genotypes_file_name = tempfile()
write_tsv(snp_tbl, genotypes_file_name)

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile( genotypes_file_name );

#extract gene positions from GTF file
#if (! file_test("-f", gene_location_file_name)) {
#  system(paste('echo "geneid\tchr\ts1\ts2" >', gene_location_file_name))
#}
#system(paste("cat ~/BTSync/FetalRNAseq/Ref/genes.gtf | awk '{if ($3 == \"gene\") print $10, $1, $4, $5}' | sed 's/[\";]//g' >>", gene_location_file_name))
#gene_location_file_name = "~/BTSync/FetalRNAseq/Github/GENEX-FB2/MatrixEQTL/geneloc.txt"

snpspos <- read_tsv(opt$snps, col_names = c("chr", "snp",	"skip",	"pos"));
snpspos <- snpspos %>% dplyr::select(snp, chr, pos) %>%
  mutate(chr = paste0("chr", chr)) %>%
  as.data.frame()


###############################################################
output_file_name_cis = opt$cis
output_file_name_tra = opt$trans

pvOutputThreshold_cis = opt$p_cis;
pvOutputThreshold_tra = opt$p_trans;

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

save.image(opt$image)
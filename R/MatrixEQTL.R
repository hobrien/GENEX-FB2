#setwd("~/BTSync/FetalRNAseq/Github/GENEX-FB2/")

library(MatrixEQTL)
library(readr)
library(dplyr)
library(tidyr)
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
  make_option(c("--cofactors"), type="character", default=NULL,
              help="Information about each sample"),
  make_option(c("--cis"), type="character", default="../MatrixEQTL/cis_eqtl.txt",
              help="table of cis eQTLs"),
  make_option(c("--trans"), type="character", default="../MatrixEQTL/trans_eqtl.txt",
              help="table of trans eQTLs"),
  make_option(c("--image"), type="character", default="../MatrixEQTL/results.RData",
              help="R image file with all outfiles"),
  make_option(c("--p_trans"), type="numeric", default=1e-8,
              help="p-value cutoff for trans eQTLs"),
  make_option(c("--p_cis"), type="numeric", default=1e-4,
              help="p-value cutoff for cis eQTLs")
 )
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#opt<- list(genotypes="Genotypes/Plink/sig_snps_gene.traw", p_trans=1e-8, p_cis=1e-4, counts=Data/expression_gene.bed.gz", cofactors="Peer/factors.txt", snps="Genotypes/Plink/sig_snps_gene.map", genes="Data/geneloc.txt", cis="MatrixEQTL/cis_eqtl_gene.txt", trans="MatrixEQTL/trans_eqtl_gene.txt", image="MatrixEQTL/results_gene.RData")
#opt<- list(genotypes="Genotypes/Plink/sig_snps_transcript.traw", p_trans=1e-8, p_cis=1e-4, counts="Data/expression_transcript.bed.gz", cofactors="Peer/factors.txt", snps="Genotypes/Plink/sig_snps_transcript.map", genes="Data/transcriptloc.txt", cis="MatrixEQTL/cis_eqtl_transcript.txt", trans="MatrixEQTL/trans_eqtltranscript.txt", image="MatrixEQTL/results_transcript.RData")


# Import normalised counts to filter genepos
print(paste("reading", opt$counts))
counts <- read_delim(opt$counts, "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  select(-one_of(c('#Chr', 'start', 'end')))

print(paste("reading", opt$genes))
genepos = read_delim(opt$genes, "\t", col_names = c("ID", "chr", "s1", "s2", "strand"), escape_double = FALSE, trim_ws = TRUE)
print(paste("finished reading", opt$genes))
genepos <- genepos %>% mutate(ID= str_replace(ID, '(ENS[GT]\\d+)\\.\\d+', '\\1')) %>%
  select(-strand) %>%
  semi_join(counts, by=c("ID")) %>%
  as.data.frame()

# import co-variates
# Seems that data are coerced into float/integer column-wiseduring import with 
# the LoadFile command even though the data are encoded row-wise. I think the 
# only solution to this is to round every thing to whole numbers.

if (!is.null(opt$cofactors)) {
  cofactors_t <- read_tsv(opt$cofactors)
  cofactors <- cofactors_t[,-1] %>% as.matrix() %>% t() %>% as.data.frame()
  colnames(cofactors) <- cofactors_t$ID
  cofactors$ID <- rownames(cofactors)
} else{
  cofactors=counts[0,]
}

print(paste("reading", opt$genotypes))
snp_tbl<- read_delim(opt$genotypes, delim='\t')
snp_tbl <- snp_tbl %>% dplyr::select(-CHR, -`(C)M`, -POS, -COUNTED, -ALT) %>%
  mutate_at(vars(-SNP), funs(str_replace(., '_.*', '')))
#snp_tbl <- snp_tbl %>%  mutate(IID = str_replace(IID, '_[ACTG]', ''))
#snp_tbl <- dplyr::rename(snp_tbl,  `18208` = `18121`)
snp_tbl <- snp_tbl[ , order(names(snp_tbl))]
snp_tbl <- rename(snp_tbl, ID=SNP)
snp_tbl <- dplyr::select(snp_tbl, ID, everything())
colnames(snp_tbl) <- str_replace(colnames(snp_tbl), '_.*', '')

print("removing samples that are missing in one or more of the files")
counts <- counts %>% dplyr::select(-one_of(setdiff(colnames(counts), colnames(cofactors))))
counts <- counts %>% dplyr::select(-one_of(setdiff(colnames(counts), colnames(snp_tbl))))
cofactors <- cofactors %>% dplyr::select(-one_of(setdiff(colnames(cofactors), colnames(counts))))
snp_tbl <- snp_tbl %>% dplyr::select(-one_of(setdiff(colnames(snp_tbl), colnames(counts))))

print("save counts file and read in as matrix EQTL object")
counts_file_name = "temp_counts.txt"
write_tsv(counts, counts_file_name)

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(counts_file_name);

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(!is.null(opt$covariates)) {
  print("save covariates file and read in as matrix EQTL object")
  covariates_file_name = tempfile()
  write_tsv(cofactors, covariates_file_name)
  cvrt$LoadFile(covariates_file_name);
}

print("save genotypes file and read in as matrix EQTL object")
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
print(paste("reading", opt$snps))
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

print("finding eQTLs")
me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold  = pvOutputThreshold_tra,
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

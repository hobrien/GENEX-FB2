library(peer)
library(readr)
library(dplyr)
library(tibble)
library(optparse)

option_list <- list(
  make_option(c("-g", "--genotypes"), type="character", default=NULL, 
              help="genotype matrix (additive)"),
  make_option(c("-c", "--counts"), type="character", default="../examples/brem_data/expression.csv", 
              help="counts matrix (normalised, vst-transformed)"),
  make_option(c("-b", "--batch"), type="character", default=NULL, 
              help="batch factors/covariates"),
  make_option(c("-p", "--pca"), type="character", default=NULL, 
              help="principle components"),
  make_option(c("-n", "--num_factors"), type="integer", default=25, 
              help="Number of PEER factors to estimate"),
  make_option(c("-r", "--residuals"), type="character", default="../examples/brem_data/residuals.txt", 
              help="Outfile for residuals of expression data"),
  make_option(c("-a", "--alpha"), type="character", default="../examples/brem_data/alpha.txt", 
              help="Outfile for alpha values of PEER factors"),
  make_option(c("-f", "--factors"), type="character", default="../examples/brem_data/factors.txt", 
              help="Outfile for PEER factors"),
  make_option(c("-e", "--exclude"), type="character", default='', 
              help="IDs of samples to exclude")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser, positional_arguments=FALSE)
exclude <- strsplit(opt$exclude, ',')[[1]]

if (!is.null(opt$genotypes)) {
  genfile <- read_csv(opt$genotypes) %>% as.data.frame()
  rownames(genfile) <- genfile[,1]
  genfile<-colnames(genfile[,-1])
}

print("reading phenotypes")
pheno <- read_tsv(opt$counts)
print("finished reading phenotypes")
pheno<-select(pheno, -one_of("#Chr", "start", "end")) # added columns in BED file for FastQTL
colnames(pheno)[1]<-"ID" #make sure ID column is consistently named
if (length(exclude) > 0) {
  pheno<-select(pheno, -one_of(exclude))
}

if (!is.null(opt$batch)) {
  print("reading covariates")
  covariates <- read.table(opt$batch, header=TRUE, stringsAsFactors = TRUE)
  colnames(covariates)[1]<-"ID"
  covariates$ID <- as.character(covariates$ID)
  covariates <- covariates[match(colnames(pheno[,-1]), covariates$ID),]
  covariates
}

if (!is.null(opt$pca)) {
  print("reading PCA")
  if (is.null(opt$batch)) {
    covariates <- pheno[1,]
  }
  pca<-read_delim(opt$pca, delim=' ', col_names=FALSE) %>% 
    select(-X1) %>%
    mutate(X2=as.character(X2))
  print("joining PCA, covariates")
  pca
  covariates<-inner_join(covariates, pca, by=c('ID' = 'X2')) %>% as.data.frame()
  print("finished joining PCA, covariates")
}


model=PEER()

PEER_setNk(model, opt$num_factors)
PEER_setPhenoMean(model, t(as.matrix(pheno[,-1])))
if (!is.null(opt$batch) | !is.null(opt$pca)) {
  print("setting covariates")
  PEER_setCovariates(model, as.matrix(mutate_all(covariates[,-1], as.numeric)))
}

print("updating model")
PEER_update(model)

save.image(file="PEER.RData")
PEER_getResiduals(model) %>% as_tibble() %>% write_tsv("residuals_temp.txt")
PEER_getX(model) %>% as_tibble() %>% write_tsv("factors_temp.txt")
PEER_getAlpha(model) %>% as_tibble() %>% write_tsv("alpha_temp.txt")

residuals <- PEER_getResiduals(model) %>% t() %>% as_tibble()
colnames(residuals) <- colnames(pheno[,-1])
residuals$ID <- pheno$ID
residuals <- dplyr::select(residuals, ID, everything())
residuals %>% write_tsv(opt$residuals)

factors <- PEER_getX(model) %>% as_tibble()
factors$ID <- colnames(pheno[,-1])
factors <- dplyr::select(factors, ID, everything())
factors %>% write_tsv(opt$factors)

alpha <- PEER_getAlpha(model) %>% as_tibble()
alpha$factor <- colnames(factors[,-1])
alpha <- dplyr::select(alpha, factor, alpha=V1) 
alpha %>% write_tsv(opt$alpha)

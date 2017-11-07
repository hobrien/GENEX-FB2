library(peer)
library(readr)
library(dplyr)
library(tibble)
library(optparse)

option_list <- list(
  make_option(c("-g", "--genotypes"), type="character", default=NULL, 
              help="genotype matrix (additive)"),
  make_option(c("-e", "--expression"), type="character", default="../examples/brem_data/expression.csv", 
              help="counts matrix (normalised, vst-transformed)"),
  make_option(c("-c", "--covariates"), type="character", default=NULL, 
              help="counts matrix (normalised, vst-transformed)"),
  make_option(c("-p", "--ploidy"), type="character", default='haploid', 
              help="haploid (genotypes=c(0,1) or diploid (genotypes=c(0,1,2))"),
  make_option(c("-n", "--num_factors"), type="integer", default=25, 
              help="Number of PEER factors to estimate"),
  make_option(c("-r", "--residuals"), type="character", default="../examples/brem_data/residuals.txt", 
              help="Outfile for residuals of expression data"),
  make_option(c("-a", "--alpha"), type="character", default="../examples/brem_data/alpha.txt", 
              help="Outfile for alpha values of PEER factors"),
  make_option(c("-o", "--out"), type="character", default="../examples/brem_data/factors.txt", 
              help="Outfile for PEER factors")
  
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser, positional_arguments=FALSE)

if (! opt$ploidy %in% c('haploid', 'diploid')) {
  stop("ploidy type not recognised. Must be one of (haploid, diploid)")
}

if (!is.null(opt$genotypes)) {
  genfile <- read_csv(opt$genotypes) %>% as.data.frame()
  rownames(genfile) <- genfile[,1]
  genfile<-colnames(genfile[,-1])
}
pheno <- read_tsv(opt$expression)
colnames(pheno)[1]<-"ID" #make sure ID column is consistently named

if (!is.null(opt$covariates)) {
  covariates <- read_csv(opt$covariates) %>% as.data.frame()
  rownames(covariates) <- covariates[,1]
  covariates<-colnames(covariates[,-1])
}


model=PEER()

PEER_setNk(model, opt$num_factors)
PEER_setPhenoMean(model, t(as.matrix(pheno[,-1])))
if (!is.null(opt$covariates)) {
    PEER_setCovariates(model, as.matrix(covariates))
}
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

factors <- PEER_getX(model) %>% t() %>% as_tibble()
factors$ID <- pheno$ID
factors <- dplyr::select(factors, ID, everything()) 
factors %>% write_tsv(opt$out)

alpha <- PEER_getAlpha(model) %>% as_tibble()
alpha$factor <- factors$ID
alpha <- dplyr::select(alpha, factor, alpha=V1) 
alpha %>% write_tsv(opt$alpha)

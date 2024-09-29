#### Colocalization Analysis
#### Using Preterm Birth and Blood eQTL as an Example



rm(list=ls())

library(remotes)
library("coloc")
library(dplyr)
library(locuscomparer)
library(stringr)
library(tibble)
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))


# 1.Data Preparation
### 1_1 Preterm birth GWAS summary
pretermgwas = fread("summary/sumstats_preterm_rs_final.txt") 
names(pretermgwas)[names(pretermgwas) == "CHR"] <- "chr"

names(pretermgwas)[names(pretermgwas) == "P"] <- "pvalue"
names(pretermgwas)[names(pretermgwas) == "BETA"] <- "beta"
names(pretermgwas)[names(pretermgwas) == "SNP"] <- "snp"
names(pretermgwas)[names(pretermgwas) == "A2"] <- "ref"
names(pretermgwas)[names(pretermgwas) == "A1"] <- "alt"

names(pretermgwas)[names(pretermgwas) == "POS(hg19)"] <- "position"
names(pretermgwas)[names(pretermgwas) == "AF1"] <- "maf"


pretermgwas$varbeta =  pretermgwas$SE ^ 2

pretermgwas$position <- as.numeric(pretermgwas$position)
head(pretermgwas)

### 1_2 Lead SNP
### The region within ± 500kb from gene start positions was examined
lead <- fread("pretermgwas/sig/coloc_gene.csv")
pretermgwas_group <- pretermgwas %>% 
  group_by(chr)

lead$chr <- as.numeric(lead$chr)
pretermgwas_group$chr <- as.numeric(pretermgwas_group$chr)

matched_rows <- lead %>% 
  left_join(pretermgwas_group, by = "chr", multiple = "all") %>% 
  filter(position.x >= position.y - 500000 & position.x <= position.y + 500000) %>% 
  dplyr::select(snp,chr,position.x,alt, ref,beta, SE,pvalue,varbeta,maf,rsid)
  head(matched_rows) 

matched_rows$pos <- sub(".*:", "", matched_rows$snp)
matched_rows$chrom <-  sub(":.*", "", matched_rows$snp)

gwas_foruse <- matched_rows[,c(12,11,3,6,9,10)]

### 1_3 Blood eQTL
eqtl <- read.table(file="phenocode-Peripheral_blood.tsv", header=T,  sep="\t")
head(eqtl)
eqtl$snp <- paste(eqtl$chrom, eqtl$pos, sep = ":")
eqtl_use <- eqtl[ ,c("snp", "geno_id", "pval", "ref", "alt")]

eqtl$ref <- tolower(eqtl$ref)
eqtl$alt<- tolower(eqtl$alt)


preterm_eqtl <- merge(eqtl, matched_rows, by="snp", all=FALSE, suffixes=c("_eqtl","gwas"))
preterm_eqtl$maf <- as.numeric(preterm_eqtl$maf)

# 2.Colocalization Analysis
all <- NULL
for (i in unique(preterm_eqtl$geno_id)){
 data2 <- preterm_eqtl[which(preterm_eqtl$geno_id==i),]
  res <- coloc.abf(dataset1 = list(snp=data2$snp, beta=data2$betagwas, varbeta=data2$varbeta, MAF=data2$maf, type='cc'),
                   dataset2 = list(snp=data2$snp, pvalues=data2$pval, type='quant',MAF=data2$maf, N=98))
  a <- data.frame(t(res$summary))
  a$geno_id <- i
  all <- rbind(all, a)
}


write.csv(all,"preterm_coloc_blood.csv")

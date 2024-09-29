##### tissus-specific enrichment 

library(dplyr)
library(plyr)
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

# tissue gene data 
ref <- fread("rna_tissue_consensus.tsv")
head(ref)
ref$nTPM <- as.numeric(ref$nTPM)

# protein-coding genes located near the lead SNPs of preterm birth
list <- readLines("sig_gene.txt")

# tissue data preparation
gene_count <- length(unique(ref$`Gene name`))
gene_count
ref2 <- ref[,2:4]
ref_use <-  reshape(ref2, idvar = "Tissue", timevar = "Gene name", direction = "wide")
ref_use1 <- ref_use[,-1]
ref_use1[is.na(ref_use1)] <- 0
zero_counts <- colSums(ref_use1 == 0, na.rm = TRUE)
zero_counts1 <- data.frame(zero_counts)
names(zero_counts1) <- "count"
zero_counts1$rownames_col <- rownames(zero_counts1)
zero_counts_45 <- zero_counts1[zero_counts1$count <45,]
zero_counts_45 <- data.frame(zero_counts_45)
genelist <- unique(zero_counts_45$rownames_col)


ref_use2 <- ref_use1[,genelist, with= FALSE]
colnames(ref_use1)
ref_rank <- apply(ref_use2, 2, function(x) rank(x, ties.method = "min"))
ref_rank <- data.frame(ref_rank)
list_use <- paste0("nTPM.", list)
ref_rank1 <- cbind(ref_use[,1],ref_rank)
ref_rank2 <-reshape2::melt(ref_rank1)
ref_rank2$group <- ifelse(ref_rank2$variable %in% list_use, "A","B")

# enrichment analysis
tissuelist <- unique(ref_rank2$Tissue)
output <- data.frame()
for(i in 1:length(tissuelist)){
  submatr <- ref_rank2[ref_rank2$Tissue == tissuelist[i],]
  submatr$value <- as.numeric(submatr$value)
  group_a <- submatr$value[submatr$group=="A"]
  group_b <- submatr$value[submatr$group=="B"]
  mean_a <- mean(group_a)
  mean_b <- mean(group_b)
  diff= mean_a -mean_b
  stat <- wilcox.test(group_a,group_b)
  result <- data.frame(tissue=tissuelist[i], p=stat$p.value, method=stat$method,diff=diff)
  output <- rbind.fill(output, result)
}

write.csv(output,file = "tissue_enrichment_res.csv")



#### Cell type enrichment analysis 


# cell type RNA data preparation
folder_path <- "E-GEAD-397.processed/count"
file_names <- list.files(folder_path, full.names = TRUE)
result_list <- list()

for (file_name in file_names) {
  data <- read.table(file_name, header = TRUE)
  data_use <- data[, -c(1,2)]
  percentages <- apply(data_use, 1, function(row) sum(row != 0) / length(row))
  new_data <- data.frame(
  gene = data[, 2],
  tissue = basename(file_name),
  percentage = percentages
  )
  result_list[[file_name]] <- new_data
}

final_data <- do.call(rbind, result_list)
rownames(final_data) <- NULL

sc_use <-  reshape(final_data, idvar = "celltype", timevar = "gene", direction = "wide")
sc_use1 <- sc_use[,-1]

sc_rank <- apply(sc_use1, 2, function(x) rank(x, ties.method = "max"))
sc_rank <- data.frame(sc_rank)
list_use_sc <- paste0("percentage.", list)
sc_rank1 <- cbind(sc_use[,1],sc_rank)
sc_rank2 <-reshape2::melt(sc_rank1)
sc_rank2$group <- ifelse(sc_rank2$variable %in% list_use_sc, "A","B")

# enrichment analysis
tissuelist <- unique(sc_rank2$`sc_use[, 1]`)
output <- data.frame()
for(i in 1:length(tissuelist)){
  submatr <- sc_rank2[sc_rank2$`sc_use[, 1]` == tissuelist[i],]
  submatr$value <- as.numeric(submatr$value)
  group_a <- submatr$value[submatr$group=="A"]
  group_b <- submatr$value[submatr$group=="B"]
  mean_a <- mean(group_a)
  mean_b <- mean(group_b)
  diff= mean_a -mean_b
  stat <- wilcox.test(group_a,group_b)
  result <- data.frame(celltype=tissuelist[i], p=stat$p.value, method=stat$method,diff=diff)
  output <- rbind.fill(output, result)
}

write.csv(output,file = "cell_enrichment_res.csv")


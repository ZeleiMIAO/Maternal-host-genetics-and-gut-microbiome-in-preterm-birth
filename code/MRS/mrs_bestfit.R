##### Develop microbial risk score


library(dplyr)
library(magrittr)
library(stringr)
library(tibble)
library(data.table)
library(openxlsx)
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
options(scipen = 999)


#### read gut microbial data
genus_preterm <- fread("data_afterclr.csv")

#### read association analysis results
output_all <- fread("sumamry.txt")

#### score construction
p.threshold <- c(0.0001,0.001,0.01,0.02,0.03,0.04,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)

score_list <- list()

for (threshold in p.threshold) {
  selected_genus_beta <- subset(output_all, p_value.x< threshold, select = c(genus.x, beta_dir, beta.x))
  selected_genus_cols <- match(selected_genus_beta$genus.x, colnames(genus_preterm_use))
  selected_genus_values <- genus_preterm_use[, selected_genus_cols,drop = FALSE]
  # weighted
  score <- apply(selected_genus_values, 1, function(row) sum(row * selected_genus_beta$beta.x* (selected_genus_beta$beta.x > 0), na.rm = TRUE))
  # unweighted
  #score <- apply(selected_genus_values, 1, function(row) sum(row * (selected_genus_beta$beta.x > 0), na.rm = TRUE))
  score_list[[as.character(threshold)]] <- score
  score_threshold_name <- paste("score_p", threshold, sep = "_")
  genus_preterm_use[[score_threshold_name]] <- score
}

#### best fit
prs.result1 <- NULL
for(i in p.threshold){
  score_threshold_name <- paste("score_p", i, sep = "_")
  formula <- paste("preterm ~", score_threshold_name)
  model <- glm(formula, data=genus_preterm_use, family = binomial(link = "logit"))
  model.r2 <- fmsb::NagelkerkeR2(model)
  prs.r2 <- model.r2$R2
  prs.coef <- summary(model)$coeff[score_threshold_name,]
  prs.beta <- as.numeric(prs.coef[1])
  prs.se <- as.numeric(prs.coef[2])
  prs.p <- as.numeric(prs.coef[4])
  prs.result1 <- rbind(prs.result1, data.frame(Threshold=i, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
}

# Best result is:
prs.result1[which.max(prs.result1$R2),]
prs.result1[which.min(prs.result1$P),]

write.table(prs.resul1, "genus_bestfit_weighted.txt", sep="\t", col.names = T, quote = F, row.names = F)


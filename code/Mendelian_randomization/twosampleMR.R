#### Two-Sample MR for preterm and gestational duration

library(devtools)
# install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)


################ crp ###########################
crp <- read_exposure_data(filename ="data_bbj/crp_sig",sep ="\t",snp_col = "SNP",beta_col = "BETA",se_col = "SE", eaf_col = "FRQ",
                              effect_allele_col = "A1",other_allele_col = "A2",pval_col = "P")

crp_lead <- fread("data_bbj/crp_lead.txt.clumped")
crp_lead <- crp_lead[,3]
crp_use <- merge(crp, crp_lead, by="SNP")

outcome_dat1 <- read_outcome_data(
  snps = crp_use$SNP,
  filename = "summary/gw_gsmr.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "freq",
  pval_col = "p",
  samplesize_col = "N"
)

harm_crp1 <- harmonise_data(exposure_dat = crp, outcome_dat = outcome_dat1) 
res_crp1 <- mr(harm_crp1)
res_crp1$type <- "crp_gw"

outcome_dat2 <- read_outcome_data(
  snps = crp_use$SNP,
  filename = "summary/preterm_gsmr.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "freq",
  pval_col = "p",
  samplesize_col = "N"
)

harm_crp2 <- harmonise_data(exposure_dat = crp, outcome_dat = outcome_dat2) 
res_crp2 <- mr(harm_crp2)
res_crp2$type <- "crp_preterm"

mr_scatter_plot(res_crp2,harm_crp2)

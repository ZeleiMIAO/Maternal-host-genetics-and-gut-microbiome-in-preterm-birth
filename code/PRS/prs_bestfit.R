#### Best P values for PRS construction


suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
options (scipen = 200)
p.threshold <- c(0.00000001,0.00000005,0.0000001, 0.000001, 0.00001, 0.0001,0.001,0.05,0.1,0.2,0.3,0.4,0.5)

# Read in the phenotype file 
phenotype <- read.table("webirth_preterm.txt")

names(phenotype) <-
  c("FID",
    "IID",
    "preterm")

# Read in the covariates 
pcs <- read.table("webirth_covarall.txt", header=F)
colnames(pcs) <- c("FID", "IID", paste0("PC",1:10),"age","parity","sexbaby") 
pheno <- merge(phenotype, pcs, by=c("FID", "IID"))

# Null model
null.model <- glm(preterm~., data=pheno[,!colnames(pheno)%in%c("FID","IID")], family = binomial(link = "logit"))
summary(null.model)
library(fmsb)
null.r2 <- fmsb::NagelkerkeR2(null.model)
null.r2$R2

# Go through the p values
prs.result <- NULL

for(i in p.threshold){
prs <- read.table(paste0("webirthpreterm_finalscore.",i,".profile"), header=T)
pheno.prs <- merge(pheno, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
pheno.prs$z <- (pheno.prs$SCORE-mean(pheno.prs$SCORE))/sd(pheno.prs$SCORE)
pheno_use <- pheno.prs[, -which(names(pheno.prs)=="SCORE")]
model <- glm(preterm~., data=pheno_use[,!colnames(pheno_use)%in%c("FID","IID")], family = binomial(link = "logit"))
model.r2 <- fmsb::NagelkerkeR2(model)
prs.r2 <- model.r2$R2-null.r2$R2
prs.coef <- summary(model)$coeff["z",]
prs.beta <- as.numeric(prs.coef[1])
prs.se <- as.numeric(prs.coef[2])
prs.p <- as.numeric(prs.coef[4])
prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
}

# Best result is:
prs.result[which.max(prs.result$R2),]
prs.result[which.min(prs.result$P),]

write.table(prs.result, "webirth_preterm_finalscorefit.txt", sep="\t", col.names = T, quote = F, row.names = F)


#!/usr/bin/env Rscript --vanilla

numberPrincipalComponents <- 100

PCAEVEC <- read.table("${cohortEvec}", header=TRUE)
PHENOTYPE <- read.table("${cohortFam}", header=TRUE)

colnames(PCAEVEC) <- c("FID","IID",paste("PC", 1:numberPrincipalComponents, sep=""),"Pheno")

PCAPHENO <- merge(PCAEVEC,PHENOTYPE)

sink("${output}")

for (i in 1:numberPrincipalComponents) {
    DATA<-as.data.frame(PCAPHENO[,c(3:(i+2))])
    print(summary(lm(PCAPHENO[,104] ~ ., data=DATA)))
}

sink()

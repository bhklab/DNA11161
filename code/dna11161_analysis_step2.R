########################
## Benjamin Haibe-Kains
##
## Nov 23, 2012
########################

rm(list=ls(all=TRUE))

source("ufoo.R")

saveres <- file.path("..", "saveres")
if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE) }

########################
## step 2: Compare Affymetrix-based breast cancer gene expression signatures computed from Affymetrix and Illumina RNA-seq data
########################

## load data
#myfns <- c("jetset"=file.path(saveres, "dna11161_jetset_common.RData"), "bestpgene"=file.path(saveres, "dna11161_bestpgene_common.RData"), "bestptranscript"=file.path(saveres, "dna11161_bestptranscript_common.RData"))

## gene signatures
library(genefu)
library(vcd)
library(epibasix)
library(Hmisc)

## tumors in common
nn <- intersect(rownames(data.affy), rownames(data.rnaseq))

##############
## GENIUS
##############
sig.affy <- genius(data=data.affy, annot=annot.affy, do.mapping=FALSE)
sig.affy$risk <- as.numeric(cut2(x=sig.affy$score, cuts=median(sig.affy$score))) - 1
names(sig.affy$risk) <- names(sig.affy$score)
sig.rnaseq <- genius(data=data.rnaseq, annot=annot.rnaseq, do.mapping=TRUE)
sig.rnaseq$risk <- as.numeric(cut2(x=sig.rnaseq$score, cuts=median(sig.rnaseq$score))) - 1
names(sig.rnaseq$risk) <- names(sig.rnaseq$score)

tt <- table("AFFY"=sig.affy$risk[nn], "RNASEQ"=sig.rnaseq$risk[nn])
attributes(tt)$dimnames$AFFY <- c("Low-risk", "High-risk")
attributes(tt)$dimnames$RNASEQ <- c("Low-risk", "High-risk")
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "genius_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for GENIUS risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

pdf(file.path(saveres, "genius_scores.pdf"))
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col="darkgrey", main="GENIUS\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="GENIUS scores (AFFY)", ylab="GENIUS scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
dev.off()

pdf(file.path(saveres, "genius_scores_risk.pdf"))
mycol <- rep("red", length(nn))
cc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[cc == 0] <- "blue"
mycol[cc == 2] <- "yellow"
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col=mycol, main="GENIUS\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="GENIUS scores (AFFY)", ylab="GENIUS scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
legend("topleft", col=c("blue", "yellow", "red"), legend=c("Low-risk", "High-risk", "Discordance"), bty="n", pch=16)
dev.off()

##############
## GGI
##############
sig.affy <- ggi(data=data.affy, annot=annot.affy, do.mapping=FALSE, hg=demo[rownames(data.affy),"grade"])
sig.rnaseq <- ggi(data=data.rnaseq, annot=annot.rnaseq, do.mapping=TRUE, hg=demo[rownames(data.rnaseq),"grade"])

tt <- table("AFFY"=sig.affy$risk[nn], "RNASEQ"=sig.rnaseq$risk[nn])
attributes(tt)$dimnames$AFFY <- c("Low-risk", "High-risk")
attributes(tt)$dimnames$RNASEQ <- c("Low-risk", "High-risk")
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "ggi_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for GGI risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

pdf(file.path(saveres, "ggi_scores.pdf"))
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col="darkgrey", main="GGI\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="GGI scores (AFFY)", ylab="GGI scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
dev.off()

pdf(file.path(saveres, "ggi_scores_risk.pdf"))
mycol <- rep("red", length(nn))
cc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[cc == 0] <- "blue"
mycol[cc == 2] <- "yellow"
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col=mycol, main="GGI\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="GGI scores (AFFY)", ylab="GGI scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
legend("topleft", col=c("blue", "yellow", "red"), legend=c("Low-risk", "High-risk", "Discordance"), bty="n", pch=16)
dev.off()

##############
## TAMR13
##############
sig.affy <- tamr13(data=data.affy, annot=annot.affy, do.mapping=TRUE)
sig.affy$risk <- as.numeric(cut2(x=sig.affy$score, cuts=median(sig.affy$score))) - 1
names(sig.affy$risk) <- names(sig.affy$score)
sig.rnaseq <- tamr13(data=data.rnaseq, annot=annot.rnaseq, do.mapping=TRUE)
sig.rnaseq$risk <- as.numeric(cut2(x=sig.rnaseq$score, cuts=median(sig.rnaseq$score))) - 1
names(sig.rnaseq$risk) <- names(sig.rnaseq$score)

tt <- table("AFFY"=sig.affy$risk[nn], "RNASEQ"=sig.rnaseq$risk[nn])
attributes(tt)$dimnames$AFFY <- c("Low-risk", "High-risk")
attributes(tt)$dimnames$RNASEQ <- c("Low-risk", "High-risk")
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "tamr13_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for TAMR13 risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

pdf(file.path(saveres, "tamr13_scores.pdf"))
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col="darkgrey", main="TAMR13\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="TAMR13 scores (AFFY)", ylab="TAMR13 scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
dev.off()

pdf(file.path(saveres, "tamr13_scores_risk.pdf"))
mycol <- rep("red", length(nn))
cc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[cc == 0] <- "blue"
mycol[cc == 2] <- "yellow"
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col=mycol, main="TAMR13\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="TAMR13 scores (AFFY)", ylab="TAMR13 scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
legend("topleft", col=c("blue", "yellow", "red"), legend=c("Low-risk", "High-risk", "Discordance"), bty="n", pch=16)
dev.off()

##############
## PIK3CAGS
##############
sig.affy <- list("score"=pik3cags(data=data.affy, annot=annot.affy, do.mapping=FALSE))
sig.affy$risk <- as.numeric(cut2(x=sig.affy$score, cuts=median(sig.affy$score))) - 1
names(sig.affy$risk) <- names(sig.affy$score)
sig.rnaseq <- list("score"=pik3cags(data=data.rnaseq, annot=annot.rnaseq, do.mapping=TRUE))
sig.rnaseq$risk <- as.numeric(cut2(x=sig.rnaseq$score, cuts=median(sig.rnaseq$score))) - 1
names(sig.rnaseq$risk) <- names(sig.rnaseq$score)

tt <- table("AFFY"=sig.affy$risk[nn], "RNASEQ"=sig.rnaseq$risk[nn])
attributes(tt)$dimnames$AFFY <- c("Low-risk", "High-risk")
attributes(tt)$dimnames$RNASEQ <- c("Low-risk", "High-risk")
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "pik3cags_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for PIK3CAGS risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

pdf(file.path(saveres, "pik3cags_scores.pdf"))
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col="darkgrey", main="PIK3CAGS\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="PIK3CAGS scores (AFFY)", ylab="PIK3CAGS scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
dev.off()

pdf(file.path(saveres, "pik3cags_scores_risk.pdf"))
mycol <- rep("red", length(nn))
cc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[cc == 0] <- "blue"
mycol[cc == 2] <- "yellow"
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col=mycol, main="PIK3CAGS\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="PIK3CAGS scores (AFFY)", ylab="PIK3CAGS scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
legend("topleft", col=c("blue", "yellow", "red"), legend=c("Low-risk", "High-risk", "Discordance"), bty="n", pch=16)
dev.off()



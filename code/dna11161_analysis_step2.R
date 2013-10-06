########################
## Benjamin Haibe-Kains
##
## Nov 23, 2012
########################

rm(list=ls(all=TRUE))

source(file.path("code", "ufoo.R"))

saveres <- file.path("saveres")
if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE) }
  
colo <- c("darkblue", "darkorange", "darkred")

########################
## step 2a: Compare breast cancer gene expression signatures computed from Affymetrix and Illumina RNA-seq data
########################

## load data
#myfns <- c("jetset"=file.path(saveres, "dna11161_jetset_common.RData"), "bestpgene"=file.path(saveres, "dna11161_bestpgene_common.RData"), "bestptranscript"=file.path(saveres, "dna11161_bestptranscript_common.RData"))

library(genefu)
library(vcd)
library(epibasix)
library(Hmisc)
library(plotrix)
library(gdata)

## read data
load(file.path(saveres, "DNA11161_demo.RData"))
load(file.path(saveres, "dna11161_affy_frma.RData"))
load(file.path(saveres, "dna11161_gene_rnaseq.RData"))

## rnaseq: selection of unique Entrez genid, keep the most variant ENSG in case of ambiguities
gid.rnaseq <- as.character(annot.gene.rnaseq[ ,"EntrezGene.ID"])
names(gid.rnaseq) <- rownames(annot.gene.rnaseq)
data.rnaseq <- data.gene.rnaseq[ , !is.na(gid.rnaseq) & annot.gene.rnaseq[ ,"best"]]
annot.rnaseq <- annot.gene.rnaseq[colnames(data.rnaseq), ,drop=FALSE]
colnames(data.rnaseq) <- rownames(annot.rnaseq) <- paste("geneid", gid.rnaseq[colnames(data.rnaseq)], sep=".")

## tumors in common
nn <- intersect(rownames(data.affy), rownames(data.rnaseq))

data.affy <- data.affy[nn, , drop=FALSE]
data.rnaseq <- data.rnaseq[nn, , drop=FALSE]
demo <- demo[nn, , drop=FALSE]

res.scores <- res.risks <- res.tabs <- NULL

sig.affy <- NULL
sig.rnaseq <- NULL

## read signatures from ignatiadis et al, JCO, 2012
## signatures from Ignatiadis et al. 2012
dwl.status <- download.file(url="http://jco.ascopubs.org/content/suppl/2012/04/16/JCO.2011.39.5624.DC1/Appendix_module_composition_file.xls", destfile=file.path(saveres, "gene_modules.xls"))
if(dwl.status != 0) { stop("Download failed, please rerun the script!") }
## read xls file
mfile <- gdata::read.xls(file.path(saveres, "gene_modules.xls"), skip=1)

	ix.delim <- c(which(mfile[ ,1] != "")[-1]-1, nrow(mfile) + 1)
	ix.f <- ix.l <- 1
	groups <- NULL
	npp <- np <- NULL
	for (i in 1:length(ix.delim)) {
		ix.l <- ix.delim[i]
		np <- c(np, as.character(mfile[ix.f,1]))
		groups <- c(groups, rep(i, ix.l - ix.f + 1))
		npp <- rbind(npp, mfile[ix.f:ix.l,2:ncol(mfile)])
		ix.f <- ix.l + 1
	}
	ugroups <- unique(groups)
	obj <- NULL
	for (j in 1:length(ugroups)) {
		obj <- c(obj, list(npp[groups == ugroups[j], ]))
	}
names(obj) <- np
## clean modules
gene.modules <- lapply(obj, function(x) {
  x <- x[complete.cases(x), , drop=FALSE]
  x <- x[!duplicated(x[ , "EntrezGene.ID"]), , drop=FALSE]
  rownames(x) <- paste("geneid", as.character(x[ , "EntrezGene.ID"]), sep=".")
  x <- cbind("probe"=rownames(x), x)
  return(x)
})

## select "CIN70", "IRMODULE", "RAS", "MAPK", "PTEN", "AKTmTOR", "IGF1", "SRC", "MYC", "E2F3", and "BetaCatenin"
gene.modules <- gene.modules[c("CIN70", "Immune1", "RAS", "MAPK", "PTEN", "AKTmTOR", "IGF1", "SRC", "MYC", "E2F3", "BetaCatenin")]
names(gene.modules)[names(gene.modules) %in% "Immune1"] <- "IRMODULE"

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
plot(t(tt), main="Contingency table for GGI risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.3g, 95%%CI [%.3g,%.3g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

res.tabs <- c(res.tabs, list("GGI"=tt))
res.risks <- c(res.risks, list("GGI"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))

cc <- cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")
ttc <- spearmanCI(x=cc, n=sum(complete.cases(sig.affy$score[nn], sig.rnaseq$score[nn])), alpha=0.05)
pdf(file.path(saveres, "ggi_scores.pdf"))
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col="darkgrey", main="GGI\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="GGI scores (AFFY)", ylab="GGI scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

ttt <- c("rho"=cc, "lower"=ttc["lower"], "upper"=ttc["upper"], "p"=ttc["p.value"])
names(ttt) <- c("rho", "lower", "upper", "p")
res.scores <- c(res.scores, list("GGI"=ttt))

pdf(file.path(saveres, "ggi_scores_risk.pdf"))
mycol <- rep("darkgrey", length(nn))
ccc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[ccc == 0] <- colo[1]
mycol[ccc == 2] <- colo[3]
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col=mycol, main="GGI\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="GGI scores (AFFY)", ylab="GGI scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
legend("topleft", col=c(colo[1], colo[3], "darkgrey"), legend=c("Low-risk", "High-risk", "Discordance"), bty="n", pch=16)
dev.off()

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
plot(t(tt), main="Contingency table for GENIUS risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.3g, 95%%CI [%.3g,%.3g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

res.tabs <- c(res.tabs, list("GENIUS"=tt))
res.risks <- c(res.risks, list("GENIUS"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))

cc <- cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")
ttc <- spearmanCI(x=cc, n=sum(complete.cases(sig.affy$score[nn], sig.rnaseq$score[nn])), alpha=0.05)
pdf(file.path(saveres, "genius_scores.pdf"))
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col="darkgrey", main="GENIUS\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="GENIUS scores (AFFY)", ylab="GENIUS scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

ttt <- c("rho"=cc, "lower"=ttc["lower"], "upper"=ttc["upper"], "p"=ttc["p.value"])
names(ttt) <- c("rho", "lower", "upper", "p")
res.scores <- c(res.scores, list("GENIUS"=ttt))

pdf(file.path(saveres, "genius_scores_risk.pdf"))
mycol <- rep("darkgrey", length(nn))
ccc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[ccc == 0] <- colo[1]
mycol[ccc == 2] <- colo[3]
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col=mycol, main="GENIUS\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="GENIUS scores (AFFY)", ylab="GENIUS scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
legend("topleft", col=c(colo[1], colo[3], "darkgrey"), legend=c("Low-risk", "High-risk", "Discordance"), bty="n", pch=16)
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

##############
## TAMR13
##############
#sig.affy <- tamr13(data=data.affy, annot=annot.affy, do.mapping=TRUE)
#sig.affy$risk <- as.numeric(cut2(x=sig.affy$score, cuts=median(sig.affy$score))) - 1
#names(sig.affy$risk) <- names(sig.affy$score)
#sig.rnaseq <- tamr13(data=data.rnaseq, annot=annot.rnaseq, do.mapping=TRUE)
#sig.rnaseq$risk <- as.numeric(cut2(x=sig.rnaseq$score, cuts=median(sig.rnaseq$score))) - 1
#names(sig.rnaseq$risk) <- names(sig.rnaseq$score)

#tt <- table("AFFY"=sig.affy$risk[nn], "RNASEQ"=sig.rnaseq$risk[nn])
#attributes(tt)$dimnames$AFFY <- c("Low-risk", "High-risk")
#attributes(tt)$dimnames$RNASEQ <- c("Low-risk", "High-risk")
#print(tt)
#tts <- assocstats(tt)
#print(tts)
#kk <- epiKappa(tt, k0=0)
#print(kk)
#pdf(file.path(saveres, "tamr13_risk_contingency_table.pdf"))
#plot(t(tt), main="Contingency table for TAMR13 risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.3g, 95%%CI [%.3g,%.3g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
#dev.off()

#res.tabs <- c(res.tabs, list("TAMR13"=tt))
#res.risks <- c(res.risks, list("TAMR13"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))

#cc <- cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")
#ttc <- spearmanCI(x=cc, n=sum(complete.cases(sig.affy$score[nn], sig.rnaseq$score[nn])), alpha=0.05)
#pdf(file.path(saveres, "tamr13_scores.pdf"))
#plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col="darkgrey", main="TAMR13\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="TAMR13 scores (AFFY)", ylab="TAMR13 scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
#mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
#abline(mm)
#dev.off()

#ttt <- c("rho"=cc, "lower"=ttc["lower"], "upper"=ttc["upper"], "p"=ttc["p.value"])
#names(ttt) <- c("rho", "lower", "upper", "p")
#res.scores <- c(res.scores, list("TAMR13"=ttt))

#pdf(file.path(saveres, "tamr13_scores_risk.pdf"))
#mycol <- rep("darkgrey", length(nn))
#ccc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
#mycol[ccc == 0] <- colo[1]
#mycol[ccc == 2] <- colo[3]
#plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col=mycol, main="TAMR13\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="TAMR13 scores (AFFY)", ylab="TAMR13 scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
#legend("topleft", col=c(colo[1], colo[3], "darkgrey"), legend=c("Low-risk", "High-risk", "Discordance"), bty="n", pch=16)
#mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
#abline(mm)
#dev.off()

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
plot(t(tt), main="Contingency table for PIK3CAGS risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.3g, 95%%CI [%.3g,%.3g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

res.tabs <- c(res.tabs, list("PIK3CAGS"=tt))
res.risks <- c(res.risks, list("PIK3CAGS"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))

cc <- cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")
ttc <- spearmanCI(x=cc, n=sum(complete.cases(sig.affy$score[nn], sig.rnaseq$score[nn])), alpha=0.05)
pdf(file.path(saveres, "pik3cags_scores.pdf"))
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col="darkgrey", main="PIK3CAGS\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="PIK3CAGS scores (AFFY)", ylab="PIK3CAGS scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

ttt <- c("rho"=cc, "lower"=ttc["lower"], "upper"=ttc["upper"], "p"=ttc["p.value"])
names(ttt) <- c("rho", "lower", "upper", "p")
res.scores <- c(res.scores, list("PIK3CAGS"=ttt))

pdf(file.path(saveres, "pik3cags_scores_risk.pdf"))
mycol <- rep("darkgrey", length(nn))
ccc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[ccc == 0] <- colo[1]
mycol[ccc == 2] <- colo[3]
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col=mycol, main="PIK3CAGS\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="PIK3CAGS scores (AFFY)", ylab="PIK3CAGS scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
legend("topleft", col=c(colo[1], colo[3], "darkgrey"), legend=c("Low-risk", "High-risk", "Discordance"), bty="n", pch=16)
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

##############
## PLAUMODULE
##############

ss <- genefu::sig.score(x=mod1$PLAU, data=data.affy, annot=annot.affy, do.mapping=FALSE)$score
rr <- as.numeric(cut2(x=ss, cuts=median(ss))) - 1
names(rr) <- names(ss)
sig.affy <- list("score"=ss, "risk"=rr)
ss <- genefu::sig.score(x=mod1$PLAU, data=data.rnaseq, annot=annot.rnaseq, do.mapping=TRUE)$score
rr <- as.numeric(cut2(x=ss, cuts=median(ss))) - 1
names(rr) <- names(ss)
sig.rnaseq <- list("score"=ss, "risk"=rr)

tt <- table("AFFY"=sig.affy$risk[nn], "RNASEQ"=sig.rnaseq$risk[nn])
attributes(tt)$dimnames$AFFY <- c("Low-risk", "High-risk")
attributes(tt)$dimnames$RNASEQ <- c("Low-risk", "High-risk")
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "plaumodule_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for PLAUMODULE risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.3g, 95%%CI [%.3g,%.3g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

res.tabs <- c(res.tabs, list("PLAUMODULE"=tt))
res.risks <- c(res.risks, list("PLAUMODULE"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))

cc <- cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")
ttc <- spearmanCI(x=cc, n=sum(complete.cases(sig.affy$score[nn], sig.rnaseq$score[nn])), alpha=0.05)
pdf(file.path(saveres, "plaumodule_scores.pdf"))
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col="darkgrey", main="PIK3CAGS\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="PIK3CAGS scores (AFFY)", ylab="PLAUMODULE scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

ttt <- c("rho"=cc, "lower"=ttc["lower"], "upper"=ttc["upper"], "p"=ttc["p.value"])
names(ttt) <- c("rho", "lower", "upper", "p")
res.scores <- c(res.scores, list("PLAUMODULE"=ttt))

pdf(file.path(saveres, "plaumodule_scores_risk.pdf"))
mycol <- rep("darkgrey", length(nn))
ccc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[ccc == 0] <- colo[1]
mycol[ccc == 2] <- colo[3]
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col=mycol, main="PIK3CAGS\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="PLAUMODULE scores (AFFY)", ylab="PIK3CAGS scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
legend("topleft", col=c(colo[1], colo[3], "darkgrey"), legend=c("Low-risk", "High-risk", "Discordance"), bty="n", pch=16)
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

##############
## STAT1MODULE
##############

ss <- genefu::sig.score(x=mod1$STAT1, data=data.affy, annot=annot.affy, do.mapping=FALSE)$score
rr <- as.numeric(cut2(x=ss, cuts=median(ss))) - 1
names(rr) <- names(ss)
sig.affy <- list("score"=ss, "risk"=rr)
ss <- genefu::sig.score(x=mod1$STAT1, data=data.rnaseq, annot=annot.rnaseq, do.mapping=TRUE)$score
rr <- as.numeric(cut2(x=ss, cuts=median(ss))) - 1
names(rr) <- names(ss)
sig.rnaseq <- list("score"=ss, "risk"=rr)

tt <- table("AFFY"=sig.affy$risk[nn], "RNASEQ"=sig.rnaseq$risk[nn])
attributes(tt)$dimnames$AFFY <- c("Low-risk", "High-risk")
attributes(tt)$dimnames$RNASEQ <- c("Low-risk", "High-risk")
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "stat1module_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for STAT1MODULE risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.3g, 95%%CI [%.3g,%.3g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

res.tabs <- c(res.tabs, list("STAT1MODULE"=tt))
res.risks <- c(res.risks, list("STAT1MODULE"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))

cc <- cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")
ttc <- spearmanCI(x=cc, n=sum(complete.cases(sig.affy$score[nn], sig.rnaseq$score[nn])), alpha=0.05)
pdf(file.path(saveres, "stat1module_scores.pdf"))
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col="darkgrey", main="PIK3CAGS\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="PIK3CAGS scores (AFFY)", ylab="STAT1MODULE scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

ttt <- c("rho"=cc, "lower"=ttc["lower"], "upper"=ttc["upper"], "p"=ttc["p.value"])
names(ttt) <- c("rho", "lower", "upper", "p")
res.scores <- c(res.scores, list("STAT1MODULE"=ttt))

pdf(file.path(saveres, "stat1module_scores_risk.pdf"))
mycol <- rep("darkgrey", length(nn))
ccc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[ccc == 0] <- colo[1]
mycol[ccc == 2] <- colo[3]
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col=mycol, main="PIK3CAGS\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="STAT1MODULE scores (AFFY)", ylab="PIK3CAGS scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
legend("topleft", col=c(colo[1], colo[3], "darkgrey"), legend=c("Low-risk", "High-risk", "Discordance"), bty="n", pch=16)
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()


########################
## step 2b: Compare NON Affymetrix-based breast cancer gene expression signatures computed from Affymetrix and Illumina RNA-seq data
########################

## load data
load(file.path(saveres, "dna11161_common_jetset.RData"))
annotc.affy <- annotc
annotc.rnaseq <- annotc

## tumors in common
nn <- intersect(rownames(datac.affy), rownames(datac.rnaseq))

##############
## ER
##############
sig.affy$score <- datac.affy[ , "geneid.2099"]
sig.rnaseq$score <- datac.rnaseq[ , "geneid.2099"]

cc <- cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")
ttc <- spearmanCI(x=cc, n=sum(complete.cases(sig.affy$score[nn], sig.rnaseq$score[nn])), alpha=0.05)
pdf(file.path(saveres, "er_expression_ihc.pdf"))
mycol <- rep("grey", length(nn))
names(mycol) <- nn
mycol[!is.na(demo[nn, "er_status"]) & demo[nn, "er_status"] == 1] <- colo[3]
mycol[!is.na(demo[nn, "er_status"]) & demo[nn, "er_status"] == 0] <- colo[1]
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col=mycol, main="ER\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="ER gene expression (AFFY)", ylab="ER gene expression (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
legend("topleft", col=c(colo[1], colo[3]), title="ER status (IHC)", legend=c("negative", "positive"), bty="n", pch=16)
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

##############
## PGR
##############
sig.affy$score <- datac.affy[ , "geneid.5241"]
sig.rnaseq$score <- datac.rnaseq[ , "geneid.5241"]

cc <- cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")
ttc <- spearmanCI(x=cc, n=sum(complete.cases(sig.affy$score[nn], sig.rnaseq$score[nn])), alpha=0.05)
pdf(file.path(saveres, "pgr_expression_ihc.pdf"))
mycol <- rep("grey", length(nn))
names(mycol) <- nn
mycol[!is.na(demo[nn, "pgr_status"]) & demo[nn, "pgr_status"] == 1] <- colo[3]
mycol[!is.na(demo[nn, "pgr_status"]) & demo[nn, "pgr_status"] == 0] <- colo[1]
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col=mycol, main="PGR\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="PGR gene expression (AFFY)", ylab="PGR gene expression (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
legend("topleft", col=c(colo[1], colo[3]), title="PGR status (IHC)", legend=c("negative", "positive"), bty="n", pch=16)
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

##############
## HER2
##############
sig.affy$score <- datac.affy[ , "geneid.2064"]
sig.rnaseq$score <- datac.rnaseq[ , "geneid.2064"]

her2status <- rep(NA, length(nn))
names(her2status) <- nn
her2status[!is.na(demo[nn, "HER2_status"]) & is.element(demo[nn, "HER2_status"], c(0, 1))] <- 0
her2status[!is.na(demo[nn, "HER2_status"]) & is.element(demo[nn, "HER2_status"], c(3))] <- 1

cc <- cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")
ttc <- spearmanCI(x=cc, n=sum(complete.cases(sig.affy$score[nn], sig.rnaseq$score[nn])), alpha=0.05)
pdf(file.path(saveres, "her2_expression_ihc.pdf"))
mycol <- rep("grey", length(nn))
names(mycol) <- nn
mycol[!is.na(her2status) & her2status == 1] <- colo[3]
mycol[!is.na(her2status) & her2status == 0] <- colo[1]
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col=mycol, main="HER2\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="HER2 gene expression (AFFY)", ylab="HER2 gene expression (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
legend("topleft", col=c(colo[1], colo[3]), title="HER2 status (IHC)", legend=c("negative", "positive"), bty="n", pch=16)
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

##############
## GENE70
##############
sig.affy <- gene70(data=datac.affy, annot=annotc.affy, do.mapping=TRUE, std="robust")
sig.rnaseq <- gene70(data=datac.rnaseq, annot=annotc.rnaseq, do.mapping=TRUE, std="robust")

tt <- table("AFFY"=sig.affy$risk[nn], "RNASEQ"=sig.rnaseq$risk[nn])
attributes(tt)$dimnames$AFFY <- c("Low-risk", "High-risk")
attributes(tt)$dimnames$RNASEQ <- c("Low-risk", "High-risk")
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "gene70_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for GENE70 risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.3g, 95%%CI [%.3g,%.3g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

res.tabs <- c(res.tabs, list("GENE70"=tt))
res.risks <- c(res.risks, list("GENE70"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))

cc <- cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")
ttc <- spearmanCI(x=cc, n=sum(complete.cases(sig.affy$score[nn], sig.rnaseq$score[nn])), alpha=0.05)
pdf(file.path(saveres, "gene70_scores.pdf"))
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col="darkgrey", main="GENE70\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="GENE70 scores (AFFY)", ylab="GENE70 scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

res.scores <- c(res.scores, list("GENE70"=c("rho"=cc, "lower"=ttc["lower"], "upper"=ttc["upper"], "p"=ttc["p.value"])))

pdf(file.path(saveres, "gene70_scores_risk.pdf"))
mycol <- rep("darkgrey", length(nn))
ccc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[ccc == 0] <- colo[1]
mycol[ccc == 2] <- colo[3]
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col=mycol, main="GENE70\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="GENE70 scores (AFFY)", ylab="GENE70 scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
legend("topleft", col=c(colo[1], colo[3], "darkgrey"), legend=c("Low-risk", "High-risk", "Discordance"), bty="n", pch=16)
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

##############
## GENE21
##############
sig.oncotypedx.save <- sig.oncotypedx
sig.oncotypedx["GSTM1", "EntrezGene.ID"] <- "2948"
## the affymetrix probes for GSTM1 are ambiguious with GSTM2, GSTM3, and GSTM4; GSTM3 is present in the data
## http://xavierlab2.mgh.harvard.edu/EnrichmentProfiler/primary/Enrichment/204550_x_at.html
sig.affy <- oncotypedx(data=datac.affy, annot=annotc.affy, do.mapping=TRUE)
sig.rnaseq <- oncotypedx(data=datac.rnaseq, annot=annotc.rnaseq, do.mapping=TRUE)
sig.oncotypedx <- sig.oncotypedx.save

tt <- table("AFFY"=sig.affy$risk[nn], "RNASEQ"=sig.rnaseq$risk[nn])
attributes(tt)$dimnames$AFFY <- c("Low-risk", "Intermediate-risk", "High-risk")
attributes(tt)$dimnames$RNASEQ <- c("Low-risk", "Intermediate-risk", "High-risk")
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "gene21_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for GENE21 risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.3g, 95%%CI [%.3g,%.3g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

res.tabs <- c(res.tabs, list("GENE21"=tt))
res.risks <- c(res.risks, list("GENE21"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))

cc <- cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")
ttc <- spearmanCI(x=cc, n=sum(complete.cases(sig.affy$score[nn], sig.rnaseq$score[nn])), alpha=0.05)
pdf(file.path(saveres, "gene21_scores.pdf"))
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col="darkgrey", main="GENE21\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="GENE21 scores (AFFY)", ylab="GENE21 scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

res.scores <- c(res.scores, list("GENE21"=c("rho"=cc, "lower"=ttc["lower"], "upper"=ttc["upper"], "p"=ttc["p.value"])))

pdf(file.path(saveres, "gene21_scores_risk.pdf"))
mycol <- rep("darkgrey", length(nn))
mycol[complete.cases(sig.affy$risk[nn], sig.rnaseq$risk[nn]) & sig.affy$risk[nn] == 0 & sig.rnaseq$risk[nn] == 0] <- colo[1]
mycol[complete.cases(sig.affy$risk[nn], sig.rnaseq$risk[nn]) & sig.affy$risk[nn] == 0.5 & sig.rnaseq$risk[nn] == 0.5] <- colo[2]
mycol[complete.cases(sig.affy$risk[nn], sig.rnaseq$risk[nn]) & sig.affy$risk[nn] == 1 & sig.rnaseq$risk[nn] == 1] <- colo[3]
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col=mycol, main="GENE21\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="GENE21 scores (AFFY)", ylab="GENE21 scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
legend("topleft", col=c(colo[1], colo[2], colo[3], "darkgrey"), legend=c("Low-risk", "Intermediate-risk", "High-risk", "Discordance"), bty="n", pch=16)
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

##############
## RORS
##############
sig.affy <- rorS(data=datac.affy, annot=annotc.affy, do.mapping=TRUE)
sig.rnaseq <- rorS(data=datac.rnaseq, annot=annotc.rnaseq, do.mapping=TRUE)

tt <- table("AFFY"=sig.affy$risk[nn], "RNASEQ"=sig.rnaseq$risk[nn])
attributes(tt)$dimnames$AFFY <- c("Low-risk", "Intermediate-risk", "High-risk")
attributes(tt)$dimnames$RNASEQ <- c("Low-risk", "Intermediate-risk", "High-risk")
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "rors_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for RORS risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.3g, 95%%CI [%.3g,%.3g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

res.tabs <- c(res.tabs, list("RORS"=tt))
res.risks <- c(res.risks, list("RORS"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))

cc <- cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")
ttc <- spearmanCI(x=cc, n=sum(complete.cases(sig.affy$score[nn], sig.rnaseq$score[nn])), alpha=0.05)
pdf(file.path(saveres, "rors_scores.pdf"))
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col="darkgrey", main="RORS\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="RORS scores (AFFY)", ylab="RORS scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

res.scores <- c(res.scores, list("RORS"=c("rho"=cc, "lower"=ttc["lower"], "upper"=ttc["upper"], "p"=ttc["p.value"])))

pdf(file.path(saveres, "rors_scores_risk.pdf"))
mycol <- rep("darkgrey", length(nn))
mycol[complete.cases(sig.affy$risk[nn], sig.rnaseq$risk[nn]) & sig.affy$risk[nn] == "Low" & sig.rnaseq$risk[nn] == "Low"] <- colo[1]
mycol[complete.cases(sig.affy$risk[nn], sig.rnaseq$risk[nn]) & sig.affy$risk[nn] == "Intermediate" & sig.rnaseq$risk[nn] == "Intermediate"] <- colo[2]
mycol[complete.cases(sig.affy$risk[nn], sig.rnaseq$risk[nn]) & sig.affy$risk[nn] == "High" & sig.rnaseq$risk[nn] == "High"] <- colo[3]
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col=mycol, main="RORS\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="RORS scores (AFFY)", ylab="RORS scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
legend("topleft", col=c(colo[1], colo[2], colo[3], "darkgrey"), legend=c("Low-risk", "Intermediate-risk", "High-risk", "Discordance"), bty="n", pch=16)
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

##############
## EndoPredict
##############
sig.affy <- endoPredict(data=datac.affy, annot=annotc.affy, do.mapping=TRUE)
sig.rnaseq <- endoPredict(data=datac.rnaseq, annot=annotc.rnaseq, do.mapping=TRUE)

tt <- table("AFFY"=sig.affy$risk[nn], "RNASEQ"=sig.rnaseq$risk[nn])
attributes(tt)$dimnames$AFFY <- c("Low-risk", "High-risk")
attributes(tt)$dimnames$RNASEQ <- c("Low-risk", "High-risk")
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "endopredict_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for ENDOPREDICT DX risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.3g, 95%%CI [%.3g,%.3g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

res.tabs <- c(res.tabs, list("ENDOPREDICT"=tt))
res.risks <- c(res.risks, list("ENDOPREDICT"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))

cc <- cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")
ttc <- spearmanCI(x=cc, n=sum(complete.cases(sig.affy$score[nn], sig.rnaseq$score[nn])), alpha=0.05)
pdf(file.path(saveres, "endopredict_scores.pdf"))
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col="darkgrey", main="ENDOPREDICT\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="ENDOPREDICT scores (AFFY)", ylab="ENDOPREDICT DX scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

res.scores <- c(res.scores, list("ENDOPREDICT"=c("rho"=cc, "lower"=ttc["lower"], "upper"=ttc["upper"], "p"=ttc["p.value"])))

pdf(file.path(saveres, "endopredict_scores_risk.pdf"))
mycol <- rep("darkgrey", length(nn))
mycol[complete.cases(sig.affy$risk[nn], sig.rnaseq$risk[nn]) & sig.affy$risk[nn] == 0 & sig.rnaseq$risk[nn] == 0] <- colo[1]
mycol[complete.cases(sig.affy$risk[nn], sig.rnaseq$risk[nn]) & sig.affy$risk[nn] == 0.5 & sig.rnaseq$risk[nn] == 0.5] <- colo[2]
mycol[complete.cases(sig.affy$risk[nn], sig.rnaseq$risk[nn]) & sig.affy$risk[nn] == 1 & sig.rnaseq$risk[nn] == 1] <- colo[3]
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col=mycol, main="ENDOPREDICT DX\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="ENDOPREDICT scores (AFFY)", ylab="ENDOPREDICT scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
legend("topleft", col=c(colo[1], colo[3], "darkgrey"), legend=c("Low-risk", "High-risk", "Discordance"), bty="n", pch=16)
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

##############
## DCN
##############

## read published signature
sigg <- read.csv(file.path("sigs", "farmer2007_dcn_signature_50.csv"), stringsAsFactor=FALSE)

ss <- genefu::sig.score(x=sigg, data=datac.affy, annot=annotc.affy, do.mapping=TRUE)$score
rr <- as.numeric(cut2(x=ss, cuts=median(ss))) - 1
names(rr) <- names(ss)
sig.affy <- list("score"=ss, "risk"=rr)
ss <- genefu::sig.score(x=sigg, data=datac.rnaseq, annot=annotc.rnaseq, do.mapping=TRUE)$score
rr <- as.numeric(cut2(x=ss, cuts=median(ss))) - 1
names(rr) <- names(ss)
sig.rnaseq <- list("score"=ss, "risk"=rr)

tt <- table("AFFY"=sig.affy$risk[nn], "RNASEQ"=sig.rnaseq$risk[nn])
attributes(tt)$dimnames$AFFY <- c("Low-risk", "High-risk")
attributes(tt)$dimnames$RNASEQ <- c("Low-risk", "High-risk")
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "dcn_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for DCN risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.3g, 95%%CI [%.3g,%.3g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

res.tabs <- c(res.tabs, list("DCN"=tt))
res.risks <- c(res.risks, list("DCN"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))

cc <- cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")
ttc <- spearmanCI(x=cc, n=sum(complete.cases(sig.affy$score[nn], sig.rnaseq$score[nn])), alpha=0.05)
pdf(file.path(saveres, "dcn_scores.pdf"))
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col="darkgrey", main="DCN\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="DCN scores (AFFY)", ylab="DCN scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

ttt <- c("rho"=cc, "lower"=ttc["lower"], "upper"=ttc["upper"], "p"=ttc["p.value"])
names(ttt) <- c("rho", "lower", "upper", "p")
res.scores <- c(res.scores, list("DCN"=ttt))

pdf(file.path(saveres, "dcn_scores_risk.pdf"))
mycol <- rep("darkgrey", length(nn))
ccc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[ccc == 0] <- colo[1]
mycol[ccc == 2] <- colo[3]
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col=mycol, main="DCN\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="DCN scores (AFFY)", ylab="DCN scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
legend("topleft", col=c(colo[1], colo[3], "darkgrey"), legend=c("Low-risk", "High-risk", "Discordance"), bty="n", pch=16)
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

##############
## STROMACD10
##############

## read published signature
sigg <- read.csv(file.path("sigs", "desmedt2012_stroma_cd10_signature_12.csv"), stringsAsFactor=FALSE)

ss <- genefu::sig.score(x=sigg, data=datac.affy, annot=annotc.affy, do.mapping=TRUE)$score
rr <- as.numeric(cut2(x=ss, cuts=median(ss))) - 1
names(rr) <- names(ss)
sig.affy <- list("score"=ss, "risk"=rr)
ss <- genefu::sig.score(x=sigg, data=datac.rnaseq, annot=annotc.rnaseq, do.mapping=TRUE)$score
rr <- as.numeric(cut2(x=ss, cuts=median(ss))) - 1
names(rr) <- names(ss)
sig.rnaseq <- list("score"=ss, "risk"=rr)

tt <- table("AFFY"=sig.affy$risk[nn], "RNASEQ"=sig.rnaseq$risk[nn])
attributes(tt)$dimnames$AFFY <- c("Low-risk", "High-risk")
attributes(tt)$dimnames$RNASEQ <- c("Low-risk", "High-risk")
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "stromacd10_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for STROMACD10 risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.3g, 95%%CI [%.3g,%.3g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

res.tabs <- c(res.tabs, list("STROMACD10"=tt))
res.risks <- c(res.risks, list("STROMACD10"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))

cc <- cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")
ttc <- spearmanCI(x=cc, n=sum(complete.cases(sig.affy$score[nn], sig.rnaseq$score[nn])), alpha=0.05)
pdf(file.path(saveres, "stromacd10_scores.pdf"))
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col="darkgrey", main="STROMACD10\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="STROMACD10 scores (AFFY)", ylab="STROMACD10 scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

ttt <- c("rho"=cc, "lower"=ttc["lower"], "upper"=ttc["upper"], "p"=ttc["p.value"])
names(ttt) <- c("rho", "lower", "upper", "p")
res.scores <- c(res.scores, list("STROMACD10"=ttt))

pdf(file.path(saveres, "stromacd10_scores_risk.pdf"))
mycol <- rep("darkgrey", length(nn))
ccc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[ccc == 0] <- colo[1]
mycol[ccc == 2] <- colo[3]
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col=mycol, main="STROMACD10\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="STROMACD10 scores (AFFY)", ylab="STROMACD10 scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
legend("topleft", col=c(colo[1], colo[3], "darkgrey"), legend=c("Low-risk", "High-risk", "Discordance"), bty="n", pch=16)
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

##############
## TOP2AINDEX
##############

## read published signature
sigg <- read.csv(file.path("sigs", "desmedt2011_top2aindex_signature.csv"), stringsAsFactor=FALSE)

ss <- genefu::sig.score(x=sigg, data=datac.affy, annot=annotc.affy, do.mapping=TRUE)$score
rr <- as.numeric(cut2(x=ss, cuts=median(ss))) - 1
names(rr) <- names(ss)
sig.affy <- list("score"=ss, "risk"=rr)
ss <- genefu::sig.score(x=sigg, data=datac.rnaseq, annot=annotc.rnaseq, do.mapping=TRUE)$score
rr <- as.numeric(cut2(x=ss, cuts=median(ss))) - 1
names(rr) <- names(ss)
sig.rnaseq <- list("score"=ss, "risk"=rr)

tt <- table("AFFY"=sig.affy$risk[nn], "RNASEQ"=sig.rnaseq$risk[nn])
attributes(tt)$dimnames$AFFY <- c("Low-risk", "High-risk")
attributes(tt)$dimnames$RNASEQ <- c("Low-risk", "High-risk")
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "top2aindex_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for TOP2AINDEX risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.3g, 95%%CI [%.3g,%.3g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

# res.tabs <- c(res.tabs, list("TOP2AINDEX"=tt))
# res.risks <- c(res.risks, list("TOP2AINDEX"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))

cc <- cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")
ttc <- spearmanCI(x=cc, n=sum(complete.cases(sig.affy$score[nn], sig.rnaseq$score[nn])), alpha=0.05)
pdf(file.path(saveres, "top2aindex_scores.pdf"))
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col="darkgrey", main="TOP2AINDEX\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="TOP2AINDEX scores (AFFY)", ylab="TOP2AINDEX scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

ttt <- c("rho"=cc, "lower"=ttc["lower"], "upper"=ttc["upper"], "p"=ttc["p.value"])
names(ttt) <- c("rho", "lower", "upper", "p")
# res.scores <- c(res.scores, list("TOP2AINDEX"=ttt))

pdf(file.path(saveres, "top2aindex_scores_risk.pdf"))
mycol <- rep("darkgrey", length(nn))
ccc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[ccc == 0] <- colo[1]
mycol[ccc == 2] <- colo[3]
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col=mycol, main="TOP2AINDEX\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="TOP2AINDEX scores (AFFY)", ylab="TOP2AINDEX scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
legend("topleft", col=c(colo[1], colo[3], "darkgrey"), legend=c("Low-risk", "High-risk", "Discordance"), bty="n", pch=16)
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

##############
## ASCORE
##############

## ASCORE
top2aindex <- read.csv(file.path("sigs", "desmedt2011_top2aindex_signature.csv"), stringsAsFactor=FALSE)
## compute the A-SCORE in affy
ss.top2aindex <- genefu::sig.score(x=top2aindex, data=datac.affy, annot=annotc.affy, do.mapping=TRUE)$score
ss.plau <- genefu::sig.score(x=mod1$PLAU, data=data.affy, annot=annot.affy, do.mapping=FALSE)$score
ss.stat1 <- genefu::sig.score(x=mod1$STAT1, data=data.affy, annot=annot.affy, do.mapping=FALSE)$score
erbb2p.proba <- genefu::bimod(x=mod1$ERBB2[1, , drop=FALSE], data=data.affy, annot=annot.affy, do.mapping=FALSE, model="V")$status1.proba
## -> HER2- : STAT1MODULE - PLAUMODULE
## -> HER2+ : ASCORE + STAT1MODULE - PLAUMODULE
scoretn <- (genefu::rescale(ss.stat1 - ss.plau, q = 0.05, na.rm = TRUE) - 0.5) * 2
scoretp <- (genefu::rescale(ss.top2aindex + ss.stat1 - ss.plau, q = 0.05, na.rm = TRUE) - 0.5) * 2
ss <- (1 - erbb2p.proba) * scoretn + erbb2p.proba * scoretp
rr <- as.numeric(cut2(x=ss, cuts=median(ss))) - 1
names(rr) <- names(ss)
sig.affy <- list("score"=ss, "risk"=rr)
## compute the A-SCORE in rnaseq
ss.top2aindex <- genefu::sig.score(x=top2aindex, data=datac.rnaseq, annot=annotc.rnaseq, do.mapping=TRUE)$score
ss.plau <- genefu::sig.score(x=mod1$PLAU, data=datac.rnaseq, annot=annotc.rnaseq, do.mapping=TRUE)$score
ss.stat1 <- genefu::sig.score(x=mod1$STAT1, data=datac.rnaseq, annot=annotc.rnaseq, do.mapping=TRUE)$score
erbb2p.proba <- genefu::bimod(x=mod1$ERBB2[1, , drop=FALSE], data=datac.rnaseq, annot=annotc.rnaseq, do.mapping=TRUE, model="V")$status1.proba
## -> HER2- : STAT1MODULE - PLAUMODULE
## -> HER2+ : ASCORE + STAT1MODULE - PLAUMODULE
scoretn <- (genefu::rescale(ss.stat1 - ss.plau, q = 0.05, na.rm = TRUE) - 0.5) * 2
scoretp <- (genefu::rescale(ss.top2aindex + ss.stat1 - ss.plau, q = 0.05, na.rm = TRUE) - 0.5) * 2
ss <- (1 - erbb2p.proba) * scoretn + erbb2p.proba * scoretp
rr <- as.numeric(cut2(x=ss, cuts=median(ss))) - 1
names(rr) <- names(ss)
sig.rnaseq <- list("score"=ss, "risk"=rr)

tt <- table("AFFY"=sig.affy$risk[nn], "RNASEQ"=sig.rnaseq$risk[nn])
attributes(tt)$dimnames$AFFY <- c("Low-risk", "High-risk")
attributes(tt)$dimnames$RNASEQ <- c("Low-risk", "High-risk")
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "ascore_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for ASCORE risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.3g, 95%%CI [%.3g,%.3g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

res.tabs <- c(res.tabs, list("ASCORE"=tt))
res.risks <- c(res.risks, list("ASCORE"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))

cc <- cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")
ttc <- spearmanCI(x=cc, n=sum(complete.cases(sig.affy$score[nn], sig.rnaseq$score[nn])), alpha=0.05)
pdf(file.path(saveres, "ascore_scores.pdf"))
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col="darkgrey", main="ASCORE\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="ASCORE scores (AFFY)", ylab="ASCORE scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

ttt <- c("rho"=cc, "lower"=ttc["lower"], "upper"=ttc["upper"], "p"=ttc["p.value"])
names(ttt) <- c("rho", "lower", "upper", "p")
res.scores <- c(res.scores, list("ASCORE"=ttt))

pdf(file.path(saveres, "ascore_scores_risk.pdf"))
mycol <- rep("darkgrey", length(nn))
ccc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[ccc == 0] <- colo[1]
mycol[ccc == 2] <- colo[3]
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col=mycol, main="ASCORE\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="ASCORE scores (AFFY)", ylab="ASCORE scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
legend("topleft", col=c(colo[1], colo[3], "darkgrey"), legend=c("Low-risk", "High-risk", "Discordance"), bty="n", pch=16)
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

##############
## Signatures in Ignatiadis et al., JCO, 2012
##############

for(i in 1:length(gene.modules)) {

## read published signature
sigg <- gene.modules[[i]]

ss <- genefu::sig.score(x=sigg, data=datac.affy, annot=annotc.affy, do.mapping=TRUE)$score
rr <- as.numeric(cut2(x=ss, cuts=median(ss))) - 1
names(rr) <- names(ss)
sig.affy <- list("score"=ss, "risk"=rr)
ss <- genefu::sig.score(x=sigg, data=datac.rnaseq, annot=annotc.rnaseq, do.mapping=TRUE)$score
rr <- as.numeric(cut2(x=ss, cuts=median(ss))) - 1
names(rr) <- names(ss)
sig.rnaseq <- list("score"=ss, "risk"=rr)

tt <- table("AFFY"=sig.affy$risk[nn], "RNASEQ"=sig.rnaseq$risk[nn])
attributes(tt)$dimnames$AFFY <- c("Low-risk", "High-risk")
attributes(tt)$dimnames$RNASEQ <- c("Low-risk", "High-risk")
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "stromacd10_risk_contingency_table.pdf"))
plot(t(tt), main=sprintf("Contingency table for %s risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", names(gene.modules)[i]), sub=sprintf("Kappa coefficient=%.3g, 95%%CI [%.3g,%.3g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

res.tabs <- c(res.tabs, list(tt))
res.risks <- c(res.risks, list(c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))
names(res.tabs)[length(res.tabs)] <- names(res.risks)[length(res.risks)] <- names(gene.modules)[i]

cc <- cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")
ttc <- spearmanCI(x=cc, n=sum(complete.cases(sig.affy$score[nn], sig.rnaseq$score[nn])), alpha=0.05)
pdf(file.path(saveres, "stromacd10_scores.pdf"))
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col="darkgrey", main=sprintf("%s\nAFFYMETRIX vs ILLUMINA RNA-seq", names(gene.modules)[i]), xlab=sprintf("%s scores (AFFY)", names(gene.modules)[i]), ylab=sprintf("%s scores (ILLUMINA RNA-seq)", names(gene.modules)[i]), sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

ttt <- c("rho"=cc, "lower"=ttc["lower"], "upper"=ttc["upper"], "p"=ttc["p.value"])
names(ttt) <- c("rho", "lower", "upper", "p")
res.scores <- c(res.scores, list(ttt))
names(res.scores)[length(res.scores)] <- names(gene.modules)[i]

pdf(file.path(saveres, sprintf("%s_scores_risk.pdf", tolower(names(gene.modules)))))
mycol <- rep("darkgrey", length(nn))
ccc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[ccc == 0] <- colo[1]
mycol[ccc == 2] <- colo[3]
plot(x=sig.affy$score[nn], y=sig.rnaseq$score[nn], pch=16, col=mycol, main=sprintf("%s\nAFFYMETRIX vs ILLUMINA RNA-seq", names(gene.modules)[i]), xlab=sprintf("%s scores (AFFY)", names(gene.modules)[i]), ylab=sprintf("%s scores (ILLUMINA RNA-seq)", names(gene.modules)[i]), sub=sprintf("Spearman correlation: %.3g, 95%%CI [%.3g,%.3g], p=%.1E", cc, ttc["lower"], ttc["upper"], ttc["p.value"]))
legend("topleft", col=c(colo[1], colo[3], "darkgrey"), legend=c("Low-risk", "High-risk", "Discordance"), bty="n", pch=16)
mm <- lm(sig.rnaseq$score[nn] ~ sig.affy$score[nn])
abline(mm)
dev.off()

}

##############
## barplot
##############

## reorder signatures
# res.scores <- res.scores[c(7, 8, 1, 5, 6, 3, 4, 9, 10, 2, 11, 12)]

## barplot for correlations for each signature
pdf(file.path(saveres, "barplot_scores_sigs.pdf"), height=7, width=7)
par(las=3, mar=c(8, 5, 5, 2))
xx <- sapply(res.scores, function(x) { return(x["rho"]) })
ll <- sapply(res.scores, function(x) { return(x["lower"]) })
ll[!is.na(ll) & ll < -1] <- -1
uu <- sapply(res.scores, function(x) { return(x["upper"]) })
uu[!is.na(uu) & uu > 1] <- 1
names(xx) <- names(ll) <- names(uu) <- names(res.scores)
oo <- order(xx, decreasing=TRUE)
mycol <- rainbow(length(xx), v=0.9)
names(mycol) <- names(xx)[oo]
co <- barplot(height=xx[oo], space=0.3, col=mycol, ylab="Spearman rho", ylim=c(0,1))
plotCI(x=co, y=xx[oo], li=ll[oo], ui=uu[oo], err="y", pch=".", add=TRUE)
dev.off()


## barplot for kappa coefficient for each signature
pdf(file.path(saveres, "barplot_risks_sigs.pdf"), height=7, width=7)
par(las=3, mar=c(8, 5, 5, 2))
xx <- sapply(res.risks, function(x) { return(x["kappa"]) })
ll <- sapply(res.risks, function(x) { return(x["lower"]) })
ll[!is.na(ll) & ll < -1] <- -1
uu <- sapply(res.risks, function(x) { return(x["upper"]) })
uu[!is.na(uu) & uu > 1] <- 1
names(xx) <- names(ll) <- names(uu) <- names(res.risks)
oo <- order(xx, decreasing=TRUE)
co <- barplot(height=xx[oo], space=0.3, col=mycol[names(xx)[oo]], ylab="Kappa coefficient", ylim=c(0,1))
plotCI(x=co, y=xx[oo], li=ll[oo], ui=uu[oo], err="y", pch=".", add=TRUE)
dev.off()


tt <- lapply(res.tabs, function(x) {
  class(x) <- "matrix"
  x <- data.frame(x)
  return(x)
  })

tt <- matrix("", nrow=sum(sapply(res.tabs, nrow)) + 4 * length(res.tabs), ncol=max(sapply(res.tabs, ncol)) + 3)
ii <- 1
for (i in 1:length(res.tabs)) {
  tt[ii + 1, 1] <- names(res.tabs)[i]
  tt[ii + 3, 1] <- names(dimnames(res.tabs[[i]]))[1]
  tt[ii + 1, 3] <- names(dimnames(res.tabs[[i]]))[2]
  tt[(ii + 3):(ii + nrow(res.tabs[[i]]) + 2), 3:(ncol(res.tabs[[i]]) + 2)] <- res.tabs[[i]]
  tt[(ii + 3):(ii + nrow(res.tabs[[i]]) + 2), 2] <- rownames(res.tabs[[i]])
  tt[ii + 2, 3:(ncol(res.tabs[[i]]) + 2)] <- colnames(res.tabs[[i]])
  ii <- ii + nrow(res.tabs[[i]]) + 4
}
colnames(tt) <- rep("", ncol(tt))
write.csv(tt, file=file.path(saveres, "sigs_risk_tables.csv"), row.names=FALSE)
  
  
save(list=c("res.scores", "res.risks"), compress=TRUE, file=file.path(saveres, "risk_score_sigs.RData"))







## end

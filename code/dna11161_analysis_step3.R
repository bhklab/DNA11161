########################
## Benjamin Haibe-Kains
##
## Nov 23, 2012
########################

rm(list=ls(all=TRUE))

source(file.path("code", "ufoo.R"))

saveres <- file.path("saveres")
if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE) }

library(genefu)
library(vcd)
library(epibasix)
library(Hmisc)
  
########################
## step 3: Compare NON Affymetrix-based breast cancer gene expression signatures computed from Affymetrix and Illumina RNA-seq data
########################

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

res.risks <- res.tabs <- NULL

#################################################
## subtype classification
#################################################

sbts <- NULL

sbtn <- c("ER-/HER2-", "HER2+", "ER+/HER2- High Prolif", "ER+/HER2- Low Prolif")
sbtnn <- c("Basal", "Her2", "LumB", "LumA", "Normal")

##############
## SCMGENE
##############

pdf(file.path(saveres, "scmgene_sbt_plot_affy.pdf"))
sig.affy <- subtype.cluster.predict(sbt.model=scmgene.robust, data=data.affy, annot=annot.affy, do.mapping=FALSE, plot=TRUE)
title(main="SCMGENE on AFFYMETRIX")
sig.affy$subtype2 <- factor(sig.affy$subtype2, levels=sbtn)
dev.off()
pdf(file.path(saveres, "scmgene_sbt_plot_rnaseq.pdf"))
sig.rnaseq <- subtype.cluster.predict(sbt.model=scmgene.robust, data=data.rnaseq, annot=annot.rnaseq, do.mapping=TRUE, plot=TRUE)
title(main="SCMGENE on ILLUMINA RNA-seq")
sig.rnaseq$subtype2 <- factor(sig.rnaseq$subtype2, levels=sbtn)
dev.off()

sbts <- cbind(sbts, "SCMGENE.AFFY"=as.character(sig.affy$subtype2[nn]), "SCMGENE.RNASEQ"=as.character(sig.rnaseq$subtype2[nn]))

tt <- table("AFFY"=sig.affy$subtype2[nn], "RNASEQ"=sig.rnaseq$subtype2[nn])
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "scmgene_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for SCMGENE subtype classification\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g, 95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

res.tabs <- c(res.tabs, list("SCMGENE"=tt))
res.risks <- c(res.risks, list("SCMGENE"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))

##############
## SCMOD1
##############

pdf(file.path(saveres, "scmod1_sbt_plot_affy.pdf"))
sig.affy <- subtype.cluster.predict(sbt.model=scmod1.robust, data=data.affy, annot=annot.affy, do.mapping=FALSE, plot=TRUE)
title(main="SCMOD1 on AFFYMETRIX")
sig.affy$subtype2 <- factor(sig.affy$subtype2, levels=sbtn)
dev.off()
pdf(file.path(saveres, "scmod1_sbt_plot_rnaseq.pdf"))
sig.rnaseq <- subtype.cluster.predict(sbt.model=scmod1.robust, data=data.rnaseq, annot=annot.rnaseq, do.mapping=TRUE, plot=TRUE)
title(main="SCMOD1 on ILLUMINA RNA-seq")
sig.rnaseq$subtype2 <- factor(sig.rnaseq$subtype2, levels=sbtn)
dev.off()

sbts <- cbind(sbts, "SCMOD1.AFFY"=as.character(sig.affy$subtype2[nn]), "SCMOD1.RNASEQ"=as.character(sig.rnaseq$subtype2[nn]))

tt <- table("AFFY"=sig.affy$subtype2[nn], "RNASEQ"=sig.rnaseq$subtype2[nn])
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "scmod1_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for SCMOD1 subtype classification\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g, 95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

res.tabs <- c(res.tabs, list("SCMOD1"=tt))
res.risks <- c(res.risks, list("SCMOD1"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))

##############
## SCMOD2
##############

pdf(file.path(saveres, "scmod2_sbt_plot_affy.pdf"))
sig.affy <- subtype.cluster.predict(sbt.model=scmod2.robust, data=data.affy, annot=annot.affy, do.mapping=FALSE, plot=TRUE)
title(main="SCMOD2 on AFFYMETRIX")
sig.affy$subtype2 <- factor(sig.affy$subtype2, levels=sbtn)
dev.off()
pdf(file.path(saveres, "scmod2_sbt_plot_rnaseq.pdf"))
sig.rnaseq <- subtype.cluster.predict(sbt.model=scmod2.robust, data=data.rnaseq, annot=annot.rnaseq, do.mapping=TRUE, plot=TRUE)
title(main="SCMOD2 on ILLUMINA RNA-seq")
sig.rnaseq$subtype2 <- factor(sig.rnaseq$subtype2, levels=sbtn)
dev.off()

sbts <- cbind(sbts, "SCMOD2.AFFY"=as.character(sig.affy$subtype2[nn]), "SCMOD2.RNASEQ"=as.character(sig.rnaseq$subtype2[nn]))

tt <- table("AFFY"=sig.affy$subtype2[nn], "RNASEQ"=sig.rnaseq$subtype2[nn])
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "scmod2_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for SCMOD2 subtype classification\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g, 95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

res.tabs <- c(res.tabs, list("SCMOD2"=tt))
res.risks <- c(res.risks, list("SCMOD2"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))


##############
## non affymetrix models
##############

## load data
load(file.path(saveres, "dna11161_common_jetset.RData"))
annotc.affy <- annotc
annotc.rnaseq <- annotc

## tumors in common
nn <- intersect(rownames(datac.affy), rownames(datac.rnaseq))

##############
## PAM50
##############
sig.affy <- intrinsic.cluster.predict(sbt.model=pam50.robust, data=datac.affy, annot=annotc.affy, do.mapping=TRUE)
sig.affy$subtype <- factor(sig.affy$subtype, levels=sbtnn)
sig.rnaseq <- intrinsic.cluster.predict(sbt.model=pam50.robust, data=datac.rnaseq, annot=annotc.rnaseq, do.mapping=TRUE)
sig.rnaseq$subtype <- factor(sig.rnaseq$subtype, levels=sbtnn)

sbts <- cbind(sbts, "PAM50.AFFY"=as.character(sig.affy$subtype[nn]), "PAM50.RNASEQ"=as.character(sig.rnaseq$subtype[nn]))

tt <- table("AFFY"=sig.affy$subtype[nn], "RNASEQ"=sig.rnaseq$subtype[nn])
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "pam50_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for PAM50 subtype classification\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g, 95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

res.tabs <- c(res.tabs, list("PAM50"=tt))
res.risks <- c(res.risks, list("PAM50"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))

##############
## SSP2006
##############
sig.affy <- intrinsic.cluster.predict(sbt.model=ssp2006.robust, data=datac.affy, annot=annotc.affy, do.mapping=TRUE)
sig.affy$subtype <- factor(sig.affy$subtype, levels=sbtnn)
sig.rnaseq <- intrinsic.cluster.predict(sbt.model=ssp2006.robust, data=datac.rnaseq, annot=annotc.rnaseq, do.mapping=TRUE)
sig.rnaseq$subtype <- factor(sig.rnaseq$subtype, levels=sbtnn)

sbts <- cbind(sbts, "SSP2006.AFFY"=as.character(sig.affy$subtype[nn]), "SSP2006.RNASEQ"=as.character(sig.rnaseq$subtype[nn]))

tt <- table("AFFY"=sig.affy$subtype[nn], "RNASEQ"=sig.rnaseq$subtype[nn])
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "ssp2006_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for SSP2006 subtype classification\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g, 95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

res.tabs <- c(res.tabs, list("SSP2006"=tt))
res.risks <- c(res.risks, list("SSP2006"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))

##############
## SSP2003
##############
sig.affy <- intrinsic.cluster.predict(sbt.model=ssp2003.robust, data=datac.affy, annot=annotc.affy, do.mapping=TRUE)
sig.affy$subtype <- factor(sig.affy$subtype, levels=sbtnn)
sig.rnaseq <- intrinsic.cluster.predict(sbt.model=ssp2003.robust, data=datac.rnaseq, annot=annotc.rnaseq, do.mapping=TRUE)
sig.rnaseq$subtype <- factor(sig.rnaseq$subtype, levels=sbtnn)

sbts <- cbind(sbts, "SSP2003.AFFY"=as.character(sig.affy$subtype[nn]), "SSP2003.RNASEQ"=as.character(sig.rnaseq$subtype[nn]))

tt <- table("AFFY"=sig.affy$subtype[nn], "RNASEQ"=sig.rnaseq$subtype[nn])
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "ssp2003_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for SSP2003 subtype classification\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g, 95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

res.tabs <- c(res.tabs, list("SSP2003"=tt))
res.risks <- c(res.risks, list("SSP2003"=c("kappa"=kk$kappa, "lower"=kk$CIL, "upper"=kk$CIU, "p"=tts$chisq_tests[1,3])))

write.csv(sbts, file=file.path(saveres, "subtype_classifs_affy_rnaseq.csv"))


##############
## barplot
##############

## barplot for kappa coefficient for each signature
pdf(file.path(saveres, "barplot_subtypes_sigs.pdf"), height=7, width=4)
par(las=3, mar=c(8, 5, 5, 2))
xx <- sapply(res.risks, function(x) { return(x["kappa"]) })
ll <- sapply(res.risks, function(x) { return(x["lower"]) })
ll[!is.na(ll) & ll < -1] <- -1
uu <- sapply(res.risks, function(x) { return(x["upper"]) })
uu[!is.na(uu) & uu > 1] <- 1
names(xx) <- names(ll) <- names(uu) <- names(res.risks)
co <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Kappa coefficient", ylim=c(0,1))
plotCI(x=co, y=xx, li=ll, ui=uu, err="y", pch=".", add=TRUE)
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
write.csv(tt, file=file.path(saveres, "sigs_subtypes_tables.csv"), row.names=FALSE)
  




## end
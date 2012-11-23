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
## step 3: Compare NON Affymetrix-based breast cancer gene expression signatures computed from Affymetrix and Illumina RNA-seq data
########################

## load data
#myfns <- c("jetset"=file.path(saveres, "dna11161_jetset_common.RData"), "bestpgene"=file.path(saveres, "dna11161_bestpgene_common.RData"), "bestptranscript"=file.path(saveres, "dna11161_bestptranscript_common.RData"))

##############
## MAMMAPRINT
##############
sig.affy <- gene70(data=data.affy, annot=annot.affy, do.mapping=TRUE)
sig.rnaseq <- gene70(data=data.rnaseq, annot=annot.rnaseq, do.mapping=TRUE)

tt <- table("AFFY"=sig.affy$risk[nn], "RNASEQ"=sig.rnaseq$risk[nn])
attributes(tt)$dimnames$AFFY <- c("Low-risk", "High-risk")
attributes(tt)$dimnames$RNASEQ <- c("Low-risk", "High-risk")
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "mammaprint_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for MAMMAPRINT risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

pdf(file.path(saveres, "mammaprint_scores.pdf"))
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col="darkgrey", main="MAMMAPRINT\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="MAMMAPRINT scores (AFFY)", ylab="MAMMAPRINT scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
dev.off()

pdf(file.path(saveres, "mammaprint_scores_risk.pdf"))
mycol <- rep("red", length(nn))
cc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[cc == 0] <- "blue"
mycol[cc == 2] <- "yellow"
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col=mycol, main="MAMMAPRINT\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="MAMMAPRINT scores (AFFY)", ylab="MAMMAPRINT scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
legend("topleft", col=c("blue", "yellow", "red"), legend=c("Low-risk", "High-risk", "Discordance"), bty="n", pch=16)
dev.off()

##############
## ONCOTYPE DX
##############
sig.oncotypedx.save <- sig.oncotypedx
sig.oncotypedx["GSTM1", "EntrezGene.ID"] <- "2948"
## the affymetrix probes for GSTM1 are ambiguious with GSTM2, GSTM3, and GSTM4; GSTM3 is present in the data
## http://xavierlab2.mgh.harvard.edu/EnrichmentProfiler/primary/Enrichment/204550_x_at.html
sig.affy <- oncotypedx(data=data.affy, annot=annot.affy, do.mapping=TRUE)
sig.rnaseq <- oncotypedx(data=data.rnaseq, annot=annot.rnaseq, do.mapping=TRUE)
sig.oncotypedx <- sig.oncotypedx.save

tt <- table("AFFY"=sig.affy$risk[nn], "RNASEQ"=sig.rnaseq$risk[nn])
attributes(tt)$dimnames$AFFY <- c("Low-risk", "Intermediate-risk", "High-risk")
attributes(tt)$dimnames$RNASEQ <- c("Low-risk", "Intermediate-risk", "High-risk")
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "oncotypedx_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for ONCOTYPE DX risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

pdf(file.path(saveres, "oncotypedx_scores.pdf"))
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col="darkgrey", main="ONCOTYPE DX\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="ONCOTYPE DX scores (AFFY)", ylab="ONCOTYPE DX scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
dev.off()

pdf(file.path(saveres, "oncotypedx_scores_risk.pdf"))
mycol <- rep("red", length(nn))
cc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[complete.cases(sig.affy$risk[nn], sig.rnaseq$risk[nn]) & sig.affy$risk[nn] == 0 & sig.rnaseq$risk[nn] == 0] <- "blue"
mycol[complete.cases(sig.affy$risk[nn], sig.rnaseq$risk[nn]) & sig.affy$risk[nn] == 0.5 & sig.rnaseq$risk[nn] == 0.5] <- "orange"
mycol[complete.cases(sig.affy$risk[nn], sig.rnaseq$risk[nn]) & sig.affy$risk[nn] == 1 & sig.rnaseq$risk[nn] == 1] <- "yellow"

plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col=mycol, main="ONCOTYPE DX\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="ONCOTYPE DX scores (AFFY)", ylab="ONCOTYPE DX scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
legend("topleft", col=c("blue", "orange", "yellow", "red"), legend=c("Low-risk", "Intermediate-risk", "High-risk", "Discordance"), bty="n", pch=16)
dev.off()

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
plot(t(tt), main="Contingency table for SCMGENE subtype classification\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

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
plot(t(tt), main="Contingency table for SCMOD1 subtype classification\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

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
plot(t(tt), main="Contingency table for SCMOD2 subtype classification\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

##############
## PAM50
##############
sig.affy <- intrinsic.cluster.predict(sbt.model=pam50.robust, data=data.affy, annot=annot.affy, do.mapping=TRUE)
sig.affy$subtype <- factor(sig.affy$subtype, levels=sbtnn)
sig.rnaseq <- intrinsic.cluster.predict(sbt.model=pam50.robust, data=data.rnaseq, annot=annot.rnaseq, do.mapping=TRUE)
sig.rnaseq$subtype <- factor(sig.rnaseq$subtype, levels=sbtnn)

sbts <- cbind(sbts, "PAM50.AFFY"=as.character(sig.affy$subtype[nn]), "PAM50.RNASEQ"=as.character(sig.rnaseq$subtype[nn]))

tt <- table("AFFY"=sig.affy$subtype[nn], "RNASEQ"=sig.rnaseq$subtype[nn])
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "pam50_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for PAM50 subtype classification\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

##############
## SSP2006
##############
sig.affy <- intrinsic.cluster.predict(sbt.model=ssp2006.robust, data=data.affy, annot=annot.affy, do.mapping=TRUE)
sig.affy$subtype <- factor(sig.affy$subtype, levels=sbtnn)
sig.rnaseq <- intrinsic.cluster.predict(sbt.model=ssp2006.robust, data=data.rnaseq, annot=annot.rnaseq, do.mapping=TRUE)
sig.rnaseq$subtype <- factor(sig.rnaseq$subtype, levels=sbtnn)

sbts <- cbind(sbts, "SSP2006.AFFY"=as.character(sig.affy$subtype[nn]), "SSP2006.RNASEQ"=as.character(sig.rnaseq$subtype[nn]))

tt <- table("AFFY"=sig.affy$subtype[nn], "RNASEQ"=sig.rnaseq$subtype[nn])
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "ssp2006_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for SSP2006 subtype classification\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

##############
## SSP2003
##############
sig.affy <- intrinsic.cluster.predict(sbt.model=ssp2003.robust, data=data.affy, annot=annot.affy, do.mapping=TRUE)
sig.affy$subtype <- factor(sig.affy$subtype, levels=sbtnn)
sig.rnaseq <- intrinsic.cluster.predict(sbt.model=ssp2003.robust, data=data.rnaseq, annot=annot.rnaseq, do.mapping=TRUE)
sig.rnaseq$subtype <- factor(sig.rnaseq$subtype, levels=sbtnn)

sbts <- cbind(sbts, "SSP2003.AFFY"=as.character(sig.affy$subtype[nn]), "SSP2003.RNASEQ"=as.character(sig.rnaseq$subtype[nn]))

tt <- table("AFFY"=sig.affy$subtype[nn], "RNASEQ"=sig.rnaseq$subtype[nn])
print(tt)
tts <- assocstats(tt)
print(tts)
kk <- epiKappa(tt, k0=0)
print(kk)
pdf(file.path(saveres, "ssp2003_risk_contingency_table.pdf"))
plot(t(tt), main="Contingency table for SSP2003 subtype classification\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

write.csv(sbts, file=file.path(saveres, "subtype_classifs_affy_rnaseq.csv"))
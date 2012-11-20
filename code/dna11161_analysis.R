########################
## Benjamin Haibe-Kains
##
## April 21, 2012
########################

rm(list=ls())

source("ufoo.R")

library(genefu)

saveres <- file.path("..", "saveres")
if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE) }

#################################################
## microarray data

myfn <- file.path(saveres, "dna11161_affy_frma.RData")
if(!file.exists(myfn)) {
library(affy)
library(frma)
library(frmaTools)
## CEL file names
rdp <- "../affy"
celfn <- list.celfiles(rdp, full.names=TRUE)
celfns <- list.celfiles(rdp, full.names=FALSE)
## experiments' names
names(celfn) <- names(celfns) <- paste(gsub("_[(]HG-U133_PLUS_2[)].CEL.GZ", "", toupper(list.celfiles(rdp))), "T", sep="_")
## chip type and date
chipt <- sapply(celfn, celfileChip)
chipd <- t(sapply(celfn, celfileDateHour))
sbt <- sapply(strsplit(names(celfn), "_"), function(x) { return(x[[1]]) })
sampleinfo <- data.frame("samplename"=names(celfn), "subtype"=sbt, "chiptype"=chipt, "hybridization.date"=chipd[ ,"day"], "hybridization.hour"=chipd[ ,"hour"])
rownames(sampleinfo) <- names(celfn)
## frma normalization
celfnt <- celfn
ns <- 60
nsu <- ceiling(length(celfnt) / ns)
data <- NULL
for(i in 1:nsu) {
message(sprintf("\nProcessing set %i/%i ...", i, nsu))
if(i < nsu) { iix <- (((i-1) * ns)+1):(i*ns) 
} else { iix <- (((i-1) * ns)+1):length(celfnt) }
## fRMA
tt <- celfnt[iix]
names(tt) <- NULL
abatch <- read.affybatch(filenames=tt)
rr <- frma(object=abatch, summarize="robust_weighted_average", verbose=TRUE)
rr2 <- exprs(rr)
colnames(rr2) <- names(celfnt)[iix]
data <- rbind(data, t(rr2))
}
## build annotation matrix
library(biomaRt)
ensembl.db <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
#write.csv(listAttributes(mart=ensembl.db), row.names=FALSE, file="mart_attributes.csv")
# gene.an <- getBM(attributes=c("affy_hg_u133_plus_2", "ensembl_gene_id", "hgnc_symbol", "entrezgene", "unigene", "description", "chromosome_name", "start_position", "end_position", "strand", "band"), filters="affy_hg_u133_plus_2", values=colnames(data), mart=ensembl.db)
# gene.an[gene.an == "" | gene.an == " "] <- NA
# gene.an <- gene.an[!is.na(gene.an[ , "affy_hg_u133_plus_2"]) & !duplicated(gene.an[ , "affy_hg_u133_plus_2"]), , drop=FALSE]
# annot <- data.frame(matrix(NA, nrow=ncol(data), ncol=ncol(gene.an)+2, dimnames=list(colnames(data), c("probe", "EntrezGene.ID", colnames(gene.an)))))
# annot[match(gene.an[ , "affy_hg_u133_plus_2"], colnames(data)), colnames(gene.an)] <- gene.an
# annot[ ,"probe"] <- colnames(data)
# annot[ , "EntrezGene.ID"] <- annot[ ,"entrezgene"]
## select the best probe for a single gene
js <- jscores(chip="hgu133plus2", probeset=colnames(data))
js <- js[colnames(data), , drop=FALSE]
## identify the best probeset for each entrez gene id
geneid1 <- as.character(js[ ,"EntrezID"])
names(geneid1) <- rownames(js)
geneid2 <- sort(unique(geneid1))
names(geneid2) <- paste("geneid", geneid2, sep=".")
gix1 <- !is.na(geneid1)
gix2 <- !is.na(geneid2)
geneid.common <- intersect(geneid1[gix1], geneid2[gix2])
## probes corresponding to common gene ids
gg <- names(geneid1)[is.element(geneid1, geneid.common)]
gid <- geneid1[is.element(geneid1, geneid.common)]
## duplicated gene ids
gid.dupl <- unique(gid[duplicated(gid)])
gg.dupl <- names(geneid1)[is.element(geneid1, gid.dupl)]
## unique gene ids
gid.uniq <- gid[!is.element(gid, gid.dupl)]
gg.uniq <- names(geneid1)[is.element(geneid1, gid.uniq)]
## data for duplicated gene ids
if(length(gid.dupl) > 0) {
	## us jetset oevrall score to select the best probeset
	myscore <- js[gg.dupl,"overall"]
	myscore <- cbind("probe"=gg.dupl, "gid"=geneid1[gg.dupl], "score"=myscore)
	myscore <- myscore[order(as.numeric(myscore[ , "score"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
	myscore <- myscore[!duplicated(myscore[ , "gid"]), , drop=FALSE]
}
js <- data.frame(js, "best"=FALSE)
js[gg.uniq, "best"] <- TRUE
js[myscore[ ,"probe"], "best"] <- TRUE
## more annotations from biomart
ugid <- sort(unique(js[ ,"EntrezID"]))
gene.an <- getBM(attributes=c("entrezgene", "ensembl_gene_id", "hgnc_symbol", "unigene", "description", "chromosome_name", "start_position", "end_position", "strand", "band"), filters="entrezgene", values=ugid, mart=ensembl.db)
gene.an[gene.an == "" | gene.an == " "] <- NA
gene.an <- gene.an[!is.na(gene.an[ , "entrezgene"]) & !duplicated(gene.an[ , "entrezgene"]), , drop=FALSE]
annot <- data.frame(matrix(NA, nrow=ncol(data), ncol=ncol(gene.an)+1, dimnames=list(colnames(data), c("probe", colnames(gene.an)))))
annot[match(gene.an[ , "entrezgene"], js[ ,"EntrezID"]), colnames(gene.an)] <- gene.an
annot[ ,"probe"] <- colnames(data)
colnames(js)[colnames(js) != "best"] <- paste("jetset", colnames(js)[colnames(js) != "best"], sep=".")
annot <- data.frame(annot, "EntrezGene.ID"=js[ ,"jetset.EntrezID"], js)
## rename objects
data.affy <- data
annot.affy <- annot[colnames(data), ,drop=FALSE]
sampleinfo.affy <- sampleinfo[rownames(data), ,drop=FALSE]
rm(list=c("data", "annot", "sampleinfo", "rr2", "abatch"))
gc()
save(list=c("data.affy", "annot.affy", "sampleinfo.affy"), compress=TRUE, file=myfn)
} else { load(myfn) }



#################################################
## RNA-seq data
myfn <- file.path(saveres, "dna11161_rnaseq.RData")
if(!file.exists(myfn)) {
## sample info
sampleinfo <- read.csv("../rnaseq/DNA111061_dnavision_ids.csv", row.names=NULL, stringsAsFactors=FALSE)
rownames(sampleinfo) <- gsub("BIS", "", toupper(as.character(sampleinfo[ ,2])))
sampleinfo <- data.frame("samplename"=rownames(sampleinfo), sampleinfo)
## alignment with tophat and expression estimation (FPKM) with cufflink
dd <- read.csv("../rnaseq/consolidated/allFpkmGenesEnsembl.txt", sep="\t", row.names=NULL, stringsAsFactors=FALSE)
#ddf <- read.csv("../rnaseq/consolidated/allFpkmGenesEnsemblStatut.txt", sep="\t", row.names=NULL, stringsAsFactors=FALSE)
#ddf <- ddf[match(dd[ ,1], ddf[ ,1]), match(colnames(dd), colnames(ddf))]
annott <- dd[ ,1:2]
ddf <- dd[ , grep("Status", colnames(dd)), drop=FALSE]
dd <- data.matrix(dd[ , grep("FPKM", colnames(dd)), drop=FALSE])
rownames(dd) <- rownames(ddf) <- as.character(annott[ ,1])
data <- t(dd)
dataf <- t(ddf)
## remove FPKMs of bad quality
data[dataf != "OK"] <- NA
## log2 transformation of the FPKM for better comparability with microarray
data <- log2(data+1)
nn <- gsub("[.]", "-", substr(colnames(dd), 1, nchar(colnames(dd)) - 5))
#sampleinfo <- sampleinfo[sample(1:nrow(sampleinfo)), ]
sampleinfo <- sampleinfo[match(nn, sampleinfo[ ,"DNAvision.ID"]), , drop=FALSE]
rownames(data) <- rownames(dataf) <- rownames(sampleinfo)
## annotations
library(biomaRt)
ensembl.db <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gid <- colnames(data)
gene.an <- getBM(attributes=c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "unigene", "description", "chromosome_name", "start_position", "end_position", "strand", "band"), filters="ensembl_gene_id", values=gid, mart=ensembl.db)
gene.an[gene.an == "" | gene.an == " "] <- NA
gene.an <- gene.an[!is.na(gene.an[ , "ensembl_gene_id"]) & !duplicated(gene.an[ , "ensembl_gene_id"]), , drop=FALSE]
annot <- data.frame(matrix(NA, nrow=length(gid), ncol=ncol(gene.an), dimnames=list(gid, colnames(gene.an))))
annot[match(gene.an[ , "ensembl_gene_id"], gid), colnames(gene.an)] <- gene.an
annot <- data.frame("probe"=rownames(annot), "EntrezGene.ID"=annot[ ,"entrezgene"], annot)
## identify the best enesemble gene id for each entrez gene id
gid1 <- as.character(annot[ ,"EntrezGene.ID"])
names(gid1) <- rownames(annot)
gid2 <- sort(unique(gid1))
names(gid2) <- paste("geneid", gid2, sep=".")
rr <- geneid.map(geneid1=gid1, data1=data, geneid2=gid2)
annot <- data.frame(annot, "best"=FALSE)
annot[names(rr$geneid1), "best"] <- TRUE
## rename objects
data.rnaseq <- data
annot.rnaseq <- annot[colnames(data), ,drop=FALSE]
sampleinfo.rnaseq <- sampleinfo[rownames(data), ,drop=FALSE]
rm(list=c("data", "annot", "sampleinfo", "dataf", "annott", "dd", "ddf"))
gc()
save(list=c("data.rnaseq", "annot.rnaseq", "sampleinfo.rnaseq"), compress=TRUE, file=myfn)
} else { load(myfn) }

#################################################
## clinical information
dd <- read.csv("../DNA11161_clinical_info_201202.csv", row.names=NULL, stringsAsFactors=FALSE)
rownames(dd) <- paste(gsub(" ", "_", gsub(" $", "", dd[ ,1])), "T", sep="_")
dd <- data.frame("samplename"=rownames(dd), dd)
nn <- sort(unique(c(row.names(sampleinfo.affy), row.names(sampleinfo.rnaseq), row.names(dd))))
demo <- data.frame(matrix(NA, nrow=length(nn), ncol=ncol(dd), dimnames=list(nn, colnames(dd))), stringsAsFactors=FALSE)
demo <- setcolclass.df(df=demo, colclass=sapply(dd, class), factor.levels=sapply(dd, levels))
demo[rownames(dd), ] <- dd
write.csv(demo, file=sprintf("%s/DNA11161_demo.csv", saveres))

#################################################
## retrict data to common tumors and genes
myfn <- sprintf("%s/dna11161_common.RData", saveres)
if(!file.exists(myfn)) {
## common tumors
nn <- intersect(rownames(sampleinfo.affy), rownames(sampleinfo.rnaseq))
## common genes
gid.affy <- annot.affy[ ,"EntrezGene.ID"]
names(gid.affy) <- rownames(annot.affy)
gid.rnaseq <- annot.rnaseq[ ,"EntrezGene.ID"]
names(gid.rnaseq) <- rownames(annot.rnaseq)
ng <- sort(unique(intersect(gid.affy[annot.affy[ ,"best"]], gid.rnaseq[annot.rnaseq[ ,"best"]])))
datac.affy <- data.affy[nn, !is.na(gid.affy) & annot.affy[ ,"best"] & (gid.affy %in% ng)]
colnames(datac.affy) <- paste("geneid", gid.affy[colnames(datac.affy)], sep=".")
datac.rnaseq <- data.rnaseq[nn, !is.na(gid.rnaseq) & annot.rnaseq[ ,"best"] & (gid.rnaseq %in% ng)]
annotc <- annot.rnaseq[colnames(datac.rnaseq), ,drop=FALSE]
colnames(datac.rnaseq) <- rownames(annotc) <- paste("geneid", gid.rnaseq[colnames(datac.rnaseq)], sep=".")
datac.rnaseq <- datac.rnaseq[ ,colnames(datac.affy),drop=FALSE]
democ <- demo[nn, ,drop=FALSE]
## remove genes with greater or equal than 50% of missing values in either datasets
myx <- apply(datac.affy, 2, function(x, y) { return(sum(is.na(x)) < (length(x) * y)) }, y=0.5) & apply(datac.rnaseq, 2, function(x, y) { return(sum(is.na(x)) < (length(x) * y)) }, y=0.5)
datac.affy <- datac.affy[ ,myx,drop=FALSE]
datac.rnaseq <- datac.rnaseq[ ,myx,drop=FALSE]
annotc <- annotc[myx, ,drop=FALSE]
save(list=c("datac.rnaseq", "datac.affy", "annotc", "democ"), compress=TRUE, file=myfn)
} else { load(myfn) }


#################################################
## correlations
nc <- ncol(datac.affy)
myfn <- sprintf("%s/cores.RData", saveres)
if(!file.exists(myfn)) {
cores <- sapply(1:nc, function(x, y, z) {
	# ct <- cor.test(x=y[ ,x], y=z[ ,x], method="spearman", use="complete.obs")
	# return(ct$estimate)
	return(cor(x=y[ ,x], y=z[ ,x], method="spearman", use="complete.obs"))
}, y=datac.affy, z=datac.rnaseq)
names(cores) <- colnames(datac.affy)
save(list=c("cores"), compress=TRUE, file=myfn)
} else { load(myfn) }

pdf(sprintf("%s/correlation_allgenes.pdf", saveres))
hist(cores, breaks=100, xlim=c(-1,1), freq=FALSE, main="Correlation for all genes\nAFFY vs ILLUMINA RNA-seq", xlab="Spearman correlation", sub=sprintf("Quantiles\t5%%: %.2g; 25%%: %.2g; 50%%: %.2g; 75%%: %.2g; 95%%: %.2g", quantile(cores, probs=0.05, na.rm=TRUE), quantile(cores, probs=0.25, na.rm=TRUE), quantile(cores, probs=0.5, na.rm=TRUE), quantile(cores, probs=0.75, na.rm=TRUE), quantile(cores, probs=0.95, na.rm=TRUE)))
dev.off()

# myfn <- sprintf("%s/pcores.RData", saveres)
# if(!file.exists(myfn)) {
# pcores.affy <- cor(datac.affy, method="spearman", use="pairwise.complete.obs")
# pcores.rnaseq <- cor(datac.rnaseq, method="spearman", use="pairwise.complete.obs")
# pcores <- sapply(1:nc, function(x, y, z) {
# 	# ct <- cor.test(x=y[ ,x], y=z[ ,x], method="spearman", use="complete.obs")
# 	# return(ct$estimate)
# 	return(cor(x=y[ ,x], y=z[ ,x], method="spearman", use="complete.obs"))
# }, y=pcores.affy, z=pcores.rnaseq)
# names(cores) <- colnames(datac.affy)
# save(list=c("pcores.affy", "pcores.rnaseq", "pcores"), compress=TRUE, file=myfn)
# } else { load(myfn) }
# 
# pdf(sprintf("%s/pairw_correlation_allgenes.pdf", saveres))
# hist(pcores, breaks=100, xlim=c(-1,1), freq=FALSE, main="Correlation of pairwise correlations for all genes\nAFFY vs ILLUMINA RNA-seq", xlab="Spearman correlation", sub=sprintf("Quantiles\t5%%: %.2g; 25%%: %.2g; 50%%: %.2g; 75%%: %.2g; 95%%: %.2g", quantile(cores, probs=0.05, na.rm=TRUE), quantile(cores, probs=0.25, na.rm=TRUE), quantile(cores, probs=0.5, na.rm=TRUE), quantile(cores, probs=0.75, na.rm=TRUE), quantile(cores, probs=0.95, na.rm=TRUE)))
# dev.off()

## test the corrgram package for visualizing pairwise correlation


#################################################
## gene signatures

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
pdf(sprintf("%s/genius_risk_contingency_table.pdf", saveres))
plot(t(tt), main="Contingency table for GENIUS risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

pdf(sprintf("%s/genius_scores.pdf", saveres))
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col="darkgrey", main="GENIUS\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="GENIUS scores (AFFY)", ylab="GENIUS scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
dev.off()

pdf(sprintf("%s/genius_scores_risk.pdf", saveres))
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
pdf(sprintf("%s/ggi_risk_contingency_table.pdf", saveres))
plot(t(tt), main="Contingency table for GGI risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

pdf(sprintf("%s/ggi_scores.pdf", saveres))
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col="darkgrey", main="GGI\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="GGI scores (AFFY)", ylab="GGI scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
dev.off()

pdf(sprintf("%s/ggi_scores_risk.pdf", saveres))
mycol <- rep("red", length(nn))
cc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[cc == 0] <- "blue"
mycol[cc == 2] <- "yellow"
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col=mycol, main="GGI\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="GGI scores (AFFY)", ylab="GGI scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
legend("topleft", col=c("blue", "yellow", "red"), legend=c("Low-risk", "High-risk", "Discordance"), bty="n", pch=16)
dev.off()

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
pdf(sprintf("%s/mammaprint_risk_contingency_table.pdf", saveres))
plot(t(tt), main="Contingency table for MAMMAPRINT risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

pdf(sprintf("%s/mammaprint_scores.pdf", saveres))
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col="darkgrey", main="MAMMAPRINT\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="MAMMAPRINT scores (AFFY)", ylab="MAMMAPRINT scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
dev.off()

pdf(sprintf("%s/mammaprint_scores_risk.pdf", saveres))
mycol <- rep("red", length(nn))
cc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[cc == 0] <- "blue"
mycol[cc == 2] <- "yellow"
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col=mycol, main="MAMMAPRINT\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="MAMMAPRINT scores (AFFY)", ylab="MAMMAPRINT scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
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
pdf(sprintf("%s/tamr13_risk_contingency_table.pdf", saveres))
plot(t(tt), main="Contingency table for TAMR13 risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

pdf(sprintf("%s/tamr13_scores.pdf", saveres))
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col="darkgrey", main="TAMR13\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="TAMR13 scores (AFFY)", ylab="TAMR13 scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
dev.off()

pdf(sprintf("%s/tamr13_scores_risk.pdf", saveres))
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
pdf(sprintf("%s/pik3cags_risk_contingency_table.pdf", saveres))
plot(t(tt), main="Contingency table for PIK3CAGS risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

pdf(sprintf("%s/pik3cags_scores.pdf", saveres))
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col="darkgrey", main="PIK3CAGS\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="PIK3CAGS scores (AFFY)", ylab="PIK3CAGS scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
dev.off()

pdf(sprintf("%s/pik3cags_scores_risk.pdf", saveres))
mycol <- rep("red", length(nn))
cc <- sig.affy$risk[nn] + sig.rnaseq$risk[nn]
mycol[cc == 0] <- "blue"
mycol[cc == 2] <- "yellow"
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col=mycol, main="PIK3CAGS\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="PIK3CAGS scores (AFFY)", ylab="PIK3CAGS scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
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
pdf(sprintf("%s/oncotypedx_risk_contingency_table.pdf", saveres))
plot(t(tt), main="Contingency table for ONCOTYPE DX risk predictions\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

pdf(sprintf("%s/oncotypedx_scores.pdf", saveres))
plot(sig.affy$score[nn], sig.rnaseq$score[nn], pch=16, col="darkgrey", main="ONCOTYPE DX\nAFFYMETRIX vs ILLUMINA RNA-seq", xlab="ONCOTYPE DX scores (AFFY)", ylab="ONCOTYPE DX scores (ILLUMINA RNA-seq)", sub=sprintf("Spearman correlation: %.2g", cor(sig.affy$score[nn], sig.rnaseq$score[nn], method="spearman", use="complete.obs")))
dev.off()

pdf(sprintf("%s/oncotypedx_scores_risk.pdf", saveres))
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
pdf(sprintf("%s/scmgene_sbt_plot_affy.pdf", saveres))
sig.affy <- subtype.cluster.predict(sbt.model=scmgene.robust, data=data.affy, annot=annot.affy, do.mapping=FALSE, plot=TRUE)
title(main="SCMGENE on AFFYMETRIX")
sig.affy$subtype2 <- factor(sig.affy$subtype2, levels=sbtn)
dev.off()
pdf(sprintf("%s/scmgene_sbt_plot_rnaseq.pdf", saveres))
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
pdf(sprintf("%s/scmgene_risk_contingency_table.pdf", saveres))
plot(t(tt), main="Contingency table for SCMGENE subtype classification\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

##############
## SCMOD1
##############
pdf(sprintf("%s/scmod1_sbt_plot_affy.pdf", saveres))
sig.affy <- subtype.cluster.predict(sbt.model=scmod1.robust, data=data.affy, annot=annot.affy, do.mapping=FALSE, plot=TRUE)
title(main="SCMOD1 on AFFYMETRIX")
sig.affy$subtype2 <- factor(sig.affy$subtype2, levels=sbtn)
dev.off()
pdf(sprintf("%s/scmod1_sbt_plot_rnaseq.pdf", saveres))
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
pdf(sprintf("%s/scmod1_risk_contingency_table.pdf", saveres))
plot(t(tt), main="Contingency table for SCMOD1 subtype classification\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

##############
## SCMOD2
##############
pdf(sprintf("%s/scmod2_sbt_plot_affy.pdf", saveres))
sig.affy <- subtype.cluster.predict(sbt.model=scmod2.robust, data=data.affy, annot=annot.affy, do.mapping=FALSE, plot=TRUE)
title(main="SCMOD2 on AFFYMETRIX")
sig.affy$subtype2 <- factor(sig.affy$subtype2, levels=sbtn)
dev.off()
pdf(sprintf("%s/scmod2_sbt_plot_rnaseq.pdf", saveres))
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
pdf(sprintf("%s/scmod2_risk_contingency_table.pdf", saveres))
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
pdf(sprintf("%s/pam50_risk_contingency_table.pdf", saveres))
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
pdf(sprintf("%s/ssp2006_risk_contingency_table.pdf", saveres))
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
pdf(sprintf("%s/ssp2003_risk_contingency_table.pdf", saveres))
plot(t(tt), main="Contingency table for SSP2003 subtype classification\nAFFYMETRIX vs ILLUMINA RNA-seq", sub=sprintf("Kappa coefficient=%.2g,  95%%CI [%.2g,%.2g], p=%.1E", kk$kappa, kk$CIL, kk$CIU, tts$chisq_tests[1,3]))
dev.off()

write.csv(sbts, file=sprintf("%s/subtype_classifs_affy_rnaseq.csv", saveres))
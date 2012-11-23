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
## step 1: Compute pairwise correlation between Affy and illumina RNA-seq
########################

########################
## retrict data to common tumors and genes
## step 1a: select best affy probeset using jetset
myfn <- file.path(saveres, "dna11161_jetset_common.RData")
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

########################
## retrict data to common tumors and genes
## step 1b: select the most coprrelated affy probesets for each gene
myfn <- file.path(saveres, "dna11161_bestpgene_common.RData")

########################
## retrict data to common tumors and genes
## step 1c: select the most coprrelated affy probesets for each trsnacript
myfn <- file.path(saveres, "dna11161_bestptranscript_common.RData")

########################
## for each dataset compute correlations

## load data
myfns <- c("jetset"=file.path(saveres, "dna11161_jetset_common.RData"), "bestpgene"=file.path(saveres, "dna11161_bestpgene_common.RData"), "bestptranscript"=file.path(saveres, "dna11161_bestptranscript_common.RData"))



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

pdf(file.path(saveres, "correlation_allgenes.pdf"))
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
# pdf(file.path(saveres, "pairw_correlation_allgenes.pdf"))
# hist(pcores, breaks=100, xlim=c(-1,1), freq=FALSE, main="Correlation of pairwise correlations for all genes\nAFFY vs ILLUMINA RNA-seq", xlab="Spearman correlation", sub=sprintf("Quantiles\t5%%: %.2g; 25%%: %.2g; 50%%: %.2g; 75%%: %.2g; 95%%: %.2g", quantile(cores, probs=0.05, na.rm=TRUE), quantile(cores, probs=0.25, na.rm=TRUE), quantile(cores, probs=0.5, na.rm=TRUE), quantile(cores, probs=0.75, na.rm=TRUE), quantile(cores, probs=0.95, na.rm=TRUE)))
# dev.off()

## test the corrgram package for visualizing pairwise correlation

########################
## Benjamin Haibe-Kains
##
## Nov 23, 2012
########################

rm(list=ls(all=TRUE))

library(mclust)

source(file.path("code", "ufoo.R"))

saveres <- file.path("saveres")
if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE) }

mycol <- c("darkblue", "darkorange", "darkred")

########################
## step 1: Compute pairwise correlation between Affy and illumina RNA-seq
########################

########################
## retrict data to common tumors and genes
## step 1a: select best affy probeset using jetset
myfn <- file.path(saveres, "dna11161_common_jetset.RData")
if(!file.exists(myfn)) {
  message("Mapping using jetset")
  ## load rnaseq gene data
  load(file.path(saveres, "dna11161_gene_rnaseq.RData"))
  ## load affy data
  load(file.path(saveres, "dna11161_affy_frma.RData"))
  ## common tumors
  nn <- intersect(rownames(sampleinfo.affy), rownames(sampleinfo.rnaseq))
  sampleinfoc <- sampleinfo.affy[nn, , drop=FALSE]
  ## common genes
  gid.affy <- as.character(annot.affy[ ,"EntrezGene.ID"])
  names(gid.affy) <- rownames(annot.affy)
  gid.rnaseq <- as.character(annot.gene.rnaseq[ ,"EntrezGene.ID"])
  names(gid.rnaseq) <- rownames(annot.gene.rnaseq)
  ng <- sort(unique(intersect(gid.affy[annot.affy[ ,"best"]], gid.rnaseq[annot.gene.rnaseq[ ,"best"]])))
  ## affy: selection with jetset
  datac.affy <- data.affy[nn, !is.na(gid.affy) & annot.affy[ ,"best"] & (gid.affy %in% ng)]
  pp <- colnames(datac.affy)
  colnames(datac.affy) <- names(pp) <- paste("geneid", gid.affy[colnames(datac.affy)], sep=".")
  ## rnaseq: selection of the most variant
  datac.rnaseq <- data.gene.rnaseq[nn, !is.na(gid.rnaseq) & annot.gene.rnaseq[ ,"best"] & (gid.rnaseq %in% ng)]
  pp2 <- colnames(datac.rnaseq)
  colnames(datac.rnaseq) <- names(pp2) <- paste("geneid", gid.rnaseq[colnames(datac.rnaseq)], sep=".")
  pp <- pp[paste("geneid", ng, sep=".")]
  pp2 <- pp2[paste("geneid", ng, sep=".")]
  datac.affy <- datac.affy[ , paste("geneid", ng, sep="."), drop=FALSE]
  datac.rnaseq <- datac.rnaseq[ , paste("geneid", ng, sep="."), drop=FALSE]
  ## annotations
  annotc <- annot.gene.rnaseq[pp2, , drop=FALSE]
  rownames(annotc) <- names(pp2)
  annotc <- annotc[paste("geneid", ng, sep="."), , drop=FALSE]
  annotc[ ,"probe"] <- rownames(annotc)
  annotc <- data.frame("probeset.affy"=pp, annotc)
  ## remove genes with greater or equal than 50% of missing values in either datasets
  myx <- apply(datac.affy, 2, function(x, y) { return(sum(is.na(x)) < (length(x) * y)) }, y=0.5) & apply(datac.rnaseq, 2, function(x, y) { return(sum(is.na(x)) < (length(x) * y)) }, y=0.5)
  datac.affy <- datac.affy[ ,myx,drop=FALSE]
  datac.rnaseq <- datac.rnaseq[ ,myx,drop=FALSE]
  annotc <- annotc[myx, ,drop=FALSE]
  save(list=c("datac.rnaseq", "datac.affy", "annotc", "sampleinfoc"), compress=TRUE, file=myfn)
} else { load(myfn) }

########################
## restrict data to common tumors and genes
## step 1b: select the most correlated affy probesets for each gene
myfn <- file.path(saveres, "dna11161_common_bestpgene.RData")
if(!file.exists(myfn)) {
  message("Mapping using the most correlated probeset/gene per entrezgene id")
  ## load rnaseq gene data
  load(file.path(saveres, "dna11161_gene_rnaseq.RData"))
  ## load affy data
  load(file.path(saveres, "dna11161_affy_frma.RData"))
  ## common tumors
  nn <- intersect(rownames(sampleinfo.affy), rownames(sampleinfo.rnaseq))
  sampleinfoc <- sampleinfo.affy[nn, , drop=FALSE]
  data.affy <- data.affy[nn, , drop=FALSE]
  data.gene.rnaseq <- data.gene.rnaseq[nn, , drop=FALSE]
  dataf.gene.rnaseq <- dataf.gene.rnaseq[nn, , drop=FALSE]
  ## common genes
  gid.affy <- annot.affy[ ,"EntrezGene.ID"]
  names(gid.affy) <- rownames(annot.affy)
  gid.rnaseq <- annot.gene.rnaseq[ ,"EntrezGene.ID"]
  names(gid.rnaseq) <- rownames(annot.gene.rnaseq) 
  
  ## for each gene in common that are not unique in both platforms, look for the best correlated pairs
  gid.common <- sort(intersect(annot.affy[ ,"EntrezGene.ID"], annot.gene.rnaseq[ ,"EntrezGene.ID"]))
  gid.affy <- as.character(annot.affy[ ,"EntrezGene.ID"])
  names(gid.affy) <- rownames(annot.affy)
  gid.rnaseq <- as.character(annot.gene.rnaseq[ ,"EntrezGene.ID"])
  names(gid.rnaseq) <- rownames(annot.gene.rnaseq)
  gid.dupl.affy <- sort(intersect(gid.affy[duplicated(gid.affy)], gid.common))
  gid.dupl.rnaseq <- sort(intersect(gid.rnaseq[duplicated(gid.rnaseq)], gid.common))
  ## ambiguous genes
  gid.ambig <- intersect(unique(c(gid.dupl.affy, gid.dupl.rnaseq)), gid.common)
  rr <- t(sapply(gid.ambig, function(x, data1, gid1, data2, gid2) {
    cc <- cor(data1[ , !is.na(gid1) & gid1 == x, drop=FALSE], data2[ , !is.na(gid2) & gid2 == x, drop=FALSE], method="spearman", use="pairwise.complete.obs")
    mm <- which.max(cc)
    idx1 <- mm %% nrow(cc)
    if(idx1 == 0) { idx1 <- nrow(cc) }
    idx2 <- ceiling(mm / nrow(cc))
    return(c(rownames(cc)[idx1], colnames(cc)[idx2]))
  }, data1=data.affy, gid1=gid.affy, data2=data.gene.rnaseq, gid2=gid.rnaseq))
  dimnames(rr) <- list(paste("geneid", gid.ambig, sep="."), c("affy", "rnaseq"))
  ## unambiguous genes and merge
  gid.uniq <- setdiff(gid.common, gid.ambig)
  rr2 <- cbind(names(gid.affy)[match(gid.uniq, gid.affy)], names(gid.rnaseq)[match(gid.uniq, gid.rnaseq)])
  dimnames(rr2) <- list(paste("geneid", gid.uniq, sep="."), c("affy", "rnaseq"))
  mapp <- rbind(rr, rr2)
  mapp <- mapp[complete.cases(mapp), , drop=FALSE]
    
  datac.affy <- data.affy[ , mapp[ ,"affy"], drop=FALSE]
  datac.rnaseq <- data.gene.rnaseq[ , mapp[ ,"rnaseq"], drop=FALSE]
  datafc.rnaseq <- dataf.gene.rnaseq[ , mapp[ ,"rnaseq"], drop=FALSE]
  annotc <- data.frame("probeset.affy"=mapp[ , "affy"], annot.gene.rnaseq[mapp[ , "rnaseq"], , drop=FALSE])
  colnames(datac.affy) <- colnames(datac.rnaseq) <- colnames(datafc.rnaseq) <- rownames(annotc) <- rownames(mapp)

  save(list=c("datac.rnaseq", "datac.rnaseq", "datac.affy", "annotc", "sampleinfoc"), compress=TRUE, file=myfn)
} else { load(myfn) }

########################
## retrict data to common tumors and genes
## step 1c: select the most coprrelated affy probesets for each trsancript
myfn <- file.path(saveres, "dna11161_common_bestptranscript.RData")
if(!file.exists(myfn)) {
  message("Mapping using the most correlated probeset/transcript per entrezgene id")
  ## load rnaseq gene data
  load(file.path(saveres, "dna11161_transcript_rnaseq.RData"))
  ## load affy data
  load(file.path(saveres, "dna11161_affy_frma.RData"))
  ## common tumors
  nn <- intersect(rownames(sampleinfo.affy), rownames(sampleinfo.rnaseq))
  sampleinfoc <- sampleinfo.affy[nn, , drop=FALSE]
  data.affy <- data.affy[nn, , drop=FALSE]
  data.transcript.rnaseq <- data.transcript.rnaseq[nn, , drop=FALSE]
  dataf.transcript.rnaseq <- dataf.transcript.rnaseq[nn, , drop=FALSE]
  ## common genes
  gid.affy <- annot.affy[ ,"EntrezGene.ID"]
  names(gid.affy) <- rownames(annot.affy)
  gid.rnaseq <- annot.transcript.rnaseq[ ,"EntrezGene.ID"]
  names(gid.rnaseq) <- rownames(annot.transcript.rnaseq) 
  
  ## for each gene in common that are not unique in both platforms, look for the best correlated pairs
  gid.common <- sort(intersect(annot.affy[ ,"EntrezGene.ID"], annot.transcript.rnaseq[ ,"EntrezGene.ID"]))
  gid.affy <- as.character(annot.affy[ ,"EntrezGene.ID"])
  names(gid.affy) <- rownames(annot.affy)
  gid.rnaseq <- as.character(annot.transcript.rnaseq[ ,"EntrezGene.ID"])
  names(gid.rnaseq) <- rownames(annot.transcript.rnaseq)
  gid.dupl.affy <- sort(intersect(gid.affy[duplicated(gid.affy)], gid.common))
  gid.dupl.rnaseq <- sort(intersect(gid.rnaseq[duplicated(gid.rnaseq)], gid.common))
  ## ambiguous genes
  gid.ambig <- intersect(unique(c(gid.dupl.affy, gid.dupl.rnaseq)), gid.common)
  rr <- t(sapply(gid.ambig, function(x, data1, gid1, data2, gid2) {
    cc <- cor(data1[ , !is.na(gid1) & gid1 == x, drop=FALSE], data2[ , !is.na(gid2) & gid2 == x, drop=FALSE], method="spearman", use="pairwise.complete.obs")
    if(!all(is.na(cc))) {
      mm <- which.max(cc)
      idx1 <- mm %% nrow(cc)
      if(idx1 == 0) { idx1 <- nrow(cc) }
      idx2 <- ceiling(mm / nrow(cc))
      return(c(rownames(cc)[idx1], colnames(cc)[idx2]))
      } else { return(c(NA, NA)) }
  }, data1=data.affy, gid1=gid.affy, data2=data.transcript.rnaseq, gid2=gid.rnaseq))
  dimnames(rr) <- list(paste("geneid", gid.ambig, sep="."), c("affy", "rnaseq"))
  ## unambiguous genes and merge
  gid.uniq <- setdiff(gid.common, gid.ambig)
  rr2 <- cbind(names(gid.affy)[match(gid.uniq, gid.affy)], names(gid.rnaseq)[match(gid.uniq, gid.rnaseq)])
  dimnames(rr2) <- list(paste("geneid", gid.uniq, sep="."), c("affy", "rnaseq"))
  mapp <- rbind(rr, rr2)
  mapp <- mapp[complete.cases(mapp), , drop=FALSE]
    
  datac.affy <- data.affy[ , mapp[ ,"affy"], drop=FALSE]
  datac.rnaseq <- data.transcript.rnaseq[ , mapp[ ,"rnaseq"], drop=FALSE]
  datafc.rnaseq <- dataf.transcript.rnaseq[ , mapp[ ,"rnaseq"], drop=FALSE]
  annotc <- data.frame("probeset.affy"=mapp[ , "affy"], annot.transcript.rnaseq[mapp[ , "rnaseq"], , drop=FALSE])
  colnames(datac.affy) <- colnames(datac.rnaseq) <- colnames(datafc.rnaseq) <- rownames(annotc) <- rownames(mapp)

  save(list=c("datac.rnaseq", "datac.rnaseq", "datac.affy", "annotc", "sampleinfoc"), compress=TRUE, file=myfn)
} else { load(myfn) }

########################
## for each dataset compute correlations

  
## load data
myfns <- c("jetset"=file.path(saveres, "dna11161_common_jetset.RData"), "bestpgene"=file.path(saveres, "dna11161_common_bestpgene.RData"), "bestptranscript"=file.path(saveres, "dna11161_common_bestptranscript.RData"))

for(i in 1:length(myfns)) {
  
  ## load datasets
  load(myfns[i])
  
  ## create directory for results
  saveres2 <- file.path("saveres", names(myfns)[i])
  if(!file.exists(saveres2)) { dir.create(saveres2, showWarnings=FALSE) }
  
  nc <- ncol(datac.affy)
  myfn <- file.path(saveres2, sprintf("cores_%s.RData", names(myfns)[i]))
  if(!file.exists(myfn)) {
  cores <- sapply(1:nc, function(x, y, z) {
    if(sum(complete.cases(y[ ,x], z[ ,x])) > 3) { return(cor(x=y[ ,x], y=z[ ,x], method="spearman", use="complete.obs")) } else { return(NA) }
  }, y=datac.affy, z=datac.rnaseq)
  names(cores) <- colnames(datac.affy)
  save(list=c("cores"), compress=TRUE, file=myfn)
  } else { load(myfn) }

  ## histogram of correlation for each gene
  pdf(file.path(saveres2, sprintf("correlation_%s_allgenes.pdf", names(myfns)[i])))
  hist(cores, breaks=100, xlim=c(-1,1), ylim=c(0, 800), freq=TRUE, main=sprintf("Correlation for all genes [%s]\nAFFY vs ILLUMINA RNA-seq", names(myfns)[i]), xlab="Spearman correlation", sub=sprintf("Quantiles\t5%%: %.2g; 25%%: %.2g; 50%%: %.2g; 75%%: %.2g; 95%%: %.2g", quantile(cores, probs=0.05, na.rm=TRUE), quantile(cores, probs=0.25, na.rm=TRUE), quantile(cores, probs=0.5, na.rm=TRUE), quantile(cores, probs=0.75, na.rm=TRUE), quantile(cores, probs=0.95, na.rm=TRUE)))
  dev.off()
  
  ## use an arbitrary cutoff of 0.7 and classify genes into low and high cor classes
  rrc <- 0.7
  rrg <- ifelse(cores >= rrc, "high.cor", "low.cor")
  rrg <- factor(rrg, levels=c("low.cor", "high.cor"))
  tt <- data.frame("cor"=cores, "cor.low.vs.high"=rrg, annotc[names(cores), , drop=FALSE])
  write.csv(tt, file=file.path(saveres2, sprintf("correlation_%s_allgenes.csv", names(myfns)[i])))
  
  ## boxplot for expression avlues of low vs high correlated genes
  ## affy
  pdf(file.path(saveres2, sprintf("boxplot_cor_vs_median_expr_%s_allgenes.pdf", names(myfns)[i])), width=14, height=7)
  par(mfrow=c(1, 2))
  ll <- list("low.cor"=apply(datac.affy[ , names(cores)[!is.na(rrg) & rrg == "low.cor"], drop=FALSE], 2, median, na.rm=TRUE), "high.cor"=apply(datac.affy[ , names(cores)[!is.na(rrg) & rrg == "high.cor"], drop=FALSE], 2, median, na.rm=TRUE))
  wt <- wilcox.test(x=ll[["low.cor"]], y=ll[["high.cor"]])
  boxplot(ll, outline=FALSE, col="lightgrey", main="Median expressions on Affymetrix", xlab=sprintf("Wilcooxn rank sum test p-value = %.1E", wt$p.value), ylab="Gene median expression")
  lapply(ll, median)
  ll <- list("low.cor"=apply(datac.rnaseq[ , names(cores)[!is.na(rrg) & rrg == "low.cor"], drop=FALSE], 2, median, na.rm=TRUE), "high.cor"=apply(datac.rnaseq[ , names(cores)[!is.na(rrg) & rrg == "high.cor"], drop=FALSE], 2, median, na.rm=TRUE))
  wt <- wilcox.test(x=ll[["low.cor"]], y=ll[["high.cor"]])
  boxplot(ll, outline=FALSE, col="lightgrey", main="Median expressions on ILLUMINA RNA-seq", xlab=sprintf("Wilcooxn rank sum test p-value = %.1E", wt$p.value), ylab="Gene median expression")
  lapply(ll, median)
  dev.off()
  
  ## same plot with mixture of 2 gaussians
  rr <- mclust::Mclust(data=cores[!is.na(cores)], G=2, modelNames="V")
  ## cutoff
  rrc <- (max(cores[!is.na(cores)][rr$classification == 1]) + min(cores[!is.na(cores)][rr$classification == 2]))/2
  xx1 <- seq(-1, 1, 0.01)
  yy1 <- dnorm(x=xx1, mean=rr$parameters$mean[1], sd=sqrt(rr$parameters$variance$sigmasq[1]))
  xx2 <- seq(-1, 1, 0.01)
  yy2 <- dnorm(x=xx2, mean=rr$parameters$mean[2], sd=sqrt(rr$parameters$variance$sigmasq[2]))
  pdf(file.path(saveres2, sprintf("correlation_mclust_%s_allgenes.pdf", names(myfns)[i])))
  hist(cores, breaks=100, xlim=c(-1,1), ylim=c(0, 800), freq=TRUE, main=sprintf("Correlation for all genes [%s]\nAFFY vs ILLUMINA RNA-seq", names(myfns)[i]), xlab="Spearman correlation", sub=sprintf("Cutoff = %.3g", rrc))
  lines(x=xx1, y=yy1*100, type="l", lty=1, lwd=2, col=mycol[3])
  lines(x=xx2, y=yy2*200, type="l", lty=1, lwd=2, col=mycol[1])
  abline(v=rrc, col="black", lty=2, lwd=2)
  legend("topleft", col=mycol[c(1, 3)], lwd=c(2,2), legend=c(sprintf("Low cor (mean=%.2g, sd=%.2g)", rr$parameters$mean[1], rr$parameters$variance$sigmasq[1]), sprintf("High cor (mean=%.2g, sd=%.2g)", rr$parameters$mean[2], rr$parameters$variance$sigmasq[2])), bty="n")
  dev.off()
  
  ## histogram of correlation for each pair of sample
  myfn <- file.path(saveres2, sprintf("corespairw_%s.RData", names(myfns)[i]))
  if(!file.exists(myfn)) {
  corespairw <- cor(t(datac.affy), t(datac.rnaseq), method="spearman", use="pairwise.complete.obs")
  save(list=c("corespairw"), compress=TRUE, file=myfn)
  } else { load(myfn) }
  
  pdf(file.path(saveres2, sprintf("correlation_pairw_heatmap_%s_allgenes.pdf", names(myfns)[i])))
  hist(diag(corespairw), breaks=10, xlim=c(0.5,1), freq=FALSE, main=sprintf("Pairwise correlation using all genes [%s]\nAFFY vs ILLUMINA RNA-seq", names(myfns)[i]), xlab="Spearman correlation", sub=sprintf("Quantiles\t5%%: %.2g; 25%%: %.2g; 50%%: %.2g; 75%%: %.2g; 95%%: %.2g", quantile(diag(corespairw), probs=0.05, na.rm=TRUE), quantile(diag(corespairw), probs=0.25, na.rm=TRUE), quantile(diag(corespairw), probs=0.5, na.rm=TRUE), quantile(diag(corespairw), probs=0.75, na.rm=TRUE), quantile(diag(corespairw), probs=0.95, na.rm=TRUE)))
  dev.off()
  
  ## scatterplot for each pair of samples
  pdf(file.path(saveres2, sprintf("scatterplot_pairw_%s_allgenes.pdf", names(myfns)[i])))
  for(j in 1:nrow(datac.affy)) {
    cc <- cor(x=datac.affy[j, ], y=datac.rnaseq[j, ], method="spearman", use="complete.obs")
    cci <- spearmanCI(cc, sum(complete.cases(datac.affy[j, ], datac.rnaseq[j, ])))
    smoothScatter(x=datac.affy[j, ], y=datac.rnaseq[j, ], main=sprintf("%s, all genes [%s]\nAFFY vs ILLUMINA RNA-seq", rownames(datac.affy)[j], names(myfns)[i]), xlab="Gene expression (AFFY)", ylab="Gene expression (ILLUMINA RNA-seq)", sub=sprintf("Sperman correlation = %.2g [%.2g,%.2g], p=%.1E", cc, cci[1], cci[2], cci[3]))
    ## add regression line
    reg1 <- lm(datac.rnaseq[j, ] ~ datac.affy[j, ])
    abline(reg1, col="black")
  }
  dev.off()
}

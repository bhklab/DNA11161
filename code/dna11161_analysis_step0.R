########################
## Benjamin Haibe-Kains
##
## Nov 23, 2012
########################

rm(list=ls(all=TRUE))

source(file.path("code", "ufoo.R"))

saveres <- file.path("saveres")
if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE) }

########################
## step 0: data normalization and formatting
########################

## microarray data
myfn <- file.path(saveres, "dna11161_affy_frma.RData")
if(!file.exists(myfn)) {
library(affy)
library(frma)
library(frmaTools)
## CEL file names
rdp <- file.path("affy")
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


## RNA-seq gene-centric and transcript data
myfn <- file.path(saveres, "dna11161_gene_rnaseq.RData")
if(!file.exists(myfn)) {
  
## sample info
message("Read sample info")
sampleinfo <- read.csv(file.path("rnaseq" ,"DNA111061_dnavision_ids.csv"), row.names=NULL, stringsAsFactors=FALSE)
rownames(sampleinfo) <- gsub("BIS", "", toupper(as.character(sampleinfo[ ,2])))
sampleinfo <- data.frame("samplename"=rownames(sampleinfo), sampleinfo)


## alignment with tophat and expression estimation (FPKM) with cufflink

## gene expressions
message("Read gene expression values generated with RNA-seq")
dd <- read.csv(file.path("rnaseq", "consolidated", "allFpkmGenesEnsembl.txt"), sep="\t", row.names=NULL, stringsAsFactors=FALSE)
dd[dd == "" | dd == " " | dd == "----"] <- NA
annott <- dd[ ,1:2]
ddf <- dd[ , grep("Status", colnames(dd)), drop=FALSE]
dd <- dd[ , grep("FPKM", colnames(dd)), drop=FALSE]
dd <- data.matrix(dd[apply(dd, 1, function (x) { return(!all(is.na(x))) }), , drop=FALSE])
rownames(dd) <- rownames(ddf) <- as.character(annott[ ,1])
data <- t(dd)
dataf <- t(ddf)
## remove FPKMs of bad quality
data[dataf != "OK"] <- NA
## log2 transformation of the FPKM for better comparability with microarray
data <- log2(data+1)
nn <- gsub("[.]", "-", substr(colnames(dd), 1, nchar(colnames(dd)) - 5))
rownames(data) <- rownames(dataf) <- nn
## rename objects
data.gene <- data
dataf.gene <- dataf
annot.gene <- annott

## transcript expressions
message("Read transcript expression values generated with RNA-seq")
dd <- read.csv(file.path("rnaseq", "consolidated", "allFpkmTranscriptsEnsembl.txt"), sep="\t", row.names=NULL, stringsAsFactors=FALSE)
annott <- dd[ ,1:3]
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
rownames(data) <- rownames(dataf) <- nn
## rename objects
data.transcript <- data
dataf.transcript <- dataf
annot.transcript <- annott

## keep only ensemble gene ids in common between FPKMs at the gene and transcript levels
gn <- intersect(as.character(annot.gene[ ,"gene_id"]), as.character(annot.transcript[ ,"gene_id"]))
gx <- is.element(as.character(annot.gene[ ,"gene_id"]), gn)
data.gene <- data.gene[ , gx, drop=FALSE]
dataf.gene <- dataf.gene[ , gx, drop=FALSE]
annot.gene <- annot.gene[gx, , drop=FALSE]
gx <- is.element(as.character(annot.transcript[ ,"gene_id"]), gn)
data.transcript <- data.transcript[ , gx, drop=FALSE]
dataf.transcript <- dataf.transcript[ , gx, drop=FALSE]
annot.transcript <- annot.transcript[gx, , drop=FALSE]

## keep the samples for which we have all data
nn <- fold(intersect, sampleinfo[ ,"DNAvision.ID"], rownames(data,gene), rownames(data.transcript))
sampleinfo <- sampleinfo[match(nn, sampleinfo[ ,"DNAvision.ID"]), , drop=FALSE]
data.gene <- data.gene[nn, , drop=FALSE]
dataf.gene <- dataf.gene[nn, , drop=FALSE]
data.transcript <- data.transcript[nn, , drop=FALSE]
dataf.transcript <- dataf.transcript[nn, , drop=FALSE]
rownames(data.gene) <- rownames(data.transcript) <- rownames(dataf.gene) <- rownames(dataf.transcript) <- rownames(sampleinfo)

## annotations
message("Annotate genes and transcripts with biomart")

## for genes
library(biomaRt)
ensembl.db <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gid <- colnames(data.gene)
filt <- "ensembl_gene_id"
gene.an <- getBM(attributes=c(filt, "entrezgene", "hgnc_symbol", "unigene", "description", "chromosome_name", "start_position", "end_position", "strand", "band"), filters=filt, values=gid, mart=ensembl.db)
gene.an[gene.an == "" | gene.an == " "] <- NA
gene.an <- gene.an[!is.na(gene.an[ , filt]) & !duplicated(gene.an[ , filt]), , drop=FALSE]
annot <- data.frame(matrix(NA, nrow=length(gid), ncol=ncol(gene.an), dimnames=list(gid, colnames(gene.an))))
annot[match(gene.an[ , filt], gid), colnames(gene.an)] <- gene.an
annot.gene  <- data.frame("probe"=rownames(annot), "EntrezGene.ID"=annot[ ,"entrezgene"], annot.gene, annot)
rownames(annot.gene) <- colnames(data.gene)

## identify the most variant enesemble gene id for each entrez gene id
library(genefu)
gid1 <- as.character(annot.gene[ ,"EntrezGene.ID"])
names(gid1) <- rownames(annot.gene)
gid2 <- sort(unique(gid1))
names(gid2) <- paste("geneid", gid2, sep=".")
rr <- genefu::geneid.map(geneid1=gid1, data1=data.gene, geneid2=gid2)
annot.gene <- data.frame(annot.gene, "best"=FALSE)
annot.gene[names(rr$geneid1), "best"] <- TRUE

## for transcripts
gid <- colnames(data.transcript)
filt <- "ensembl_transcript_id"
gene.an <- getBM(attributes=c(filt, "ensembl_gene_id", "entrezgene", "hgnc_symbol", "unigene", "description", "chromosome_name", "start_position", "end_position", "strand", "band"), filters=filt, values=gid, mart=ensembl.db)
gene.an[gene.an == "" | gene.an == " "] <- NA
gene.an <- gene.an[!is.na(gene.an[ , filt]) & !duplicated(gene.an[ , filt]), , drop=FALSE]
annot <- data.frame(matrix(NA, nrow=length(gid), ncol=ncol(gene.an), dimnames=list(gid, colnames(gene.an))))
annot[match(gene.an[ , filt], gid), colnames(gene.an)] <- gene.an
annot.transcript <- data.frame("probe"=rownames(annot), "EntrezGene.ID"=annot[ ,"entrezgene"], annot.transcript, annot)
rownames(annot.transcript) <- colnames(data.transcript)

## identify the most variant transcript id for each entrez gene id
gid1 <- as.character(annot.transcript[ ,"EntrezGene.ID"])
names(gid1) <- rownames(annot.transcript)
gid2 <- sort(unique(gid1))
names(gid2) <- paste("geneid", gid2, sep=".")
rr <- genefu::geneid.map(geneid1=gid1, data1=data.transcript, geneid2=gid2)
annot.transcript <- data.frame(annot.transcript, "best"=FALSE)
annot.transcript[names(rr$geneid1), "best"] <- TRUE

## rename objects
data.gene.rnaseq <- data.gene
dataf.gene.rnaseq <- dataf.gene
annot.gene.rnaseq <- annot.gene[colnames(data.gene), ,drop=FALSE]
data.transcript.rnaseq <- data.transcript
dataf.transcript.rnaseq <- dataf.transcript
annot.transcript.rnaseq <- annot.transcript[colnames(data.transcript), ,drop=FALSE]
sampleinfo.rnaseq <- sampleinfo
rm(list=c("data", "annot", "sampleinfo", "dataf", "annott", "dd", "ddf"))
gc()
## gene
save(list=c("data.gene.rnaseq", "dataf.gene.rnaseq", "annot.gene.rnaseq", "sampleinfo.rnaseq"), compress=TRUE, file=myfn)
## transcript
save(list=c("data.transcript.rnaseq", "dataf.transcript.rnaseq", "annot.transcript.rnaseq", "sampleinfo.rnaseq"), compress=TRUE, file=file.path(saveres, "dna11161_transcript_rnaseq.RData"))

} else { load(myfn) }

## clinical information
dd <- read.csv(file.path("clinic_info", "Clinical Data PNC Seq with survival data_FINAL.csv"), row.names=NULL, stringsAsFactors=FALSE)
dd[dd == "UK" | dd == "" | dd == " "] <- NA
dd[ , "sample_id"] <- gsub("^LA", "LUMA", dd[ , "sample_id"])
dd[ , "sample_id"] <- gsub("^LB", "LUMB", dd[ , "sample_id"])
nn <- paste(gsub("[ ]|-", "_", gsub(" $", "", dd[ ,1])), "T", sep="_")
rownames(dd) <- nn
dd <- data.frame("samplename"=rownames(dd), dd)
nn <- sort(unique(c(rownames(sampleinfo.affy), rownames(sampleinfo.rnaseq), rownames(dd))))
demo <- dd

# demo <- data.frame(matrix(NA, nrow=length(nn), ncol=ncol(dd), dimnames=list(nn, colnames(dd))), stringsAsFactors=FALSE)
# demo <- setcolclass.df(df=demo, colclass=sapply(dd, class), factor.levels=sapply(dd, levels))
# demo[rownames(dd), ] <- dd

write.csv(demo, file=file.path(saveres, "DNA11161_demo.csv"), row.names=FALSE)
save(list=c("demo"), compress=TRUE, file=file.path(saveres, "DNA11161_demo.RData"))


## save all data for GEO
load(file.path(saveres, "dna11161_affy_frma.RData"))
load(file.path(saveres, "dna11161_gene_rnaseq.RData"))
load(file.path(saveres, "dna11161_transcript_rnaseq.RData"))
save(list=c("data.affy", "annot.affy", "data.gene.rnaseq", "annot.gene.rnaseq", "data.transcript.rnaseq", "annot.transcript.rnaseq", "dataf.gene.rnaseq", "dataf.transcript.rnaseq", "sampleinfo.affy", "sampleinfo.rnaseq"), compress=TRUE, file=file.path(saveres, "dna11161_data_all_full.RData"))
## keep only the common 57 patients
nn <- fold(intersect, rownames(data.affy), rownames(data.gene.rnaseq), rownames(data.transcript.rnaseq), rownames(dataf.gene.rnaseq), rownames(dataf.transcript.rnaseq), rownames(demo))
data.affy <- data.affy[nn, , drop=FALSE]
data.gene.rnaseq <- data.gene.rnaseq[nn, , drop=FALSE]
data.transcript.rnaseq <- data.transcript.rnaseq[nn, , drop=FALSE]
dataf.gene.rnaseq <- dataf.gene.rnaseq[nn, , drop=FALSE]
dataf.transcript.rnaseq <- dataf.transcript.rnaseq[nn, , drop=FALSE]
sampleinfo.affy <- sampleinfo.affy[nn, , drop=FALSE]
sampleinfo.rnaseq <- sampleinfo.rnaseq[nn, , drop=FALSE]
demo <- demo[nn, , drop=FALSE]
save(list=c("demo", "data.affy", "annot.affy", "data.gene.rnaseq", "annot.gene.rnaseq", "data.transcript.rnaseq", "dataf.gene.rnaseq", "dataf.transcript.rnaseq", "annot.transcript.rnaseq", "sampleinfo.affy", "sampleinfo.rnaseq"), compress=TRUE, file=file.path(saveres, "dna11161_data_all.RData"))



## end

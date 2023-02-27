library(bsseq)
library(DSS)
library(tidyverse)
makeBSseqData <- function (dat, sampleNames) {
  n0 <- length(dat)
  if (missing(sampleNames)) 
    sampleNames <- paste("sample", 1:n0, sep = "")
  alldat <- dat[[1]]
  if (any(alldat[, "N"] < alldat[, "X"], na.rm = TRUE)) 
    stop("Some methylation counts are greater than coverage.\n")
  ix.X <- which(colnames(alldat) == "X")
  ix.N <- which(colnames(alldat) == "N")
  colnames(alldat)[ix.X] <- "X1"
  colnames(alldat)[ix.N] <- "N1"
  if (n0 > 1) {
    for (i in 2:n0) {
      thisdat <- dat[[i]]
      if (any(thisdat[, "N"] < thisdat[, "X"], na.rm = TRUE)) 
        stop("Some methylation counts are greater than coverage.\n")
      ix.X <- which(colnames(thisdat) == "X")
      ix.N <- which(colnames(thisdat) == "N")
      colnames(thisdat)[c(ix.X, ix.N)] <- paste(c("X", 
                                                  "N"), i, sep = "")
      alldat <- merge(alldat, thisdat, all = TRUE, by = c("chr", "pos"))
    }
  }
  alldat <- alldat[order(alldat$chr, alldat$pos), ]
  ix.X <- grep("X", colnames(alldat))
  ix.N <- grep("N", colnames(alldat))
  alldat[is.na(alldat)] <- 0
  M <- as.matrix(alldat[, ix.X, drop = FALSE])
  Cov <- as.matrix(alldat[, ix.N, drop = FALSE])
  colnames(M) <- colnames(Cov) <- sampleNames
  idx <- split(1:length(alldat$chr), alldat$chr)
  M.ordered <- M
  Cov.ordered <- Cov
  pos.ordered <- alldat$pos
  for (i in seq(along = idx)) {
    thisidx = idx[[i]]
    thispos = alldat$pos[thisidx]
    dd = diff(thispos)
    if (min(dd) < 0) {
      warning(paste0("CG positions in chromosome ", names(idx)[i], 
                     " is not ordered. Reorder CG sites.\n"))
      iii = order(thispos)
      M.ordered[thisidx, ] <- M[thisidx, ][iii, ]
      Cov.ordered[thisidx, ] <- Cov[thisidx, ][iii, ]
      pos.ordered[thisidx] <- alldat$pos[thisidx][iii]
    }
  }
  result <- BSseq(chr = alldat$chr, pos = pos.ordered, M = M.ordered, 
                  Cov = Cov.ordered)
  result
}

prepare_for_DSS<-function(query){
  colnames(query) <-c("chr", "pos", "N", "X")
  query$chr <- paste(rep("chr",nrow(query)), query$chr, sep = "")
  query$N<-as.numeric(paste(query$N))
  query$X<-as.numeric(paste(query$X))
  return(query)
}

#mC
C0<-read.table(file="~/C0/mCmethylation.bsseq.tsv,  header=T, sep="\t")

C1<-read.table(file="~/C1/mCmethylation.bsseq.tsv",  header=T, sep="\t")

C2<-read.table(file="~/C2/mCmethylation.bsseq.tsv", header=T, sep="\t")
  
C3<-read.table(file="~/C3/mCmethylation.bsseq.tsv",  header=T, sep="\t")

CORT12<-read.table(file="~/CORT12/mCmethylation.bsseq.tsv",  header=T, sep="\t")

CORT13<-read.table(file="~/CORT13/mCmethylation.bsseq.tsv", header=T, sep="\t")

CORT14<-read.table(file="~/CORT14/mCmethylation.bsseq.tsv",  header=T, sep="\t")

CORT15<-read.table(file="~/CORT15/mCmethylation.bsseq.tsv",  header=T, sep="\t")

C0_DSS<-prepare_for_DSS(C0)
C1_DSS<-prepare_for_DSS(C1)
C2_DSS<-prepare_for_DSS(C2)
C3_DSS<-prepare_for_DSS(C3)
CORT12_DSS<-prepare_for_DSS(CORT12)
CORT13_DSS<-prepare_for_DSS(CORT13)
CORT14_DSS<-prepare_for_DSS(CORT14)
CORT15_DSS<-prepare_for_DSS(CORT15)

rm(C0, C1, C2, C3, CORT12, CORT13, CORT14, CORT15)

BSobjmc <- makeBSseqData( list(C0_DSS, C1_DSS, C2_DSS, C3_DSS, CORT12_DSS, CORT13_DSS, CORT14_DSS, CORT15_DSS),
                        c("C0","C1", "C2", "C3", "CORT12", "CORT13", "CORT14", "CORT15") )

rm(C0_DSS, C1_DSS, C2_DSS, C3_DSS, CORT12_DSS, CORT13_DSS, CORT14_DSS, CORT15_DSS)

dmlTest= DMLtest(BSobjhmc, group1=c("CORT12", "CORT13", "CORT14", "CORT15"), group2=c("C0","C1", "C2", "C3"), smoothing=FALSE)

dmrs = callDMR(dmlTest, p.threshold=0.05)

write.table(dmlTest, file= "dmlTest_mCG.tsv", append = FALSE, sep = "\t", dec = ".", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(dmrs, file= "/dmrs.tsv", append = FALSE, sep = "\t", dec = ".", quote = FALSE, row.names = FALSE, col.names = TRUE)

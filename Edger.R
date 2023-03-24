library(EdgeR) # This script comes from Edger user guide and has been adapted for my data.

Group <- c(rep("Control", 4), rep( "CORT", 4))
Sample<-c("C0", "C1", "C2", "C3","CORT12", "CORT13", "CORT14", "CORT15")
File<-c("~/C0/methylation.edgeR.tsv","~/C1/methylation.edgeR.tsv","~/C2/methylation.edgeR.tsv","~/C3/methylation.edgeR.tsv","~/CORT12/methylation.edgeR.tsv","~/CORT13/methylation.edgeR.tsv","~/CORT14/methylation.edgeR.
tsv","~/CORT15/methylation.edgeR.tsv")
MedianQuality<-c(12.8, 14.0 ,11.7 ,12.7 ,7.3 ,11.6 ,13.6 , 13.3)
Medianreadlength<-c(4980, 5870, 5763, 5888, 5696, 7956, 6099, 6277)     
batch<-c(2,2,1,0,1,0,2,2)
designSL <- model.matrix(~0+Group+MedianQuality+Medianreadlength+batch)
design <- modelMatrixMeth(designSL)
yall <- readBismark2DGE(File, sample.names=Sample)
TSS <- nearestTSS(yall$genes$Chr, yall$genes$Locus, species="Mm")
yall$genes$EntrezID <- TSS$gene_id
yall$genes$Symbol <- TSS$symbol
yall$genes$Strand <- TSS$strand
yall$genes$Distance <- TSS$distance
yall$genes$Width <- TSS$width
Methylation <- gl(2,1,ncol(yall), labels=c("Me","Un"))
Me <- yall$counts[, Methylation=="Me"]
Un <- yall$counts[, Methylation=="Un"]
Coverage <- Me + Un

#KEEP only positions where all teh samples(8) have at least 5 coverage.
HasCoverage <- rowSums(Coverage >= 5) ==8
HasBoth <- rowSums(Me) > 0 & rowSums(Un) > 0
table(HasCoverage, HasBoth)

y <- yall[HasCoverage & HasBoth,, keep.lib.sizes=FALSE]

#now we have to adjust the library sizes so all samples are the same for me and un.
TotalLibSize <- y$samples$lib.size[Methylation=="Me"] + y$samples$lib.size[Methylation=="Un"]
y$samples$lib.size <- rep(TotalLibSize, each=2)

y <- estimateDisp(y, design=design, trend="none")
contr <- makeContrasts(GroupCORTvsControl = GroupCORT - GroupControl, levels=design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, contrast=contr)
summary(decideTests(lrt))
tags<-topTags(lrt, n=nrow(lrt2$table))

write.table(tags$table, "~/genesltrF0DE_edger.tsv")



#Promoter analysis
InPromoter <- yall$genes$Distance >= -1000 & yall$genes$Distance <= 2000
yIP <- yall[InPromoter,,keep.lib.sizes=FALSE]
yIP$samples$group<-c(rep(1,8), rep(2,8))
yIP$samples$group<-as.factor(yIP$samples$group)
yIP<-yIP[!is.na(rownames(yIP)),] #clean NA  
ypr <- rowsum(yIP, yIP$genes$EntrezID , reorder=FALSE) 

#filter coverage 
Mepr <- ypr$counts[,Methylation=="Me"]
Unpr <- ypr$counts[,Methylation=="Un"]
Coveragepr <- Mepr + Unpr
HasCoveragepr <- rowSums(Coveragepr >= 8) == 12
HasCoveragepr <- rowSums(Coveragepr >= 8) == 8
HasBothpr <- rowSums(Mepr) > 0 & rowSums(Unpr) > 0

ypr <- ypr[HasCoveragepr & HasBothpr,,keep.lib.sizes=FALSE]

TotalLibSizepr <- 0.5*ypr$samples$lib.size[Methylation=="Me"] +
  + 0.5*ypr$samples$lib.size[Methylation=="Un"]
ypr$samples$lib.size <- rep(TotalLibSizepr, each=2)

ypr <- estimateDisp(ypr, design, trend="none")

fitpr <- glmFit(ypr, design)
lrtpr <- glmLRT(fitpr, contrast=contr)

#get a table with promoter and methylation change
tags<-topTags(lrtpr, n=nrow(lrtpr$table))

#get table of significantly change promoters
write.table(tags$table, file= "~/promoters_mc_norm2_edger.tsv", append = FALSE, sep = "\t", dec = ".", quote = FALSE, row.names = FALSE, col.names = T)
sigpr<-tags[tags$table$FDR<0.05,]
write.table(sigpr$table, file= "~/sig_promoters_mc_norm2_edger.tsv", append = FALSE, sep = "\t", dec = ".", quote = FALSE, row.names = FALSE, col.names = T)

#plot and ttest from a count matrix.
#normdata would be a count matrix with samples as columns and CpGs as lines.
#grtabla is a genomic ranges object with all the CpG positions from normdata. Ibtained from the BSobject.

extractregion<-function(matrix, chr, start, end){
  gr<- GRanges(seqnames = chr,
                   ranges = IRanges(start = start, end=end))
  overlaps<-matrix[queryHits(findOverlaps(grtabla, gr)),]
  means<-colMeans(overlaps)
  data<-data.frame(means, c(rep("Control",4),rep("CORT",4)))
  Control<-data[data$c.rep..Control...4...rep..CORT...4..=="Control",]
  CORT<-data[data$c.rep..Control...4...rep..CORT...4..=="CORT",]
  df2<-data.frame(Control, CORT)
  ttest<-t.test(df2$means,df2$means.1)

  plot<-ggplot(data, aes(x=c.rep..Control...4...rep..CORT...4.., y=means, fill=c.rep..Control...4...rep..CORT...4..)) +
    geom_boxplot() + geom_point()
  return(list(ttest, plot))
}
#Igf2
extractregion(normdata, "chr7", 142575503, 142582140)
extractregion(normdata[3:10], "chr7", 142575503, 142582140)
extractregion(normdata, "chrchr7", 142575503, 142582140)
extractregion(df[3:14], "chrchr7", 142575503, 142582140)

#GR1 chr18 39490671 39490731
extractregion(normdata, "chr18", 39490671, 39490731)
extractregion(normdata_01, "chrchr18", 39490671, 39490731)
#GR2 #chr18 39489794 39489854
extractregion(normdata_01, "chrchr18", 39489794, 39489854)
extractregion(df[3:14], "chrchr18", 39489794, 39489854)

#RASGRF1
matrix<-df[3:14]
extractregion(normdata, "chr9", 89879568, 89880045)
extractregion(normdata, "chrchr9", 89879568, 89880045)
extractregion(df[3:14], "chrchr9", 89879568, 89880045)
extractregion(matrix, "chrchr9", 89879568, 89880045)

#GTL2/DLK1
extractregion(normdata, "chrchr12", 109526740, 109528845)
extractregion(normdata, "chr12", 109526740, 109528845)

#BDNF
#chr2 109675944 109676004 +      promoters       Bdnf_3  .
extractregion(normdata_01[3:14], "chrchr2", 109675944, 109676004)
extractregion(normdata_01[3:10], "chr2", 109675944, 109676004)
#chr2 109677035 109677095 +      promoters       Bdnf_2  .
extractregion(normdata_01[3:14], "chrchr2", 109677035, 109677095)
extractregion(normdata_01[3:10], "chr2", 109677035, 109677095)
#chr2 109692418 109692478 +      promoters       Bdnf_4  .
extractregion(normdata_01[3:14], "chrchr2", 109692418, 109692478)
extractregion(normdata_01[3:10], "chr2", 109692418, 109692478)
#chr2 109693566 109693626 +      promoters       Bdnf_1  .
extractregion(normdata_01[3:14], "chrchr2", 109693566, 109693626)
extractregion(normdata_01[3:10], "chr2", 109693566, 109693626)
#chr2 109720465 109720525 +      promoters       Bdnf_5  .
extractregion(normdata_01[3:14], "chrchr2", 109720465, 109720525)
extractregion(normdata_01[3:10], "chr2", 109720465, 109720525)


#spaguetti plots

dmr<-genetable[queryHits(findOverlaps(genegr, dmrgr)),]
Position<-(rep(gene[,3], 8))
gene<-gene[4:11]
gene<-dmr[4:11]
gene<-dmr[4:15]
colnames(gene)<-c("C0","C1", "C2", "C3", "CORT12", "CORT13", "CORT14", "CORT15")
colnames(gene)<-c("C4413","C4414", "C4423", "C4424", "C453", "C454", "PatCORT171", "PatCORT172", "PatCORT181", "PatCORT182", "PatCORT191", "PatCORT192")
Treatment = c(rep("Control",11600), rep( "CORT",11600))
Treatment = c(rep("Control",24), rep( "CORT",24))
Treatment = c(rep("Control",36), rep( "CORT",36))
gene<-gene %>% 
  gather(key="Sample", value="Normalized frequency", )
gene$Treatment<-Treatment
Position<-c(rep(gene$V2,8))
Position<-c(rep(dmr$V2,8))
Position<-c(rep(dmr$V2,12))
gene$Position<-Position
gene$Position<-as.numeric(gene$Position)
#gene$Position<-rep(1:2900, 8)


# Make the plot
ggplot(data=gene, aes(x=Position, y=`Normalized frequency`, color=Treatment)) + 
  geom_line() 

ggplot(data=gene, aes(x=Position, y=`Normalized frequency`, color=Sample, linetype=Treatment)) + 
  geom_line()

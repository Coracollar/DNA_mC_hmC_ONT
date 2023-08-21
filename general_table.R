
library(tidyverse)
library(plyranges)
#all_suptable_simple_intronsintergenic.tsv files are obtained from overlapping our methylation frequency matrix by CpG with a table containing info for our regions of interest. 
#The ID column is obtained by colapsing the first three columns of the resulting bed file to obtain something like: chr:start-end.
met<-read.table("met_all_suptable_simple_intronsintergenic.tsv")
hmet<-read.table("hmet_all_suptable_simple_intronsintergenic.tsv")

colnames(met)<-c("ID", "type", "subtype", "name", "C0","C1", "C2", "C3", "CORT12", "CORT13", "CORT14", "CORT15", "C4413","C4414", "C4423", "C4424", "C453", "C454", "PatCORT171", "PatCORT172", "PatCORT181", "PatCORT182", "PatCORT191", "PatCORT192")
colnames(hmet)<-c("ID", "type", "subtype", "name", "C0","C1", "C2", "C3", "CORT12", "CORT13", "CORT14", "CORT15", "C4413","C4414", "C4423", "C4424", "C453", "C454", "PatCORT171", "PatCORT172", "PatCORT181", "PatCORT182", "PatCORT191", "PatCORT192")

#collapse CpGs for a particular ID, which corresponds to individual genes or CpG islands, or line element, etc.
df1<-met %>% group_by(ID)%>% summarise(across(where(is.numeric), mean), type) %>% 
  distinct()
df1h<-hmet %>% group_by(ID)%>% summarise(across(where(is.numeric), mean), type) %>% 
  distinct()

df1<-as.data.frame(df1)
dfh<-as.data.frame(df1h)
#choose genomic regions
df3<-df1[df1$type=="gene" | df1$type=="Intergenic" | df1$type=="promoters" | df1$type=="5UTR"| df1$type=="CDS" | df1$type=="intron" | df1$type=="3UTR" | df1$type=="enhancer"  | df1$type=="CpGisland" | df1$type=="ZFP57_peaks" | df1$type=="CTCFBS" | df1$type=="ICR" | df1$type=="LINE"|
             df1$type=="LTR" | df1$type=="Simple_repeat" | df1$type=="tRNA" | df1$type=="SINE",]

df3h<-df1h[df1h$type=="gene" | df1h$type=="Intergenic" | df1h$type=="promoters" | df1h$type=="5UTR"| df1h$type=="CDS" | df1h$type=="intron" | df1h$type=="3UTR" | df1h$type=="enhancer"  | df1h$type=="CpGisland" | df1h$type=="ZFP57_peaks" | df1h$type=="CTCFBS" | df1h$type=="ICR" | df1h$type=="LINE"|
          df1h$type=="LTR" | df1h$type=="Simple_repeat" | df1h$type=="tRNA" | df1h$type=="SINE",]

#covert type into factor
df3h$type<-factor(df3h$type, levels=c("gene","5UTR", "CDS", "3UTR", "Intergenic", "promoters", "enhancer", "CpGisland", "ICR", "tRNA",
                                      "CTCFBS","ZFP57_peaks","LINE", "SINE", "LTR", "Simple_repeat" ))
df3$type<-factor(df3$type, levels=c("gene","5UTR", "CDS", "3UTR", "Intergenic", "promoters", "enhancer", "CpGisland", "ICR", "tRNA",
                 "CTCFBS","ZFP57_peaks","LINE", "SINE", "LTR", "Simple_repeat" ))


#get mean methylation per regions of our experimental groups
df3h_means<-data.frame(df3h$type, rowMeans(df3h[2:5]), rowMeans(df3h[6:9]), rowMeans(df3h[10:15]), rowMeans(df3h[16:21]))
df3_means<-data.frame(df3$type, rowMeans(df3[2:5]), rowMeans(df3[6:9]), rowMeans(df3[10:15]), rowMeans(df3[16:21]) )
colnames(df3_means)<-c("Genomic region", "F0 Control", "F0 CORT", "F1 Control", "F1 PatCORT")
colnames(df3h_means)<-c("Genomic region", "F0 Control", "F0 CORT", "F1 Control", "F1 PatCORT")

#plot
#transform data fram so it is suitable for ggplot
densitydf<-df3_means %>% gather(key="Treatment", value="Methylation frequency", `F0 Control`, `F0 CORT`, `F1 Control`, `F1 PatCORT`)
densitydfh<-df3h_means %>% gather(key="Treatment", value="Methylation frequency", `F0 Control`, `F0 CORT`, `F1 Control`, `F1 PatCORT`)

#actual plot
ggplot(densitydf) + geom_density(aes(x=`Methylation frequency`, col=Treatment))  +
+   facet_wrap(~`Genomic region`,scales = "free") 
ggplot(densitydfh) + geom_density(aes(x=`Methylation frequency`, col=Treatment))  +
+   facet_wrap(~`Genomic region`,scales = "free") 

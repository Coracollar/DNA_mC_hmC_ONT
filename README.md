# DNA_mC_hmC_ONT

## Data preparation and methylation calls. 

Fast5s from each sample need to be located in a separate folder. For each sample we run megalodon_mC_hmC.slurm.

We obtain methylation and hydroxymethylation frequency files from megalodon in the form of bed files.
```bash
head modified_bases.5mC.bed
chr10	3125029	3125030	.	11	+	3125029	3125030	0,0,0	11	90.9
chr10	3125030	3125031	.	8	-	3125030	3125031	0,0,0	8	100.0
chr10	3125246	3125247	.	14	+	3125246	3125247	0,0,0	14	85.7
chr10	3125247	3125248	.	8	-	3125247	3125248	0,0,0	8	100.0
chr10	3125654	3125655	.	10	+	3125654	3125655	0,0,0	10	80.0
chr10	3125655	3125656	.	9	-	3125655	3125656	0,0,0	9	88.9
chr10	3125928	3125929	.	15	+	3125928	3125929	0,0,0	15	80.0
chr10	3125929	3125930	.	10	-	3125929	3125930	0,0,0	10	100.0
chr10	3125940	3125941	.	14	+	3125940	3125941	0,0,0	14	92.9
chr10	3125941	3125942	.	9	-	3125941	3125942	0,0,0	9	100.0
```
Sequencing data stats were obtained using NanoStat.

## PCA analysis.

For PCA analysis we created a matrix of counts per sample per CpG. This can be obtained directly from a bsseq object (BSobj) like the one created in DSS.R.

```R
library(factoextra)
methylation_data<-getMeth(BSobj, type="raw")
methylation_data<-na.omit(methylation_data)
PCA<-prcomp(t(methylation_data), scale = FALSE)

#visualize PCA
fviz_eig(PCA) # inspect how many components there are and their contribution to the total.
fviz_pca_ind(PCA,
             axes = c(1, 2), #Number of the components to be plotted
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
```
For ONT data, read length and quality have been known to affect methylation calls. Moreover, batch effect or other features can be hiding the effect of the particular variable we are studying.
In order to find which variables might be affecting PCA, we can obtain an R squared value in two different ways, with a correlation test for numeric variables or analysis of variance summary for non-numeric variables.

```R
#create a design. Example from F0 data.
Treatment = c(rep("Control",4), rep( "CORT",4))
Quality = c(12.0 , 13.0 , 10.9 , 11.5 , 8.2 ,10.9 , 12.6 , 12.3)
sample<-c("C0","C1", "C2", "C3", "CORT12", "CORT13", "CORT14", "CORT15")
readlenghth<-c(7141,	8066,	7430,	7907,	7461,	11384,	8137,	8806)
batch<-c(2,2,1,0,1,0,2,2)
design = data.frame(sample, Treatment, Quality, readlenghth, batch)
PCA_df <- as.data.frame(PCA$x)

PCA_metadata <- cbind(design,PCA_df)

#get R
cor.test(pca_adj_metadata$PC1, PCA_metadata$batch, method="spearman") # for numeric values
summary.lm(aov(PC1~ Treatment, PCA_metadata))$adj.r.squared # For non numeric values

biol <- model.matrix(~ design$Treatment)
covar <- model.matrix(~ design$MedianQuality + design$Medianreadlength + design$batch)

#Transform data  
adjusted_counts <- limma::removeBatchEffect(methylation_data, covariates=qual, design=biol)
```
PCA of adjusted counts will show us whether our transformation has improved our separation due to Group/Treatment. Once we know the covariates that are important from our data we can use them to improve our discovery power or our variable of interest.


## DMR analysis using DSS.
In order to transform the data for statistical analysis using bsseq we create two new files for each sample, one for mC and other for hmC.
```bash
awk '{print $1, $2, $10, $4=int($11*$10/100+0.5) }' modified_bases.5mC.bed > mCmethylation.bsseq.tsv
awk '{print $1, $2, $10, $4=int($11*$10/100+0.5) }' modified_bases.5hmC.bed > hmCmethylation.bsseq.tsv

head mCmethylation.bsseq.tsv
chr10 3125029 11 9
chr10 3125030 8 8
chr10 3125246 14 11
chr10 3125247 8 8
chr10 3125654 10 8
chr10 3125655 9 8
chr10 3125928 15 12
chr10 3125929 10 10
chr10 3125940 14 13
chr10 3125941 9 9
```
Once we have created all necessary .bsseq.tsv we can run DSS.R to find Differentially methylated regions. Two scripts are provided: DSS.R DSS_design.R. DSS.R does not take into account any known covariates, is the most basic test; DSS_design.R. Same scripts were used for mC and hmC analysis. Scripts need to be adapted depending on the number of samples and the path where they are stored. The data exemplified in our scripts in our F0 sperm data (4vs4). R scripts can be run as bash jobs for large datasets using the script DSS.slurm.


## Promoter and whole gene analysis using EdgeR.

Again, we transformed our .bed files from megalodon to fit EdgeR model. Once the methylation.edgeR.tsv file are created the script Edger.R was run ( the script provided was the one run for methylation of F0, it has to be modified for each set of samples).
```bash
awk '{print $1, $2, $3, $11, int($10*$11/100 + 0.5), int($10*(1-$11/100) + 0.5)}' modified_bases.5mC.bed | sed 's/^/chr/' | sed 's/ /\t/g' > methylation.edgeR.tsv

awk '{print $1, $2, $3, $11, int($10*$11/100 + 0.5), int($10*(1-$11/100) + 0.5)}' modified_bases.5hmC.bed | sed 's/^/chr/' | sed 's/ /\t/g' > hmethylation.edgeR.tsv
```

The different plots for methylation data can be found in the script plots.R.

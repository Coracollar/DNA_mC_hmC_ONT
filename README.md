# DNA_mC_hmC_ONT

## Data preparation and methylation 

Fast5s from each sample need to be located in a separate folder. For each sample we run megalodon_mC_hmC.slurm.

We obtain methylation and hydroxymethylation frequency files from megalodon in the form of be files.
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
Once we have created all necessary .bsseq.tsv we can run DSS.R to find Differentially methylated regions. Two scripts are provided: DSS.R DSS_design.R. DSS.R does not take into account any known covariates, is the most basic test; DSS_design.R. Same scripts were used for mC and hmC analysis. Scripts need to be adapted depending on the number of samples and the path where they are stored. The data exemplified in our scripts in our F0 sperm data (4vs4).



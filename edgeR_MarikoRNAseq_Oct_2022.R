# Analysis of Mariko RNA_seq data of Daphnia from infected = treated and uninfected = control at different time points


library(ggplot2)
library(edgeR)
#source("http://bioconductor.org/biocLite.R")
#biocLite()
# biocLite("edgeR")  https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
library(data.table)
library(tidyverse)

#####
#
# User-defined variables
#
#####

# set working directory; must only contain directories containing eXpress output, no other files or directories 
setwd("/Users/j.stillman/switchdrive/Jonathon_Analysis_of_Mariko_RNAseq/count/")
# set output directory
out_dir <- "/Users/j.stillman/switchdrive/Jonathon_Analysis_of_Mariko_RNAseq/"
# get annotation file
# "trinotate_annotation_report.txt" is a Trinotate output excel file, exported as tab-delimited
transcripts <- read.delim("/Users/jhstillman/switchdrive/Mariko_RNAseq/functional_annotation_link.for_Jonathon_26072022.txt",stringsAsFactors=FALSE,header=FALSE) # annotations sent from Peter.  
#Need to modify a column to match to gene### in the mapped reads.  According to Peter: "Here gene and rna should be the same (I think), so if you see gene0 it is the same as rna0. So if you were to do sed 's/rna/gene/g' file1.txt > file2.txt I think that should match them up with the count table you have. The APZ42_#### is the Uniprot ID for the protein."
library(stringr)
transcripts$V1 = str_replace(transcripts$V1,"rna","gene")
names(transcripts)=c("gene","description","UniprotID")
## note, there are fewer transcripts than there are genes (26525 transcripts annotated vs. 26646 genes that data were mapped to)

#group
##########
#
# Import expression estimates, total counts
#
##########

# creates a vector of the sequencing library directory names, "libs"
libs <- list.files()

# load in the files
counts = fread(file[1],skip=4,header=FALSE,col.names=c("gene",file[1]),colClasses = c("character","numeric","NULL","NULL"))
counter=0

for (file in libs) {
  cat('reading in', file, '\n')
  tempframe <- fread(file,skip=4,header=FALSE,colClasses = c("numeric"),col.names=c(file),select=c("V2"))
  assign(file, tempframe[order(tempframe$gene),]) 
  counts=cbind(counts, tempframe) # generate dataframe, organizing rows by target_id
}
counts=counts[,-2]

# check sums of counts
sums=colSums(counts[,-1])
sums2=sums/10^6

# Make an FPKM dataframe for heatmap and plots, "all_fpkm"
fpkm=cbind(as.data.frame(mapply('/', counts[,-1], sums2)))
fpkmsums=colSums(fpkm[,-1]) # sanity check to be sure that each column sums to 1e+06 - they do.


# create row names from genes
  all_count=counts
  rownames(all_count)=all_count[,1]
  rownames(all_count) # confirm that they are genes
  all_count=all_count[,-1]


# now the data is ready.  Confirm similar heterogeneity and distribution with boxplot:

boxplot(log2(all_count[,-1]+1), las=2) # before filter

# define groups
samples = colnames(all_count)
time = c(rep(0.125,9),rep(0.25,9),rep(0.5,9),rep(1,9),rep(2,9),rep(4,9),rep(6,9),rep(8,9),rep(16,9))
treatment=c(rep("uninfected",3),rep("infected",6),rep("uninfected",3),rep("infected",6),rep("uninfected",3),rep("infected",6),rep("uninfected",3),rep("infected",6),rep("uninfected",3),rep("infected",6),rep("uninfected",3),rep("infected",6),rep("uninfected",3),rep("infected",6),rep("uninfected",3),rep("infected",6),rep("uninfected",3),rep("infected",6))
group = data.frame(samples,as.factor(time),as.factor(treatment))
names(group)=c("samples","time","treatment")
combtreat=paste(time,substr(treatment,1,1),sep = "_")

group2 = data.frame(rbind(as.factor(time),as.factor(treatment)))
colnames(group2)=samples
rownames(group2)=c("time","treatment")

##########
#
# Set up differential gene expression (DGE) analysis, generate summary plots:
#   define a model matrix *
#   estimate dispersion
#   make a plot of biological coefficient of variation vs. log(CPM)
#   make a multidimensional scaling plot of total counts
#   fit GLM for each transcript
#     * Note, how to define the statistical model is completely dependent on your experimental design and what question(s) you want to address.
#
##########
all_count_m = as.matrix(all_count)
group_m = as.matrix(group)
group2_m = as.matrix(group2)
# DGEList makes an EdgeR object
y <- DGEList(counts=all_count, group=group2) 


# define the statistical model, a design matrix using the model.matrix function
design <- model.matrix(~group, data=y$samples)
colnames(design) <- levels(y$samples$group)
colnames(design)
design
# filter out transcripts that do not have at least two reads out of 1,000,000 reads mapped in at least 4 samples
y <- y[rowSums(1e+06 * y$counts/expandAsMatrix(y$samples$lib.size, dim(y)) > 2) >= 4, ]
# reset the library sizes after filtering
y$samples$lib.size <- colSums(y$counts)
# TMM normalization compensates not just for library size but also the relative expression level among transcripts
y <- calcNormFactors(y, method="upperquartile")
# estimate dispersion for GLM fit, common, trended, and tagwise
# see edgeR docs for references on dispersion estimates, e.g., "?estimateGLMCommonDisp"
y <- estimateGLMCommonDisp(y,design) 
# estimate trended dispersion for use in tagwise dispersion estimate
y <- estimateGLMTrendedDisp(y,design) 
# estimate tagwise dispersion to be used in glmFit()
y <- estimateGLMTagwiseDisp(y,design)
y$samples
# plot genewise biological coefficient of variation against gene abundance
plotBCV(y)
# make plot as a pdf
pdf(file=paste(out_dir, "BCV_plot.pdf", sep=""), height=6, width=6)
plotBCV(y, main = "Biological Coefficient of Variation")
dev.off()
# MDS plot
plotMDS(y , main = "MDS Plot for Count Data", labels = combtreat, cex=0.7)
# make a pdf
pdf(file=paste(out_dir, "MDS_plot.pdf", sep=""), height=6, width=6)
plotMDS(y , main = "MDS Plot for Count Data", labels = colnames( y$counts ), cex=0.7)
plotMDS(y , main = "MDS Plot for Count Data", labels = group, cex=0.7)
dev.off()

# GLM
# fit the negative binomial generalized linear model (GLM) for each tag, creating a new fit object
fit <- glmFit(y, design=design)
colnames(fit)


####
# 
# identify DE Genes
#
####

# first - within each treatment by time:
T1v2_u <- glmLRT(fit, contrast=c(-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) # Time 1 vs 2 uninfected
T2v3_u <- glmLRT(fit, contrast=c(0,0,-1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0)) # Time 2 vs 3 uninfected
T3v4_u <- glmLRT(fit, contrast=c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,-1,0,0,0)) # Time 3 vs 4 uninfected
T4v5_u <- glmLRT(fit, contrast=c(0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,1,0,0)) # Time 4 vs 5 uninfected
T5v6_u <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,-1,0,0)) # Time 5 vs 6 uninfected
T6v7_u <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,1,0)) # Time 6 vs 7 uninfected
T7v8_u <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,-1,0)) # Time 7 vs 8 uninfected
T8v9_u <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,1)) # Time 8 vs 9 uninfected

summary(decideTestsDGE(T1v2_u, p=0.05, adjust="BH"))
summary(decideTestsDGE(T2v3_u, p=0.05, adjust="BH"))
summary(decideTestsDGE(T3v4_u, p=0.05, adjust="BH"))
summary(decideTestsDGE(T4v5_u, p=0.05, adjust="BH"))
summary(decideTestsDGE(T5v6_u, p=0.05, adjust="BH"))
summary(decideTestsDGE(T6v7_u, p=0.05, adjust="BH"))
summary(decideTestsDGE(T7v8_u, p=0.05, adjust="BH"))
summary(decideTestsDGE(T8v9_u, p=0.05, adjust="BH"))

T1v2_i <- glmLRT(fit, contrast=c(0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) # Time 1 vs 2 infected
T2v3_i <- glmLRT(fit, contrast=c(0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0)) # Time 2 vs 3 infected
T3v4_i <- glmLRT(fit, contrast=c(0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0)) # Time 3 vs 4 infected
T4v5_i <- glmLRT(fit, contrast=c(0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0)) # Time 4 vs 5 infected
T5v6_i <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0)) # Time 5 vs 6 infected
T6v7_i <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0)) # Time 6 vs 7 infected
T7v8_i <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0)) # Time 7 vs 8 infected
T8v9_i <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0)) # Time 8 vs 9 infected

summary(decideTestsDGE(T1v2_i, p=0.05, adjust="BH"))
summary(decideTestsDGE(T2v3_i, p=0.05, adjust="BH"))
summary(decideTestsDGE(T3v4_i, p=0.05, adjust="BH"))
summary(decideTestsDGE(T4v5_i, p=0.05, adjust="BH"))
summary(decideTestsDGE(T5v6_i, p=0.05, adjust="BH"))
summary(decideTestsDGE(T6v7_i, p=0.05, adjust="BH"))
summary(decideTestsDGE(T7v8_i, p=0.05, adjust="BH"))
summary(decideTestsDGE(T8v9_i, p=0.05, adjust="BH"))

# Next compare uninfected vs. infected at each time point
T1u_i <- glmLRT(fit, contrast=c(-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) # Time 1 uninfected vs infected
T2u_i <- glmLRT(fit, contrast=c(0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) # Time 2 uninfected vs infected
T3u_i <- glmLRT(fit, contrast=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,-1,0,0,0)) # Time 3 uninfected vs infected
T4u_i <- glmLRT(fit, contrast=c(0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0)) # Time 4 uninfected vs infected
T5u_i <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,-1,0,0)) # Time 5 uninfected vs infected
T6u_i <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0)) # Time 6 uninfected vs infected
T7u_i <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,-1,0)) # Time 7 uninfected vs infected
T8u_i <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0)) # Time 8 uninfected vs infected
T9u_i <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,-1)) # Time 9 uninfected vs infected


summary(decideTestsDGE(T1u_i, p=0.05, adjust="BH"))
summary(decideTestsDGE(T2u_i, p=0.05, adjust="BH"))
summary(decideTestsDGE(T3u_i, p=0.05, adjust="BH"))
summary(decideTestsDGE(T4u_i, p=0.05, adjust="BH"))
summary(decideTestsDGE(T5u_i, p=0.05, adjust="BH"))
summary(decideTestsDGE(T6u_i, p=0.05, adjust="BH"))
summary(decideTestsDGE(T7u_i, p=0.05, adjust="BH"))
summary(decideTestsDGE(T8u_i, p=0.05, adjust="BH"))
summary(decideTestsDGE(T9u_i, p=0.05, adjust="BH"))

# Tablulate expression data for each comparison
aa=data.frame(topTags(T1v2_u,n=14256))
bb=data.frame(topTags(T2v3_u,n=14256))
cc=data.frame(topTags(T3v4_u,n=14256))
dd=data.frame(topTags(T4v5_u,n=14256))
ee=data.frame(topTags(T5v6_u,n=14256))
ff=data.frame(topTags(T6v7_u,n=14256))
gg=data.frame(topTags(T7v8_u,n=14256))
hh=data.frame(topTags(T8v9_u,n=14256))

ii=data.frame(topTags(T1v2_i,n=14256))
jj=data.frame(topTags(T2v3_i,n=14256))
kk=data.frame(topTags(T3v4_i,n=14256))
ll=data.frame(topTags(T4v5_i,n=14256))
mm=data.frame(topTags(T5v6_i,n=14256))
nn=data.frame(topTags(T6v7_i,n=14256))
oo=data.frame(topTags(T7v8_i,n=14256))
pp=data.frame(topTags(T8v9_i,n=14256))

qq=data.frame(topTags(T1u_i,n=14256))
rr=data.frame(topTags(T2u_i,n=14256))
ss=data.frame(topTags(T3u_i,n=14256))
tt=data.frame(topTags(T4u_i,n=14256))
uu=data.frame(topTags(T5u_i,n=14256))
vv=data.frame(topTags(T6u_i,n=14256))
ww=data.frame(topTags(T7u_i,n=14256))
xx=data.frame(topTags(T8u_i,n=14256))
yy=data.frame(topTags(T9u_i,n=14256))

aa = tibble::rownames_to_column(aa, "Gene") # convert row names to columns
bb = tibble::rownames_to_column(bb, "Gene") # convert row names to columns
cc = tibble::rownames_to_column(cc, "Gene") # convert row names to columns
dd = tibble::rownames_to_column(dd, "Gene") # convert row names to columns
ee = tibble::rownames_to_column(ee, "Gene") # convert row names to columns
ff = tibble::rownames_to_column(ff, "Gene") # convert row names to columns
gg = tibble::rownames_to_column(gg, "Gene") # convert row names to columns
hh = tibble::rownames_to_column(hh, "Gene") # convert row names to columns

ii = tibble::rownames_to_column(ii, "Gene") # convert row names to columns
jj = tibble::rownames_to_column(jj, "Gene") # convert row names to columns
kk = tibble::rownames_to_column(kk, "Gene") # convert row names to columns
ll = tibble::rownames_to_column(ll, "Gene") # convert row names to columns
mm = tibble::rownames_to_column(mm, "Gene") # convert row names to columns
nn = tibble::rownames_to_column(nn, "Gene") # convert row names to columns
oo = tibble::rownames_to_column(oo, "Gene") # convert row names to columns
pp = tibble::rownames_to_column(pp, "Gene") # convert row names to columns

qq = tibble::rownames_to_column(qq, "Gene") # convert row names to columns
rr = tibble::rownames_to_column(rr, "Gene") # convert row names to columns
ss = tibble::rownames_to_column(ss, "Gene") # convert row names to columns
tt = tibble::rownames_to_column(tt, "Gene") # convert row names to columns
uu = tibble::rownames_to_column(uu, "Gene") # convert row names to columns
vv = tibble::rownames_to_column(vv, "Gene") # convert row names to columns
ww = tibble::rownames_to_column(ww, "Gene") # convert row names to columns
xx = tibble::rownames_to_column(xx, "Gene") # convert row names to columns
yy = tibble::rownames_to_column(yy, "Gene") # convert row names to columns




Uninfect_time_DEGlist <- rbind(subset(bb,bb$FDR<=1e-5),subset(cc,cc$FDR<=1e-5),subset(dd,dd$FDR<=1e-5),subset(ee,ee$FDR<=1e-5),subset(ff,ff$FDR<=1e-5),subset(gg,gg$FDR<=1e-5),subset(hh,hh$FDR<=1e-5))

unique_Uninfect_time_DEGlist <- unique(Uninfect_time_DEGlist$Gene)
length(unique_Uninfect_time_DEGlist) # uninfected number of uniquely differentially expressed transcripts across time comparisons = 2699


Infect_time_DEGlist <- rbind(subset(ii,ii$FDR<=1e-5),subset(jj,jj$FDR<=1e-5),subset(kk,kk$FDR<=1e-5),subset(ll,ll$FDR<=1e-5),subset(mm,mm$FDR<=1e-5),subset(nn,nn$FDR<=1e-5),subset(oo,oo$FDR<=1e-5),subset(pp,pp$FDR<=1e-5))

unique_Infect_time_DEGlist <- unique(Infect_time_DEGlist$Gene)
length(unique_Infect_time_DEGlist) # Infected number of uniquely differentially expressed transcripts across time comparisons = 7030

UvsI_DEGlist = rbind(subset(rr,rr$FDR<=1e-5),subset(ss,ss$FDR<=1e-5),subset(tt,tt$FDR<=1e-5),subset(uu,uu$FDR<=1e-5),subset(vv,vv$FDR<=1e-5),subset(ww,ww$FDR<=1e-5),subset(xx,xx$FDR<=1e-5),subset(yy,yy$FDR<=1e-5))
unique_UvsI_DEGlist <- unique(UvsI_DEGlist$Gene)
length(unique_UvsI_DEGlist) # Infected number of uniquely differentially expressed transcripts across time comparisons = 3238

All_three_Combined_DEGlist = c(unique_UvsI_DEGlist,unique_Infect_time_DEGlist,unique_Uninfect_time_DEGlist)
Unique_All_three = unique(All_three_Combined_DEGlist) # 8263 Genes.

# make a new data frame that has normalized counts data for each of the genes in the Unique_All_three list
UA3 = Unique_All_three
UA3names = paste("gene",UA3,sep="")

Counts_UA3 = data.matrix(all_count[counts$gene %in% UA3names,]) 
Log_Counts_UA3 = log2(Counts_UA3+1)
centered_data = as.data.frame(t(scale(t(Log_Counts_UA3), scale=F)))
centered_data$gene=UA3names
write.table(centered_data, file="UA3_Centered_Log2_counts", quote=FALSE, row.names=FALSE, sep="\t")

# make a new data frame that has normalized counts data for each of the genes in the uninfected list
UInf = unique_Uninfect_time_DEGlist
UInfnames = paste("gene",UInf,sep="")

Counts_UInf = data.matrix(all_count[counts$gene %in% UInfnames,]) 
Log_Counts_UInf = log2(Counts_UInf+1)
UInf_centered_data = as.data.frame(t(scale(t(Log_Counts_UInf), scale=F)))
UInf_centered_data$gene=UInfnames # add gene names
Uninfected_data = UInf_centered_data[,c(82,1:3,10:12,19:21,28:30,37:39,46:48,55:57,64:66,73:75)] # save only columns for uninfected data
write.table(Uninfected_data, file="UInfected_Log2Counts.txt", quote=FALSE, row.names=FALSE, sep="\t")

# make a new data frame that has normalized counts data for each of the genes in the Infected list
Infe = unique_Infect_time_DEGlist
Infnames = paste("gene",Infe,sep="")

Counts_Inf = data.matrix(all_count[counts$gene %in% Infnames,]) 
Log_Counts_Inf = log2(Counts_Inf+1)
Inf_centered_data = as.data.frame(t(scale(t(Log_Counts_Inf), scale=F)))
Inf_centered_data$gene=Infnames
Infected_data = Inf_centered_data[,c(82,4:9,13:18,22:27,31:36,40:45,49:54,58:63,67:72,76:81)] # save only columns for uninfected data
write.table(Infected_data, file="Infected_counts.txt", quote=FALSE, row.names=FALSE, sep="\t")

# make a new data frame that has normalized counts data for each of the genes in the Uninfected vs Infected list
Un_vs_I = unique_UvsI_DEGlist
Un_vs_I_names = paste("gene",Un_vs_I,sep="")

Counts_Un_vs_I = data.matrix(all_count[counts$gene %in% Un_vs_I_names,]) 
Log_Counts_Un_vs_I = log2(Counts_Un_vs_I+1)
Un_vs_I_centered_data = as.data.frame(t(scale(t(Log_Counts_Un_vs_I), scale=F)))
Un_vs_I_centered_data$gene=Un_vs_I_names
#reorder columns
Un_vs_I_centered_data = Un_vs_I_centered_data[,c(82,1)]
write.table(Un_vs_I_centered_data, file="Un_vs_I_centered_data.txt", quote=FALSE, row.names=FALSE, sep="\t")

# Add annotation information
# make a table of transcript annotations that matches the data.
Transcripts_Un_vs_I_centered_data_data = transcripts[transcripts$gene %in% Un_vs_I_centered_data$gene,] # this yielded a set of transcript annotations that exist in the data, but the number of genes in "Infected_data" is 3238 and the number of transcripts is 3159, so cannot merge these records directly.

# try using InnerJoin to do this:
Un_vs_I_centered_data_annotations = inner_join(Un_vs_I_centered_data,Transcripts_Un_vs_I_centered_data_data,by="gene")
# This seems to have worked and produced a dataframe with 3159 entries (those that were annotated our of the 3238 - use these.)
write.table(Un_vs_I_centered_data_annotations, file="Uninf_vs_Inf_annot.txt", quote=FALSE, row.names=FALSE, sep="\t")


# Add annotation information
# make a table of transcript annotations that matches the data.
Transcripts_Infected_data = transcripts[transcripts$gene %in% Infected_data$gene,] # this yielded a set of transcript annotations that exist in the data, but the number of genes in "Infected_data" is 7030 and the number of transcripts is 6875, so cannot merge these records directly.

# try using InnerJoin to do this:
library(dplyr)
Infected_data_annotations = inner_join(Infected_data,Transcripts_Infected_data,by="gene")
# This seems to have worked and produced a dataframe with 6875 entries (those that were annotated - use these.)
write.table(Infected_data_annotations, file="Infected_counts_annot.txt", quote=FALSE, row.names=FALSE, sep="\t")

# for Uninfected data
Transcripts_Uninfected_data = transcripts[transcripts$gene %in% Uninfected_data$gene,] # this yielded a set of transcript annotations that exist in the data, but the number of genes in "Uninfected_data" is 2699 and the number of transcripts is 2643, so cannot merge these records directly.

# try using InnerJoin to do this:
library(dplyr)
Uninfected_data_annotations = inner_join(Uninfected_data,Transcripts_Uninfected_data,by="gene")
# This seems to have worked and produced a dataframe with 2643 entries (those that were annotated - use these.)
write.table(Uninfected_data_annotations, file="Uninfected_counts_annot.txt", quote=FALSE, row.names=FALSE, sep="\t")


### Stopped at this point Oct 9 2022 and opened the output files in Cluster 3.0 for Mac for clustering.  Reopen clustered data files for plotting. 

#####
#
# Clusters of Uninfected Daphnia
#
####


# Filter: MaxVal - MinVal ≥ 1.0 (2416 out of 2643 passed)
# Clustering: Examined clusters with 6, 8, 10 and 12 k-means clusters.  Decided that n=6 clusters was good (visually no two clusters looked the same, although two clusters had the same patterns but different intensities of expression)
# manually enter cluster information in excel according to the visualized data in Treeview (Java)
# Read in clustered datafile 
Uninfected_Clustered_Data=read.csv(file.choose())
# Add annotations to the data
Uninfected_Clustered_Data_annot = inner_join(Uninfected_Clustered_Data,Transcripts_Uninfected_data,by="gene")
write.csv(Uninfected_Clustered_Data_annot[,c(1,3,33,34)],"Uninfected_Clustered_Data_annot.csv")

# prepare data for plotting using reshape
library(reshape2)
clusterdata2=melt(Uninfected_Clustered_Data_annot[,c(3:30)],id.var=c("cluster"))
library(plyr)
library(ggplot2)
head(clusterdata2)
clusterdata2$group = substr(clusterdata2$variable,4,6) # get just the hour info

# make means and variances for each group
new_plotting_data_by_cluster_avgd<-ddply(clusterdata2,.(cluster,group), function(d) mean(d$value)) 

sd_by_cluster_avgd<-ddply(clusterdata2,.(cluster,group), function(d) sd(d$value)) 
library(plotrix)
se_by_cluster_avgd<-ddply(clusterdata2,.(cluster,group), function(d) std.error(d$value)) 
head(new_plotting_data_by_cluster_avgd)
names(new_plotting_data_by_cluster_avgd)=c("Cluster","Group","Mean") # rename
names(sd_by_cluster_avgd)=c("Cluster","Group","SD") # rename
names(se_by_cluster_avgd)=c("Cluster","Group","SE")
plotting_data_by_cluster_avgd=merge(new_plotting_data_by_cluster_avgd,sd_by_cluster_avgd)
plotting_data_by_cluster_avgd=merge(plotting_data_by_cluster_avgd,se_by_cluster_avgd)
plotting_data_by_cluster_avgd$pct=plotting_data_by_cluster_avgd$SE*1.98 # 95% confidence interval
head(plotting_data_by_cluster_avgd)

# make time column by only hours for plotting purposes
plotting_data_by_cluster_avgd$hours = rep(c(3,96,144,6,192,12,384,24,48),6)

#plot with 95% conf intervals
c <- ggplot(plotting_data_by_cluster_avgd, aes(x = hours, y = Mean, ymax = Mean + pct, ymin = Mean - pct))
c + geom_errorbar(color="black", width=0.1) + geom_point(size = 2) + geom_line( size=0.5, linetype=1) + labs(title="Uninfected",x = "Time(hours)", y = "Median centered log2 FPKM") + facet_wrap(~Cluster, ncol=1,scales="free")+theme_bw()+scale_x_continuous(breaks=c(3,6,12,24,48,96,144,192,384))

# Saved this plot as "Uninfected Clusters.pdf" (8X10 size)

# Add data from infected individuals for these genes and then make plots so we can see how these genes vary between infected and uninfected individuals.
# Add annotations to the data
Uninfected_Clustered_Data_add_infected = inner_join(Uninfected_Clustered_Data,Infected_data,by="gene")
# doing this reduced teh gene set to 2326 genes

# remove unneeded columns
Uninfected_Clustered_Data_add_infected$description=NULL
Uninfected_Clustered_Data_add_infected$UniprotID=NULL

clusterdata3=melt(Uninfected_Clustered_Data_add_infected[,c(2:83)],id.var=c("cluster"))
head(clusterdata3)
clusterdata3$group = substr(clusterdata3$variable,4,6) # get just the hour info
clusterdata3$treatment = substr(clusterdata3$variable,8,12)
head(clusterdata3)

new_plotting_data_by_cluster_avgd3<-ddply(clusterdata3,.(cluster,group,treatment), function(d) mean(d$value)) 
names(new_plotting_data_by_cluster_avgd3)=c("Cluster","Group","Treatment","Mean") # rename
head(new_plotting_data_by_cluster_avgd3)

sd_by_cluster_avgd3<-ddply(clusterdata3,.(cluster,group,treatment), function(d) sd(d$value)) 
names(sd_by_cluster_avgd3)=c("Cluster","Group","Treatment","SD") # rename
head(sd_by_cluster_avgd3)

se_by_cluster_avgd3<-ddply(clusterdata3,.(cluster,group,treatment), function(d) std.error(d$value)) 
names(se_by_cluster_avgd3)=c("Cluster","Group","Treatment","SE")
head(se_by_cluster_avgd3)

plotting_data_by_cluster_avgd3=merge(new_plotting_data_by_cluster_avgd3,sd_by_cluster_avgd3)
plotting_data_by_cluster_avgd3=merge(plotting_data_by_cluster_avgd3,se_by_cluster_avgd3)
head(plotting_data_by_cluster_avgd3)
plotting_data_by_cluster_avgd3$pct=plotting_data_by_cluster_avgd3$SE*1.98 # 95% confidence interval
head(plotting_data_by_cluster_avgd3)

# make time column by only hours for plotting purposes
plotting_data_by_cluster_avgd3$hours = rep(c(3,3,96,96,144,144,6,6,192,192,12,12,384,384,24,24,48,48),6)
head(plotting_data_by_cluster_avgd3)

#plot with 95% conf intervals
#plot with 95% conf intervals
d <- ggplot(plotting_data_by_cluster_avgd3, aes(x = hours, y = Mean, ymax = Mean + pct, ymin = Mean - pct, color=Treatment, fill=Treatment, group=Treatment))
d + geom_errorbar(color="black", width=0.1) + geom_point(aes(shape=Treatment, fill=Treatment),size = 2) + geom_line(aes(group=Treatment),size=0.5, linetype=1) + labs(title="Uninfected_identified_genes + data for those genes from Infected Individuals",x = "Time (hours)", y = "Median centered log2 FPKM") + facet_wrap(~Cluster, ncol=1,scales="free")+theme_bw()+scale_x_continuous(breaks=c(3,6,12,24,48,96,144,192,384))

# Analysis of graph: Gene expression patterns of Uninfected (=contr) and Infected (=treat) Daphnia are only similar in cluster 4.  In general, the large changes in gene expression observed over time in the Uninfected Daphnia were not observed for the same genes in Infected Daphnia.

#####
#
# Clusters of Infected Daphnia
#
####

# Filter: MaxVal - MinVal ≥ 1.0 (6372 out of 6875 passed)
# Clustering: Examined clusters with 6, 8, 10 and 12 k-means clusters.  Decided that n=6 clusters was good (visually no two clusters looked the same)
# manually enter cluster information in excel according to the visualized data in Treeview (Java)

Infected_Clustered_Data=read.csv(file.choose()) 
# Add annotations to the data
Infected_Clustered_Data_annot = inner_join(Infected_Clustered_Data,Transcripts_Infected_data,by="gene")
write.csv(Infected_Clustered_Data_annot[,c(1,2,59,60)],"Infected_Clustered_Data_annot.csv")

# prepare data for plotting using reshape
clusterdata4=melt(Infected_Clustered_Data_annot[,c(2:57)],id.var=c("cluster"))
head(clusterdata4)
clusterdata4$group = substr(clusterdata4$variable,4,6) # get just the hour info
# make means and variances for each group
new_plotting_data_by_cluster_avgd4<-ddply(clusterdata4,.(cluster,group), function(d) mean(d$value)) 

sd_by_cluster_avgd4<-ddply(clusterdata4,.(cluster,group), function(d) sd(d$value)) 
se_by_cluster_avgd4<-ddply(clusterdata4,.(cluster,group), function(d) std.error(d$value)) 
head(new_plotting_data_by_cluster_avgd4)
names(new_plotting_data_by_cluster_avgd4)=c("Cluster","Group","Mean") # rename
names(sd_by_cluster_avgd4)=c("Cluster","Group","SD") # rename
names(se_by_cluster_avgd4)=c("Cluster","Group","SE")
plotting_data_by_cluster_avgd4=merge(new_plotting_data_by_cluster_avgd4,sd_by_cluster_avgd4)
plotting_data_by_cluster_avgd4=merge(plotting_data_by_cluster_avgd4,se_by_cluster_avgd4)
plotting_data_by_cluster_avgd4$pct=plotting_data_by_cluster_avgd4$SE*1.98 # 95% confidence interval
head(plotting_data_by_cluster_avgd4,12)
# for some reason there is a group called "cri" in each 10th place - I have no idea where this came from.  Delete it.
plotting_data_by_cluster_avgd4=plotting_data_by_cluster_avgd4[-c(10,20,30,40,50,60),]

# make time column by only hours for plotting purposes
plotting_data_by_cluster_avgd4$hours = rep(c(3,96,144,6,192,12,384,24,48),6)

#plot with 95% conf intervals
c <- ggplot(plotting_data_by_cluster_avgd4, aes(x = hours, y = Mean, ymax = Mean + pct, ymin = Mean - pct))
c + geom_errorbar(color="black", width=0.1) + geom_point(size = 2) + geom_line( size=0.5, linetype=1) + labs(title="Infected",x = "Time (hours)", y = "Median centered log2 FPKM") + facet_wrap(~Cluster, ncol=1,scales="free")+theme_bw()+scale_x_continuous(breaks=c(3,6,12,24,48,96,144,192,384))

# Saved this plot as "Infected.pdf" (8X10 size)

# Add data from infected individuals for these genes and then make plots so we can see how these genes vary between infected and uninfected individuals.
# Add annotations to the data
Infected_Clustered_Data_add_uninfected = inner_join(Infected_Clustered_Data,Uninfected_data,by="gene")
# doing this reduced the gene set to 2390 genes

# remove unneeded columns
Infected_Clustered_Data_add_uninfected$description=NULL
Infected_Clustered_Data_add_uninfected$UniprotID=NULL

clusterdata5=melt(Infected_Clustered_Data_add_uninfected[,c(2:83)],id.var=c("cluster"))
head(clusterdata5)
clusterdata5$group = substr(clusterdata5$variable,4,6) # get just the hour info
clusterdata5$treatment = substr(clusterdata5$variable,8,12)
head(clusterdata5)

new_plotting_data_by_cluster_avgd5<-ddply(clusterdata5,.(cluster,group,treatment), function(d) mean(d$value)) 
names(new_plotting_data_by_cluster_avgd5)=c("Cluster","Group","Treatment","Mean") # rename
head(new_plotting_data_by_cluster_avgd5)

sd_by_cluster_avgd5<-ddply(clusterdata5,.(cluster,group,treatment), function(d) sd(d$value)) 
names(sd_by_cluster_avgd5)=c("Cluster","Group","Treatment","SD") # rename
head(sd_by_cluster_avgd5)

se_by_cluster_avgd5<-ddply(clusterdata5,.(cluster,group,treatment), function(d) std.error(d$value)) 
names(se_by_cluster_avgd5)=c("Cluster","Group","Treatment","SE")
head(se_by_cluster_avgd5)

plotting_data_by_cluster_avgd5=merge(new_plotting_data_by_cluster_avgd5,sd_by_cluster_avgd5)
plotting_data_by_cluster_avgd5=merge(plotting_data_by_cluster_avgd5,se_by_cluster_avgd5)
head(plotting_data_by_cluster_avgd5)
plotting_data_by_cluster_avgd5$pct=plotting_data_by_cluster_avgd5$SE*1.98 # 95% confidence interval
head(plotting_data_by_cluster_avgd5)

# make time column by only hours for plotting purposes
plotting_data_by_cluster_avgd5$hours = rep(c(3,3,96,96,144,144,6,6,192,192,12,12,384,384,24,24,48,48),6)
head(plotting_data_by_cluster_avgd5)

#plot with 95% conf intervals
e <- ggplot(plotting_data_by_cluster_avgd5, aes(x = hours, y = Mean, ymax = Mean + pct, ymin = Mean - pct, color=Treatment, fill=Treatment, group=Treatment))
e + geom_errorbar(color="black", width=0.1) + geom_point(aes(shape=Treatment, fill=Treatment),size = 2) + geom_line(aes(group=Treatment),size=0.5, linetype=1) + labs(title="Infected_identified_genes + data for those genes from Uninfected Individuals",x = "Time (hours)", y = "Median centered log2 FPKM") + facet_wrap(~Cluster, ncol=1,scales="free")+theme_bw()+scale_x_continuous(breaks=c(3,6,12,24,48,96,144,192,384))

# Analysis of graph: Gene expression patterns of Uninfected (=contr) and Infected (=treat) Daphnia are only similar in cluster 1 and maybe 2.  In general, the large changes in gene expression observed over time in the Infected Daphnia were not observed for the same genes in Uninfected Daphnia.


####
#
#  Uninfected vs. Infected Genes
#
####

# Filter: MaxVal - MinVal ≥ 1.0 (2975 out of 3159 passed)
# Clustering: Examined clusters with 6, 8, 10 and 12 k-means clusters.  Decided that n=8 clusters was good (visually no two clusters looked the same)
# manually enter cluster information in excel according to the visualized data in Treeview (Java)

UvI_Clustered_Data=read.csv(file.choose()) 
# Add annotations to the data
UvI_Clustered_Data_annot = inner_join(UvI_Clustered_Data,Transcripts_Un_vs_I_centered_data_data,by="gene")
write.csv(UvI_Clustered_Data_annot[,c(1,2,86,87)],"UvI_Clustered_Data_annot.csv")

# prepare data for plotting using reshape
clusterdata6=melt(UvI_Clustered_Data_annot[,c(2:83)],id.var=c("cluster"))
head(clusterdata6)
clusterdata6$group = substr(clusterdata6$variable,4,6) # get just the hour info
clusterdata6$treatment = substr(clusterdata6$variable,8,12)
# make means and variances for each group
new_plotting_data_by_cluster_avgd6<-ddply(clusterdata6,.(cluster,group,treatment), function(d) mean(d$value)) 
 
names(new_plotting_data_by_cluster_avgd6)=c("Cluster","Group","Treatment","Mean") # rename
head(new_plotting_data_by_cluster_avgd6)

sd_by_cluster_avgd6<-ddply(clusterdata6,.(cluster,group,treatment), function(d) sd(d$value)) 
names(sd_by_cluster_avgd6)=c("Cluster","Group","Treatment","SD") # rename
head(sd_by_cluster_avgd6)

se_by_cluster_avgd6<-ddply(clusterdata6,.(cluster,group,treatment), function(d) std.error(d$value)) 
names(se_by_cluster_avgd6)=c("Cluster","Group","Treatment","SE")
head(se_by_cluster_avgd6)

plotting_data_by_cluster_avgd6=merge(new_plotting_data_by_cluster_avgd6,sd_by_cluster_avgd6)
plotting_data_by_cluster_avgd6=merge(plotting_data_by_cluster_avgd6,se_by_cluster_avgd6)
head(plotting_data_by_cluster_avgd6)
plotting_data_by_cluster_avgd6$pct=plotting_data_by_cluster_avgd6$SE*1.98 # 95% confidence interval
head(plotting_data_by_cluster_avgd6)

# make time column by only hours for plotting purposes
plotting_data_by_cluster_avgd6$hours = rep(c(3,3,96,96,144,144,6,6,192,192,12,12,384,384,24,24,48,48),8)
head(plotting_data_by_cluster_avgd6)

#plot with 95% conf intervals
f <- ggplot(plotting_data_by_cluster_avgd6, aes(x = hours, y = Mean, ymax = Mean + pct, ymin = Mean - pct, color=Treatment, fill=Treatment, group=Treatment))
f + geom_errorbar(color="black", width=0.1) + geom_point(aes(shape=Treatment, fill=Treatment),size = 2) + geom_line(aes(group=Treatment),size=0.5, linetype=1) + labs(title="Genes that differed between Infected and Uninfected Daphnia in at least one timepoint",x = "Time (hours)", y = "Median centered log2 FPKM") + facet_wrap(~Cluster, ncol=1,scales="free")+theme_bw()+scale_x_continuous(breaks=c(3,6,12,24,48,96,144,192,384))

# Analysis of graph: Gene expression patterns of Uninfected (=contr) and Infected (=treat) Daphnia very similar across timepoints in general. Interpretation is that at most timepoints there were not big differences between Infected and Uninfected Individuals, though at some timepoints gene expression varied between the two groups, especially at the longer timepoints.

######
#
# Pick genes of interest for molt cycle based on the Mykles Frontiers review paper: 
# https://www.frontiersin.org/articles/10.3389/fendo.2021.674711/full
#
######


MoltGeneList = c('ecdys','mTOR','mtor','rapamycin','activin','myostat','TGF','Mstn','CHH','RXR','MIH','MAP','AMP','ErbB','Hedgehog','HIF','Jak','Hippo','NF-Kappa','Notch','TNF','Wnt','secretin','gpcr','metabotropic glutamate','gonadotropin-releasing','tachykinin','relaxin','rhodopsin G0','methusela','dopamine D2','opsin','serotonin','neuropeptide f','bursicon','itpr-like','moody','neuroparsin','inotocin','vasopressin','eclosion','itp','gonad inhibiting','mandibular inhibiting','corazonin','ptth','prothororacic','dilp8','disembodied','phantom','spookier','phospholipase C','PLC','akt','insulin-like peptide','ilp','insulin receptor','insR','egfr','Ras KEGG','Pi3K','VEGF','EGFR','Sp-vtg','sp-cyclin','tgfb','tgfr','camkii','nodal','smad','babo','torso','daw','dawdle','mstn','bambi','tsc','rheb','phm','dib','gi-ef2','gi elongation factor','bombyxin','FXPRLamide','Bommo-myosuppressin','FMRFamide','CCHamide','crustacean cardioactive peptide','ccap','diuretic hormones DH31')

MoltGenes=transcripts[grepl(paste(MoltGeneList, collapse = '|'), transcripts$description, ignore.case=T),] # syntax from https://stackoverflow.com/questions/26347947/subset-a-dataframe-based-on-a-list-of-words

MoltGenesData = inner_join(MoltGenes,counts,by="gene")
Log_MoltGenesData = log2(MoltGenesData[,4:84]+1)
MoltGenesData_centered_data = as.data.frame(t(scale(t(Log_MoltGenesData), scale=F)))
MoltGenesData_centered_data_annot = cbind(MoltGenes,MoltGenesData_centered_data)
write.table(MoltGenesData_centered_data_annot,file="MoltGenesData.txt", quote=FALSE, row.names=FALSE, sep="\t")

# clustered into 9 clusters in Cluster3.0.  There are 3 similar clusters that have wnt genes, but kept them separate as expression patterns not identical.  

library(plyr)
library(plotrix)

# Read in clustered file:
Molt_Clustered_data=read.csv(file.choose()) 
colnames(Molt_Clustered_data)[3]=c("cluster")

# prepare data for plotting using reshape
clusterdata7=melt(Molt_Clustered_data[,c(3:84)],id.var=c("cluster"))
head(clusterdata7)
clusterdata7$group = substr(clusterdata7$variable,4,6) # get just the hour info
clusterdata7$treatment = substr(clusterdata7$variable,8,12)
# make means and variances for each group
new_plotting_data_by_cluster_avgd7<-ddply(clusterdata7,.(cluster,group,treatment), function(d) mean(d$value)) 

names(new_plotting_data_by_cluster_avgd7)=c("Cluster","Group","Treatment","Mean") # rename
head(new_plotting_data_by_cluster_avgd7)

sd_by_cluster_avgd7<-ddply(clusterdata7,.(cluster,group,treatment), function(d) sd(d$value)) 
names(sd_by_cluster_avgd7)=c("Cluster","Group","Treatment","SD") # rename
head(sd_by_cluster_avgd7)

se_by_cluster_avgd7<-ddply(clusterdata7,.(cluster,group,treatment), function(d) std.error(d$value)) 
names(se_by_cluster_avgd7)=c("Cluster","Group","Treatment","SE")
head(se_by_cluster_avgd7)

plotting_data_by_cluster_avgd7=merge(new_plotting_data_by_cluster_avgd7,sd_by_cluster_avgd6)
plotting_data_by_cluster_avgd7=merge(plotting_data_by_cluster_avgd7,se_by_cluster_avgd6)
head(plotting_data_by_cluster_avgd7)
plotting_data_by_cluster_avgd7$pct=plotting_data_by_cluster_avgd7$SE*1.98 # 95% confidence interval
head(plotting_data_by_cluster_avgd7)

# make time column by only hours for plotting purposes
plotting_data_by_cluster_avgd7$hours = rep(c(3,3,96,96,144,144,6,6,192,192,12,12,384,384,24,24,48,48),8)
head(plotting_data_by_cluster_avgd7)

#plot with 95% conf intervals
f <- ggplot(plotting_data_by_cluster_avgd7, aes(x = hours, y = Mean, ymax = Mean + pct, ymin = Mean - pct, color=Treatment, fill=Treatment, group=Treatment))
f + geom_errorbar(color="black", width=0.1) + geom_point(aes(shape=Treatment, fill=Treatment),size = 2) + geom_line(aes(group=Treatment),size=0.5, linetype=1) + labs(title="Molt Related Genes",x = "Time (hours)", y = "Median centered log2 Counts") + facet_wrap(~Cluster, ncol=1,scales="free")+theme_bw()+scale_x_continuous(breaks=c(3,6,12,24,48,96,144,192,384))

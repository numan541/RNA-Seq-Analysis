
#---------------------------------#
#RNASeq data analysis in R   #
#---------------------------------#
# means comments
# After installing R from, get the full latest version of bioconductor using following commands or Load the bioconductor (www.bioconductor.org) packages for working with RNASeq data:
#Bioconductor provides tools for the analysis and comprehension of high-throughput genomic data. 
#Bioconductor uses the R statistical programming language, and is open source and open development. 
#It has two releases each year, 1024 software packages, and an active user community.
# Go to R console and type following command
library(edgeR) # will upload edgeR, if not you will use BiocManager() as under;
install.packages("BiocManager")
BiocManager::install("edgeR") # will upload edgeR directly from server, once downloaded, it will 
library(edgeR) # load in youer R Session 
# BiocManager::install# will get all the updated packages in it takes longer, we can do them one by one 
#### Same way load others also as;
library(affy) # if not loading, use biocLite as  
BiocManager::install("affy") 
# if not loading, use biocLite as above
library(affycoretools) # if not loading, use biocLite as above
BiocManager::install("affycoretools") 
library(limma) # plotDensities
# during installation if it asks for go with yes or all or yes to all
dir()# will tell what files and directories/folders you have in your current folder
library() # will tell you what libraries you have already installed
ls() # will tell what variables you have so far, variable is anything that can store value/s, type of a variable will be specified by the value it stores.
#### Check the current working directory, which is where R reads from and writes to unlessfull path names are given. Change the working directory to the one where the data is.
getwd()  ### will tell you where exactly you are right now
setwd("F:/MS/Semester2/DGE_assignment/Modified/") # will change to where we want
# Need help? want to know about what any commond does, just type 
help() # will open help server for you, you can specify your question also as help(anything) e.g 
help("ls") 
help("objects")
help("class") # or
?class# will help if its parent package is already loaded
??class# will search for what package it exists in R
x# will give an error since there is no variable defined as x
x<-1# will create a variable x (a numeric vector) and stores a value of 1 in it. <- is assignment operator, similar to =
x# now will return var x and it's value
y<-1:10# will store 1 to 10 in y
y
z<-1,3,4# will give an error, the right way is to use c(), concatenate function
z<-c(1,2,3,3)
z
# both x and y are also numeric vectors
v<-c ('a', 'b', 'c') # we can store characters in v. V will become a character vector. Let's check it's class
class(v)
rm() # you can remove any variable
rm(list=ls()) # you can remove all variables
length(z) # will return you how many items are there in a vector
unique(z) # will return a non redundant list
w<-matrix(1:10, nrow=2, ncol=5, byrow=T)# will create a matrix type object & fill it row-wise try byrow=F also. Matrix stores numeric vectors
w
dim(w) # will tell how many rows and columns we have in w
data.frame(1,10,c('a', 'b')) # will create a dataframe that can have different types of vectors like integer, character etc
factor() # will create a factor class object that have different levels
z<-c(1,2,2,3,3,5)
factor(z)
as.factor() # 'as' will change the data type to as what specified
as.numeric() # will change the data type to numeric
is.factor() # will check whether the data type is factor
is.numeric # will check whether the data type is numeric 
apply() # will apply any function on any data
rep(1,3) # will print 1 three times
rep(c("A","B"), 3) # will print A,B three times alternatively
rep(c("A","B"),each=3) # will print A, 3 times and then B, 3 times 
ls() #### Show the objects in workspace; if there are any, remove them.
rm(list=ls())
###################### RNASeq Data Analysis ###########################
#The data is a subset from an in-house RNASeq project, comparing wild type with a mutant,
#using 3 independent replicates of each. The data file has gene ids, their transcript lengths and counts of reads mapping to these 100 genes 
#for each of the six samples obtained after aligning the reads with A.thaliana transcripts 
#obtained from TAIR10 using bowtie aligner and reads mapped to each gene/transcript are counted 
#using in-house python script getCountMatrix.py. Lets input this file as csv##
raw.data<- read.csv("AT_data_modified.csv",header=T, row.names=1)
head(raw.data)
# will read a .csv file, header=T will take first line as column names
# We can read also a text file (.txt) as under
#raw.data <- read.delim("RNASeq_data.txt",stringsAsFactors=F) 
# default T, takes numeric values as factors. Once the pro mpt returns, verify that you have created a new object in the workspace, and see its kind.
ls()
class(raw.data)# Because it's a data frame, we can ask what the dimensions are:
dim(raw.data)  # how many rows and columns we have
head(raw.data) # Look at the first few (6) rows of the data
raw.data[1:5,] # Look at the first 5 rows of the data
raw.data[1:5,3:14] # Look at the first 5 rows and 2nd to 7th coulmns of the data
tail(raw.data) # Look at the last few lines of the data
raw.data$X35STFL2   # will print this column or create a new one if it does not exist
raw.data$TFL1
#### We need to tell R the overall experimental design, i.e what treatments are there in which column/s, how many total treatments are there and their replicates etc. This experimental design is simple, we tell R as:
expt<-rep(c("X35STFL","X152","TFL", "WT"),each=3)
expt
expt<- factor(expt)  # ensure the order
expt
expt<-factor(expt, levels=c("X35STFL","X152","TFL", "WT"))#to set the levels and their order manually 
expt
#### Lets take a break, save the workspace and 	history, then practicequitting R and restarting your analysis
save.image("RNASeqDemo.RData")
savehistory("RNASeqDemo.Rhistory")
q()  #say no; if say yes, will save files with default names of ".RData" and ".Rhistory"
#### Start R again, change working directory, load .RData and .Rhistory
setwd("C:/Myworkshop")
load("RNASeqDemo.RData")
loadhistory("RNASeqDemo.Rhistory")
ls()  #see that objects have been loaded
#### Use up and down arrows to see through history (.Rhistory)
#### Must load libraries again, every time you start R!
BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("affy")
BiocManager::install("affycoretools")
library(edgeR)
library(limma)
library(affy)
library(affycoretools)

##########################Data Exploration############################
#### Because RNA-Seq count data is so very different from continuous microarray data, it's a good idea to do some basic exploration of the data. Let's start by looking at the range of counts per sample:
summary(raw.data)       # function from affy package
summary(raw.data[,3:14])   
boxplot(raw.data[,3:14]) # from graphics package, boxplots here are weird because of much variation in the data. We might transform data
boxplot(log2(raw.data[,3:14])) # Give warnings as log0 is -Infboxplot(log2(raw.data[,2:7]+0.01)) # small constant added to avoid 0s, lets make it colorful
boxplot(log2(raw.data[,3:14]+0.01), col=rep(c("green","Red", "Yellow","pink"), each=3), xlab= "Samples", ylab="Raw Expression in log2 scale")
## Boxplot does not show means in R. We might ask R to show the means 
points(log2(colMeans(raw.data[,3:14]+0.01)), col="blue", pch=18)
#### All 6 samples have some 0 counts, which is typical for RNA-Seq data. We can also examine the overall distributions of counts to see if any of the samples are different. 
#Because of the extreme range in the data, converting to log2 scale helps in plotting; however, you can't take the log of 0, so we need to add a small constant before taking the logs and plotting the distributions. Since the smallest count value is 1, pick a constant smaller than that, like 0.01
##All graphics can be saved using the menu File -> Save as
plotDensities(log2(raw.data[,3:14]), legend="topright") # Warnings
plotDensities(log2(raw.data[,3:14]+0.01), legend="topright") # add constant if don't like to be warned    
#### The shapes are similar, except mu3 is slightly different. It has many more low values (likely zeros), so maybe it has a smaller number of total counts. Let's check the library size for each sample by summing all the counts:
library.sizes <- colSums(raw.data[,3:14])
library.sizes
#x11()   #This will open a new graphing window; without it, previous graph will be replaced. But it slows the system.
barplot(library.sizes, col=rep(c("green","Red", "Yellow","pink"), each=3), xlab= "Samples", ylab="Total RNA/sample", width=1000000, ylim=c(0, 12000000))
plot(library.sizes, col=rep(c("green","Red", "Yellow","pink"), each=3), xlab= "Samples", ylab="Total RNA/sample")
heatmap(as.matrix(raw.data[1:10,3:14])) # Heat map of counts of first 10 genes
heatmap(as.matrix(raw.data[1:10,3:14]),col = cm.colors(256)) # better colors
heatmap(as.matrix(raw.data[1:10,3:14]), Rowv=NA, Colv=NA) # Heat map without clutering
heatmap(as.matrix(raw.data[1:10,3:14]), Colv=NA) # Heat map without sample clutering
###########################RPKM Normalization#######################
##We can calculate RPKM from the gene lengths and the library.sizes. We need to do this separately for each sample; there are many ways this can be done, but here is one. First, set up a matrix with the appropriate number of rows and columns to hold the RPKM data:
rpkm.data <- matrix(NA,nrow=nrow(raw.data),ncol=12,dimnames=list(row.names(raw.data),colnames(raw.data[,3:14])))
# look at what you created:
rpkm.data[1:10,]
# Now, compute the RPKM values and put in the rpkm.data object:
for (i in 1:12) {
rpkm.data[,i] <- raw.data[,2+i] / (raw.data$length/1000) / (library.sizes[i]/1000000)
    }
rpkm.data[1:10,]
write.csv(rpkm.data, "RPKMS.csv")
# We can also compute the RPM values and put in the rpkm.data object:
rpm.data <- matrix(NA,nrow=nrow(raw.data),ncol=12,dimnames=list(row.names(raw.data),colnames(raw.data[,3:14])))
# look at what you created:
rpm.data[1:10,]
for (i in 1:12) {
    rpm.data[,i] <- raw.data[,1+i] / (library.sizes[i]/1000000)
    }
rpm.data[1:5,]
raw.data[1:5,]
#### RPKM values appear to be continuous, but they are actually based on count data. Examine the range of values per sample:
summary(rpkm.data)
#you can also save the graph as graphics 
pdf("boxplotsRPKM.pdf")
plot_rpkm<-boxplot(log2(rpkm.data+0.01), col=rep(c("green","Red", "Yellow","pink"), each=3), xlab= "Samples", ylab="RPKM Expression in log2 scale")
plot_rpkm
## Boxplot does not show means in R. We might ask R to show the means 
plot_rpkm<-points(log2(colMeans(rpkm.data+0.01)),col="blue",pch=18)
plot_rpkm
dev.off()#
barplot(colSums(rpkm.data), col=rep(c("green","Red", "Yellow","pink"), each=3), xlab= "Samples", ylab="Total RNA in RPKM/sample")
heatmap(rpkm.data[1:50,])
########################Clustering and Separation######################
#### Do a quick-and-dirty cluster analysis to see how similar the samples are. 
#We will use a fast hierarchical clustering from the WGCNA package, which is useful for seeing outliers and if there are any major groupings of samples. 
#We will use some of WGCNA functions must load the WGCNA package. We will do the clustering on the RPKM values, but again they need to have a small constant added and the log2 taken. 
#Then, since we want to cluster the samples, the data matrix must be transposed so that the rows X column are samples X Genes. 
#Finally, we calculate a distance metric between the samples and perform the 
#hierarchical cluster. All this can be done in:
hc.rpkm <- hclust(dist(t(log2(rpkm.data+0.01))),method="average")
# Now plot it, and change the size of the graphing window:
#sizeGrWindow(10,6)
plot(hc.rpkm,hang = -1, main = "Hierarchical Clustering", sub = "", xlab = "",cex=0.9)
#### PCA plot. This will do a principle components analysis, 
#which compresses the the information into just a few 'principle components', and plots the deg. 
#The plot can either be a scree plot (which shows how much variation each principle component explains) or a plot of 2 PCs.
#First, check the screeplot to see how many of the PCs explain much variation. By definition, PC1 > PC2 > PC3, etc.
plotPCA(log2(rpkm.data+0.01), screeplot=T) # Affycoretools
plotPCA(log2(rpkm.data+0.01))
# Adding names on the sample
PCA_names<-plotPCA(log2(rpkm.data+0.01),addtext=colnames(raw.data[,3:14]), legend=F)
PCA_names
# You can see that even using RPKM values, the mu3 sample is very different. To output RPKM data to use in other programs:
write.csv(rpkm.data,file="rpkm_values.csv")
save.image("RNASeq2.RData")  #don't forget to save often!!
savehistory("RNASeqDemo.Rhistory")
gc()    #this helps to reclaim used memory
#########Differential Gene Analyses using package edgeR ###############
####edgeR requires 3 pieces of data:
# "1. counts: a matrix of counts where each row represents a gene/exon (or whatever genomic feature is being tracked) and each column is a different sample. The row names are transcript IDs. We have these in raw.data, but are the row names the gene IDs? Check the 1st few:
head(raw.data)
# 2. group: a factor with length equal to the number of columns denoting the experimental group. We have this in the expt object:
expt
# 3. lib.size: vector of the same length as group giving the total numberof reads aligned/mapped within each sample. Its library.sizes:
library.sizes
#### We need to put this information into an object of a specific class,DGEList. There is a special function to create this object; access the help file to see how to use it:
?DGEList
# put all 3 pieces into a DGEList object:
d<- DGEList(counts=as.matrix(raw.data[,3:14]), lib.size=library.sizes, group=expt)
# What kind of object is it?
class(d)
# Since it's a list, we can find out the names of the items in it
names(d)
# The counts are stored in the $counts:
d$counts[1:5,]
# The group info and library sizes are in:
d$samples$group
# Notice that the norm.factors column are all 1s. In addition to library size,
# read counts can be affected by a few genes with very high expression. For example, here are the proportions of total counts for the top 10 genes in the mu1 sample:
sort(d$counts[,1]/library.sizes[1],decreasing=T)[1:10]
# The top 2 genes each make up ~15% of the total reads! You can imagine that any changes in these high abundance genes could affect the total pool of RNAs and hence the counts of other genes. 
#edgeR suggest an additional normalization factor using a TMM method; see ?calcNormFactors for details:
d<- calcNormFactors(d)
d$samples
#### edgeR has a function to cluster the samples based on multidimensional scaling,
#### which is sort of like PCA, but has a distance measure appropriate for count data:
plotMDS.DGEList(d,main="MDS plot",col=rep(1:4,each=3))
# In contrast to the PCA on RPKM values, using tools appropriate for count data makes it so the mu3 sample is no longer an outlier!
save.image("RNASeqDemo.RData")  #don't forget to save often!!
savehistory("RNASeqDemo.Rhistory")
########## Now do DE testing using edgeR #############.
#Accounting for the variance among replicates is termed "estimating dispersions" in edgeR. You can read about the different options in the edgeR vignette, but as this experiment has a single factor, we will use the qCML method, common dispersion estimate.
d<- estimateCommonDisp(d)
#### Expression differences can be tested using an exact test with a negativebinomial distribution:
nbt1<- exactTest(d, pair = c("X35STFL","WT"))
nbt2<- exactTest(d, pair = c("X152","WT"))
nbt3<- exactTest(d, pair = c("TFL","WT"))
# See what is in the de.com object:
names(nbt1)
names(nbt2)
names(nbt3)
#the deg we want are in the $table
nbt1$table[1:6,]
nbt2$table[1:6,]
nbt3$table[1:6,]
# The p.values have not been adjusted for multiple hypothesis testing; This can be done with the topTags function, similar to the topTable limma function. It will perform FDR (False Discovery Rate correction)
nbt.corrected1<- topTags(nbt1,n=Inf)
nbt.corrected1[1:5,] #see deg for top 5 genes
nbt.corrected2<- topTags(nbt2,n=Inf)
nbt.corrected2[1:5,] #see deg for top 5 genes
nbt.corrected3<- topTags(nbt3,n=Inf)
nbt.corrected3[1:5,] #see deg for top 5 genes
#### This result can be output by
write.csv(nbt.corrected1$table, file="NBT1.csv")
write.csv(nbt.corrected2$table, file="NBT2.csv")
write.csv(nbt.corrected3$table, file="NBT3.csv")
# How many genes are significant at FDR p < 0.05?
sum(nbt.corrected1$table$FDR<0.05)
sum(nbt.corrected2$table$FDR<0.05)
sum(nbt.corrected3$table$FDR<0.05)
nbttable1<- nbt.corrected1$table
nbttable2<- nbt.corrected2$table
nbttable3<- nbt.corrected3$table
deg1<- nbttable1[nbttable1$FDR<0.05,]
head(deg1)
dim(deg1)
deg2<- nbttable2[nbttable2$FDR<0.05,]
head(deg2)
dim(deg2)
deg3<- nbttable3[nbttable3$FDR<0.05,]
head(deg3)
dim(deg3)
deg1$ID<-rownames(deg1)
head(deg1)
deg2$ID<-rownames(deg2)
head(deg2)
deg3$ID<-rownames(deg3)
head(deg3)
# separate up-regulated and down-regulated genes
deg_up1<- deg1[(deg1$FDR<0.05&deg1$logFC>0),]
deg_up1<- deg_up1[order(deg_up1$logFC, decreasing=T),]
dim(deg_up1)
head(deg_up1)
deg_up2<- deg2[(deg2$FDR<0.05&deg2$logFC>0),]
deg_up2<- deg_up2[order(deg_up2$logFC, decreasing=T),]
dim(deg_up2)
head(deg_up2)
deg_up3<- deg3[(deg3$FDR<0.05&deg3$logFC>0),]
deg_up3<- deg_up3[order(deg_up3$logFC, decreasing=T),]
dim(deg_up3)
head(deg_up3)
#Down_regulated
deg_down1<- deg1[(deg1$FDR<0.05&deg1$logFC<0),]
deg_down1<- deg_down1[order(deg_down1$logFC),]
dim(deg_down1)
head(deg_down1)
deg_down2<- deg2[(deg2$FDR<0.05&deg2$logFC<0),]
deg_down2<- deg_down2[order(deg_down2$logFC),]
dim(deg_down2)
head(deg_down2)
deg_down3<- deg3[(deg3$FDR<0.05&deg3$logFC<0),]
deg_down3<- deg_down3[order(deg_down3$logFC),]
dim(deg_down3)
head(deg_down3)
#Output your deg
write.csv(deg_up1,file="UPRegulated1.csv")
write.csv(deg_up2,file="UPRegulated2.csv")
write.csv(deg_up3,file="UPRegulated3.csv")
write.csv(deg_down1,file="DownRegulated1.csv")
write.csv(deg_down2,file="DownRegulated2.csv")
write.csv(deg_down3,file="DownRegulated3.csv")
# combine deg with raw counts
raw.data$ID<-rownames(raw.data)
head(raw.data)
combined1<-merge(raw.data, deg1, by='ID', all.ID=T)
head(combined1)
combined2<-merge(raw.data, deg2, by='ID', all.ID=T)
head(combined2)
combined3<-merge(raw.data, deg3, by='ID', all.ID=T)
head(combined3)
deg_up1<- deg_up1[order(deg_up1$logFC, decreasing=T),]
deg_up2<- deg_up2[order(deg_up2$logFC, decreasing=T),]
deg_up3<- deg_up3[order(deg_up3$logFC, decreasing=T),]
combined1<-combined1[order(combined1$logFC, decreasing=T),]
head(combined1)
combined2<-combined2[order(combined2$logFC, decreasing=T),]
head(combined2)
combined3<-combined3[order(combined3$logFC, decreasing=T),]
head(combined3)
heatmap(as.matrix(combined1[1:10,4:15]), Colv=NA, Rowv=NA)
heatmap(as.matrix(combined2[1:10,4:15]), Colv=NA, Rowv=NA)
heatmap(as.matrix(combined3[1:10,4:15]), Colv=NA, Rowv=NA)

write.csv(combined1, "DEG1.csv")
write.csv(combined2, "DEG2.csv")
write.csv(combined3, "DEG3.csv")


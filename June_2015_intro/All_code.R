# Gene Expression analysis in R

## Section 0 - packages and libraries
# later on we might need some extra tools and systems
# so we'll download them now - select the lines below, and use Edit->Execute to run them
# remember, everything prefaced by the symbol '#' is a comment, and is purely for making notes in the code

# if any of them ask whether you want to update other packages, press 'a' for all
# and if any ask you if they should be put into other directories, press 'y' for yes

source("http://bioconductor.org/biocLite.R")

biocLite("Biobase")
biocLite("marray")
biocLite("limma")
biocLite("survival")
biocLite("oligo")
biocLite("affy")
biocLite("annotate")
biocLite("gplots")

library(Biobase)
library(marray)
library(limma)
library(survival)
library(oligo)
library(affy)
library(annotate)
library(gplots)

########################################################
## Section 1 - we have to do a little bit of R basics, and then do some plotting
# we will access the data from http://www.nejm.org/doi/full/10.1056/NEJMoa021967
# it is a microarray gene expression dataset of 295 breast cancers
# first we will read the data into R from a file
# we will use an appropriate function to read this data - readExpressionSet
NKI295 <- readExpressionSet(exprsFile="data/NKI295.exprs.txt", sep="\t", header=T, stringsAsFactors=F) 

# now, the object 'NKI295' is a data structure holding the expression data

# we'll include some extra data 
fData(NKI295) <- read.delim("data/NKI295.fdata.txt", sep="\t", header=T, row.names=1, stringsAsFactors=F) 
pData(NKI295) <- read.table("data/NKI295.pdata.txt", sep="\t", header=T, row.names=1) 
# fData is information about the probes
# pData is information about the samples
# these are now associated with the expression set

# we will first look at the values for ESR1
exprs(NKI295)[which(fData(NKI295)$symbol == "ESR1"),]

# lets make some plots!
# a barplot is a solid way to start visualising this
barplot(exprs(NKI295)[which(fData(NKI295)$symbol == "ESR1"),])

# we can clean it up a little with some extra parameters
barplot(exprs(NKI295)[which(fData(NKI295)$symbol == "ESR1"),], cex.names=0.2, las=2)

# and make it more interesting with some colours
subtype.colours <- c("red", "deeppink","dark blue", "light blue","orange" )
barplot(exprs(NKI295)[which(fData(NKI295)$symbol == "ESR1"),], cex.names=0.1, las=2, col=subtype.colours)

# it's a little confusing, so lets sort the bars to make a waterfall plot
NKI295 <- NKI295[,order(exprs(NKI295)[which(fData(NKI295)$symbol == "ESR1"),])]
barplot(exprs(NKI295)[which(fData(NKI295)$symbol == "ESR1"),], cex.names=0.2, las=2, col=subtype.colours[NKI295$subtype], main="ER expression by sample")

# lets save this plot into the data folder, it will be named 'subtype.ER.waterfall.pdf'
pdf("subtype.ER.waterfall.pdf", width=12, height=6)
barplot(exprs(NKI295)[which(fData(NKI295)$symbol == "ESR1"),], cex.names=0.2, las=2, col=subtype.colours[NKI295$subtype], main="ER expression by sample")
legend("bottomright", legend=levels(NKI295$subtype), fill=subtype.colours, bty="n")
dev.off()

# let's do something similar but use a boxplot instead
pdf("subtype.ER.boxplot.pdf")
boxplot(exprs(NKI295)[which(fData(NKI295)$symbol == "ESR1"),]~NKI295$subtype, col=subtype.colours, main="ER expression by subtype")
dev.off()

## End of section 1
# Optional tasks
# Try plotting a different gene - not ESR1!
# Try some different colours.
# Try ordering the bars backwards 


########################################################
## Section 2 - more complicated R, looping, and differential expression
# we will use a library called limma for this section
# it will help us conduct some differential expression from the dataset we were using earlier
	
# first we set up some parameters
limma.parameters <- cbind(intercept=rep(1, ncol(NKI295)), NKI295$subtype)

# fit linear models for the data set
limma.fit <- lmFit(exprs(NKI295), limma.parameters)

# apply Bayes smoothing
limma.fit <- eBayes(limma.fit)

# the top ten genes are collated
limma.table <- topTable(limma.fit, coef=2, number=nrow(NKI295), sort.by="none")

# check with Al about this bit 
means <- matrix(0, nrow=nrow(NKI295), ncol=length(levels(NKI295$subtype)))
for(i in 1:length(levels(NKI295$subtype))){
	means[,i] <- apply(exprs(NKI295)[,which(NKI295$subtype ==levels(NKI295$subtype)[i])], 1, mean)
	}
means.df <- as.data.frame(means)
names(means.df) <- levels(NKI295$subtype)
limma.table <- cbind.data.frame(limma.table, means.df, fData(NKI295))
rownames(limma.table) <- featureNames(NKI295)
	
# access the significantly affected genes
subtype.significant <- limma.table[which(limma.table$adj.P.Val < 1e-10),]
write.table(subtype.significant, "NKI295.subtype.significant.txt", sep="\t", row.names=F, na="")
NKI295_subtype.significant <- NKI295[rownames(subtype.significant),]
	
#	then make a heatmap
multiHeatmap(NKI295_subtype.significant, phenotypes="subtype", pheno.colours= subtype.colours, device="PDF", project="NKI295.subtype", width=6, height=9, plot.sample.names=F)

#	then read in the 70 gene list
signature_70genes <- read.table("data/vantVeer.70genes.txt", sep="\t", header=T)
signature_metastasis <- read.table("data/vantVeer.metastasis.txt", sep="\t", header=T)

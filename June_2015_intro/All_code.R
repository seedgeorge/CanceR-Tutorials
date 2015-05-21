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
## Section 1 - we have to do a little bit of R basics
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
# Optional tasks!
# Try plotting a different gene - not ESR1 - 


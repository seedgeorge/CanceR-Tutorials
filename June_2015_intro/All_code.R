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

library(Biobase)
library(marray)
library(limma)
library(survival)
library(oligo)
library(affy)
library(annotate)


## Section 1 - we have to do a little bit of R basics
# we will access the data from http://www.nejm.org/doi/full/10.1056/NEJMoa021967
# it is a microarray gene expression dataset of 295 breast cancers
# first we will read the data into R from a file
# this command looks for a specific file and makes it into a table that R can deal with
NKI295_data <- read.table("data/NKI295.expression.data.txt", sep="\t", header=T)

# what kind of data is this?
# the following commands will give you clues about the data without opening it to read
class(NKI295_data)
summary(NKI295_data)
nrow(NKI295_data)
ncol(NKI295_data)

# what is the difference between a data frame and a matrix?
# what kind of data is in each column?
head(NKI295_data)
apply(NKI295_data, 2, class)

# lets make a quick plot of all the scores 
# we will use the boxplot command, and select all the columns except the first six
boxplot(NKI295_data[,-1:-6], pch=16, cex=0.1)

#	but we also want to know which sample is which and what their phenotypes are
#	so we will read in some more data
NKI295_samples <- read.table("data/NKI295.samples.txt", sep="\t", header=T)
class(NKI295_samples)
head(NKI295_samples)

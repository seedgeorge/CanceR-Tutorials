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

#	what kind of data is this?
class(NKI295_data)
summary(NKI295_data)

#	what is the difference between a data frame and a matrix?
#	what kind of data is in each column?
apply(NKI295_data, 2, class)




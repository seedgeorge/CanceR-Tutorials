# Gene Expression analysis in R

# Section 0 - packages and libraries
# later on we might need some extra tools and systems
# so we'll download them now - select the lines below, and use Edit->Execute to run them
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

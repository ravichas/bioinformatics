# High-throughput Experiment

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(devtools)
library(rafalib)
BiocManager::install("Biobase")
BiocManager::install("geneplotter")
library(geneplotter)
library(Biobase)

# Example1: Let us explore a sample Expression Set data

library(Biobase)
# class??
class(4)
class(pi); pi
# Other way to get help
# ?ExpressionSet  # use tab to get more help
# ?"ExpressionSet"
data(sample.ExpressionSet)
sample.ExpressionSet
head(exprs(sample.ExpressionSet))

# Help in Bioconductor
# vignette("ExpressionSetIntroduction")
# vignette(package="Biobase")
# HTML Help
# browseVignettes(package="Biobase")
# help.start()


#---------------------------------------------------
# Microarray 
## MICROARRAY DATA
##--------------------------------------------------
# HIGH-THROUGHPUT EXPERIMENT 
#  How Bioconductor works? Example1 (based on Dr. Izarry/Mike Love'snotes)
# High-throughput technologies 
# features could be genes, locations of genome (single base)
#          genomic regions etc
# Samples (not related to statistical term, derived from population)
#         : Samples are also called Experimental units 
#         : Samples could be from diff. parts of a tumor


library(devtools)
install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
library(Biobase)
mypar(1,1)

.oldls = ls()
data(GSE5859Subset)
.newls = ls()
newstuff = setdiff(.newls, .oldls)
cls <- sapply( newstuff,function(x){class(get(x))}) 
cls
dim(geneExpression) # 8793 x 24 

# What are the experimental units?
# In this case 24 is the experimental unit 24 samples of blood
# How many features?
# 8793 genes dervied from blood
# What is the experimental type?
# Gene expression measurement experiment
# GSE5859 full data set was derived for testing expression
# pattern for differnt ethnic groups
# GSE5859Subset was restricted just one ethnic group
# but divided into two subgroups

boxplot(date~factor(ethnicity), data=sampleInfo, ylab="chip date")

# to find information about the gene BRCA2

all.equal(colnames(geneExpression), sampleInfo$filename)
all.equal(rownames(geneExpression), geneAnnotation$PROBEID)

# explore sapleInfo
head(sampleInfo)
g <- sampleInfo$group
g
# explore geneAnnotation
head(geneAnnotation)

#Method 1 (convoluted)
geneExpression[grep("BRCA2",geneAnnotation$SYMBOL),]

# Method 2 (Taken From Dr. Mike Love's notes)
ind = which(geneAnnotation$SYMBOL=="BRCA2")
boxplot(geneExpression[ind,]~sampleInfo$ethnicity, ylab="BRCA2 by hgfocus")

# It works in a complicated way. Bioconductor has simplified 
# the commonds in the following way

es1 <- ExpressionSet(geneExpression)
es1
pData(es1) <- sampleInfo     # phenotypic data
fData(es1) <- geneAnnotation # feature Data
es1
# Note we can use the same pData(es1) and fData(es1)
# to explore how es1 is constructed
class(pData(es1))

#boxplot(date~factor(ethnicity), data=sampleInfo, ylab="chip date")
boxplot(es1$date ~ es1$ethnicity, xlab = "Ethnicity", ylab="chip date")

# How many samples were created on 2005-10-28?
library(dplyr)
filter(sampleInfo, date == "2005-10-28")

# How many genes are part of CHR 2 ?
tbl <- table(geneAnnotation$CHR)
tbl
tbl[names(tbl)== "chr2"]

barplot(tbl, col = "green")


# log expresson of the gene BRCA1
id <- which(geneAnnotation$SYMBOL == "ERH")
geneExpression[id,]
dplyr::filter(geneAnnotation, SYMBOL == "ERH")


##****************************************
# Hypothesis Test                        *
# t-test, errors, alpha and beta         *
# Type-1 and Type-II                     *
##****************************************


# t-test and p-values (random variable)

population <- matrix(rnorm(211032), 8793,24)
nttest <- vector(mode = "numeric", nrow(population))
for (i in 1:nrow(population)) {
  nttest[i] <- t.test(population[i,1:12],population[i,13:24])$p.value
}
hist(nttest)

# How many genes are significant at 0.05 level 
sum(nttest < 0.05)/length(nttest)
hist(nttest)

#### Let us talk about t-test and how to carry out a t-test
e0 <- geneExpression[10,g == 0]
e1 <- geneExpression[10,g == 1]
e0
e1
hist(e0)
hist(e1)
t.test(e0,e1)

# let us do a q-q plot
library(rafalib)
mypar(1,2)
qqnorm(e0)
qqline(e0)

qqnorm(e1)
qqline(e1)

# random variables 

#### Let us talk about t-test and how to carry out a t-test
e0 <- geneExpression[1000,g == 0]
e1 <- geneExpression[1000,g == 1]
e0
e1
hist(e0)
hist(e1)
t.test(e0,e1)

# thousands of t-test
ttest <- vector(mode = "numeric", nrow(geneExpression))
g <- sampleInfo$group
for (i in 1:nrow(geneExpression)) {
   ttest[i] <- t.test(geneExpression[i,g == 0],geneExpression[i,g == 1])$p.value
}

## When we do the t-test? 
## We make mistakes 
## the 

##                                  OUR-CALL     
##                         Called Sig    Called insig   Total
##        Null is True     V            m0-V           m0
## TRUTH   Alt is True     S            m1-S           m1
##        Total True       R            m-R            m

# For any given set, we should be able to identify 
# the total TRUE correctly
# Type 1 error: False positive, V (number); alpha = prob = V/m0
# Type-II: error: False negative: S (number); beta = prob = S/m1

# In reality, we only know R and we do not know 
# What we know is the significant genes or units are 
# much smaller than the non-significant ones

# Note V and S are random numbers 

# Before genomics, in most cases, we do only one test
# In genomics, we do many many hundreds of tests
# so, accidently finding a test with < 0.05 is very possible.
# remember in a case where you know it is a null and 
# and carry out tests, the p-value distribution is uniform
# 

# error rates

# Probability of making at least one type-1 error among 
#      the total number of tests 

#---------------FWER-----
# FWER Family-wise Error Rate (FWER)
# Probability of making at least one error 
#  one or more false discoveries ot Type-1 errors when 
# mulforming multiple hypothesis tests
#   FWER = Pr(V >= 1)
#   FWER = 1 - Pr(V = 0 )

# Bonferroni Correction 

# Banajmini Hochberg correction
#  FDR = V/R
# two ways to estimate FDR  often called expected FDR
# Benjamini Hochberg is on eof them

# compute the p-values 
# sort the p-values smallest to largest
# select the p(i) <=  (i/m) * alpha 
# alpha = 0.05 ; m is th eotal number of tests
# i is just the index of p values 
# the number of genes that correspond to the 
# the index will have FDR of < 0.05 (in this case, alpha)

# Second method to esimate the FDR is BH metho
# BH method (Storey 2002)
# less conservative approach then Benjamini Hochberg 

# Let k be the largest i such that 
# pi_0 * P(i) <= (i/m) * alpha
# Note pi_0 is the proportion of genes for which H0 is true
# this is an unknown quantity and often estimated
# then reject H_i for i = 1...k 
# when pi_0 is one, we fall back to BH method

#EDA in high throughout analysis

# let us compute p-values, adj. pvalues, FDR 

library(genefilter)
library(GSE5859Subset)
data(GSE5859Subset)
g <- factor(sampleInfo$group)
tt_results <- rowttests(geneExpression,g)
p_values <- tt_results$p.value < 0.05
num_p_values <- length(p_values) 
num_p_values # 1383 genes 
# Apply Benajimi correction to achieve a FWER of 0.05
# let us find out how many genes are now significant
k <- 0.05/num_p_values
sum(tt_results$p.value < k)



library(genefilter)
library(GSE5859Subset)
data(GSE5859Subset)
g <- factor(sampleInfo$group)
tt_results <- rowttests(geneExpression,g)
hist(tt_results$p.value)

# let us set up a null

NR <- nrow(geneExpression)
NC <- ncol(geneExpression)
rData <- matrix(rnorm(NR*NC),NR,NC)
ntt_results <- rowttests(rData,g)
hist(ntt_results$p.value)

# Reporting only P-values is wrong
# we can have a small effect size and also small p-values
# report both The plot that accomplishes this is called 
# Volcano Plot
mypar(1,1)
plot(tt_results$dm, -log10(tt_results$p.value), 
       xlab = "Effect Size", ylab = "-log10(P-value)")

# EDA 2 tool boxplot can show the gene expression is 
# comparable sample distribution in the dataset?
# any abnormality in the data 
dim(geneExpression)
e <- geneExpression
# create an error 
e[,12] <- e[,12]/log2(exp(0.9)) 
boxplot(e, names=1:ncol(e),col=rainbow(24))
# go back to the original data
e[,12] <- e[,12]/log2(exp(0.9)) 
boxplot(e, names=1:ncol(e),col=rainbow(24))

## MA Plot in Microarrays

#Bland-Altman plot or in genomics MA plot
# The name MA comes from plots of 
# red log intensity minus (M) green intensities 
# versus average (A) log intensities

x <- e[,1]
y <- e[,2]
mypar(1,2)
plot(x,y)
# most data points are towards the bottom-left. 
# let us view them using smoothscatter plot
smoothScatter(x,y, col = "blue")
A <- (x+y)/2; M <- x-y 
plot(A,M)
smoothScatter(A,M, col = "blue")
# most points are twoards the bottom-left

# GSE5859 
# Let us load the full GSE5859 set

gse <- getGEO("GSE5859")
e <- exprs(gse[[1]])

# This time you get a single expression set constructed and 
# loaded for you to work with 
#  hgfocus array
dim(e)


mypar(1,1)
boxplot(e[,1:20], col = rainbow(20))

#boxplot(ge, names=1:ncol(e),col=rainbow(208))

#### ENDS MICROARRAY DATA


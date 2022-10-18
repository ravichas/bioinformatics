# based on bioconductor workshop

# libraries used in this tutorial
# limma
# Glimma
# edgeR
# Mus.musculus
# RColorBrewer
# gplots

##   SessionInfo 
# R version 4.1.1 (2021-08-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.6
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] gplots_3.1.3                              RColorBrewer_1.1-3                       
# [3] Mus.musculus_1.3.1                        TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0
# [5] org.Mm.eg.db_3.14.0                       GO.db_3.14.0                             
# [7] OrganismDbi_1.36.0                        GenomicFeatures_1.46.5                   
# [9] GenomicRanges_1.46.1                      GenomeInfoDb_1.30.1                      
# [11] AnnotationDbi_1.56.2                      IRanges_2.28.0                           
# [13] S4Vectors_0.32.4                          edgeR_3.36.0                             
# [15] Glimma_2.4.0                              limma_3.50.3                             
# [17] BiocManager_1.30.18                       Biobase_2.54.0                           
# [19] BiocGenerics_0.40.0                      
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-7                matrixStats_0.62.0          bit64_4.0.5                
# [4] filelock_1.0.2              progress_1.2.2              httr_1.4.3                 
# [7] tools_4.1.1                 utf8_1.2.2                  R6_2.5.1                   
# [10] KernSmooth_2.23-20          DBI_1.1.3                   colorspace_2.0-3           
# [13] tidyselect_1.1.2            prettyunits_1.1.1           DESeq2_1.34.0              
# [16] bit_4.0.4                   curl_4.3.2                  compiler_4.1.1             
# [19] graph_1.72.0                cli_3.3.0                   xml2_1.3.3                 
# [22] DelayedArray_0.20.0         rtracklayer_1.54.0          caTools_1.18.2             
# [25] scales_1.2.0                genefilter_1.76.0           RBGL_1.70.0                
# [28] rappdirs_0.3.3              Rsamtools_2.10.0            stringr_1.4.0              
# [31] digest_0.6.29               R.utils_2.12.0              XVector_0.34.0             
# [34] pkgconfig_2.0.3             htmltools_0.5.3             MatrixGenerics_1.6.0       
# [37] dbplyr_2.2.1                fastmap_1.1.0               htmlwidgets_1.5.4          
# [40] rlang_1.0.4                 rstudioapi_0.13             RSQLite_2.2.15             
# [43] BiocIO_1.4.0                generics_0.1.3              jsonlite_1.8.0             
# [46] gtools_3.9.3                BiocParallel_1.28.3         dplyr_1.0.9                
# [49] R.oo_1.25.0                 RCurl_1.98-1.8              magrittr_2.0.3             
# [52] GenomeInfoDbData_1.2.7      Matrix_1.4-1                Rcpp_1.0.9                 
# [55] munsell_0.5.0               fansi_1.0.3                 lifecycle_1.0.1            
# [58] R.methodsS3_1.8.2           yaml_2.3.5                  stringi_1.7.8              
# [61] SummarizedExperiment_1.24.0 zlibbioc_1.40.0             BiocFileCache_2.2.1        
# [64] grid_4.1.1                  blob_1.2.3                  parallel_4.1.1             
# [67] crayon_1.5.1                lattice_0.20-45             Biostrings_2.62.0          
# [70] splines_4.1.1               annotate_1.72.0             hms_1.1.1                  
# [73] KEGGREST_1.34.0             locfit_1.5-9.6              pillar_1.8.0               
# [76] rjson_0.2.21                geneplotter_1.72.0          biomaRt_2.50.3             
# [79] XML_3.99-0.10               glue_1.6.2                  png_0.1-7                  
# [82] vctrs_0.4.1                 gtable_0.3.0                purrr_0.3.4                
# [85] assertthat_0.2.1            cachem_1.0.6                ggplot2_3.3.6              
# [88] xtable_1.8-4                restfulr_0.0.15             survival_3.3-1             
# [91] tibble_3.1.8                GenomicAlignments_1.30.0    memoise_2.0.1              
# [94] ellipsis_0.3.2             



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RNAseq123", ask=FALSE)
BiocManager::install("Glimma", ask = FALSE)


url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"  
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb")   
utils::untar("GSE63310_RAW.tar", exdir = ".")  
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
           "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")

for(i in paste(files, ".gz", sep=""))  
  R.utils::gunzip(i, overwrite=TRUE) 

# Goals
# learn how to analyse RNA-seq data
# identify methods for pre-processing data
# understand linear models used in differential expression analysis
# examine plots for data exploration and result representation
# 
# Objectives

# read in count data and format as a DGEList-object
# annotate Entrez gene identifiers with gene information
# filter out lowly expressed genes
# normalise gene expression values
# unsupervised clustering of samples (standard and interactive plots)
# linear modelling for comparisons of interest
# remove heteroscedascity
# examine the number of differentially expressed genes
# mean-difference plots (standard and interactive plots)
# heatmaps


suppressPackageStartupMessages({
  library(limma)
  library(Glimma)
  library(edgeR)
  library(Mus.musculus)
})


getwd()

setwd("~/Downloads")
# 
getwd()


# STOP
dir.create("Law_RNAseq123")
#> Warning in dir.create("Law_RNAseq123"): 'Law_RNAseq123' already exists
setwd("Law_RNAseq123")

url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb")
utils::untar("GSE63310_RAW.tar", exdir = ".")

files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt",
           "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt",
           "GSM1545541_JMS8-4.txt", "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt",
           "GSM1545545_JMS9-P8c.txt")

for(i in paste(files, ".gz", sep=""))
  R.utils::gunzip(i, overwrite=TRUE)

setwd("./..")
read.delim(file.path("Law_RNAseq123", files[1]), nrow=5)


x <- readDGE(file.path("Law_RNAseq123", files), columns=c(1,3))
class(x)

dim(x)


samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
samplenames


colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", 
                     "Basal", "ML", "LP"))
x$samples$group <- group
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x$samples$lane <- lane
x$samples

# Organizing gene annotations

geneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENTREZID")
head(genes)


genes <- genes[!duplicated(genes$ENTREZID),]


x$genes <- genes
x


# Data Transform

cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

# Removing low expressing genes (esp. look for zero counts across 9 samples)
table(rowSums(x$counts==0)==9)


keep.exprs <- rowSums(cpm>1)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)


# plot the data

library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")


# normalize gene expression distributions

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors


# To give a better visual representation of the 
# effects of normalisation, the data was duplicated then 
# adjusted so that the counts of the first sample are 
# reduced to 5% of their original values, and in the 
# second sample they are inflated to be 5-times larger.

x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5


# plot the results

par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors
#> [1] 0.05472223 6.13059440 1.22927355 1.17051887 1.21487709 1.05622968
#> [7] 1.14587663 1.26129350 1.11702264
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")


# unsupervised clustering of samples

lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")

#plot MDS

glMDSPlot(lcpm, labels=paste(group, lane, sep="_"), 
          groups=x$samples[,c(2,5)], launch=FALSE)

# create design matrix

design <- model.matrix(~0+group+lane)
colnames(design) <- gsub("group", "", colnames(design))
design


contr.matrix <- makeContrasts(
  BasalvsLP = Basal-LP, 
  BasalvsML = Basal - ML, 
  LPvsML = LP - ML, 
  levels = colnames(design))
contr.matrix

# removing heteroscadescity from the data

par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v


vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

# Est the number of DE genes

summary(decideTests(efit))


tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)


de.common <- which(dt[,1] != 0 & dt[,2] != 0)
length(de.common)

head(tfit$genes$SYMBOL[de.common], n = 20)

# VD showing the # of genes DE in the comparison 
# between basal vs LP only, basal vs ML only and the
# common ones. 
# The # that are not in DE are shown in bottom right

vennDiagram(dt[,1:2], circle.col = c("turquoise","salmon"))

write.fit(tfit, dt, file="results.txt")

# Examining DE genes from top to bottom

basal.vs.lp <- topTreat(tfit, coef=1, n=Inf)
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf)
head(basal.vs.lp)

head(basal.vs.ml)



# graphical rep of DE results

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))


glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         side.main="ENTREZID", counts=x$counts, groups=group, launch=FALSE)


# Heat map for the top 50 genes

# Heatmap of log-CPM values for top 100 genes DE in basal versus LP. Expression across each gene (or row) have been scaled so that mean expression is zero and standard deviation is one. Samples with relatively high expression of a given gene are marked in red and samples with relatively low expression are marked in blue. Lighter shades and white represent genes with intermediate expression levels. Samples and genes have been reordered by the method of hierarchical clustering. A dendrogram is shown for the sample clustering.


library(gplots)
basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(v$E[i,], scale="row",
          labRow=v$genes$SYMBOL[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")

sessionInfo()

## Import 5 different packages for Single-Cell RNA-Seq differential analysis

## BASiCS

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BASiCS", version = "3.8")

## DEsingle

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DEsingle", version = "3.8")

## Linnorm

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Linnorm", version = "3.8")

##monocle

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("monocle", version = "3.8")

## MAST

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MAST", version = "3.8")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocParallel", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rhdf5lib", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocParallel", version = "3.8")


install.packages("MASS")

browseVignettes("BASiCS")

library("BASiCS")

Data <- makeExampleBASiCS_Data()
Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)

Chain <- BASiCS_MCMC(Data = Data, N = 1000, Thin = 10, Burn = 500, 
                     PrintProgress = FALSE, Regression = TRUE)


Chain <- BASiCS_MCMC(Data, N = 1000, Thin = 10, Burn = 500,
                     Regression = FALSE, PrintProgress = FALSE, WithSpikes = TRUE)
DC <- BASiCS_DenoisedCounts(Data, Chain)
DR <- BASiCS_DenoisedRates(Data, Chain)

plot(Chain, Param = 'mu', Gene = 1)
plot(Chain, Param = 'phi', Cell = 1)
plot(Chain, Param = 'theta', Batch = 1)

ChainSummary <- Summary(Chain)

plot(ChainSummary, Param = 'mu', main = 'All genes')
plot(ChainSummary, Param = 'mu', Genes = 1:10, main = 'First 10 genes')
plot(ChainSummary, Param = 'phi', main = 'All cells')
plot(ChainSummary, Param = 'phi', Cells = 1:5, main = 'First 5 cells')
plot(ChainSummary, Param = 'theta')


DataHVG <- BASiCS_DetectHVG(Chain, VarThreshold = 0.6, ProbThreshold = NULL, EFDR = 0.05,
                 OrderVariable = "Prob", Plot = TRUE)

DataLVG <- BASiCS_DetectLVG(Chain, VarThreshold = 0.40,
                              EFDR = 0.05, Plot = TRUE)



data(ChainSC)
data(ChainRNA)
Test <- BASiCS_TestDE(Chain1 = ChainSC, Chain2 = ChainRNA,
                      GroupLabel1 = 'SC', GroupLabel2 = 'P&S',
                      EpsilonM = log2(1.5), EpsilonD = log2(1.5),
                      OffSet = TRUE)
head(Test$TableMean)


## Data with spike-ins
# Expression counts
set.seed(1)
Counts <- matrix(rpois(200*20, 2), ncol = 20)
rownames(Counts) <- c(paste0('Gene', 1:40), paste0('ERCC', 1:10))
# Technical information
Tech <- grepl("ERCC", rownames(Counts))
# Spikes input number of molecules
set.seed(2)
SpikeInfo <- data.frame(gene = rownames(Counts)[Tech],
                        amount = rgamma(10, 1, 1))
# Creating a BASiCS_Data object (no batch effect)
DataExample <- newBASiCS_Data(Counts, Tech = Tech, SpikeInfo = SpikeInfo)
# Creating a BASiCS_Data object (with batch effect)
BatchInfo <- c(rep(1, 5), rep(2, 5))
DataExample <- newBASiCS_Data(Counts, Tech = Tech,
                              SpikeInfo = SpikeInfo, BatchInfo = BatchInfo)

Chain <- BASiCS_MCMC(DataExample, N = 1000, Thin = 10, Burn = 500,
                     Regression = TRUE, PrintProgress = FALSE, WithSpikes = TRUE)

## Data without spike-ins (BatchInfo is required)
# Expression counts
set.seed(1)
Counts <- matrix(rpois(70*30, 2), ncol = 30)
rownames(Counts) <- paste0('Gene', 1:70)
BatchInfo <- c(rep(1, 15), rep(2, 15))
# Creating a BASiCS_Data object (with batch effect)
DataExample <- newBASiCS_Data(Counts, BatchInfo = BatchInfo)

Chain <- BASiCS_MCMC(DataExample, N = 1000, Thin = 10, Burn = 500,
                     Regression = TRUE, PrintProgress = FALSE, WithSpikes = FALSE)

Counts2 <- matrix(rpois(70*30, 2), ncol = 30)
rownames(Counts2) <- paste0('Gene', 1:70)
BatchInfo2 <- c(rep(1, 15), rep(2, 15))
# Creating a BASiCS_Data object (with batch effect)
DataExample2 <- newBASiCS_Data(Counts, BatchInfo = BatchInfo)

Chain2 <- BASiCS_MCMC(DataExample2, N = 1000, Thin = 10, Burn = 500,
                     Regression = TRUE, PrintProgress = FALSE, WithSpikes = FALSE)


Test <- BASiCS_TestDE(Chain1 = Chain, Chain2 = Chain2,
                      GroupLabel1 = 'Group1', GroupLabel2 = 'Group2',
                      EpsilonM = log2(1.5), EpsilonD = log2(1.5),
                      OffSet = TRUE)

label <- c(rep(1, 5), rep(0, 5))

simulateddata <- RnaXSim(expMatrix, distribution="NB", NumRep=10, NumDiff = 100, NumFea = 1000)
BatchInfo <- c(rep(1, 5), rep(2, 5))
trt <- newBASiCS_Data(simulateddata$data[, 1:10], BatchInfo = BatchInfo)
grp1 <- BASiCS_MCMC(trt, N = 100, Thin = 10, Burn = 50,
                     Regression = TRUE, PrintProgress = FALSE, WithSpikes = FALSE)
ctl <- 


library("DEsingle")

data(TestData)

group <- factor(c(rep(1,50), rep(2,100)))
results <- DEsingle(counts = counts, group = group)
results.classified <- DEtype(results = results, threshold = 0.05)


group <- factor(c(rep(1,10), rep(2,10)))
results <- DEsingle(counts = simulateddata$data, group = group)
results.classified <- DEtype(results = results, threshold = 0.05)
summary(results.classified)
DEsingle_sig <- rownames(results[which(results$pvalue.adj.FDR <= 0.05),])
DEsingle.s <- as.numeric(gsub("\\D", "", DEsingle_sig)) 


truev <- c(rep(0, 1000))
truev[simulateddata$DiffList] <- 1

DEsingle.v <- c(rep(0, 1000))
DEsingle.v[DEsingle.s] <- 1

install.packages("e1071")
library(e1071)

install.packages("caret")
library(caret)

install.packages ('ROCR')
library(ROCR)

cf6=confusionMatrix(as.factor(DEsingle.v), as.factor(truev))
cf6$positive
cf7=confusionMatrix(as.factor(truev), as.factor(DEsingle.v))
cf7

a <- table(DEsingle.v, truev)

FDR <- a[2, 1]/(a[2,1]+a[2,2])
sensitivity <- a[2, 2]/(a[1,2] + a[2, 2])
specificity <- a[1, 1]/(a[1,1]+a[2,1])
accuracy <- (a[1, 1] + a[2, 2])/ngene

install.packages("pROC")
library(pROC)
roc_obj <- roc(truev, DEsingle.v)
aucv <- auc(roc_obj)

pred=prediction(DEsingle.v, truev)
roc.perf=performance(pred, measure="tpr", x.measure="fpr")
plot(roc.perf, colorize=TRUE)
abline(0,1, lty=2)
auc=performance(pred, "auc")
b <- auc@y.values

simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

simp_roc <- simple_roc(DEsingle.v, truev)
with(simp_roc, lines(1 - FPR, TPR, col="blue", lty=2))




library("Linnorm")
data(Islam2011)
HClust.results <- Linnorm.HClust(Islam2011, Group=c(rep("ESC",48), rep("EF",44), rep("NegCtrl",4)))
results <- Linnorm.HVar(Islam2011)
PCA.results <- Linnorm.PCA(Islam2011)
StableGenes <- Linnorm.SGenes(Islam2011)
tSNE.results <- Linnorm.tSNE(Islam2011)



#Obtain example matrix:
data(LIHC)
#Transformation:
transformedExp <- Linnorm(LIHC)
normalizedExp <- Linnorm(LIHC)



designmatrix <- c(rep(1,5),rep(2,5))
designmatrix <- model.matrix(~ 0+factor(designmatrix))
colnames(designmatrix) <- c("group1", "group2")
rownames(designmatrix) <- colnames(LIHC)
#DEG analysis
DEGResults <- Linnorm.limma(LIHC, designmatrix)
sum(DEGResults$adj.P.Val <= 0.05)



designmatrix <- model.matrix(~ 0+factor(group))
colnames(designmatrix) <- c("group1", "group2")
DEGResults <- Linnorm.limma(simulateddata$data, designmatrix)
dim(DEGResults)
sum(DEGResults$adj.P.Val <= 0.05)


#Obtain example matrix:
data(SEQC)
expMatrix <- SEQC
#Example for Negative Binomial distribution
simulateddata <- RnaXSim(expMatrix, distribution="NB", NumRep=5, NumDiff = 200, NumFea = 2000)
#Example for Poisson distribution
simulateddata <- RnaXSim(expMatrix, distribution="Poisson", NumRep=5, NumDiff = 200, NumFea = 2000)
#Example for Log Normal distribution
simulateddata <- RnaXSim(expMatrix, distribution="LogNorm", NumRep=5, NumDiff = 200, NumFea = 2000)
#Example for Gamma distribution
simulateddata <- RnaXSim(expMatrix, distribution="Gamma", NumRep=5, NumDiff = 200, NumFea = 2000)

summary(simulateddata$data)
length(simulateddata$DiffList)
length(simulateddata$CorGroup1)
length(simulateddata$CorGroup2)
head(simulateddata$data)

library(monocle)

# This example performs a greatly simplified version of the single-cell RNA-Seq
# analysis of skeletal myoblast differentiation
# described in Trapnell, Cacchiarelli et al (Nature Biotechnology, 2014).


label2 <- data.frame(c(rep(1,5),rep(2,5))) 
rownames(label2) <- colnames(simulateddata$data)
colnames(label2) <- "Group"
pd <- new("AnnotatedDataFrame", data = label2)
mydata <- newCellDataSet(simulateddata$data, phenoData = pd)

# Count how many cells each gene is expressed in, and how many
# genes are expressed in each cell
HSMM <- detectGenes(mydata, min_expr = 0.1)

# Get a list of genes expressed in at least 50 cells
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 0))
# Test the above genes for differential expression in response from switch from GM to DM
# Note: this step can take several hours on a single core, so you might want to parallelize it
# with the 'cores' argument 
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,], fullModelFormulaStr="~sm.ns(Group, df=1)", relative_expr = FALSE, cores=24)


diff_test_res <- differentialGeneTest(mydata, fullModelFormulaStr="~sm.ns(Group, df=1)", relative_expr = FALSE, cores=24)

monoResults <- differentialGeneTest(mydata2[expressed_genes,], fullModelFormulaStr="~1 + Group", relative_expr = FALSE, cores=24)

sum(diff_test_res$qval <= 0.05)

monocle_diff <- which(diff_test_res$qval <= 0.05)

Siggenelist <- simulateddata$DiffList

intersect(monocle_diff, Siggenelist)

# Use the differentially expressed genes as the basis for ordering the cells
# by progress through differentiation

library(MAST)

data(vbetaFA)
svbeta <- subset(vbetaFA, ncells==1)
svbeta <- svbeta[freq(svbeta)>.4,]
window <- function(x1) lapply(assays(x1), function(x2) x2[1:3, 1:6])
#total residuals of the response
z1 <- zlm(~ Stim.Condition, svbeta, hook=discrete_residuals_hook)
window(collectResiduals(z1, svbeta))
z2 <- zlm(~ Stim.Condition, svbeta, hook=continuous_residuals_hook)
window(collectResiduals(z2, svbeta))
z3 <- zlm(~ Stim.Condition, svbeta, hook=combined_residuals_hook)
window(collectResiduals(z3, svbeta))
#total deviance residuals
z4 <- zlm(~ Stim.Condition, svbeta, hook=deviance_residuals_hook)
window(collectResiduals(z4, svbeta))
#partial residuals
colData(svbeta)$ngeneson <- colMeans(assay(svbeta)>0)
z5 <- zlm(~ Stim.Condition + ngeneson, svbeta)
partialScore(z5, 'Stim.Condition')


data(vbetaFA)
zlmVbeta <- zlm(~ Stim.Condition+Population, subset(vbetaFA, ncells==1)[1:10,])
#Coefficients and standard errors
coef(zlmVbeta, 'D')
coef(zlmVbeta, 'C')
se.coef(zlmVbeta, 'C')
#Test for a Population effect by dropping the whole term (a 5 degree of freedom test)
lrTest(zlmVbeta, 'Population')
#Test only if the VbetaResponsive cells differ from the baseline group
lrTest(zlmVbeta, CoefficientHypothesis('PopulationVbetaResponsive'))
# Test if there is a difference between CD154+/Unresponsive and CD154-/Unresponsive.
# Note that because we parse the expression
# the columns must be enclosed in backquotes
# to protect the \quote{+} and \quote{-} characters.
lrTest(zlmVbeta, Hypothesis('`PopulationCD154+VbetaUnresponsive` -
                            `PopulationCD154-VbetaUnresponsive`'))
waldTest(zlmVbeta, Hypothesis('`PopulationCD154+VbetaUnresponsive` -
                              `PopulationCD154-VbetaUnresponsive`'))

suppressPackageStartupMessages({
  library(ggplot2)
  library(GGally)
  library(GSEABase)
  library(limma)
  library(reshape2)
  library(data.table)
  library(knitr)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(stringr)
  library(NMF)
  library(rsvd)
  library(RColorBrewer)
  library(MAST)
})

install.packages("GGally")
install.packages("NMF")
install.packages("rsvd")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GSEABase", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", version = "3.8")

options(mc.cores = 1)
knitr::opts_chunk$set(message = FALSE,error = FALSE,warning = FALSE,cache = FALSE,fig.width=8,fig.height=6, auto.dep=TRUE)

freq_expressed <- 0.2
FCTHRESHOLD <- log2(1.5)

data(maits, package='MAST')
dim(maits$expressionmat)
maits$expressionma[1:3, 1:10]
head(maits$expressionmat)
head(maits$cdat)
head(maits$fdat)


scaRaw <- FromMatrix(t(maits$expressionmat), maits$cdat, maits$fdat)
dim(scaRaw)

dev.off()
x <- as.matrix(assay(scaRaw[1:100, ]))

aheatmap(x)

aheatmap(assay(scaRaw[1:1000,]), labRow='', annCol=as.data.frame(colData(scaRaw)[,c('condition', 'ourfilter')]), distfun='spearman')


set.seed(123)
plotPCA <- function(sca_obj){
  projection <- rpca(t(assay(sca_obj)), retx=TRUE, k=4)$x
  colnames(projection)=c("PC1","PC2","PC3","PC4")
  pca <- data.table(projection,  as.data.frame(colData(sca_obj)))
  print(ggpairs(pca, columns=c('PC1', 'PC2', 'PC3', 'libSize', 'PercentToHuman', 'nGeneOn', 'exonRate'),
                mapping=aes(color=condition), upper=list(continuous='blank')))
  invisible(pca)
}

plotPCA(scaRaw)

filterCrit <- with(colData(scaRaw), pastFastqc=="PASS"& exonRate >0.3 & PercentToHuman>0.6 & nGeneOn> 4000)

sca <- subset(scaRaw,filterCrit)
eid <- select(TxDb.Hsapiens.UCSC.hg19.knownGene,keys = mcols(sca)$entrez,keytype ="GENEID",columns = c("GENEID","TXNAME"))
ueid <- unique(na.omit(eid)$GENEID)
sca <- sca[mcols(sca)$entrez %in% ueid,]
## Remove invariant genes
sca <- sca[sample(which(freq(sca)>0), 6000),]

cdr2 <-colSums(assay(sca)>0)
qplot(x=cdr2, y=colData(sca)$nGeneOn) + xlab('New CDR') + ylab('Old CDR')
colData(sca)$cngeneson <- scale(cdr2)
plotPCA(sca)


scaSample <- sca[sample(which(freq(sca)>.1), 20),]
flat <- as(scaSample, 'data.table')
ggplot(flat, aes(x=value))+geom_density() +facet_wrap(~symbolid, scale='free_y')

thres <- thresholdSCRNACountMatrix(assay(sca), nbins = 20, min_per_bin = 30)
par(mfrow=c(5,4))
plot(thres)

assays(sca) <- list(thresh=thres$counts_threshold, tpm=assay(sca))
expressed_genes <- freq(sca) > freq_expressed
sca <- sca[expressed_genes,]


label2 <- data.frame(c(rep(1, 10), rep(2, 10))) 
rownames(label2) <- colnames(simulateddata$data)
colnames(label2) <- "Group"
scaRaw2 <- FromMatrix(simulateddata$data, label2, check_sanity = FALSE)
dim(scaRaw2)

cond<-factor(colData(scaRaw2)$Group)
#cond<-relevel(cond,"Unstim")
colData(scaRaw2)$Group<-cond
zlmCond <- zlm(~Group, scaRaw2)
summaryCond <- summary(zlmCond, doLRT='Group2')
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast=='Group2' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='Group2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
MAST_sig <- fcHurdle$primerid[which(fcHurdle$fdr <= 0.05)]
MAST.s <- as.numeric(gsub("\\D", "", MAST_sig))


fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)


cond<-factor(colData(sca)$condition)
cond<-relevel(cond,"Unstim")
colData(sca)$condition<-cond
zlmCond <- zlm(~condition + cngeneson, sca)

summaryCond <- summary(zlmCond, doLRT='conditionStim')
print(summaryCond, n=4)

summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast=='conditionStim' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='conditionStim' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)

entrez_to_plot <- fcHurdleSig[1:50,primerid]
symbols_to_plot <- fcHurdleSig[1:50,symbolid]
flat_dat <- as(sca[entrez_to_plot,], 'data.table')
ggbase <- ggplot(flat_dat, aes(x=condition, y=thresh,color=condition)) + geom_jitter()+facet_wrap(~symbolid, scale='free_y')+ggtitle("DE Genes in Activated MAIT Cells")
ggbase+geom_violin() 

flat_dat[,lmPred:=lm(thresh~cngeneson + condition)$fitted, key=symbolid]
ggbase +aes(x=cngeneson) + geom_line(aes(y=lmPred), lty=1) + xlab('Standardized Cellular Detection Rate')


MM <- model.matrix(~condition,unique(colData(sca)[,c("condition"),drop=FALSE]))
rownames(MM) <- str_extract(rownames(MM), 'Stim|Unstim')
predicted <- predict(zlmCond,modelmatrix=MM)

## Avert your eyes...
predicted[, primerid:=as.character(primerid)]
predicted_sig <- merge(mcols(sca), predicted[primerid%in%entrez_to_plot], by='primerid')
predicted_sig <- as.data.table(predicted_sig)

## plot with inverse logit transformed x-axis
ggplot(predicted_sig)+aes(x=invlogit(etaD),y=muC,xse=seD,yse=seC,col=sample)+
  facet_wrap(~symbolid,scales="free_y")+theme_linedraw()+
  geom_point(size=0.5)+scale_x_continuous("Proportion expression")+
  scale_y_continuous("Estimated Mean")+
  stat_ell(aes(x=etaD,y=muC),level=0.95, invert='x')

dev.off()
mat_to_plot <- assay(sca[entrez_to_plot,])
rownames(mat_to_plot) <- symbols_to_plot
aheatmap(mat_to_plot,annCol=colData(sca)[,"condition"],main="DE genes",col=rev(colorRampPalette(colors = brewer.pal(name="PiYG",n=10))(20)))

table(colData(sca)$beta, exclude=NULL)
scaHasBeta <- subset(sca, !is.na(beta))

scaDE <- sca[entrez_to_plot,]
zlmResidDE <- zlm(~condition + cngeneson, scaDE, hook=deviance_residuals_hook)
residDE <- zlmResidDE@hookOut
residDEMatrix <- do.call(rbind, residDE)

assays(scaDE) <- c(assays(scaDE), list(resid=residDEMatrix))
scaResidFlat <- as(scaDE, 'data.table')
scaResidFlat[1:4,]

ggplot(scaResidFlat, aes(x=ngeneson, y=resid))+geom_point(aes(col=condition))+geom_smooth()+facet_wrap(~symbolid)

boots <- bootVcov1(zlmCond, R=4)

module <- "BTM"
min_gene_in_module <- 5
packageExt <- system.file("extdata", package='MAST')
module_file <- list.files(packageExt, pattern = module, full.names = TRUE)
gene_set <- getGmt(module_file)
gene_ids <- geneIds(gene_set)
gene_ids <- gene_ids[!names(gene_ids)%like%"TBA"&!names(gene_ids)%like%"B cell"]
sets_indices <- limma::ids2indices(gene_ids, mcols(sca)$symbolid)
# Only keep modules with at least min_gene_in_module
sets_indices <- sets_indices[sapply(sets_indices, length) >= min_gene_in_module]

gsea <- gseaAfterBoot(zlmCond, boots, sets_indices, CoefficientHypothesis("conditionStim")) 
z_stat_comb <- summary(gsea, testType='normal')

sigModules <- z_stat_comb[combined_adj<.01]
gseaTable <- melt(sigModules[,.(set, disc_Z, cont_Z, combined_Z)], id.vars='set')
ggplot(gseaTable, aes(y=set, x=variable, fill=value))+geom_raster() + scale_fill_distiller(palette="PiYG")



#############################################################################################################
#
#
# Single cell RNA-seq analysis examples
#
#
#############################################################################################################

data(Islam2011)
Islam2011[1:20, 1:10]
dim(Islam2011)


## remove spike in data and negative controls from the count data
example.data <- Islam2011[9:dim(Islam2011)[1], 1:92]
head(example.data)
dim(example.data)
# 14905    92

sum(example.data == 0)/(dim(example.data)[1]*dim(example.data)[2])
#[1] 0.7173986

## lable of the group 1 has 48 mouse embryonic stem cells and group 2 has 44 mouse embryonic fibroblasts
group <- factor(c(rep(1, 48), rep(2, 44)))

#DEsingle method
results <- DEsingle(counts = example.data, group = group)
head(results)
sum(results$pvalue.adj.FDR <= 0.05)
#[1] 9031
DEsingle.sig <- which(results$pvalue.adj.FDR <= 0.05)

## Linnorm method
designmatrix <- model.matrix(~ 0+factor(group))
colnames(designmatrix) <- c("group1", "group2")
DEGResults <- Linnorm.limma(example.data, designmatrix)
dim(DEGResults)
sum(DEGResults$adj.P.Val <= 0.05)
#[1] 5423
Linnorm.sig <- which(DEGResults$adj.P.Val <= 0.05)

## monocle method

label2 <- data.frame(c(rep(1,48),rep(2,44))) 
rownames(label2) <- colnames(example.data)
colnames(label2) <- "Group"
pd <- new("AnnotatedDataFrame", data = label2)
mydata <- newCellDataSet(example.data, phenoData = pd)
example.data2 <- detectGenes(mydata, min_expr = 0.1)
# Get a list of genes expressed in at least 50 cells
expressed_genes <- row.names(subset(fData(example.data2), num_cells_expressed >= 0))
 
diff_test_res <- differentialGeneTest(example.data2[expressed_genes,], fullModelFormulaStr="~sm.ns(Group, df=1)", relative_expr = FALSE, cores=24)

diff.adjp <- p.adjust(diff_test_res$pval, "BH")
sum(diff.adjp <= 0.05)
#[1] 13892
monoResults <- differentialGeneTest(example.data2[expressed_genes,], fullModelFormulaStr="~1 + Group", relative_expr = FALSE, cores=24)
sum(diff_test_res$qval <= 0.05)
#[1] 13892
sum(monoResults$qval <= 0.05)
#[1] 13892
monocle.sig <- which(diff_test_res$qval <= 0.05)


## MAST
MAST_data <- FromMatrix(example.data, label2, check_sanity = FALSE)
cond<-factor(colData(MAST_data)$Group)
#cond<-relevel(cond,"Unstim")
colData(MAST_data)$Group<-cond
zlmCond <- zlm(~Group, MAST_data)
summaryCond <- summary(zlmCond, doLRT='Group2')
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast=='Group2' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='Group2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
sum(fcHurdle$fdr <= 0.05)
#[1] 8671
MAST.sig <- which(fcHurdle$fdr <= 0.05)

## DESeq2

group.con <- data.frame(cbind(c(rep("Group1", 48), rep("Group2", 44)), c(rep("single-read", 92))))
colnames(group.con) <- c("condition", "type")
rownames(group.con) <- colnames(example.data)
DESeq2.dds <- DESeqDataSetFromMatrix(countData = example.data, colData = group.con, design = ~condition)
DESeq2.test <- DESeq(DESeq2.dds, quiet = TRUE)
DESeq2.pval <- results(DESeq2.test)$padj
sum(DESeq2.pval <= 0.05, na.rm = TRUE)
#[1] 5108
DESeq2.sig <- which(DESeq2.pval<= 0.05) 
#make signifnicant genes vector

install.packages("vennDiagram")
library(VennDiagram)
## venn diagram plot
venn.plot <- venn.diagram(x=list("DEsingle" = DEsingle.sig, "Linnorm" = Linnorm.sig, "monocle" = monocle.sig,
                                 "MAST"= MAST.sig, "DESeq2" = DESeq2.sig), filename = "/Users/dli3/Box Sync/Conferences/Huston_2018/Islam2011vennDiagram.png",
                          col = "black",
                          fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                          alpha = 0.50,
                          cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
                                  1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
                          cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                          cat.cex = 0.8,
                          cat.fontface = "bold",
                          margin = 0.05)

reject <- c(9031, 5423, 13892, 8671, 5108)
#barplot(reject)
#barplot(reject, main="Empirical Power of Different Single-cell RNA-seq Differential Analysis Methods",
#        xlab="Number of Significant Genes Detected", col=c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
#        legend = c("DEsingle", "Linnorm", "moncle", "MAST", "DESeq2"), beside=TRUE)

## Make the frequencies numbers (rather than factors)
## dat$freqs <- as.numeric(as.character(dat$freqs))
## Find a range of y's that'll leave sufficient space above the tallest bar
ylim <- c(0, 1.1*max(reject))
## Plot, and store x-coordinates of bars in xx
xx <- barplot(reject, xaxt = 'n', xlab = '', width = 0.85, ylim = ylim,
              main = "", col=c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
              ylab = "Frequency")
## Add text at top of bars
text(x = xx, y = reject, label = reject, pos = 3, cex = 1.0, col = "black")
## Add x-axis labels 
axis(1, at=xx, labels=c("DEsingle", "Linnorm", "moncle", "MAST", "DESeq2"), tick=FALSE, las=2, line=-0.5, cex.axis=1.0)



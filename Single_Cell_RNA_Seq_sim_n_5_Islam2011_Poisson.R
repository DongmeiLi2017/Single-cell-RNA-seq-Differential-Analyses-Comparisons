############################################################################
############################################################################
###                                                                      ###
### Simulation for comparing Single-Cell RNA-Seq methods                 ###
###                                                                      ###
###                                                                      ###
###                                                                      ###
###                                                                      ###
###                                                                      ###
############################################################################
############################################################################

############################################################################
### Install needed packages
###

## Import 5 different packages for Single-Cell RNA-Seq differential analysis

## BASiCS

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("BASiCS", version = "3.8")

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

install.packages("MASS")
install.packages("pROC")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.8")

install.packages("foreign")
install.packages("gamlss.dist")
install.packages("FNN")

#library(BASiCS)
library(DEsingle)
library(Linnorm)
library(monocle)
library(MAST)
library(DESeq2)
library(pROC)

## simulation set up

num<-1
total<-20
#diffpro<-c(0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50)
diffpro<-c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50)
ngene<-1000
ndiff <- ngene*diffpro
samplesize <- 10
trtsize <- samplesize/2
label <- c(rep(1, trtsize), rep(0, samplesize-trtsize))

data(Islam2011)
example.data <- Islam2011[9:dim(Islam2011)[1], 1:92]
expMatrix <- example.data
#data(SEQC)
#expMatrix <- SEQC
#resample<-1000

## Generate simulation data

#simulateddata <- RnaXSim(expMatrix, distribution="NB", NumRep=trtsize, NumDiff = ndiff, NumFea = ngene)
set.seed(123)

alpha<-0.05

FDR.DEsingle <- sen.DEsingle <- spe.DEsingle <- acu.DEsingle <- auc.DEsingle <- matrix(NA, nrow=length(diffpro), ncol=total)
FDR.Linnorm <- sen.Linnorm <- spe.Linnorm <- acu.Linnorm <- auc.Linnorm <- matrix(NA, nrow=length(diffpro), ncol=total)
FDR.monocle <- sen.monocle <- spe.monocle <- acu.monocle <- auc.monocle <- matrix(NA, nrow=length(diffpro), ncol=total)
FDR.MAST <- sen.MAST <- spe.MAST <- acu.MAST <- auc.MAST <- matrix(NA, nrow=length(diffpro), ncol=total)
FDR.DESeq2 <- sen.DESeq2 <- spe.DESeq2 <- acu.DESeq2 <- auc.DESeq2<- matrix(NA, nrow=length(diffpro), ncol=total)
propzero <- matrix(NA, nrow=length(diffpro), ncol=total)

for (j in 1:length(diffpro))
{
  num<-1
  while (num<=total)
  {
    #simulateddata <- RnaXSim(expMatrix, distribution="NB", NumRep=3, NumDiff = ndiff, NumFea = ngene)
    simulateddata <- RnaXSim(expMatrix, distribution="Poisson", NumRep=trtsize, NumDiff = ndiff[j], NumFea = ngene)
    truev <- c(rep(0, ngene))
    truev[simulateddata$DiffList] <- 1
    
    ## BASiCS
    #trt <- newBASiCS_Data(simulateddata$data[, 1:trtsize], BatchInfo = label[1:trtsize])
    
    #rownames(group) <- paste("Gene",1:ngene)
    #colnames(group)<-c(paste("Control", 1:(samplesize-trtsize)), paste("Treatment", 1:trtsize))
    #design <- cbind(Grp1=1,Grp2vs1=c(rep(0, (samplesize-trtsize)), rep(1, trtsize)))
    #options(digit=3)
    
    propzero[j, num] <- sum(simulateddata$data == 0)/(ngene*samplesize)
  
    ## DEsingle
    
    group <- factor(c(rep(1, trtsize), rep(2, samplesize-trtsize)))
    results <- DEsingle(counts = simulateddata$data, group = group)
    DEsingle_sig <- rownames(results[which(results$pvalue.adj.FDR <= alpha),])
    DEsingle.s <- as.numeric(gsub("\\D", "", DEsingle_sig)) 
    #make signifnicant genes vector
    DEsingle.v <- c(rep(0, ngene))
    DEsingle.v[DEsingle.s] <- 1
    
    a <- table(DEsingle.v, truev)
    
    if (dim(a)[1] == 2) {
      FDR.DEsingle[j, num] <- a[2, 1]/(a[2,1]+a[2,2])
      sen.DEsingle[j, num] <- a[2, 2]/ndiff[j]
      spe.DEsingle[j, num] <- a[1, 1]/(a[1,1]+a[2,1])
      acu.DEsingle[j, num] <- (a[1, 1] + a[2, 2])/ngene
      roc_DEsingle <- roc(truev, DEsingle.v)
      auc.DEsingle[j, num] <- as.numeric(auc(roc_DEsingle))
    }
     else {
      FDR.DEsingle[j, num] <- 0
      sen.DEsingle[j, num] <- 0
      spe.DEsingle[j, num] <- 1
      acu.DEsingle[j, num] <- a[1, 1]/ngene
      roc_DEsingle <- roc(truev, DEsingle.v)
      auc.DEsingle[j, num] <- as.numeric(auc(roc_DEsingle))
     }

    ## Linnorm
    
    designmatrix <- model.matrix(~0+factor(group))
    colnames(designmatrix) <- c("group1", "group2")
    LinResults <- Linnorm.limma(simulateddata$data, designmatrix)
    Linnorm_sig <- rownames(LinResults[which(LinResults$adj.P.Val <= alpha),])
    Linnorm.s <- as.numeric(gsub("\\D", "", Linnorm_sig)) 
    #make signifnicant genes vector
    Linnorm.v <- c(rep(0, ngene))
    Linnorm.v[Linnorm.s] <- 1
    
    b <- table(Linnorm.v, truev)
    
    if (dim(b)[1] == 2){
      FDR.Linnorm[j, num] <- b[2, 1]/(b[2,1]+b[2,2])
      sen.Linnorm[j, num] <- b[2, 2]/ndiff[j]
      spe.Linnorm[j, num] <- b[1, 1]/(b[1,1]+b[2,1])
      acu.Linnorm[j, num] <- (b[1, 1] + b[2, 2])/ngene
      roc_Linnorm <- roc(truev, Linnorm.v)
      auc.Linnorm[j, num] <- as.numeric(auc(roc_Linnorm))
    }
     else {
       FDR.Linnorm[j, num] <- 0
       sen.Linnorm[j, num] <- 0
       spe.Linnorm[j, num] <- 1
       acu.Linnorm[j, num] <- b[1, 1]/ngene
       roc_Linnorm <- roc(truev, Linnorm.v)
       auc.Linnorm[j, num] <- as.numeric(auc(roc_Linnorm))
     }
    
    

    ## monocle
    
    label2 <- data.frame(c(rep(1, trtsize), rep(2, samplesize-trtsize))) 
    rownames(label2) <- colnames(simulateddata$data)
    colnames(label2) <- "Group"
    pd <- new("AnnotatedDataFrame", data = label2)
    mydata <- newCellDataSet(simulateddata$data, phenoData = pd)
    mydata2 <- detectGenes(mydata, min_expr = 0)
    expressed_genes <- row.names(subset(fData(mydata2), num_cells_expressed >= 0))
    #monoResults <- differentialGeneTest(mydata2[expressed_genes,], fullModelFormulaStr="~sm.ns(Group, df=1)", relative_expr = FALSE, cores=24)
    monoResults <- differentialGeneTest(mydata2[expressed_genes,], fullModelFormulaStr="~1 + Group", relative_expr = FALSE, cores=24)
    monocle.s <- which(monoResults$qval <= 0.05)
    monocle.v <- c(rep(0, ngene))
    monocle.v[monocle.s] <- 1
    
    c <- table(monocle.v, truev)
    
    if (dim(c)[1] == 2) {
      FDR.monocle[j, num] <- c[2, 1]/(c[2,1]+c[2,2])
      sen.monocle[j, num] <- c[2, 2]/ndiff[j]
      spe.monocle[j, num] <- c[1, 1]/(c[1,1]+c[2,1])
      acu.monocle[j, num] <- (c[1, 1] + c[2, 2])/ngene
      roc_monocle <- roc(truev, monocle.v)
      auc.monocle[j, num] <- as.numeric(auc(roc_monocle))
    }
    
    else {
      FDR.monocle[j, num] <- 0
      sen.monocle[j, num] <- 0
      spe.monocle[j, num] <- 1
      acu.monocle[j, num] <- c[1, 1]/ngene
      roc_monocle <- roc(truev, monocle.v)
      auc.monocle[j, num] <- as.numeric(auc(roc_monocle))
    }
    
    ## MAST
    label2 <- data.frame(c(rep(1, trtsize), rep(2, samplesize-trtsize))) 
    rownames(label2) <- colnames(simulateddata$data)
    colnames(label2) <- "Group"
    MAST_data <- FromMatrix(simulateddata$data, label2, check_sanity = FALSE)
    cond<-factor(colData(MAST_data)$Group)
    #cond<-relevel(cond,"Unstim")
    colData(MAST_data)$Group<-cond
    zlmCond <- zlm(~Group, MAST_data)
    summaryCond <- summary(zlmCond, doLRT='Group2')
    summaryDt <- summaryCond$datatable
    fcHurdle <- merge(summaryDt[contrast=='Group2' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast=='Group2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    
    fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    MAST_sig <- fcHurdle$primerid[which(fcHurdle$fdr <= 0.05)]
    MAST.s <- as.numeric(gsub("\\D", "", MAST_sig))
    MAST.v <- c(rep(0, ngene))
    MAST.v[MAST.s] <- 1
    
    d <- table(MAST.v, truev)
    
    if (dim (d)[1] == 2){
      FDR.MAST[j, num] <- d[2, 1]/(d[2,1]+d[2,2])
      sen.MAST[j, num] <- d[2, 2]/ndiff[j]
      spe.MAST[j, num] <- d[1, 1]/(d[1,1]+d[2,1])
      acu.MAST[j, num] <- (d[1, 1] + d[2, 2])/ngene
      roc_MAST <- roc(truev, MAST.v)
      auc.MAST[j, num] <- as.numeric(auc(roc_MAST))
    }
    
    else {
      
      FDR.MAST[j, num] <- 0
      sen.MAST[j, num] <- 0
      spe.MAST[j, num] <- 1
      acu.MAST[j, num] <- d[1, 1]/ngene
      roc_MAST <- roc(truev, MAST.v)
      auc.MAST[j, num] <- as.numeric(auc(roc_MAST))
    }
    
    ## DESeq2
    
    group.con <- data.frame(cbind(c(rep("Treatment", samplesize/2), rep("Control", samplesize/2)), c(rep("single-read", samplesize))))
    colnames(group.con) <- c("condition", "type")
    rownames(group.con) <- c(paste("Sample_", 1:samplesize, sep = ""))
    DESeq2.dds <- DESeqDataSetFromMatrix(countData = simulateddata$data, colData = group.con, design = ~condition)
    DESeq2.test <- DESeq(DESeq2.dds, quiet = TRUE)
    DESeq2.pval <- results(DESeq2.test)$padj
    DESeq2.s <- which(DESeq2.pval<= alpha) 
    #make signifnicant genes vector
    DESeq2.v <- c(rep(0, ngene))
    DESeq2.v[DESeq2.s] <- 1
    
    e <- table(DESeq2.v, truev)
    
    if (dim(e)[1] == 2){
      FDR.DESeq2[j, num] <- e[2, 1]/(e[2,1]+e[2,2])
      sen.DESeq2[j, num] <- e[2, 2]/ndiff[j]
      spe.DESeq2[j, num] <- e[1, 1]/(e[1,1]+e[2,1])
      acu.DESeq2[j, num] <- (e[1, 1] + e[2, 2])/ngene
      roc_DESeq2 <- roc(truev, DESeq2.v)
      auc.DESeq2[j, num] <- as.numeric(auc(roc_DESeq2))
    }
    else {
      FDR.DESeq2[j, num] <- 0
      sen.DESeq2[j, num] <- 0
      spe.DESeq2[j, num] <- 1
      acu.DESeq2[j, num] <- e[1, 1]/ngene
      roc_DESeq2 <- roc(truev, DESeq2.v)
      auc.DESeq2[j, num] <- as.numeric(auc(roc_DESeq2))
    }
    
     
    num<-num+1
  }
}

FDRmean.DEsingle <- senmean.DEsingle <- spemean.DEsingle <- acumean.DEsingle <-  aucmean.DEsingle <- rep(NA, length(diffpro))
FDRmean.Linnorm <- senmean.Linnorm <- spemean.Linnorm <- acumean.Linnorm <-  aucmean.Linnorm <- rep(NA, length(diffpro))
FDRmean.monocle <- senmean.monocle <- spemean.monocle <- acumean.monocle <-  aucmean.monocle <- rep(NA, length(diffpro))
FDRmean.MAST <- senmean.MAST <- spemean.MAST <- acumean.MAST <- aucmean.MAST <- rep(NA, length(diffpro))
FDRmean.DESeq2 <- senmean.DESeq2 <- spemean.DESeq2 <- acumean.DESeq2 <- aucmean.DESeq2 <- rep(NA, length(diffpro))
propzeromean <- rep(NA, length(diffpro))

for (i in 1:length(diffpro))
{
  FDRmean.DEsingle[i]<-mean(FDR.DEsingle[i, ], na.rm=TRUE)
  senmean.DEsingle[i]<-mean(sen.DEsingle[i, ], na.rm=TRUE)
  spemean.DEsingle[i]<-mean(spe.DEsingle[i, ], na.rm=TRUE)
  acumean.DEsingle[i]<-mean(acu.DEsingle[i, ], na.rm=TRUE)
  aucmean.DEsingle[i]<-mean(auc.DEsingle[i, ], na.rm=TRUE)
  
  FDRmean.Linnorm[i]<-mean(FDR.Linnorm[i, ], na.rm=TRUE)
  senmean.Linnorm[i]<-mean(sen.Linnorm[i, ], na.rm=TRUE)
  spemean.Linnorm[i]<-mean(spe.Linnorm[i, ], na.rm=TRUE)
  acumean.Linnorm[i]<-mean(acu.Linnorm[i, ], na.rm=TRUE)
  aucmean.Linnorm[i]<-mean(auc.Linnorm[i, ], na.rm=TRUE)
  
  FDRmean.monocle[i]<-mean(FDR.monocle[i, ], na.rm=TRUE)
  senmean.monocle[i]<-mean(sen.monocle[i, ], na.rm=TRUE)
  spemean.monocle[i]<-mean(spe.monocle[i, ], na.rm=TRUE)
  acumean.monocle[i]<-mean(acu.monocle[i, ], na.rm=TRUE)
  aucmean.monocle[i]<-mean(auc.monocle[i, ], na.rm=TRUE)
  
  FDRmean.MAST[i]<-mean(FDR.MAST[i, ], na.rm=TRUE)
  senmean.MAST[i]<-mean(sen.MAST[i, ], na.rm=TRUE)
  spemean.MAST[i]<-mean(spe.MAST[i, ], na.rm=TRUE)
  acumean.MAST[i]<-mean(acu.MAST[i, ], na.rm=TRUE)
  aucmean.MAST[i]<-mean(auc.MAST[i, ], na.rm=TRUE)
  
  FDRmean.DESeq2[i]<-mean(FDR.DESeq2[i, ], na.rm=TRUE)
  senmean.DESeq2[i]<-mean(sen.DESeq2[i, ], na.rm=TRUE)
  spemean.DESeq2[i]<-mean(spe.DESeq2[i, ], na.rm=TRUE)
  acumean.DESeq2[i]<-mean(acu.DESeq2[i, ], na.rm=TRUE)
  aucmean.DESeq2[i]<-mean(auc.DESeq2[i, ], na.rm=TRUE)
  
  propzeromean[i]<-mean(propzero[i, ], na.rm=TRUE)
}

print(cbind(propzeromean, FDRmean.DEsingle, FDRmean.Linnorm, FDRmean.monocle, FDRmean.MAST, FDRmean.DESeq2))
print(cbind(propzeromean, senmean.DEsingle, senmean.Linnorm, senmean.monocle, senmean.MAST, senmean.DESeq2))
print(cbind(propzeromean, spemean.DEsingle, spemean.Linnorm, spemean.monocle, spemean.MAST, spemean.DESeq2))
print(cbind(propzeromean, acumean.DEsingle, acumean.Linnorm, acumean.monocle, acumean.MAST, acumean.DESeq2))
print(cbind(propzeromean, aucmean.DEsingle, aucmean.Linnorm, aucmean.monocle, aucmean.MAST, aucmean.DESeq2))


meanFDR<-cbind(propzeromean, FDRmean.DEsingle, FDRmean.Linnorm, FDRmean.monocle, FDRmean.MAST, FDRmean.DESeq2)
meansen<-cbind(propzeromean, senmean.DEsingle, senmean.Linnorm, senmean.monocle, senmean.MAST, senmean.DESeq2)
meanspe<-cbind(propzeromean, spemean.DEsingle, spemean.Linnorm, spemean.monocle, spemean.MAST, spemean.DESeq2)
meanacu<-cbind(propzeromean, acumean.DEsingle, acumean.Linnorm, acumean.monocle, acumean.MAST, acumean.DESeq2)
meanauc<-cbind(propzeromean, aucmean.DEsingle, aucmean.Linnorm, aucmean.monocle, aucmean.MAST, aucmean.DESeq2)


write.csv(meanFDR, file="D:/Box Sync/Conferences/Huston_2018/meanFDR_DESeq2_n_5_propzero_Islam2011_Poisson.csv")
write.csv(meansen, file="D:/Box Sync/Conferences/Huston_2018/meansen_DESeq2_n_5_propzero_Islam2011_Poisson.csv")
write.csv(meanspe, file="D:/Box Sync/Conferences/Huston_2018/meanspe_DESeq2_n_5_propzero_Islam2011_Poisson.csv")
write.csv(meanacu, file="D:/Box Sync/Conferences/Huston_2018/meanacu_DESeq2_n_5_propzero_Islam2011_Poisson.csv")
write.csv(meanauc, file="D:/Box Sync/Conferences/Huston_2018/meanauc_DESeq2_n_5_propzero_Islam2011_Poisson.csv")


colnames(meanFDR) <- colnames(meansen) <- colnames(meanspe) <- colnames(meanacu) <- colnames(meanauc) <- c("Proportion of Zero", "DEsingle", "Linnorm", "monocle", "MAST", "DESeq2")

#pdf(file="RNA_Seq_comparison_nel_n_24.pdf", width=4, height=4)
#par(mfrow=c(2,2), cex=0.7, mar=c(5, 5, 2.5, 2.5))
par(mfrow=c(2, 2), mar=c(8,6,2,2)+0.1, mgp=c(4, 1, 0), cex=0.7)
boxplot(meanFDR[, -1], col=c("blue","red", "darkgreen", "cyan", "orange"), xlab=NULL, ylab="FDR", las=2)
boxplot(meansen[, -1], col=c("blue","red", "darkgreen", "cyan", "orange"), xlab=NULL, ylab="Sensitivity", las=2)
boxplot(meanspe[, -1], col=c("blue","red", "darkgreen", "cyan", "orange"), xlab=NULL, ylab="Specificity", las=2)
boxplot(meanacu[, -1], col=c("blue","red", "darkgreen", "cyan", "orange"), xlab=NULL, ylab="accuracy", las=2)
dev.off()

par(mfrow=c(1, 1), mar=c(8,6,2,2)+0.1, mgp=c(4, 1, 0), cex=1.0)
boxplot(meanauc[, -1], col=c("blue","red", "darkgreen", "cyan", "orange"), xlab=NULL, ylab="AUC", las=2)
dev.off()



#pdf(file="RNA_Seq_comparison_nel_n_24_lineplot.pdf", width=4, height=4)
par(mfrow=c(2,2), cex=0.7, mgp=c(3, 1, 0), mar=c(4.5, 4.5, 1.5, 0.7))
g_range<-range(meanFDR, na.rm=TRUE)
plot(diffpro, FDRmean.DEsingle, type="l", lty=1, lwd=2, col="blue", ylim=g_range, ann=FALSE)
box()
lines(diffpro, FDRmean.Linnorm, type="l", lty=2, lwd=2, col="red")
lines(diffpro, FDRmean.monocle, type="l", lty=3, lwd=2, col="darkgreen")
lines(diffpro, FDRmean.MAST, type="l", lty=4, lwd=2, col="black")
lines(diffpro, FDRmean.DESeq2, type="l", lty=5, lwd=2, col="orange")
title(xlab=expression(pi[1]), ylab="FDR")
legend("topright", c("DEsingle","Linnorm", "monocle", "MAST", "DESeq2"), cex=0.8,
       col=c("blue","red", "darkgreen", "black", "orange"), lwd=2, lty=c(1:5, 1:5), bty="n")
#text(1, g_range[2], "A")


g_range<-range(meansen, na.rm=TRUE)
plot(diffpro, senmean.DEsingle, type="l", lty=1, lwd=2, col="blue", ylim=g_range, ann=FALSE)
box()
lines(diffpro, senmean.Linnorm, type="l", lty=2, lwd=2, col="red")
lines(diffpro, senmean.monocle, type="l", lty=3, lwd=2, col="darkgreen")
lines(diffpro, senmean.MAST, type="l", lty=4, lwd=2, col="black")
lines(diffpro, senmean.DESeq2, type="l", lty=5, lwd=2, col="orange")
title(xlab=expression(pi[1]), ylab="Sensitivity")
#legend(0.6, g_range[2], c("DEsingle","Linnorm", "monocle", "MAST"), cex=0.8,
#       col=c("blue","red", "darkgreen", "black"), lwd=2, lty=c(1:5, 1:5), bty="n")

g_range<-range(meanspe, na.rm=TRUE)
plot(diffpro, spemean.DEsingle, type="l", lty=1, lwd=2, col="blue", ylim=g_range, ann=FALSE)
box()
lines(diffpro, spemean.Linnorm, type="l", lty=2, lwd=2, col="red")
lines(diffpro, spemean.monocle, type="l", lty=3, lwd=2, col="darkgreen")
lines(diffpro, spemean.MAST, type="l", lty=4, lwd=2, col="black")
lines(diffpro, spemean.DESeq2, type="l", lty=5, lwd=2, col="orange")
title(xlab=expression(pi[1]), ylab="Specificity")
#legend(0.6, g_range[2], c("DEsingle","Linnorm", "monocle", "MAST"), cex=0.8,
#       col=c("blue","red", "darkgreen", "black"), lwd=2, lty=c(1:4, 1:4), bty="n")


g_range<-range(meanacu, na.rm=TRUE)
plot(diffpro, acumean.DEsingle, type="l", lty=1, lwd=2, col="blue", ylim=g_range, ann=FALSE)
box()
lines(diffpro, acumean.Linnorm, type="l", lty=2, lwd=2, col="red")
lines(diffpro, acumean.monocle, type="l", lty=3, lwd=2, col="darkgreen")
lines(diffpro, acumean.MAST, type="l", lty=4, lwd=2, col="black")
lines(diffpro, acumean.DESeq2, type="l", lty=5, lwd=2, col="orange")
title(xlab=expression(pi[1]), ylab="Accuracy")
#legend(0.6, g_range[2], c("DEsingle","Linnorm", "monocle", "MAST"), cex=0.8,
#       col=c("blue","red", "darkgreen", "black"), lwd=2, lty=c(1:4, 1:4), bty="n")

dev.off()


par(mfrow=c(1, 1), cex=1.0, mgp=c(3, 1, 0), mar=c(4.5, 4.5, 1.5, 0.7))
g_range<-range(meanauc, na.rm=TRUE)
plot(diffpro, aucmean.DEsingle, type="l", lty=1, lwd=2, col="blue", ylim=g_range, ann=FALSE)
box()
lines(diffpro, aucmean.Linnorm, type="l", lty=2, lwd=2, col="red")
lines(diffpro, aucmean.monocle, type="l", lty=3, lwd=2, col="darkgreen")
lines(diffpro, aucmean.MAST, type="l", lty=4, lwd=2, col="black")
lines(diffpro, aucmean.DESeq2, type="l", lty=5, lwd=2, col="orange")
title(xlab=expression(pi[1]), ylab="AUC")
legend("right", c("DEsingle","Linnorm", "monocle", "MAST", "DESeq2"), cex=1.0,
       col=c("blue","red", "darkgreen", "black", "orange"), lwd=2, lty=c(1:4, 1:4), bty="n")

dev.off()


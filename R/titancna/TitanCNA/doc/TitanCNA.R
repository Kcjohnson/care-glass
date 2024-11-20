### R code from vignette source 'TitanCNA.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: TitanCNA.Rnw:76-81
###################################################
library(TitanCNA)
infile <- system.file("extdata", "test_alleleCounts_chr2.txt", 
                      package = "TitanCNA")
data <- loadAlleleCounts(infile, genomeStyle = "NCBI")
names(data)


###################################################
### code chunk number 2: TitanCNA.Rnw:90-96
###################################################
tumWig <- system.file("extdata", "test_tum_chr2.wig", package = "TitanCNA")
normWig <- system.file("extdata", "test_norm_chr2.wig", package = "TitanCNA")
gc <- system.file("extdata", "gc_chr2.wig", package = "TitanCNA")
map <- system.file("extdata", "map_chr2.wig", package = "TitanCNA")
cnData <- correctReadDepth(tumWig, normWig, gc, map, genomeStyle = "NCBI")
head(cnData)


###################################################
### code chunk number 3: TitanCNA.Rnw:103-105
###################################################
logR <- getPositionOverlap(data$chr, data$posn, cnData)
data$logR <- log(2^logR)  #transform the log ratio to natural logs


###################################################
### code chunk number 4: TitanCNA.Rnw:115-117
###################################################
data <- filterData(data, c(1:22, "X", "Y"), minDepth = 10, maxDepth = 200, 
                   positionList = NULL, centromere = NULL, centromere.flankLength = 10000)


###################################################
### code chunk number 5: TitanCNA.Rnw:125-131
###################################################
numClusters <- 2
params <- loadDefaultParameters(copyNumber = 5, 
                                numberClonalClusters = numClusters,
                                symmetric = TRUE, hetBaselineSkew = 0, 
                                alleleEmissionModel = "binomial", data = data)
params


###################################################
### code chunk number 6: TitanCNA.Rnw:145-147 (eval = FALSE)
###################################################
## params$ploidyParams$phi_0 <- 2 # for diploid or
## params$ploidyParams$phi_0 <- 4 # for tetraploid/ployploid


###################################################
### code chunk number 7: TitanCNA.Rnw:152-155
###################################################
K <- length(params$genotypeParams$alphaKHyper)
params$genotypeParams$alphaKHyper <- rep(500, K)
params$ploidyParams$phi_0 <- 1.5


###################################################
### code chunk number 8: TitanCNA.Rnw:163-165 (eval = FALSE)
###################################################
## library(doMC)
## registerDoMC(cores = 4) #use 4 cores on a single machine


###################################################
### code chunk number 9: TitanCNA.Rnw:180-187
###################################################
convergeParams <- runEMclonalCN(data, params, 
                                maxiter = 3, maxiterUpdate = 50, 
                                useOutlierState = FALSE, txnExpLen = 1e15, 
                                txnZstrength = 5e5, 
                                normalEstimateMethod = "map", 
                                estimateS = TRUE, estimatePloidy = TRUE)
names(convergeParams)


###################################################
### code chunk number 10: TitanCNA.Rnw:207-209
###################################################
optimalPath <- viterbiClonalCN(data, convergeParams)
head(optimalPath)


###################################################
### code chunk number 11: TitanCNA.Rnw:220-230
###################################################
results <- outputTitanResults(data, convergeParams, optimalPath,
                              filename = NULL, posteriorProbs = FALSE,
                              subcloneProfiles = TRUE, correctResults = TRUE, 
                              proportionThreshold = 0.05, 
                              proportionThresholdClonal = 0.05,
                              is.haplotypeData = FALSE)
names(results)
head(results$corrResults) ## corrected results
convergeParams <- results$convergeParam ## use corrected parameters
results <- results$corrResults ## use corrected results


###################################################
### code chunk number 12: eval
###################################################
segs <- outputTitanSegments(results, id = "test", convergeParams, 
                            filename = NULL, igvfilename = NULL)
head(segs)


###################################################
### code chunk number 13: eval
###################################################
# get the estimated tumor ploidy
ploidy <- tail(convergeParams$phi, 1)
# get the estimated normal
normal <- tail(convergeParams$n, 1)
# apply tumor purity and ploidy correction 
corrIntCN.results <- correctIntegerCN(results, segs, 1 - normal, ploidy, 
                                      maxCNtoCorrect.autosomes = 8, 
                                      maxCNtoCorrect.X = NULL, 
                                      correctHOMD = FALSE, minPurityToCorrect = 0.2, 
                                      gender = "female", chrs = 2)
head(corrIntCN.results$segs)
# re-assign to results and segs objects
results <- corrIntCN.results$cn
segs <- corrIntCN.results$segs


###################################################
### code chunk number 14: TitanCNA.Rnw:340-342 (eval = FALSE)
###################################################
## outparam <- paste("test_cluster02_params.txt", sep = "")
## outputModelParameters(convergeParams, results, outparam, S_Dbw.scale = 1)


###################################################
### code chunk number 15: TitanCNA.Rnw:362-363 (eval = FALSE)
###################################################
## outputModelParameters(convergeParams, results, outparam, S_Dbw.scale = 10)


###################################################
### code chunk number 16: TitanCNA.Rnw:373-379
###################################################
ploidy <- tail(convergeParams$phi, 1)
ploidy
normal <- tail(convergeParams$n, 1)
normal
plotCNlogRByChr(results, segs = segs, chr = 2, ploidy = ploidy, normal = normal, 
		ylim = c(-2, 2), cex = 0.25, xlab = "", main = "Chr 2")


###################################################
### code chunk number 17: TitanCNA.Rnw:392-394
###################################################
plotAllelicRatio(results, chr = 2, ylim = c(0, 1), cex = 0.25, 
                 xlab = "", main = "Chr 2")


###################################################
### code chunk number 18: TitanCNA.Rnw:408-413
###################################################
norm <- tail(convergeParams$n, 1) 
norm # estimated normal contamination
1 - convergeParams$s[, ncol(convergeParams$s)] # estimated cellular prevalence 
plotClonalFrequency(results, chr = 2, normal = norm, ylim = c(0, 1), 
                    cex = 0.25, xlab = "", main = "Chr 2")


###################################################
### code chunk number 19: TitanCNA.Rnw:422-423
###################################################
plotSubcloneProfiles(results, chr = 2, cex = 1, spacing = 2, main = "Chr 2")


###################################################
### code chunk number 20: TitanCNA.Rnw:431-433
###################################################
plotSegmentMedians(segs, chr=2, resultType = "LogRatio", plotType = "CopyNumber", 
                   plot.new = TRUE, ylim = c(0, 4), main="Chr 2")


###################################################
### code chunk number 21: TitanCNA.Rnw:440-441
###################################################
toLatex(sessionInfo())



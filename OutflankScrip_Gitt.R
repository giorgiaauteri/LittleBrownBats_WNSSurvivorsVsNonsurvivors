#OUtflank
library(OutFLANK)
library(vcfR)
##############################

#Read in outflank data
##Converted from snp.vcf to Eigensoft in PGD spyder, then added " " to separate snps
o <- read.table("D:/Mirror/Mirror/Projects/PhD Thesis/SurvMortWNS/Data/SurvMort21/Pops_SMAll/OutflankOutput.geno",
                sep = "")
levels(factor(as.numeric(unlist(o))))
ot <- t(o)

#Read in list of sample IDs
p <- read.table("D:/Mirror/Mirror/Projects/PhD Thesis/SurvMortWNS/Data/SurvMort21/Pops_SMAll/popnamesOutflank.txt",
                sep = "")

#Read in list of loci names
l <- read.csv("D:/Mirror/Mirror/Projects/PhD Thesis/SurvMortWNS/Data/SurvMort21/Pops_SMAll/lociNames.csv",
              header = F)

#Convert to Outflank format
out <- MakeDiploidFSTMat(SNPmat=ot, locusNames=l$V2, popNames=p$V2)

#corrected vs uncorrectted fst
plot(out$FST, out$FSTNoCorr, 
     xlim=c(-0.01,0.3), ylim=c(-0.01,0.3),
     pch=20)
abline(0,1)

#Check if fst distribution violates chi squared
plot(out$He, out$FSTNoCorr, pch=20, col="grey") #Above He 0.1, seems mostly fine
HeCutoff <- 0.1 #Adjust to check what the He cutoff should be
abline(v = HeCutoff, col = "blue", lty= 2, lwd= 2) 

#run outflank
outRes <- OutFLANK(out,  
                   LeftTrimFraction =  0.3, #adjust to get best possible fit to histogram (below)
                   RightTrimFraction = 0.05, #adjust to get best possible fit to histogram (below)
                   NumberOfSamples = 2,
                   Hmin = HeCutoff,
                   qthreshold = 0.05) 

#Fits without outliers
OutFLANKResultsPlotter(outRes, withOutliers = FALSE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.01, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
OutFLANKResultsPlotter(outRes, withOutliers = FALSE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.01, Zoom =
                         TRUE, RightZoomFraction = 0.05, titletext = NULL)

#Fits with outliers
OutFLANKResultsPlotter(outRes, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.01, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
OutFLANKResultsPlotter(outRes, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.01, Zoom =
                         TRUE, RightZoomFraction = 0.05, titletext = NULL)

#Also check the P value histogram
hist(outRes$results$pvaluesRightTail)

outRes$numberHighFstOutliers
sum(outRes$results$qvalues<=0.01, na.rm=TRUE) #More conservative q-correction with 0.01 alpha
sum(outRes$results$qvalues <=0.05, na.rm=TRUE)

#########################################################
## Using estimated neutral mean FST and df to calculate P-values for all loci

#Now that we've estimated neutral mean FST and df to a quasi-independent set of SNPs, we can go back and calculate P-values and q-values for all the loci in our dataset. 
P1 <- pOutlierFinderChiSqNoCorr(out, Fstbar = outRes$FSTNoCorrbar, 
                                dfInferred = outRes$dfInferred, qthreshold = 0.05, Hmin=HeCutoff)

#Plot He vs Fst
plot(P1$He[P1$He >= 0.1], P1$FSTNoCorr[P1$He >= 0.1], pch=19, col=rgb(0,0,0,0.1))
points(P1$He[P1$OutlierFlag==TRUE], P1$FSTNoCorr[P1$OutlierFlag==TRUE], col="blue") #Significant above 0.05
points(P1$He[P1$qvalues <= 0.01], P1$FSTNoCorr[P1$qvalues <= 0.01], pch = 21, cex = 2, col="red")

###################################################

#write file of candidate loci
selectedLoci <- na.omit(P1[P1$qvalues <= 0.05,])
write.csv(selectedLoci, "D:/Mirror/Mirror/Projects/PhD Thesis/SurvMortWNS/Data/SurvMort21/Pops_SMAll/SelectedLoci_OutFlank.csv")

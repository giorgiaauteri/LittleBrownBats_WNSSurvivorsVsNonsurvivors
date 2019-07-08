#Bootstrapping Fst 
library("diveRsity")

set.seed(1234)

#set working directory
setwd("/scratch/lsa_flux/gauteri/Myotis1_SurvMort/ProcRad/Mapped21/Bootstrap_SMAll")

#calculate fst with CI bootstrap
calcs <- diffCalc(infile = "populations.snps.genepop",
         outfile = "BootstrapResults",
         fst = TRUE,
         bs_locus = TRUE,
         ci_type = "loci",
         boots = 1000,
         para = TRUE)

#Assign interpretable loci names
##Read in a file with loci names in order
fdat <- read.table("populations.fst_mortality-survivor.tsv",
                   header = FALSE, 
                   sep = "\t")
colnames(fdat) <- c('LocusID',	'Pop1',	'Pop2',	'CHROM',	'POS', 'Column',	'OverallPi',	'AMOVAFst',	'FishersP',	'OddsRatio',	'CILow',	'CIHigh.All',	'LOD',	'CorrectedAMOVAFst', 'SmoothedAMOVAFst',	'SmoothedAMOVAFst P-value',	'WindowSNPCount')
##Place information with corresponding rows
calcs$bs_locus$Fst$locusScaffold <- NA #Create empty columns
calcs$bs_locus$Fst$locusPosition <- NA
calcs$bs_locus$Fst$locusScaffold <- fdat$CHROM
calcs$bs_locus$Fst$locusPosition <- fdat$POS

############################
#subset outliers

#Which loci have a minimum CI that doesn't drop below the threshold of significance
##cutoff is mean + (x * 1SD)
meanFST <- mean(calcs$bs_locus$Fst$actual)
SD.1 <- sd(calcs$bs_locus$Fst$actual); SD.1 #What is one SD?
Sds <- 5 #How many standard deviations above the mean should be used for the cutoff?

cutoff <- meanFST + (Sds * SD.1) #Fst threshold
cutoff #what was the fst cutoff?

outliers <- na.omit(calcs$bs_locus$Fst[calcs$bs_locus$Fst$lower >= cutoff,])

#######################3
#Write output files (fulll output and outlier only)

write.csv(outliers, 
          "bootstrapOutliersOver5SD.csv",
          row.names = FALSE)

write.csv(calcs$bs_locus$Fst, 
          "bootstrapFST.csv",
          row.names = FALSE)

write.csv(calcs$bs_locus$gst, 
          "bootstrap_gst.csv",
          row.names = FALSE)

write.csv(calcs$bs_locus$GGst, 
          "bootstrap_GGST.csv",
          row.names = FALSE)

write.csv(calcs$bs_locus$D, 
          "bootstrap_D.csv",
          row.names = FALSE)

write.csv(calcs$bs_locus$Fis, 
          "bootstrap_Fis.csv",
          row.names = FALSE)

write.csv(calcs$bs_locus$Fit, 
          "bootstrap_Fit.csv",
          row.names = FALSE)

#Read in outliers identified via different methods and then see which SNPs are the same

#Read in outlier SNPs from STACKS
sel <- read.csv("D:/Mirror/Mirror/Projects/PhD Thesis/SurvMortWNS/Data/SurvMort21/Pops_SMAll/SelectedLoci.csv")
sel$CHROM <- as.character(sel$CHROM)
sel$POS <- as.numeric(sel$POS)

#Read in outliers from diveRsity
sel2 <- read.csv("D:/Mirror/Mirror/Projects/PhD Thesis/SurvMortWNS/Data/SurvMort21/Bootstrap/bootstrap3/bootstrapOutliersOver5SD.csv")
colnames(sel2) <- c("Locus", "diversityFST", "lower", "upper", "CHROM", "POS")
sel2$CHROM <- as.character(sel2$CHROM)
sel2$POS <- as.numeric(sel2$POS)

#Read in outliers from outFLANK
sel3 <- read.csv("D:/Mirror/Mirror/Projects/PhD Thesis/SurvMortWNS/Data/SurvMort21/Pops_SMAll/SelectedLoci_OutFlank.csv")
sel3 <- sel3[ ,c(2:5, 8, 11:15)]
colnames(sel3) <- c("CHROM", "POS", "He_OF", "FST_OF", "FSTnoCorr_OF", 
                    "meanAlleleFreq_OF", "pvalues_OF", "rightTail_OF",
                    "qvalues_OF", "OutlierFLag_OF")
sel3$CHROM <- as.character(sel3$CHROM)
sel3$POS <- as.numeric(sel3$POS)

##############################################################
#Find Sites identified by all three methods

#Combine results
final <- merge(x = sel, y = sel2, by = c("CHROM", "POS"), all = TRUE)
final <- merge(x = final, y = sel3, by = c("CHROM", "POS"), all = TRUE)
generalCandidates <- final #keep a list of all SNPs that were identified via any method
final <- na.omit(final)

#write to csv files
##Top candidates 
write.csv(final,  
           quote = FALSE,
          row.names = FALSE,
           file="D:/Mirror/Mirror/Projects/PhD Thesis/SurvMortWNS/Data/SurvMort21/Pops_SMAll/SelectedLoci_PreIDs.csv")
##General candidates 
write.csv(generalCandidates,  
          quote = FALSE,
          row.names = FALSE,
          file="D:/Mirror/Mirror/Projects/PhD Thesis/SurvMortWNS/Data/SurvMort21/Pops_SMAll/LociOfInterest_PreIDs.csv")

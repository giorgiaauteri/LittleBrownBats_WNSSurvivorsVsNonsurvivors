#Giorgia Auteri
#Subset Structure file to single SNP
#This list of single-SNPs can then be used to subset other datasets

##In this dataset column names look like this: ###_##
##the first number is the fragment number and the second is the base-pair position of the SNP

#Read structure file in as data frame
dat <- read.table("populations_33Inds_NoHead.stru",
                  sep="\t", header=TRUE) 
names(dat)[1] <- "sample"
names(dat)[2] <- "group"

#Remove multiple loci
##Create a dataframe of locus/SNP info
fullPosition <- names(dat)[3:(length(names(dat)))]  #get loc-bp info
bps <- gsub(".*_","",fullPosition)#Isolate bp by removing all before and up to "_"
fragments <- gsub("_.*","",fullPosition) #Remove bp info from fragments (just loci left)
allLocs <- as.data.frame(cbind(fullPosition, fragments, bps)) #Combine into data frame
allLocs$status <- NA #Create new column for duplicate loci status

##use duplicated.random function (found here: https://stat.ethz.ch/pipermail/r-help/2010-June/241244.html)
set.seed(28) #For reproducibility; My age when I started this project forever ago
duplicated.random = function(x, incomparables = FALSE, ...)
{ 
  if ( is.vector(x) )
  {
    permutation = sample(length(x))
    x.perm      = x[permutation]
    result.perm = duplicated(x.perm, incomparables, ...)
    result      = result.perm[order(permutation)]
    return(result)
  }
  else if ( is.matrix(x) )
  {
    permutation = sample(nrow(x))
    x.perm      = x[permutation,]
    result.perm = duplicated(x.perm, incomparables, ...)
    result      = result.perm[order(permutation)]
    return(result)
  }
  else
  {
    stop(paste("duplicated.random() only supports vectors",
               "matrices for now."))
  }
}

allLocs$status <- duplicated.random(fragments) #run the function on the data

length(allLocs$status[allLocs$status == FALSE]) #How many loci are not duplicates?
expectedSNPs <- 14345; expectedSNPs #!!! enter number of unique loci from your dataset

#Almost there... are there any SNPs in particular you want to be sure to keep?

prioritySNPs <- c("X399602_82", "X224520_42", "X224875_117", "X533095_71", "X165245_21", "X84078_35", "X399059_46",
                  "X533077_74", "X471816_80") #!!! change to reflect any SNPs that you want to prioritize
#i.e. if a priority SNP was assigned to be discarded, its place in the final SNP list
#will be restored and the other SNP(s) on that fragment will be switched to "discard"

##see how many actually have duplicates
allLocs[allLocs$fullPosition %in% prioritySNPs,]

lociToAdjust <- gsub("_.*","",prioritySNPs)

for (i in 1:(length(lociToAdjust))){ #For each locus/SNP
  tempdf <- allLocs[allLocs$fragments == lociToAdjust[i],] #subset just SNPs for that loci
  #See if it has a duplicate or not
  if ((length(tempdf$fullPosition))== 1){
    next #go to next loop iteration
  } 
  else { #There must be more than one SNP for that locus, so...
  #Turn everything in the main DF for that locus to TRUE (i.e. marked as duplicate)
  allLocs$status[allLocs$fragments %in% lociToAdjust[i]] <- TRUE
  #Edjust the priority SNP to keep (i.e. mark as FALSE; not a duplicate to discard)
  allLocs$status[allLocs$fullPosition == prioritySNPs[i]] <- FALSE
  } 
}

#Double check that you still have the right number of SNPs
length(allLocs$status[allLocs$status == FALSE]) == expectedSNPs #Should be TRUE

####################
#Write output files (after formating)

#Reduce the allLocs file to just what's needed for the new Structure file
keepSNPs <- c("sample", "group", as.character(allLocs$fullPosition[allLocs$status == FALSE])) #list of headings to keep
dat2 <- dat[,keepSNPs]
names(dat2) <- gsub("X", "", names(dat2))#remove the X in front of the name
#remove the "sample" and "group" headings manually once file is written

head(dat2[1:20, 1:20]) #Check
#write file
write.table(dat2, quote = FALSE, row.names = FALSE, 
            file = "populations_33Inds_NoHead_SingleSNP.stru")


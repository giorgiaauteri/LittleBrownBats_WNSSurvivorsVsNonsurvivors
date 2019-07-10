
# title: "IDSNPStoRef"
#author: "Giorgia Auteri"
#date: "June 27, 2018"
#output: pdf_document


########################################################
#Bring in data

#Read in reference genome info
## If using file from Auteri Github repository, you will need to unzip this input file first
ref <- read.table("ref_Myoluc2.0_top_level_reduced.txt", 
                  fill = TRUE, 
                  sep ="\t",
                  stringsAsFactors = FALSE)

#Add column names & format
colnames(ref) <- c("alignment", "prog", "type", "start", "end", c(6:(length(ref))))
ref$alignment <- as.character(ref$alignment)
ref$start <- as.numeric(ref$start)
ref$end <- as.numeric(ref$end)

#Remove uninformative rows 
##(those that are just reporting the scaffold length, not info about a particular segment)
ref <- na.omit(ref[ref$prog != "RefSeq",])

#Read in data on selected SNPs
#switch to LociOfInterest_PreIDs.csv to lookup info for SNPs identified via even just one method
sel <- read.csv("SelectedLoci_PreIDs.csv")
sel$CHROM <- as.character(sel$CHROM)
sel$POS <- as.numeric(sel$POS)

#####################################################################################################
#Match up SNPs against reference genome

#SNP loop to id what annotated loci SNPs of interest fall in (and if they are not in annotated loci,
#identify the closest anotated loci in each direction)

numSelected <- length(sel$LocusID) #set the number of selected loci
dat <- c("SNP", "Location", "Distance", colnames(ref)) #initiate what will become the data frame

for (i in 1:numSelected){
  NthAss <- sel$CHROM[i] #take the first identified scaffold
  NthPos <- sel$POS[i] #get the nth identified position
  selLociN <- ref[ref$alignment == NthAss & ref$start <= NthPos & ref$end >= NthPos,] #return all
  if (nrow(selLociN) > 0){
    for (ii in 1:nrow(selLociN)){
      toAdd <- unname(c(NthPos, "Within", 0, selLociN[ii,])) #new data to add
      dat <- as.data.frame(rbind(dat, toAdd)) #Bind to growing data frame
    }
  }
  else { #Find closest annotated region in each direction & add to data frame
    
    #Create a vector of start and end positions for the appropriate alignment
    sides <- ref[ref$alignment == NthAss, c("start","end")]
    
    #Left
    #A vevtor of the differences between the SNP and the ending locus position
    combos <- outer(NthPos, sides$end,FUN="-")
    #The closest loci to the left
    minL <- min(combos[combos > 0]) 
    
    #Right
    combos <- outer(sides$start, NthPos, FUN="-")
    minR <- min(combos[combos > 0])
    
    #Which is the smallest difference (to the Left or Right)?
    closest <- min(c(minR, minL))
    secondclosest <- max(c(minR, minL))
    
    ##Subset the closest gene from the reference genome 
    selLociClose <- ref[ref$alignment == NthAss & (ref$end == (NthPos-closest) | 
                                                     (ref$start ==   (closest + NthPos))),] 
    #return all rows and all columns (each column explicity identified)
    
    ##Second closest
    selLociClose2 <- ref[ref$alignment == NthAss & ((ref$start == (minR + NthPos))),] 
    #return all rows and all columns (each column explicity identified)
    
    for (ii in nrow(selLociClose)){
      toAdd <- unname(c(NthPos, "Outside_Closest", closest, selLociClose[ii,]))#Find 
      dat <- as.data.frame(rbind(dat, toAdd)) #Bind to growing data frame
    }
    for (ii in nrow(selLociClose2)){
      toAdd <- unname(c(NthPos, "Outside_SecondClose", secondclosest, selLociClose2[ii,])) #Find
      dat <- as.data.frame(rbind(dat, toAdd)) #Bind to growing data frame
    }
  }
}


#####################################################################################
#output results file to .csv

#Fix the data frame
colnames(dat) <- as.character(unlist(dat[1,])) # the first row will be the header
dat <- dat[-1, ]          # now remove the first row.

dat <- apply(dat,2,as.character)

write.csv(dat, 
          quote = FALSE,
          file="SelectedLociRefIDs_ThreeUnion.csv")#Write to csv

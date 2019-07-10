#Giorgia Auteri
#Script to assess relatedness of Myotis lucifugus in ddRad data
# June 10, 2019


library("related")

#read in rad-seq data
input <- readgenotypedata("relatedInput.txt")

###########################################################################
#Subsetting data so there is a lower percent missing data

#Add a new row to store the percent
input$gdata[nrow(input$gdata)+1,] <- NA
#Add count of missing data per locus
input$gdata[nrow(input$gdata),] <- colSums(na.omit(input$gdata == 0))
input$gdata[nrow(input$gdata), 1] <- "num_miss" #Add appropriate row label
tail(input$gdata[, 1:20]) #Check the name is updated and number missing are present
(input$gdata[, 20:40])

#how many loci have information for all individuals?
lastRow <- (input$gdata[nrow(input$gdata),] ) #subset the last row
length(lastRow[lastRow == 0]) #20 loci (40/2) are present for all individuals
length(lastRow[lastRow == 1]) #There are 276 loci (552/2) only missing from one individual

#subset the to select only alleles missing from 0, 1, 0r 2 individuals
tgdat <- as.data.frame(t(input$gdata))
temp <- tgdat[tgdat$`39` %in% c(0,1,2),] #in (0,..) specify how many can be missing that locus

#add sample ID data back in
temp <- as.data.frame(rbind( tgdat[1,],temp))
#Remove the last column (the summary of missing data)
temp <- temp[,1:(ncol(temp)-1)]
#re-transpose the data
filtGdat <- as.data.frame(t(temp))

#check missing data per individual
dim(filtGdat) #output/2 = # loci total
rowSums(na.omit(filtGdat == 0))/(ncol(filtGdat)-1) # % missing per individual

#subselect loci from original data frame based on loci that passed filter
LociKeep <- names(filtGdat)
gDat_filtered <- input$gdata[names(input$gdata) %in% LociKeep]
#remove last row (the number missing sum)
gDat_filtered <- gDat_filtered[1:nrow(gDat_filtered)-1,]

#Also subset the allele frequencies
LociKeep_num <- substring(LociKeep, 2) #remove the "V" from start of name
LociKeep_num <- as.numeric(LociKeep_num) #convert to numeric 
LociKeep_num <- LociKeep_num[which(LociKeep_num %% 2 == 0)] # Keep only the even numbers
LociKeep_num <- LociKeep_num / 2 #divide by two to get the locus order (as opposed to allele order)
LociKeep_num <- as.character(LociKeep_num) #Convert back to character
LociKeep_num <- paste("locus", LociKeep_num, sep = "") #Add "locus" to name to match input$freqs

#Only allow dfs in the list that are in LociKeep_num
gFreqs_filtered <- input$freqs[names(input$freqs) %in% LociKeep_num] 

#############################################################################################
## Also subset by minor allele frequency

#Create data frame to store results in
minorFreqDF <- data.frame(matrix(NA, nrow = length(names(gFreqs_filtered)), ncol = 2))
counter <- 1 #start loop counter
minorAlleleThreshold <- 0.01 #Desired minor allele frequency (minimum threshold) 

#Filter by minor allele frequency
for (i in names(gFreqs_filtered)){            #For each locus that meets the missing data requirement
  tempMin <- min(gFreqs_filtered[[i]]$frequency) #Get the minor allele frequency
  minorFreqDF[counter,1] <- i                #Store the locus name in one column
  if (tempMin >= minorAlleleThreshold){      #If the minor allele frequency is high enough
    minorFreqDF[counter,2] <- "keep"         #indicate the locus should be kept
  }
  else{
    minorFreqDF[counter, 2] <- "discard"     #Else, indicate the locus should be discarded
  }
  counter <- counter + 1                     #Add one to the loop counter
}

#Subset just the locus names that met the threshold
lociKeepNew <- minorFreqDF$X1[minorFreqDF$X2 == "keep"]
length(lociKeepNew) #how many loci were kept
#Use this to subset the frequency data frame (used only for simulations)
gFreqs_filtered <- input$freqs[names(input$freqs) %in% lociKeepNew]
##Also re-subset the gdata frame (used for simulations and empiricle estimates)
LociKeep_g <- substring(lociKeepNew, 6) #remove the "locus" from start of name
LociKeep_g <- as.numeric(LociKeep_g)#make numeric
LociKeep_g <- LociKeep_g * 2 #Multiply by two to get allele numbers
LociKeep_g_alleles <- LociKeep_g + 1 #get corresponding odd numbers for alleles
LociKeep_g <- sort(c(LociKeep_g, LociKeep_g_alleles)) #Combine even and odd (alleles 1 & 2)
LociKeep_g <- as.character(LociKeep_g) #Convert back to character
LociKeep_g <- paste("V", LociKeep_g, sep = "") #Add the "V" to start of name
LociKeep_g <- c("V1", LociKeep_g) #Add V1 to the start (so sample names are kept)
#Finalize the gdata subsetting
gDat_filtered <- input$gdata[names(input$gdata) %in% LociKeep_g]
#remove last row (it's not an individual, just the number missing sum)
gDat_filtered <- gDat_filtered[1:nrow(gDat_filtered)-1,]


#######################################################################
## Simulations to compare estimators

#If there are too many loci, randomly subset 250 of them
## otherwise the simulator will crash
if (length(gFreqs_filtered) > 250){
  set.seed(1234) #set seed to make the subset repeatable
  gFreqs_filtered <- na.omit(gFreqs_filtered[sample(length(gFreqs_filtered), 250)])
} else {
  cat("No subsetting needed")
}

#simulations to determine the degree of resolution expected in the dataset
sim <- familysim(gFreqs_filtered, 100) 

#Comparison simulations
output <- coancestry (sim, ritland = 1, trioml = 1, quellergt = 1 , wang = 1)
#Cleanup so that only individuals with relevant relatedness are included
simrel <- cleanuprvals(output$relatedness, 100)

#Write simulation results to store them
write.table(simrel, 
             "related_simResults_Miss012_MinAF_01_250Rand.txt",
             quote = FALSE,
             row.names = FALSE)

#read in previously simulated data (see above) if needed
#simrel <- read.table("related_simResults_Miss012_MinAF_01_250Rand.txt",
#                     header = TRUE)

################################################################################
#Prepping just the ritland values

#Boxplot of results
relvalues <- simrel [, 9] #ninth column is where Ritland values are

label1 <- rep ("PO", 100)
label2 <- rep (" Full ", 100)
label3 <- rep (" Half ", 100)
label4 <- rep (" Unrelated ", 100)

labels <- c(label1, label2, label3, label4)

#plot
plot (as.factor(labels), relvalues, ylab = "Relatedness Value", xlab = "Relatedness")

summary(relvalues[1:100]) #quantiles for parent-offspring pairs
summary(relvalues[101:200]) #quantiles for full siblings
summary(relvalues[201:300]) #quantiles for half siblings
summary(relvalues[301:400]) #quantiles for unrelated individuals

############################################################################
# Evaluating empiracle data

# # of loci (don't need to worry about program constraints for this)
(length(gDat_filtered) -1)/2 #number of loci (#alleles - row labels, /2)

# Calculate empirical relatedness scores with filtered data
## and 95% CI
relFilt <- coancestry(gDat_filtered, trioml = 1, wang = 1, quellergt = 1, 
                      ritland = 2, lynchli = 1, lynchrd = 1) 

#write point estimate (empiracle) results
write.table(relFilt$relatedness,
           "relRelatedness_filtered_Final_Miss012_MinAF_01.txt",
           quote = FALSE)

#write Confidence intervals
write.table(relFilt$relatedness.ci95,
           "relRelatedness_filtered_Final_Miss012_MinAF_01_95CI.txt",
           quote = FALSE)

#Basic stats
max(relFilt$relatedness[,9])
mean(relFilt$relatedness[,9])
sd(relFilt$relatedness[,9])

#Combine simulated data for plotting
backgroundSim <- as.data.frame(cbind(as.factor(labels), as.numeric(relvalues)))

#plot
ggplot(backgroundSim, aes(x= labels, y = relvalues)) +
  geom_boxplot() +
  geom_jitter(data=relFilt$relatedness, aes(x= backgroundSim$V1[316], y= ritland), 
              color = "blue", alpha = 0.2) +
    theme_bw()

#plot point estimates and confidence intervals together
## first create a dataframe with point estimates and CIs
plotDf <- merge(x = relFilt$relatedness, y = relFilt$relatedness.ci95, 
                by = "pair.no", all = TRUE)
#Sort based on point estimate of relatedness
plotDf <- plotDf[order(as.numeric(plotDf$ritland)),]
#Add and index for order
plotDf$order <- NA
plotDf$order <- c(1:length(plotDf$pair.no))

#plot with CIs
ggplot(plotDf, aes(x = order, y = ritland)) +
  geom_point() +
  geom_linerange(aes (x = order, ymin = ritland.low,
                      ymax = ritland.high)) +
  #line info for unrelated individuals
  geom_hline(yintercept = 0.003444, color = "blue") + #add the mean for unrelated estimate
  geom_hline(yintercept = -0.041075, linetype = "dashed", color = "blue") + #1st quantile
  geom_hline(yintercept = 0.027600, linetype = "dashed", color = "blue") + #third quantile
  #line infor for half-sibling
   geom_hline(yintercept = 0.2520, color = "red") + #mean
  geom_hline(yintercept = 0.1454, linetype = "dashed", color = "red") + #1st quantile
  geom_hline(yintercept = 0.3506, linetype = "dashed", color = "red") + #3rd quantile
  theme_bw()
  
#zoomed in plot with CIs
ggplot(plotDf[plotDf$order > 600,], aes(x = order, y = ritland)) +
  geom_point() +
  geom_linerange(aes (x = order, ymin = ritland.low,
                      ymax = ritland.high)) +
  #line info for unrelated individuals
  geom_hline(yintercept = 0.003444, color = "blue") + #add the mean for unrelated estimate
  geom_hline(yintercept = -0.041075, linetype = "dashed", color = "blue") + #1st quantile
  geom_hline(yintercept = 0.027600, linetype = "dashed", color = "blue") + #third quantile
  #line infor for half-sibling
  geom_hline(yintercept = 0.2520, color = "red") + #mean
  geom_hline(yintercept = 0.1454, linetype = "dashed", color = "red") + #1st quantile
  geom_hline(yintercept = 0.3506, linetype = "dashed", color = "red") + #3rd quantile
  theme_bw()


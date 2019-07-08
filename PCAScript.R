#  title: "generate PCA from structure file"
#author: "Giorgia Auteri"
#last updated "July 16, 2018"
# Auteri & Knowles, 2019
# "!!!" denote areas that will need to be changed to adapt the script to other datasets 
  
library(adegenet)
library(plyr)

#Clear working directory
rm(list = ls())

################################################################################
#Import Data
###create genind object from structure file

#!!!IMPORTANT For Plink structure file, one row per individual to TRUE, if from Stacks FALSE
# For plink structure file (ending in .strct_in), 1) save a copy with new ending (.stru) and
#    2) delete first two rows in .stru file
# For structure file from stacks, save copy with .stru ending
gind0 <- read.structure("D:/Mirror/Mirror/Projects/PhD Thesis/SurvMortWNS/Data/SurvMort21/Pops_SMAll/populations_33Inds_NoHead.stru", 
                        n.ind = 33,    #!!! Change to match the number of individuals
                        n.loc = 19797, #!!! Change to match the number of variant sites
                        col.lab = 1, 
                        row.marknames = 0, 
                        col.pop = 2, 
                        onerowperind = FALSE, 
                        NA.char=0, 
                        ask = FALSE)

gind <- scaleGen (gind0, NA.method = "mean", scale = FALSE, center = FALSE) #default center & scale = TRUE,

dfgind <- as.data.frame(as.matrix(unlist(gind), nrow=(length(gind[,1])), byrow=T),stringsAsFactors=FALSE) #convert to data frame

##########################################################################################
#filter by minor allele frequency

#Create vector of alleles
gindNames <- colnames(dfgind)
#Subset names of second allele (minor)
minorNames <- gindNames[seq(2,length(gindNames),2)]
#Subset dataframe by minor alleles
tempGind <- dfgind[,minorNames]
#Create a row on the bottom to sum up allele frequency
tempGind[nrow(tempGind)+1,] <- NA
#Add sum of column to empty row
tempGind[nrow(tempGind),] <- colSums(na.omit(tempGind))
#Set the minor allele threshold
FreqCutoff <- 0.05 #set the percent cutoff you want (e.g. 0.05 if you want threshold of 5%)
#Calculate the minimum value a locus can have in order to meet the cutoff
##0.5 is one copy, 1 is two copies
MafCutoff <-(nrow(dfgind)) * FreqCutoff #Number of individuals * % cutoff
#Only select loci that will have the minimum threshold
tt <- as.data.frame(t(tempGind)) #transpose
goodAlleles <- tt[tt$`34` >= MafCutoff,] #select data that passes threshold
tempGind <- as.data.frame(t(goodAlleles)) #transpose back
#get column names that passed the threshold
goodAlNames <- colnames(tempGind)
#Remove the last two characters (so it represents the whole locus, not just one allele)
goodAlNames <- substr(goodAlNames, 1, nchar(goodAlNames)-2)
#The full list of loci
allNames <- colnames(dfgind)
#Cycle through the acceptable loci and pull out matching alleles
for(i in 1:length(goodAlNames)){
  #Get the two matching alleles for a particular locus
  temp <- grep(goodAlNames[i], allNames, value = TRUE) 
  if (i == 1){ #if 1st time through the loop, store value to initiate the vector
    filterList <- temp
  }
  else{
    filterList <- c(filterList, temp)
  }
  
}
#Subset the original data frame by the good loci
dfgind <- dfgind[,filterList]

################################################################################
#Prepare data for calculating a PCA just based on a subset, then projecting all other points onto it

zero_varA <- which(apply(dfgind, 2, var)==0)   #make a vector with all zero-variance columns
length(zero_varA) 

#Split Data
rownames(gind) #Use this to see what order samples are in
Morts <- dfgind[1:25,]   #!!! Change these based on what groups you want to use. This will be the group which is projected onto the first PCA
Survs <- dfgind[26:33, ] #!!! This will be whatever group the pca is first calculated with

#Check for zero variance columns in 1st group (for which PCA will be calculated)
zero_var <- which(apply(Survs, 2, var)==0)             #make a vector with all zero-variance columns
length(zero_var)                                       #how many sites have no variance (x2)
SurvsVar <- Survs[,!names(Survs) %in% names(zero_var)] #make a subset of the data with only variable sites #!!! Do not execute if there are 0 non-variant sites within the group, just make survsVar <- Survs here

#Before subsetting Mortalities, how many sites are non variable
zero_varm1 <- which(apply(Morts, 2, var)==0)
length(zero_varm1)

#Prepare second set of data for pca projection
Morts <- Morts[ , (names(Morts) %in% names(SurvsVar))] #Remove sites with zero-variability in the first group from the second group !!!Do not do this if your first group (above) had 0 non-variant sites

#Of the remaining SNPs, how many have no variance among Mortalities?
zero_varm <- which(apply(Morts, 2, var)==0)             #make a vector with all zero-variance columns
length(zero_varm) # the number of non-variable allelic sites

MortsVar <- Morts[,!names(Morts) %in% names (zero_varm)]
SurvsVar <- SurvsVar[,names(SurvsVar) %in% names(MortsVar)]
length(SurvsVar) #divide by two to get number of SNPs
length(MortsVar)


#######################################  PCA ###########################################
#Using the subset data to generate a PCA

#Make calculate PCA using first group
pcas <- prcomp(SurvsVar, center = TRUE, scale = TRUE)
#Examine the PCA
summary(pcas) #Shows 1) how many axes are needed to explain all the data (# columns) and 2) what cumulative proportion of the data are explained
#Graph of proportion of data explained by each axis
plot(pcas, type = "l", main = "Proportion of data explained by each PCA axis for survivors")

#Apply center and scale from first pca to the new group
pca2 <- scale(MortsVar, pcas$center, pcas$scale) %*% pcas$rotation 
#Examine pca for second set of data (requires some calculating), figure out SDs then use these to calculate proportion of variance 
SDs <- rep(NA,length(colnames(pca2)))             #create vector for empty SD values
for (i in 1:length(colnames(pca2))){              # a loop to fill vector with SDs for each PCA axis
  SDs[i] <- sd(pca2[,i])
}
PropVars <- rep(NA,length(colnames(pca2)))        #Get the proportion of varinaces explained by each axis
for (i in 1:length(colnames(pca2))){              #A loop to fill vector with proportion explained by each axis
  PropVars[i] <- ((sd(pca2[,i]))^2)/sum(SDs^2)
}
sum(PropVars)                                     #Proportions of all variances should sum to 1
PropVars                                          #list proportion of variance explained by each axis for group 2
plot(PropVars, type = "b", main = "Proportion of data explained by each PCA axis for Mortalities")

########################################################################################################
#Plot PCA

#Combine data for plot
##Survivors
plotDat <- as.data.frame((pcas$x))
plotDat$group <- NA
plotDat$group <- "Survivor"
plotDat$side <- NA
plotDat$side <- c(rep("West", 1), rep("East", 3), "West", rep("East", 3))
##Non-survivors
plotDat2 <- as.data.frame((pca2))
plotDat2$group <- NA
plotDat2$group <- "Non-survivor"
plotDat2$side <- NA
plotDat2$side <- c(rep("East", 3), rep("West", 22))

#Bind the two data frames together
plotDat <- as.data.frame(rbind(plotDat, plotDat2))

library(ggplot2)

ggplot(plotDat, aes(x=PC1, y=PC2, group=group)) +
  geom_point(aes(shape=side, fill=group), size = 3) + 
  scale_fill_manual(values=c('#FEB701', "#127A75")) + 
  #scale_color_manual(values = c('#4562C3', "#F23C3D")) +
  scale_shape_manual(values = c(21, 23)) + #16, 18
  theme_bw()


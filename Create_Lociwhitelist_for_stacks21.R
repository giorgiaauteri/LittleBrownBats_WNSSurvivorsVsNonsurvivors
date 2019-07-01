#November 17, 2015 - Andrea Thomaz for Stacks version 1.35
#Modified by Sarp Kaya - May 20th, 2017 for Stacks version 1.45
#Modified by Giorgia Auteri - June 12, 2018 for stacks 2.0b/2.1
#read .vcf file from stacks output (WITHOUT write_random_SNP flag) for:
#plot the frequency of variable sites per position along all loci
#calculate theta based on number of segregating sites and individuals to create blacklist to delete very variable loci.

#clear workspace
rm(list=ls())

#required packages
require(plyr)
require(pegas)

#prevent scientific notation
options(scipen = 999)

#READ VCF
data <- read.table('D:/Mirror/Mirror/Projects/PhD Thesis/SurvMortWNS/Data/PostStructure/populations.snpsSurvMort.vcf', header = FALSE, sep = "\t")
head(data[1:20,1:20])

#SEQUENCE LENGTH
seq_len <- 140 #MODIFY HERE according the length of the sequence (deletes positions at the end of seq, should be 10 less than orig seq)

#creates dataframe with loci ID, the variable positions and the number of individuals in each loci
tempVec <- sub(":+", "", data$V3, fixed=TRUE)
tempVec <- sub(":-", "", tempVec, fixed=TRUE)
new_data <- data.frame(loci_ID = as.numeric(sub(":.*", "", data$V3)), #sub(pattern, replacement, x)
                       pos_vcf1 = data[,2],
                       pos = as.numeric(sub(".*:", "", tempVec)),  
                       ind = rowSums(data[,10:length(data)] != "./.:0:.,."))

head(new_data)
min(new_data$pos)#should always be position 5 (first positions are the adapters)
max(new_data$pos)#should always be position 139 !!!UNLESS!!! data is mapped to reference genome, may be longer if you allowed for some insertions
table(new_data$pos)
length(unique(new_data$loci_ID))#how many loci do I have


par(mar = rep(2, 4))
#saving graph with frequency of variable sites along the loci
#pdf("./SNPdistr_pos140bp.pdf")
hist(new_data$pos, breaks = seq_len, xlab = 'Position along the loci', main = 'The position of segregating sites'); #for unmapped loci
#play with first abline number to determine cutoff
abline(4500, 0, col = "red")#helps to find where starts to increase toward the end, last positions have strong increase
abline(v = 138, col = "red")#helps to figure out where to cut off before increase in bad calls
#move the lines around to visualize depending on my case
#dev.off()

#BASE ON THE GRAPH, CHOOSE HOW MANY POSITION TO DELETE FROM THE END
to_del <-2 #how many sites to delete in the end of the sequence
#10 is based on the 130 I chose for the ab line above
seq_len_cut <- seq_len - to_del

#create a whitelist to exclude those 10 (to_del) positions
whitelist <- subset(new_data, pos < seq_len_cut)[,c(1,3,4)]
#pdf("./SNPdistr_pos_cutto125bp.pdf")
hist(whitelist$pos, xlim = c(0,seq_len_cut), breaks = c(seq(-1, seq_len_cut -1 , by=1)), xlab = 'Position along the loci', main = 'The position of segregating sites');
#dev.off()

#calculating theta for all loci
var.sites <- count(whitelist, "loci_ID")
length(var.sites$loci_ID) #83,411
max(var.sites$freq) #max variable sites in one loci; 38
theta_calc <- merge(unique(whitelist[,-2]), var.sites, by = "loci_ID")
theta_calc$theta <- 0
head(theta_calc)
for (i in 1:length(theta_calc$theta)){
  theta_calc[i,4] <- theta.s(theta_calc$freq[i], theta_calc$ind[i])/seq_len_cut
}

#calculating the 95% quantile to exclude loci extremely variable
quant <- quantile(theta_calc$theta, probs=0.95) #set the value to be variability threshold
quant #0.02587015 for Myotis lib 1 (mapped to ref)
#pdf("./theta125bp.pdf")
hist(theta_calc$theta)
abline(v = quant, col="red")
#dev.off()

#what is the maximum number of mutations in a loci
max(theta_calc$freq) #38; max theta (variable positions before ble positions in one loci)
x <- subset(theta_calc, theta < quant)
max(x$freq) #14; max theta after, make sure is realistic for a 140 bp sequence
#think about what mutation rate the spp might have

#saving whitelist for re-run populations in stacks (write blacklist and subtract it from whitelist)
blacklist <- subset(theta_calc, theta > quant)[,1]
#write.table(blacklist, file="blacklist.txt", sep = '\n', row.names = F, col.names = F)
#above blacklist would only be for highly variable loci

#removes the blacklist from the whitelist and write off white list
whitelist$blacklist <- match(whitelist$loci_ID, blacklist, nomatch = 0)
whitelist_final <- subset(whitelist, blacklist == 0)
length(unique(whitelist_final$loci_ID)) #final number of unique loci, this is the number I need to get out with "write random loci"
length(whitelist_final$loci_ID) #final number of snps
write.table(whitelist_final[,1:2], file="whitelist_BRACHY.txt", sep = '\t', row.names = F, col.names = F)

#re-enable scientific notation
options(scipen = 0)

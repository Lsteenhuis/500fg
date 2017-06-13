library(gCMAP)
library(pheatmap)
library(reshape2)
library(ggplot2)
#####################
# Loading data 
#####################
load("/Volumes/MacOS/500fg/data/nbt/processedNetworks/D_250.Rdata")
immuneTraits <- read.table("/Volumes/MacOS/500fg/data/IRT_immuneTraits_500FG.csv", header = T, row.names = 1,stringsAsFactors = F, sep = ",")
immuneInfo <- read.table("/Volumes/MacOS/500fg/data/full_code_sampleInfo.csv", header = T, row.names = 1,stringsAsFactors = F, sep = ",")
rownames(immuneTraits) <- gsub(pattern = "H", replacement = "HV", x = rownames(immuneTraits))
immuneTraits <- as.data.frame(t(immuneTraits))
immuneTraits <- immuneTraits[c(1:41,50:57),]
immuneTraitNames <- immuneInfo[which(rownames(immuneInfo) %in%rownames(immuneTraits)),1]
immuneTraitNames[duplicated(immuneTraitNames)] <- paste(immuneTraitNames[duplicated(immuneTraitNames)],"_dup",sep="")
rownames(immuneTraits) <-immuneTraitNames 

#####################
# removing percentages and duplicates in connectivty scores matrix
scoresSub <- ds2[c(1:41,50:57),]
rownames(scoresSub) <-immuneTraitNames 
scoresSub <- as.matrix(scoresSub)
scoresSub.m <- melt(scoresSub)
# order scores from high to low
scoresSub.m <- scoresSub.m[order(abs(scoresSub.m$value),decreasing = T),]
scoresSub.m$Var2 <- gsub(pattern = ".csv",replacement = "", x = scoresSub.m$Var2)


# subsetting results from Kidd et al
nbtRes <- scoresSub.m[grep("guanfacine|trihexyphenidyl|propofol|spironolactone|clioquinol", scoresSub.m$Var2, ignore.case = T),]
nbtP <- pvalSub.m[grep("guanfacine|trihexyphenidyl|propofol|spironolactone|clioquinol", pvalSub.m$Var2, ignore.case = T),]

# subsetting positve and negative scores
negScore <- scoresSub.m[which(sign(scoresSub.m$value) == -1),]
posScore <- scoresSub.m[which(sign(scoresSub.m$value) == 1),]

# create  counts of amount of pertubations
amountOfScoresPerIT <- apply(scoresSub, 1, function(x){
  length(which(x != 0))
})

# plotting top 5 pos scores
pp<-ggplot(data=posScore[1:5,], aes(x=reorder(Var2[1:5], value), y=value)) +
  geom_bar(stat="identity") + 
  xlab("drug name") +
  ylab("connectivity score") +
  coord_flip()  

#plotting top 5 neg scors
negS <-negScore[1:5,]
pn<-ggplot(data=negS, aes(x=reorder(Var2, -value), y=value)) +
  geom_bar(stat="identity") +
  xlab("drug name") +
  ylab("connectivity score") +
  coord_flip()

# create plot of both
multiplot(pn,pp)


#
scoresSub <- as.data.frame(scoresSub)
listOfPert <- sapply(scoresSub[1:10,],function(x){
  length(which(x != 0))
  
})
names(listOfPert)  <- colnames(scoresSub)


# subsetting pvalues of p values from con scores
pvalSub <- dpp[c(1:41,50:57),]
rownames(pvalSub) <-immuneTraitNames 
pvalSub <- as.matrix(pvalSub)
pvalSub.m <- melt(pvalSub)
# adjusting pvalue
pvalSub.m$value <- p.adjust(pvalSub.m$value,method = "BH")
pvalSub.m <- pvalSub.m[order(abs(pvalSub.m$value)),]



top <- head(pvalSub.m,13)
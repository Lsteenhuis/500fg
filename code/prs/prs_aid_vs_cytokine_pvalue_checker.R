############################
# Script to calculate and plot the significance of the calculated pvalues between auto immune diseases and cytokines
############################
#loading data
###########################
# PRS data
# rows: GWAS trait
# columns: sample names
# values: Risk score
prs <- as.data.frame(read.table(file = "/Volumes/MacOS/500fg/geneticRiskScore/output.test/TT/rawScoreMatrix.txt", 
                                header = T, row.names = 1,
                                sep = "\t", stringsAsFactors = F))
grepString <- "gout|celiac_disease|Ulcerative_Colitis|Juvenile_Idiopathic_Arthritis|multiple_sclerosis|Narcolepsy|primary_biliary_cirrhosis|psoriasis|T1D|systemic_lupus_erythematosus|Rheumatoid_Arthritis|Crohns_disease|Ulcerative_colitis|Inflammatory_Bowel_Disease|Asthma"
immuneSub <- prs[grep(grepImmunoString,rownames(prs)),]

# cytokine data
# rownames: cytokine name
# colnames: sample name
# values: normalized cytokine levels
cytokine <- as.data.frame(read.table(file = "/Volumes/MacOS/500fg/data/pheno_91cytokines_4Raul.csv" ,
                                     header = T, row.names = 1,
                                     stringsAsFactors = F, sep = ","))
cytokine <- cytokine[,-1]
# create pvalue and corrlation matrices based on the genome wide diseases
prsCytoPvalMat <- getCorPvalMat(cytokine,prs)
prsCytokinesCorMatGWPrs <- getCorMat(cytokine,prs[grep("E-8",rownames(prs)),])
prsCytoPvalMatSubset <- getCorPvalMat(cytokine,immuneSub)
# set nan to 0.5 for randomness 
prsCytoPvalMat[is.nan(prsCytoPvalMat)] <- 0.5
prsCytoPvalMatSubset[is.nan(prsCytoPvalMatSubset)] <- 0.5
############################
# functions
############################
# function to validate pvalues created from getCorPvalMat function
# pvalMatrix: full pval matrix
# perm: amount of permutations there have to be used
# nEvents: amount of events of a row which are < 0.05
# significanceThreshold: pvalue cutoff
# use: creates a p value for rows to check the significance of calculated pvalues
# by sampling random p values from the full matrix 
# output: list of random pval distribution and p value 
getRandomPvalues <- function(pValMatrix, perm = 1000, nEvents , observed, significanceThreshold=0.05){
  pValMatrix[is.nan(pValMatrix)] <- 0.5
  dist <- unlist(lapply(1:perm, FUN=function(x){
    sum(sample(pValMatrix, size = nEvents, replace = FALSE) <= significanceThreshold)
  }))
  pval<- sum(dist > observed) / perm
  return(list(dist=dist, pval=pval))
}

############################
length(which(prsCytoPvalMat < 0.05))
length(which(prsCytoPvalMat > 0.05))
prsCytoPvalMat <- prsCytoPvalMat[!rownames(prsCytoPvalMat) %in% c("Celiac_disease_2010_20190752_hg19.txt.gz_P5.0E-8","Ulcerative_Colitis_2011_21297633_hg19.txt.gz_P5.0E-8") ,]

cytokineInfo_Final <- read.table("/Volumes/MacOS/500fg/data/cytokineInfo_Final.csv", 
                                 ";", stringsAsFactors = F, header =T,row.names = 1)


## check pvalues of immune subset from PRS
pval <- c()
prsCytoPvalMat[is.nan(prsDisease)] <- NA
for(i in 1:ncol(immuneSub)){
  pval <-c(pval, getRandomPvalues(pValMatrix = prsCytoPvalMat[grep("E-8",rownames(prsCytoPvalMat)),colnames(immuneSub)[i]], 
                           perm=10000, nEvents = length(immuneSub[grep("E-8",rownames(immuneSub)),i]), 
                           observed = sum(immuneSub[,i] <= 0.05),
                           significanceThreshold = 0.05)$pval)
}
qval <- p.adjust(pval, method = 'fdr')

pData <- data.frame(qValue = qval, id= colnames(immune), cytokineInfo_Final[colnames(immune),])
pData$pValue <- pData$qValue
# adjusting 0 to range
pData$pValue[pData$qValue == 0] <- 0.0001

pData$id <- factor(pData$id, levels = as.character(pData$id[order(pData$pValue)]))

bPlot_enrichment <- ggplot(pData, aes(x =id, y=-log10(pValue), fill=cellType))+
  geom_bar(position = 'dodge', stat = 'identity')+
  geom_hline(yintercept = -log10(0.05), color='red', alpha = 0.7)+
  xlab('Cytokine')+ 
  coord_flip()+
  theme_minimal()+
  scale_fill_futurama()+
  theme(legend.position = 'bottom', text = element_text(family = 'Helvetica', size=9))
pdf("/Volumes/partition/Users/umcg-lsteenhuis/Dropbox/lsteenhuis_internship/500FG/plots/enrichment.subset.pdf", height = 7, width = 8)
print(bPlot_enrichment)
dev.off()



# check the pvalues from all diseases and plot the results
pval <- c()
for(i in 1:ncol(prsCytoPvalMatSubset)){
  pval <-c(pval, getRandomPvalues(pValMatrix = prsCytoPvalMat[,colnames(prsCytoPvalMatSubset)[i]], 
                                  perm=10000, nEvents = length(prsCytoPvalMatSubset[grep("E-8",rownames(prsCytoPvalMatSubset)),i]), 
                                  observed = sum(prsCytoPvalMatSubset[,i] <= 0.05),
                                  significanceThreshold = 0.05)$pval)
}
qval <- p.adjust(pval, method = 'fdr')
pData <- data.frame(qValue = qval, id= colnames(prsCytoPvalMatSubset), cytokineInfo_Final[colnames(prsCytoPvalMatSubset),])
pData$pValue <- pData$qValue
pData$pValue[pData$qValue == 0] <- 0.0001
pData$id <- factor(pData$id, levels = as.character(pData$id[order(pData$pValue)]))

ggplot(pData, aes(x =id, y=-log10(pValue), fill=cellType))+
  geom_bar(position = 'dodge', stat = 'identity')+
  geom_hline(yintercept = -log10(0.05), color='red', alpha = 0.7)+
  coord_flip()+
  theme_minimal()+
  scale_fill_uchicago()+
  theme(legend.position = 'bottom')

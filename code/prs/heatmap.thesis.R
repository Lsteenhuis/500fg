# script to create plots for thesis
library(Hmisc)
library(ggplot2)
library(reshape2)
library(gtools)
source("/Volumes/MacOS/500fg/code/r/prs/plots.R")
##############################################################################################################
# data

# creates melted matrix of -log10 (pvalue ) * sign of correlation
# corR <- correlation matrix
# corP <- pvalue matrix
# threshold <- amount of hits there need to be in a column
getLog10 <- function(corR,corP,threshold){

  log10Val <- as.matrix(sign(corR) * -log10(corP))
  log10Val <- log10Val[order(rownames(log10Val)),]
  log10Val <- log10Val[,order(colnames(log10Val))]
  log10Val[which(abs(log10Val) < 1.3)] <- NA 
  b.test <- lapply(colnames(log10Val), function(colName) {
    column <- log10Val[,which(colnames(log10Val) == colName)]
    if ( sum( !is.na(column) ) >= threshold ) {
      return(colName)
    } else {
      return(NA)
    }
  })
  
  log10Val <-log10Val[,c(unlist(b.test)[!is.na(unlist(b.test))])]
  log10Val.m <- melt(log10Val)
  log10Val.m
}

#immuno modulator data
immunoMod <- read.table("/Volumes/MacOS/500fg/500FG/data_for_lars/20151117resultsTotalReFormat.csv", 
                      row.names = 1, header = T,sep = ",")
immunoMod <- log2(immunoMod)

#cytokine data
cytokine <- as.data.frame(read.table(file = "/Volumes/MacOS/500fg/data/pheno_91cytokines_4Raul.csv" ,
                                   header = T, row.names = 1,
                                   stringsAsFactors = F, sep = ","))
cytokine <- cytokine[,-1]

#Immune traits
immuneTraits <- read.table("/Volumes/MacOS/500fg/data/IRT_immuneTraits_500FG.csv", header = T, row.names = 1,stringsAsFactors = F, sep = ",")
immuneInfo <- read.table("/Volumes/MacOS/500fg/data/full_code_sampleInfo.csv", header = T, row.names = 1,stringsAsFactors = F, sep = ",")
rownames(immuneTraits) <- gsub(pattern = "H", replacement = "HV", x = rownames(immuneTraits))
immuneTraits <- as.data.frame(t(immuneTraits))
immuneTraits <- immuneTraits[c(1:41,50:57),]
immuneTraitNames <- immuneInfo[which(rownames(immuneInfo) %in%rownames(immuneTraits)),1]
immuneTraitNames[duplicated(immuneTraitNames)] <- paste(immuneTraitNames[duplicated(immuneTraitNames)],"_dup",sep="")
rownames(immuneTraits) <-immuneTraitNames 

microBiomPath <- as.data.frame(read.table(file = "/Volumes/MacOS/500fg/500FG/data_for_lars/500FG_normalized_raw_pathways.csv", 
                                header = T, row.names = 1,
                                sep = ",", stringsAsFactors = F))

microBiomTax <- as.data.frame(read.table(file = "/Volumes/MacOS/500fg/500FG/data_for_lars/500FG_normalized_raw_taxonomy.csv", 
                                          header = T, row.names = 1,
                                          sep = ",", stringsAsFactors = F))
#PRS AID DISEASES  
prs <- as.data.frame(read.table(file = "/Volumes/MacOS/500fg/geneticRiskScore/output.test/TT/rawScoreMatrix.txt", 
                                header = T, row.names = 1,
                                sep = "\t", stringsAsFactors = F))
prsGW <- prs[grep("P5.0E-8",rownames(prs)),]

# grabbing AID diseases
grepDis <- "Ulcerative_Colitis|coro|Inflammatory_Bowel_Disease|Crohns_disease_EUR|Coronary_artery_disease|T1D|Rheumatoid_Arthritis|Juvenile_Idiopathic_Arthritis|Asthma|Primary_biliary_cirrhosis|Psoriasis|Multiple_sclerosis|Celiac_disease|Multiple_sclerosis|Inflammatory_Bowel_Disease|Systemic_lupus_erythematosus|ALS|Major_depression|Alzheimers_disease"
#grepDis <- "Type_2_Diabetes|ALS|Alzheimer|Depression|Inflammatory_Bowel_Disease|Crohn|celiac|Ulcerative_Colitis|Arthritis|sclerosis|cirrhosis|psoriasis|T1D|systemic_lupus_erythematosus|Crohns|Ulcerative_colitis|Inflammatory|Asthma"
ibdDis <- prsGW[grep(grepDis,rownames(prsGW),ignore.case = T),]
ibdDis<- as.data.frame(getAbbr(ibdDis))

basicPheno <- read.csv("/Volumes/MacOS/500fg/500FG/data_for_lars/500FG_basicPhenos.csv", sep = ",", stringsAsFactors = F)
subs <-basicPheno[,c(1,3,4,5,6)]
subs <- t(subs)
colnames(subs) <- subs[1,]
subs <- as.data.frame(subs[-1,])

#platelet info
#plat <- read.table("/Volumes/MacOS/500fg/500FG/data_for_lars/Gro_PLT_request_160628.txt", sep = "\t", stringsAsFactors = F,header = T, row.names = 1)
#plat <- as.matrix(plat)
#plat <- log2(plat)

#metabolite info
meta <- read.table("/Volumes/MacOS/500fg/500FG/data_for_lars/500FG_normalized_raw_metabolome.csv", sep=",",header = T,row.names = 1)
meta <- meta[,-grep("perc",colnames(meta))]

# IG levels
ig <- read.table("/Volumes/MacOS/500fg/500FG/data_for_lars/igLevels500FG.csv", sep=",",header = T,row.names = 1)



#######################################################################################
# Hormone levels VS PRS                                                                   #
#######################################################################################

horm <- read.table("/Volumes/MacOS/500fg/500FG/data_for_lars/log2HormoneLevels500FG.txt", sep = " ", header = T, row.names = 1)

corResults <- getCorMat(horm,ibdDis, threshold = 1)
corP <- as.data.frame(corResults[1])
corR <- as.data.frame(corResults[2])
cor10.m <- getLog10(corP = corP, corR = corR, threshold=1)
cor10.m <- cor10.m[grep("T1D|JIA|UC|RA|IBD", cor10.m$Var1, ignore.case = T),]
p <- createHeatMap(na.omit(cor10.m), yLab = "diseases",xLab="hormones",lLab="-log10(pvalue)")

pdf("/Volumes/MacOS/500fg/plots/ownHeatmap/AID.vs.hormones.pdf")
p +   theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.y = element_text(hjust = 1, size = 10),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 11))
dev.off()
####################################################################################################################################
##  PRS VS PHENOTYPES
##
corResults <- getCorMat(subs,ibdDis)
corP <- as.data.frame(corResults[1])
corR <- as.data.frame(corResults[2])
cor10.m <- getLog10(corP = corP, corR = corR, threshold=1)
p <- createHeatMap(na.omit(cor10.m), yLab = "diseases",xLab="phenotypes",lLab="-log10(pvalue)")
pdf("/Volumes/MacOS/500fg/plots/ownHeatmap/AID.vs.pheno.pdf")
p +   theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.y = element_text(hjust = 1, size = 7),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
dev.off()
####################################################################################################################################
##  PRS VS IMMUNOMOD
##
corResults <- getCorMat(immunoMod,ibdDis)
corP <- as.data.frame(corResults[1])
corR <- as.data.frame(corResults[2])
cor10.m <- getLog10(corP = corP, corR = corR, threshold=1)
cor10.m <- cor10.m[grep("As|Celiac|Crohns|JIA", cor10.m$Var2, ignore.case = T),]
p <- createHeatMap(na.omit(cor10.m), yLab = "diseases",xLab="immunoModualtor",lLab="-log10(pvalue)")
pdf("/Volumes/MacOS/500fg/plots/ownHeatmap/AID.vs.immunomod.pdf")
p +   theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.y = element_text(hjust = 1, size = 10),
            axis.text.x = element_text(angle = 45, hjust = 1, size =15))
dev.off()
####################################################################################################################################
##  PRS VS IT
##
corResults <- getCorMat(ibdDis,immuneTraits)
corP <- as.data.frame(corResults[1])
corR <- as.data.frame(corResults[2])
cor10.m <- getLog10(corP = corP, corR = corR, threshold=1)
cor10.m <- cor10.m[grep("Crohn|IBD|UC|T1D|T2D", cor10.m$Var2, ignore.case = T),]
p <- createHeatMap(na.omit(cor10.m), yLab = "immunetraits",xLab="diseases",lLab="-log10(pvalue)")
pdf("/Volumes/MacOS/500fg/plots/ownHeatmap/AID.vs.immunotraits.pdf")
p +   theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.y = element_text(hjust = 1, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = ))
dev.off()


####################################################################################################################################
##  PRS VS CYTOKINES
##
corResults <- getCorMat(ibdDis,cytokine)
corP <- as.data.frame(corResults[1])
corR <- as.data.frame(corResults[2])
cor10.m <- getLog10(corP = corP, corR = corR, threshold=1)
cor10.m <- cor10.m[grep("Crohn|IBD|UC|RA|depression|Alz", cor10.m$Var2, ignore.case = T),]
p <- createHeatMap(na.omit(cor10.m), yLab = "Cytokines",xLab="diseases",lLab="-log10(pvalue)")

pdf("/Volumes/MacOS/500fg/plots/ownHeatmap/AID.vs.cytokine.pdf")
p +   theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.y = element_text(hjust = 1, size = 6.5),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
dev.off()





####################################################################################################################################
##  PRS VS METABIOME
##
corResults <- getCorMat(ibdDis, meta , threshold = 5)
corP <- as.data.frame(corResults[1])
corR <- as.data.frame(corResults[2])
cor10.m <- getLog10(corP = corP, corR = corR, threshold=20)
p <- createHeatMap(na.omit(cor10.m), yLab = "Metabolites",xLab="diseases",lLab="-log10(pvalue)")
pdf("/Volumes/MacOS/500fg/plots/ownHeatmap/prs.vs.meta.pdf")
p +   theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.y = element_text(hjust = 1, size = 4),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
dev.off()

####################################################################################################################################
##  CYTOKINES VS METABIOME
##
corResults <- getCorMat(meta,cytokine)
corP <- as.data.frame(corResults[1])
corR <- as.data.frame(corResults[2])
cor10.m <- getLog10(corP = corP, corR = corR, threshold=15)
# removing 
cor10.m <- cor10.m[-grep("gln|ace|crea|alb|gp",cor10.m$Var2,ignore.case = T),]
p <- createHeatMap(na.omit(cor10.m), yLab = "cytokines",xLab="metabolites",lLab="-log10(pvalue)")
pdf("/Volumes/MacOS/500fg/plots/ownHeatmap/cytokine.vs.meta.pdf")
p +   theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.y = element_text(hjust = 1, size = 6),
            axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
dev.off()

####################################################################################################################################
##  cytokine VS Ig levels
##
corResults <- getCorMat(ig,cytokine)
corP <- as.data.frame(corResults[1])
corR <- as.data.frame(corResults[2])
cor10.m <- getLog10(corP = corP, corR = corR, threshold=7)
p <- createHeatMap(na.omit(cor10.m), yLab = "Cytokines",xLab="Immunoglubolin ",lLab="-log10(pvalue)")
pdf("/Volumes/MacOS/500fg/plots/ownHeatmap/cytokine.ig.pdf")
p +   theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.y = element_text(hjust = 1, size = 9),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
dev.off()

####################################################################################################################################
##  CYTOKINE VS METABIOME
##
corResults <- getCorMat(microBiomPath, cytokine)
corP <- as.data.frame(corResults[1])
corR <- as.data.frame(corResults[2])
cor10.m <- getLog10(corP = corP, corR = corR, threshold=20)
p <- createHeatMap(na.omit(cor10.m), yLab = "Metabolites",xLab="diseases",lLab="-log10(pvalue)")
pdf("/Volumes/MacOS/500fg/plots/ownHeatmap/cyto.vs.metabPath.pdf")
p +   theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.y = element_text(hjust = 1, size = 4),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
dev.off()

corResults <- getCorMat(microBiomTax, cytokine)
corP <- as.data.frame(corResults[1])
corR <- as.data.frame(corResults[2])
cor10.m <- getLog10(corP = corP, corR = corR, threshold=10)
p <- createHeatMap(na.omit(cor10.m), yLab = "Metabolites",xLab="diseases",lLab="-log10(pvalue)")
#pdf("/Volumes/MacOS/500fg/plots/ownHeatmap/cyto.vs.metabPath.pdf")
p +   theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.y = element_text(hjust = 1, size = 4),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
dev.off()
####################################################################################################################################
##  CYTOKINE VS GENES
##
rc <- read.table("/Volumes/MacOS/500fg/500FG/data_for_lars/500FG_normalized_filtered_read_counts.csv",header = T,row.names = 1,sep=",")
corResults <- getCorMat(rc,cytokine)
corP <- as.data.frame(corResults[1])
corR <- as.data.frame(corResults[2])
#write.table(corP,file = "/Volumes/MacOS/500fg/data/prs_info/rc.cyto.P.csv",quote = F,row.names = T,col.names =T, sep = ",")
#write.table(corR,file = "/Volumes/MacOS/500fg/data/prs_info/rc.cyto.R.csv",quote = F,row.names = T,col.names =T, sep = ",")
rc.cyto.P <- read.table("/Volumes/MacOS/500fg/data/prs_info/rc.cyto.P.csv",header = T,row.names = 1)
rc.cyto.R <- read.table("/Volumes/MacOS/500fg/data/prs_info/rc.cyto.R.csv",header = T,row.names = 1)
cor10.m <- getLog10(corP = corP, corR = corR, threshold = 1 )

write.table(cor10.m,file = "/Volumes/MacOS/500fg/data/prs_info/rc.cyto.log10P.csv",quote = F,row.names = T,col.names =T, sep = ",")
rm(c(cor10.m,corP,corR,rc))

corResults <- getCorMat(rc,ibdDis)
corP <- as.data.frame(corResults[1])
corR <- as.data.frame(corResults[2])
write.table(corP,file = "/Volumes/MacOS/500fg/data/prs_info/rc.prs.P.csv",quote = F,row.names = T,col.names =T, sep = ",")
write.table(corR,file = "/Volumes/MacOS/500fg/data/prs_info/rc.prs.R.csv",quote = F,row.names = T,col.names =T, sep = ",")
cor10.m <- getLog10(corP = corP, corR = corR, threshold = 1)
write.table(cor10.m,file = "/Volumes/MacOS/500fg/data/prs_info/rc.prs.log10P.csv",quote = F,row.names = T,col.names =T, sep = ",")
rm(c(cor10.m,corP,corR,rc))


p <- createHeatMap(na.omit(cor10.m), yLab = "Metabolites",xLab="diseases",lLab="-log10(pvalue)")


pdf("/Volumes/MacOS/500fg/plots/ownHeatmap/cytokine.ig.pdf")
p +   theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.y = element_text(hjust = 1, size = 4),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
dev.off()


####################################################################################################################################
##  CYTOKINE VS IT
##
corResults <- getCorMat(cytokine,immuneTraits, threshold = 5)
corP <- as.data.frame(corResults[1])
corR <- as.data.frame(corResults[2])

p <- createHeatMap(na.omit(cor10.m), yLab = "immunetraits",xLab="diseases",lLab="-log10(pvalue)")
pdf("/Volumes/MacOS/500fg/plots/ownHeatmap/cytokine.vs.immunotraits.pdf")
p +   theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 5))
dev.off()


####################################################################################################################################
##  CYTOKINE VS IMMUNETRAIT
##
corResults <- getCorMat(cytokine,immuneTraits, threshold = 5)
corP <- as.data.frame(corResults[1])
corR <- as.data.frame(corResults[2])

p <- createHeatMap(na.omit(cor10.m), yLab = "immunetraits",xLab="diseases",lLab="-log10(pvalue)")
pdf("/Volumes/MacOS/500fg/plots/ownHeatmap/cytokine.vs.immunotraits.pdf")
p +   theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.y = element_text(hjust = 1, size = 5),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
dev.off()

####################################################################################################################################
## CYTOKINE VS IMMOD
corResults <- getCorMat(t(immunoMod), cytokine)
corP <- as.data.frame(corResults[1])
corR <- as.data.frame(corResults[2])
cor10.m <- getLog10(corP = corP, corR = corR, threshold=1)
p <- createHeatMap(na.omit(cor10.m), yLab = "immunetraits",xLab="immune modulators",lLab="-log10(pvalue)")
pdf("/Volumes/MacOS/500fg/plots/ownHeatmap/cytokine.vs.immunoModulator.pdf")
p +   theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.y = element_text(hjust = 1, size = 6),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
dev.off()
####################################################################################################################################
##  CYTOKINES VS PHENOTYPES
##
corResults <- getCorMat(subs, cytokine)
corP <- as.data.frame(corResults[1])
corR <- as.data.frame(corResults[2])
cor10.m <- getLog10(corP = corP, corR = corR, threshold=1)
p <- createHeatMap(na.omit(cor10.m), yLab = "Cytokines",xLab="phenotypes",lLab="-log10(pvalue)")
pdf("/Volumes/MacOS/500fg/plots/ownHeatmap/cytokine.vs.pheno.pdf")
p +   theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.y = element_text(hjust = 1, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
dev.off()



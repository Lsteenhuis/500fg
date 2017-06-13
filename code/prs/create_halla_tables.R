#################################
# This script creates tables for components of the immune system for HALLA
# Two tables will be intersected based on sample names and subsets will be created.
# The resulting tables will be writen to the output folder for use in halla
#################################
# loading data
#################################
# metabolite data 
# rows: metabolite data
# columns: sample names
# values: normalized matabolome data 
metaData <- read.table("/Volumes/MacOS/500fg/500FG/data_for_lars/500FG_normalized_raw_metabolome.csv", 
                       row.names = 1, header = T, 
                       sep = ",")
metaData <- as.data.frame(t(metaData))

# hormone data
# rows: hormone ID
# columns: sample names
# values normalized hormone levels
horm <- read.csv("/Volumes/MacOS/500fg/500FG/data_for_lars/log2HormoneLevels500FG.txt", 
                 row.names = 1, header = T, 
                 sep=" ")
horm <- as.data.frame(t(horm))

# platelete data
# rows: platelet ID
# columns: sample names 
# values: platelet levels (need to be normalized)
plat <- read.table("/Volumes/MacOS/500fg/500FG/data_for_lars/Gro_PLT_request_160628.txt",
                   row.names =1 , header = T,
                   sep = "\t", dec = ",")
plat <- as.data.frame(t(plat))
plat <- log2(plat)

# PRS data
# rows: GWAS trait
# columns: sample names
# values: Risk score
prs <- as.data.frame(read.table(file = "/Volumes/MacOS/500fg/geneticRiskScore/output.test/TT/rawScoreMatrix.txt", 
                                header = T, row.names = 1,
                                sep = "\t", stringsAsFactors = F))

grepImmunoString <- "gout|celiac|Ulcerative_Colitis|Juvenile_Idiopathic_Arthritis|multiple_sclerosis|Narcolepsy|primary_biliary_cirrhosis|psoriasis|T1D|systemic_lupus_erythematosus|Rheumatoid_Arthritis|Crohns_disease|Ulcerative_colitis|Inflammatory_Bowel_Disease|Asthma"
prsImmuno <- prs[grep(grepImmunoString,rownames(prs)),]


#################################
# functions
#################################

# input: two matrices
# output: list of input with intersecting columns
# use: make matrices overlap 
makeMatricesOverlap <- function(xMat,yMat){
  xSamples <- colnames(xMat)
  ySamples <- colnames(yMat)
  
  intersectSample <- intersect(xSamples,ySamples)
  xMat <- xMat[,which(colnames(xMat) %in% intersectSample)]
  yMat <- yMat[,which(colnames(yMat) %in% intersectSample)]
  xMat <- xMat[,order(colnames(xMat))]
  yMat <- yMat[,order(colnames(yMat))]
  list(xMat,yMat)
}


#################################
# logic
#################################
# creating halla tables for PRS immuno diseases from  and metabolites 
prsImmuno <- prs[grep(grepImmunoString,rownames(prs)),]
horm <- read.csv("/Volumes/MacOS/500fg/500FG/data_for_lars/log2HormoneLevels500FG.txt", 
                 row.names = 1, header = T, 
                 sep=" ")
horm <- as.data.frame(t(horm))
listOfMatrix <- makeMatricesOverlap(prsImmuno,metaData)
prsImmuno <- as.data.frame(listOfMatrix[1])
metaData <- as.data.frame(listOfMatrix[2])

write.table(prsImmuno,file = "//Volumes/MacOS/500fg/500FG/data_for_lars/halla_data/prsImmuno.vs.metabolites/prsImmuno.metabolites.tsv", 
            sep="\t", quote = F,
            row.names = T, col.names = T,
            na = "0")
write.table(metaData,file = "/Volumes/MacOS/500fg/500FG/data_for_lars/halla_data/prsImmuno.vs.metabolites/metabolites.prsImmuno.tsv", 
            sep="\t", quote = F,
            row.names = T, col.names = T,
            na = "0")


#creating halla tables for PRS immuno diseases and hormone levels
prsImmuno <- prs[grep(grepImmunoString,rownames(prs)),]
immunoMod <-  as.data.frame(t(immunoMod))
listOfMatrix <- makeMatricesOverlap(prsImmuno,immunoMod)
prsImmuno <- as.data.frame(listOfMatrix[1])
immunoMod <- as.data.frame(listOfMatrix[2])

write.table(prsImmuno,file = "/Volumes/MacOS/500fg/halla/immunoMod.prsDis/dis.tsv", 
            sep="\t", quote = F,
            row.names = T, col.names = T,
            na = "0")
write.table(immunoMod,file = "/Volumes/MacOS/500fg/halla/immunoMod.prsDis/immunoMod.tsv", 
            sep="\t", quote = F,
            row.names = T, col.names = T,
            na = "0")


# creating halla tables for PRS immuno diseases and plataletes
prsImmuno <- prs[grep(grepImmunoString,rownames(prs)),]
listOfMatrix <- makeMatricesOverlap(prsImmuno,plat)
prsImmuno <- as.data.frame(listOfMatrix[1])
plat <- as.data.frame(listOfMatrix[2])

write.table(prsImmuno,file = "/Volumes/MacOS/500fg/500FG/data_for_lars/halla_data/prsImmuno.vs.plataletes/prsImmuno.plataletes.tsv", 
            sep="\t", quote = F,
            row.names = T, col.names = T,
            na = "0")
write.table(plat,file = "/Volumes/MacOS/500fg/500FG/data_for_lars/halla_data/prsImmuno.vs.plataletes/plataletes.prsImmuno.tsv", 
            sep="\t", quote = F,
            row.names = T, col.names = T,
            na = "0")

# creating halla tables for PRS immuno diseases and microbiome
prsImmuno <- prs[grep(grepImmunoString,rownames(prs)),]
taxomy <- read.table("/Volumes/MacOS/500fg/500FG/data_for_lars/500FG_normalized_raw_taxonomy.csv", 
           row.names = 1, header = T,
           sep = ",")
taxomy <- as.data.frame(t(taxomy))
path <- read.table("/Volumes/MacOS/500fg/500FG/data_for_lars/500FG_normalized_raw_pathways.csv", 
           row.names = 1, header = T,
           sep = ",")
path <- as.data.frame(t(path))

listOfMatrix <- makeMatricesOverlap(prsImmuno,taxomy)
prsImmuno <- as.data.frame(listOfMatrix[1])
taxomy <- as.data.frame(listOfMatrix[2])

write.table(prsImmuno,file = "/Volumes/MacOS/500fg/500FG/data_for_lars/halla_data/prsImmuno.vs.microbiome/prsImmuno.microbiomeTaxomy.tsv", 
            sep="\t", quote = F,
            row.names = T, col.names = T,
            na = "0")
write.table(taxomy,file = "/Volumes/MacOS/500fg/500FG/data_for_lars/halla_data/prsImmuno.vs.microbiome/microbiomeTaxomy.prsImmuno.tsv", 
            sep="\t", quote = F,
            row.names = T, col.names = T,
            na = "0")

listOfMatrix <- makeMatricesOverlap(prsImmuno,path)
prsImmuno <- as.data.frame(listOfMatrix[1])
path <- as.data.frame(listOfMatrix[2])

write.table(prsImmuno,file = "/Volumes/MacOS/500fg/500FG/data_for_lars/halla_data/prsImmuno.vs.microbiome/prsImmuno.microbiomePathway.tsv", 
            sep="\t", quote = F,
            row.names = T, col.names = T,
            na = "0")
write.table(path,file = "/Volumes/MacOS/500fg/500FG/data_for_lars/halla_data/prsImmuno.vs.microbiome/microbiomePathway.prsImmuno.tsv", 
            sep="\t", quote = F,
            row.names = T, col.names = T,
            na = "0")


# prs vs immuno mod
prsImmuno <- prs[grep(grepImmunoString,rownames(prs)),]
immunoMod <- read.table("/Volumes/MacOS/500fg/500FG/data_for_lars/20151117resultsTotalReFormat.csv", 
                     row.names = 1, header = T,
                     sep = ",")
immunoMod <- as.data.frame(t(immunoMod))
listOfMatrix <- makeMatricesOverlap(prsImmuno,immunoMod)
prsImmuno <- as.data.frame(listOfMatrix[1])
immunoMod <- as.data.frame(listOfMatrix[2])
immunoMod <- log2(immunoMod)
write.table(prsImmuno,file = "/Volumes/MacOS/500fg/500FG/data_for_lars/halla_data/prsImmuno.vs.immuneModulator/prsImmuno.mimmuneModulator.tsv", 
            sep="\t", quote = F,
            row.names = T, col.names = T,
            na = "0")
write.table(immunoMod,file = "/Volumes/MacOS/500fg/500FG/data_for_lars/halla_data/prsImmuno.vs.immuneModulator/immuneModulator.prsImmuno.tsv", 
            sep="\t", quote = F,
            row.names = T, col.names = T,
            na = "0")

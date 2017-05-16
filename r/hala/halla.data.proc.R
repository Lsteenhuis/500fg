#######################################################################################
# Reading Data
#######################################################################################

prs <- read.table(file = "/Volumes/MacOS/500fg/geneticRiskScore/output.test/TT/rawScoreMatrix.txt", sep = "\t", header = F, stringsAsFactors = F)
meta <- read.table("/Volumes/MacOS/500fg/halla/500Fg_nmr_zeroAsmissing_log2trans_removeOutliers_std.txt", sep = "\t", stringsAsFactors = F, header = F)
cytokine <- read.table(file = "/Volumes/MacOS/500fg/data/pheno_91cytokines_4Raul.csv", stringsAsFactors = F, sep = ",", na.strings = c("NA"))
cytokine <- cytokine[,-2]

#######################################################################################
# Functions
#######################################################################################
getMatrix <- function(df1,df2) {
  df1[1,1] <- "#"
  df2[1,1] <- "#"
  df1[is.na(df1)] <- 0
  df2[is.na(df2)] <- 0
  
  # check if df is in right format 
  if(grepl(pattern="HV\\d+",x= df1[,1], perl =T)[2]){
    df1 <- as.data.frame(t(df1))
  }
  if(grepl(pattern="HV\\d+",x= df2[,1], perl =T)[2]){
    df2 <- as.data.frame(t(df2))
    }
    
    # get subset of both df's 
    df1.samples <- unname(unlist(df1[1,]))
    df2.samples <- unname(unlist(df2[1,]))
    df.intersect <- intersect(df1.samples,df2.samples)
    df1 <- df1[,which(df1.samples %in% df.intersect)]
    df2 <- df2[,which(df2.samples %in% df.intersect)]
    # order the file by HV
    df1 <- df1[,order(df1[1,])]
    df2 <- df2[,order(df2[1,])]
    
    myList <- list("a"=df1,"b"=df2)
    return(myList)
  } 

prepareCt <- function(){
  # read in ct files
  ctInfo <- read.table("/Volumes/MacOS/500fg/data/full_code_sampleInfo.csv", stringsAsFactors=F,sep = ",", header=T)
  ct <- read.table(file="/Volumes/MacOS/500fg/data/IRT_immuneTraits_500FG.csv", header = F,stringsAsFactors = F,sep=",")
  #removing all immuno celltypes which are percentages
  ct <- ct[,-which(ctInfo[,5] == "percentage")]
  # changing IT numbers with cell type names and replacing H with HV so names match with other files
  ct[1,] <- c("#",ctInfo[which(ctInfo[,1] %in% ct[1,]),2])
  ct[,1] <- sub(pattern = "H", replacement = "HV", ct[,1])
  ct <- as.data.frame(t(ct), stringsAsFactors = F)
  ct <- ct[,order(ct[1,])]
  ct[duplicated(ct[,1]),1] <- paste(ct[duplicated(ct[,1]),1],"_2",sep="") 
  ct
}

writeTable <- function(df,name) {
  write.table(df, 
              file=paste("/Volumes/MacOS/500fg/halla/",name, sep=""),
              sep="\t",
              quote = F,
              row.names = F,
              col.names = F,
              na = "0")
}

################################################################################
# Cytokine VS Immuno Cell type
#######################################################################################

ct <- prepareCt()
df.l <- getMatrix(cytokine,ct)


writeTable(df.l[1],"cyto.it.halla.tsv")
writeTable(df.l[2],"it.cyto.halla.tsv")

#######################################################################################
#Monocytes & T Cells VS Immuno Cell Type
#######################################################################################

cyto <- as.data.frame(df.l[1])
monocytes <- as.data.frame(cyto[c(1,grep(pattern = "IL1b|IL6|TNFA", x= cyto[,1], perl = T)) , ])
tcells <- as.data.frame(cyto[c(1,grep(pattern = "IL17|IFNy|IL22", x = cyto[,1], perl = T)) , ])

writeTable(monocytes,"mono.it.halla.tsv")
writeTable(tcells,"tcell.it.halla.tsv")

rm(cyto,monocytes,tcells)      

#######################################################################################
# Immuno Cell Type VS PRS results
#######################################################################################

df.l <- getMatrix(prs,ct) 

writeTable(df.l[1],"prs.it.halla.tsv")
writeTable(df.l[2],"it.prs.halla.tsv")

#######################################################################################
# Immuno Cell Type VS Metabolites
#######################################################################################

df.l <- getMatrix(meta,ct)

writeTable(df.l[1],"prs.diseases.it.halla.tsv")
writeTable(df.l[2],"it.prs.diseases.halla.tsv")

#######################################################################################
# Metabolites VS Cytokine
#######################################################################################

df.l <- getMatrix(meta,cytokine)
writeTable(df.l[1],"meta.cyto.halla.tsv")
writeTable(df.l[2],"cytokine.meta.halla.tsv")

#######################################################################################
# Cytokine VS PRS
#######################################################################################

df.l <- getMatrix(cytokine,prs)
writeTable(df.l[1], "cyto.prs.halla.tsv")
writeTable(df.l[2], "prs.diseases.cyto.halla.tsv")

#######################################################################################
# PRS diseases vs Metabolites
#######################################################################################

df.l <- getMatrix(prs,meta)
grepString <- "gout|disease|Disease|Major_depression|ALS|celiac_disease|Ulcerative_Colitis|celiac_disease|Juvenile_Idiopathic_Arthritis|multiple_sclerosis|Narcolepsy|primary_biliary_cirrhosis|psoriasis|T1D|systemic_lupus_erythematosus|Rheumatoid_Arthritis|Crohns_disease|Ulcerative_colitis|Inflammatory_Bowel_Disease|Asthma|Eczema|Type_2_Diabetes|Chronic_Kidney_Disease"
prsM <- as.data.frame(df.l[1])
prsM <- prsM[c(1,grep(pattern = grepString, x = prsM[,1])),]
writeTable(prsM, "prs.meta.halla.tsv")
writeTable(df.l[2], "meta.prs.halla.tsv")
rm(prsM)

#######################################################################################

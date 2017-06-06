#######################################################################################
# Functions                                                                           #
#######################################################################################
getMatrix <- function(df1,df2) {
  #colnames(df1) <- c("#", colnames(df1))
  #colnames(df2) <- c("#", colnames(df2))
  df1[is.na(df1)] <- 0
  df2[is.na(df2)] <- 0
  
  # check if df is in right format 
  if(grepl(pattern="HV\\d+",x= rownames(df1), perl =T)[2]){
    df1 <- as.data.frame(t(df1))
  }
  if(grepl(pattern="HV\\d+",x= rownames(df2), perl =T)[2]){
    df2 <- as.data.frame(t(df2))
  }
  
  # get subset of both df's 
  df1.samples <- colnames(df1)
  df2.samples <- colnames(df2)
  df.intersect <- intersect(df1.samples,df2.samples)
  df1 <- df1[,which(df1.samples %in% df.intersect)]
  df2 <- df2[,which(df2.samples %in% df.intersect)]
  # order the file by HV
  df1 <- df1[,order(colnames(df1))]
  df2 <- df2[,order(colnames(df2))]
  
  myList <- list(df1,df2)
  return(myList)
} 

prepareCt <- function(){
  # read in ct files
  ctInfo <- read.table("/Volumes/MacOS/500fg/data/full_code_sampleInfo.csv", stringsAsFactors=F,sep = ",", header=T)
  ct <- read.table(file="/Volumes/MacOS/500fg/data/IRT_immuneTraits_500FG.csv", 
                   header = T, row.names = 1,
                   stringsAsFactors = F,sep=",")
  #removing all immuno celltypes which are percentages
  ct <- ct[,-which(ctInfo[,5] == "percentage")]
  # changing IT numbers with cell type names and replacing H with HV so names match with other files
  rownames(ct) <- c("#",ctInfo[which(ctInfo[,1] %in% rownames(ct)),2])
  rownames(ct) <- sub(pattern = "H", replacement = "HV", rownames(ct))
  ct <- as.data.frame(t(ct), stringsAsFactors = F)
  ct <- ct[,order(colnames(ct))]
  rownames(ct)[duplicated(rownames(ct))] <- paste(rownames(ct)[duplicated(rownames(ct))],"_2",sep="")
  ct[duplicated(rownames(ct)),1] <- paste(ct[duplicated(ct[,1]),1],"_2",sep="") 
  ct
}

writeTable <- function(df,name) {
  write.table(df, 
              file=paste("/Volumes/MacOS/500fg/halla/",name, sep=""),
              sep="\t",
              quote = F,
              row.names = T,
              col.names = T,
              na = "0")
}

#######################################################################################
# Reading Data                                                                        #  
#######################################################################################
prs <- as.data.frame(read.table(file = "/Volumes/MacOS/500fg/geneticRiskScore/output.test/TT/rawScoreMatrix.txt", 
                                header = T, row.names = 1,
                                sep = "\t", stringsAsFactors = F))
prsAd <- prs[grep("Alzheimers",rownames(prs),ignore.case = T), ]
prsRa <- prs[grep("Rheumatoid",rownames(prs),ignore.case = T), ]

cytokine <- as.data.frame(read.table(file = "/Volumes/MacOS/500fg/data/pheno_91cytokines_4Raul.csv" ,
                                     header = T, row.names = 1,
                                     stringsAsFactors = F, sep = ","))
cytokine <- cytokine[,-1]

meta <- read.table("/Volumes/MacOS/500fg/halla/500Fg_nmr_zeroAsmissing_log2trans_removeOutliers_std.txt",
                   header = T, row.names = 1,
                   sep = "\t", stringsAsFactors = F)

geneCounts <- as.data.frame(read.table(file = "/Volumes/MacOS/500fg/norm_gene_counts/geneExpression/500FG_normalized_filtered_read_counts.csv",header = T, row.names = 1,sep=",",stringsAsFactors = F))
geneCounts <- as.data.frame(t(geneCounts))

#######################################################################################
# Cytokine VS Immuno Cell type                                                        #
#######################################################################################

ct <- prepareCt()
df.l <- getMatrix(cytokine,ct)


writeTable(df.l[1],"cyto.it.halla.tsv")
writeTable(df.l[2],"it.cyto.halla.tsv")

#######################################################################################
#Monocytes & T Cells VS Immuno Cell Type                                              #
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
# PRS diseases vs Metabolites                                                         #
#######################################################################################

df.l <- getMatrix(prs,meta)
grepString <- "gout|disease|Disease|Major_depression|ALS|celiac_disease|Ulcerative_Colitis|celiac_disease|Juvenile_Idiopathic_Arthritis|multiple_sclerosis|Narcolepsy|primary_biliary_cirrhosis|psoriasis|T1D|systemic_lupus_erythematosus|Rheumatoid_Arthritis|Crohns_disease|Ulcerative_colitis|Inflammatory_Bowel_Disease|Asthma|Eczema|Type_2_Diabetes|Chronic_Kidney_Disease"
prsM <- as.data.frame(df.l[1])
prsM <- prsM[c(1,grep(pattern = grepString, x = prsM[,1])),]
writeTable(prsM, "prs.meta.halla.tsv")
writeTable(df.l[2], "meta.prs.halla.tsv")
rm(prsM)

#######################################################################################
# Cytokine VS Alzheimers Disease                                                      #
#######################################################################################

df.l <- getMatrix(prsAd,cytokine)
writeTable(df.l[1], "prsAd.cytokine.tsv")
writeTable(df.l[2], "cytokine.prsAd.tsv")

#######################################################################################
# Cytokine VS Rheumatoid Arthritis                                                    #
#######################################################################################
dim(prsRa)
dim(cytokine)
df.l <- getMatrix(prsRa,cytokine)
writeTable(df.l[1], "prsRa.cytokine.tsv")
writeTable(df.l[2], "cytokine.prsRa.tsv")

#######################################################################################
# GeneCounts VS PRS                                                                   #
#######################################################################################
grepString <- "gout|disease|Disease|Major_depression|ALS|celiac_disease|Ulcerative_Colitis|celiac_disease|Juvenile_Idiopathic_Arthritis|multiple_sclerosis|Narcolepsy|primary_biliary_cirrhosis|psoriasis|T1D|systemic_lupus_erythematosus|Rheumatoid_Arthritis|Crohns_disease|Ulcerative_colitis|Inflammatory_Bowel_Disease|Asthma|Eczema|Type_2_Diabetes|Chronic_Kidney_Disease"
prs <- prs[grep(grepString, rownames(prs), ignore.case = T), ]
df.l <- getMatrix(prs,geneCounts)
writeTable(df.l[1], "prs.genecounts.tsv")
writeTable(df.l[2], "genecounts.prs.tsv")

################################
# This script is used for the comparison of HALLA figures vs own calculations
################################

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

# function to get correlation from 2 unhomogenous matrices
# xMat: matrix A
# yMat: matrix B
# output: correlation matrix
getCorMat <- function(xMat,yMat){

  corResMat <- apply(xMat,1,function(x){
    corRes <- apply(yMat,1,function(y){
      res <- rcorr(as.matrix(x), as.matrix(y),type = "spearman")
    })
  })
  #colnames(corResMat) <- rownames(xMat)
  #rownames(corResMat) <- rownames(yMat)
  corResMat
  
}

# function to get correlation from 2 unhomogenous matrices
# xMat: matrix A
# yMat: matrix B
# output: p value matrix
getCorPvalMat <- function(xMat,yMat){
  corResMat <- apply(xMat,1,function(x){
    corRes <- apply(yMat,1,function(y){
      res <- rcorr(x, y,
                   type = "spearman")
      res$P[2,1]
    })
  })
  colnames(corResMat) <- rownames(xMat)
  rownames(corResMat) <- rownames(yMat)
  corResMat
}

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

# immuno modulator data
# rows: immuno modulator names
# columns: sample names
# values: immuno modulator levels (need to be normalized)
immunoMod <- read.table("/Volumes/MacOS/500fg/500FG/data_for_lars/20151117resultsTotalReFormat.csv", 
                        row.names = 1, header = T,sep = ",")
immunoMod <- as.data.frame(t(immunoMod))
immunoMod <- log2(immunoMod)

taxomy <- read.table("/Volumes/MacOS/500fg/500FG/data_for_lars/500FG_normalized_raw_taxonomy.csv", 
                     row.names = 1, header = T,
                     sep = ",")
taxomy <- as.data.frame(t(taxomy))

# PRS data
# rows: GWAS trait
# columns: sample names
# values: Risk score
prs <- as.data.frame(read.table(file = "/Volumes/MacOS/500fg/geneticRiskScore/output.test/TT/rawScoreMatrix.txt", 
                                header = T, row.names = 1,
                                sep = "\t", stringsAsFactors = F))
grepImmunoString <- "celiac|Ulcerative_Colitis|Arthritis|sclerosis|cirrhosis|psoriasis|T1D|systemic_lupus_erythematosus|Crohns|Ulcerative_colitis|Inflammatory|Asthma"
prsImmuno <- prs[grep(grepImmunoString,rownames(prs), ignore.case = T),]
####################################################################
# Creates HALLA tables per threshold of each data set to see if correlations hold up between thresholds
# pvalue thresholds
thresholds <- c("P0.01","P0.001","P1.0E-4","P1.0E-5","P5.0E-8")
# vector of dataset names as string
datas <- c("metaData","horm","plat","immunoMod","taxomy")

# for each threshold
lapply(thresholds,function(pval){
  #  set mainDir and subDir for creation of folders
  mainDir ="/Volumes/MacOS/500fg/500FG/data_for_lars/halla_data/"
  subDir = pval
  dir.create(file.path(mainDir, subDir),showWarnings = F)
  # subset the auto immune diseases based on threshold
  pvalSub <- prsImmuno[grep(pval,rownames(prsImmuno)),]
  # for each data set
  lapply(datas,function(dataS){
    listOfMatrices <- makeMatricesOverlap(pvalSub,get(dataS))
    pvalSub <- as.data.frame(listOfMatrices[1])
    dataDf <- as.data.frame(listOfMatrices[2])

    mainDir = paste(mainDir,subDir,"",sep="/")
    print(mainDir)
    dir.create(file.path(mainDir, dataS),showWarnings = F)
    
    
    
    # write the tables in HALLA format
    pvalSub <- cbind(rownames(pvalSub),pvalSub)
    colnames(pvalSub)[1] <- "#"
    write.table(pvalSub,file =paste(mainDir,dataS,"/prs.",pval,".",dataS,".tsv",sep = ""),
                sep = "\t", quote = F,
                row.names = F,col.names = T,
                na="0")
    
    dataDf <- cbind(rownames(dataDf),dataDf)
    colnames(dataDf)[1] <- "#"
    
    write.table(dataDf,file =paste(mainDir,dataS,"/prs.",dataS,".",pval,".tsv",sep = ""),
                sep = "\t", quote = F,
                row.names = F,col.names = T,
                na = "0")
  })
})

####################################################################
# create  heatmaps based on the HALLA input tables to compare them with HALLA result
lapply(thresholds,function(pval){
  lapply(datas,function(dataS){
    # list halla input files based on threshold and data set
    listFile <- list.files(paste("/Volumes/MacOS/500fg/500FG/data_for_lars/halla_data/",pval,"/",dataS,"/", sep = ""), pattern = "tsv$", full.names = T)
    fileX <- read.table(listFile[1], sep="\t",header = T, row.names = 1, comment.char = "$")
    fileY <- read.table(listFile[2], sep="\t",header = T, row.names = 1, comment.char = "$")
    
    cormat <- getCorMat(fileX,fileY)
    dir.create(file.path("/Volumes/MacOS/500fg/500FG/data_for_lars/halla_data/own_heatmap_plots/", pval),showWarnings = F)
    if (dataS == "plat") {
      cormat <- as.data.frame(t(cormat))
    }
    
    #write.table(cormat,paste("/Volumes/MacOS/500fg/500FG/data_for_lars/halla_data/own_heatmap_plots/",pval,"/",dataS,".correlations.tsv",sep=""),quote = F,row.names = T,col.names = T)

    pMat <- getCorPvalMat(fileX,fileY)
    #write.table(pMat,paste("/Volumes/MacOS/500fg/500FG/data_for_lars/halla_data/own_heatmap_plots/",pval,"/",dataS,".pvalues.tsv",sep=""),quote = F,row.names = T,col.names = T)
    # calculate significance of pvalues 
    testP <- sapply(pMat,function(x){
      testP <- getRandomPvalues(pValMatrix = pMat, 
                       nEvents = length(x), 
                       observed = sum(x < 0.05)
                       
        )$pval
    })
    
    # see the direction of -log10 pvalue, so heatmaps are easier to read 
    log10Pval = (sign(cormat)) * -log10(pMat)
    log10Pval <- as.matrix(log10Pval)
    # set insignificant values to 0 so they dont show up in the heatmap (-log(0.05) == 1.3)
    log10Pval[which(abs(log10Pval) < 1.3)] <- 0
    
    # plotting logic
    paletteLength <- 50
    
    #create equal breaks for legend
    if (abs(min(log10Pval)) > max(log10Pval)){
      myBreaks <- c(seq(min(log10Pval), 0, length.out=ceiling(paletteLength/2) + 1), 
                    seq(max(log10Pval)/paletteLength, abs(min(log10Pval)), length.out=floor(paletteLength/2)))
    } else {
      myBreaks <- c(seq(-max(log10Pval), 0, length.out=ceiling(paletteLength/2) + 1), 
                    seq(max(log10Pval)/paletteLength, max(log10Pval), length.out=floor(paletteLength/2)))
    }
    # get unique breaks and color gradient
    myBreaks = unique(myBreaks)
    myColor <- colorRampPalette(c("blue", "white", "red"))(length(myBreaks))
    
    pdf(file = paste("/Volumes/MacOS/500fg/500FG/data_for_lars/halla_data/own_heatmap_plots/",pval,"/",dataS,".pdf",sep=""),onefile = F)
    log10Pval <- log10Pval[, colSums(abs(log10Pval)) != 0]
    log10Pval <- log10Pval[rowSums(abs(log10Pval)) != 0, ]
    if (dataS == "taxomy"){
      log10Pval <- as.data.frame(t(log10Pval))
      pheatmap(t(log10Pval), show_colnames = T, show_rownames = F,
               main = paste(pval, dataS, sep=" "),
               fontsize_row = 2.5 , fontsize_col = 5,
               cluster_cols = T,cluster_rows = T,
               breaks = myBreaks, color = myColor)
    } else if (dataS == "metaData") {
      pheatmap(t(log10Pval), show_colnames = T, show_rownames = T,
               main = paste(pval, dataS, sep=" "),
               fontsize_row = 2.5, fontsize_col = 3,
               cluster_cols = T,cluster_rows = T,
               breaks = myBreaks, color = myColor,
               cellheight = 2.5,cellwidth =3 )
    } else {
      pheatmap(log10Pval, show_colnames = T, show_rownames = T,
               main = paste(pval, dataS, sep=" "),
               fontsize_row = 5, fontsize_col = 5,
               cluster_cols = T,cluster_rows = T,
               breaks = myBreaks, color = myColor)
    }
    dev.off()
  })
})
  ####################################################################
#sim = paste("/Volumes/MacOS/500fg/500FG/data_for_lars/halla_data/",pval,"/",dataS,"/out/similarity_table.txt", sep = "")
# if(file.exists(sim)){
#   simT <- read.table(sim,sep="\t", header=T,row.names = 1, comment.char = "@",stringsAsFactors = F)
#   print(dim(simT))
#   if (nrow(simT) >=2 & ncol(simT) >=2){
#     tryCatch(
#       {
#         pheatmap(simT, show_colnames = T, show_rownames = T,
#                  main = paste(pval, dataS, sep=" "),
#                  fontsize_row = 5, fontsize_col = 5,
#                  cluster_cols = T,cluster_rows = T)
#  
#       }, error = function(e){
#         message("error with a file")
#         message(e)
#       }, finally={
#         message(paste("processed",sim, sep = ""))
#       }
#      )
#   }
# }
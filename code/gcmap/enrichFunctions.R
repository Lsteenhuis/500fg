checkDir <- function(fp){
  if (!file.exists(fp)) {
    dir.create(fp)
  }
}

# creates countSet for gCMAP based on the immuno type which has been supplied
# it: it subset
getItExprMatrix <- function(IT){
  # takes n% of high and low expressed genes from it subset
  percentage= ceiling(nrow(immuneTrait.subset) / 100 * 15)
  itValues <- immuneTrait.subset[,IT]
  immuneTrait.subset.order <- order(itValues,decreasing = T)
  
  lowValues <- rnaSeqSamples[rev(immuneTrait.subset.order)[1:percentage]]
  highValues <- rnaSeqSamples[immuneTrait.subset.order[1:percentage]]
  
  # create subset of the extremes 
  extremeValues <- c(lowValues,highValues)
  ItMatrix <- as.data.frame(geneLevelExpression[,extremeValues])
  varMeta <- matrix(c(rep("Control",15),rep("Case",15)),
                    byrow = T,
                    nrow = nrow(ItMatrix), 
                    ncol= ncol(ItMatrix))
  colnames(varMeta) <- c(paste("control",1:15,sep="_"),paste("case",1:15,sep="_"))
  varMeta <- as.data.frame(varMeta)
  
  annotatedIt <- AnnotatedDataFrame(as.data.frame(t(ItMatrix)), varMeta, dimLabels=c("gene id's","samples"))
  # create a count set with the top n% of gene expression sets
  countSet <- newCountDataSet(ItMatrix, 
                             conditions = c(rep("low",15),rep("high",15)),
                             phenoData = annotatedIt
  )
  countSet
}

# generates a drugset collection based on the column of the drug Table
# and amount of genes which are looked at.
# will be used as a database for gCMAP
generateDrugGeneSetCollection <- function(drugColIndex,perc){
  subsetLength <- nrow(drugTable) * perc
  drugGeneNames <- unique(drugTable[,2])
  drugCol <- drugTable[,drugColIndex]
  drugName <- colnames(drugTable)[drugColIndex]
  return(c(drugGeneSetLo,drugGeneSetHi))
}

# calculates enrichment between high/low expressed genes and 
# high/low drug ranks. Different modes can be used by mode.
# geneDirection: high/low expressed genes
# mode : the mode which will be used for the calculation of enrichment.
useGcmapFunction <- function(geneDirection, mode){
  # checks if cde (nchannelSet) has been loaded, loads if it isn't
  if (!exists("cde")) {
    load("data/gcMAP/nchannelSet")
  }
  perc = ceiling(nrow(cde) * 0.05)
  
  # For each geneDirection each both directions (low and high expressed)
  sapply(1:2,function(setIndex){
    if (setIndex == 1) {
      set = drugGeneLow
      setString = "drugGeneLow"
    } else if (setIndex == 2) {
      set = drugGeneHigh
      setString = "drugGeneHigh"
    }
    print(paste("Current object: ", geneDirection, setString, sep = " "))
    
    # runs the chosen function for each sample in cde
    exp <- sapply(1:114, function(i){
      # Retrieving log_fc changes from cde data. infinite cases and inf values are set to 0
      profile <- assayDataElement(cde[,i], "log_fc")
      profile[!complete.cases(profile)] <- 0
      profile[!is.finite(profile)] <- 0
      
      # switch case for deciding between high expressed DE genes or Low expressed DE genes
      if (geneDirection == "highExpr") {
        prof.ordered <- profile[order(profile,decreasing = T)]
        names(prof.ordered) <- rownames(profile)[order(profile,decreasing = T)]
      } else {
        prof.ordered <- profile[order(profile,decreasing = F)]
        names(prof.ordered) <- rownames(profile)[order(profile,decreasing = F)]
      } 
      # creating query from log_fc changes in prof.ordered
      query <- prof.ordered[1:perc]
      
      # switch case for gcMAP fucntion
      if (mode == "wilcox"){
        exp <- wilcox_score(query, set)
      } else if (mode == "probability") {
        query <- GeneSet(names(query))
        exp <- fisher_score(query = query, sets = set, universe = universe)
      }
      
      # retrieving results from exp object results
      return(exp)
    })
    
    # saving list of GcMapRestult objects
    fileLocation = paste("data/gcMAP/",mode,"/", sep="")
    fileName = paste(mode,geneDirection,setString, sep=".")
    print(paste("saving file: ",fileLocation, fileName, sep = ""))
    save(exp, file= paste(fileLocation, fileName, sep = ""))
    return(exp)
  })
  return(exp)
}

# creates heatmaps from the result files created by useGcmapFunction
createHeatMap <- function(resultFile){
  # resultFile = result file from useGcMapFunction
  load(resultFile)
  
  # creates one matrix from all gCMAP result objects
  resultMatrix <- sapply(1:114,function(drugResult){
    drugResult = exp[[drugResult]]
    resultMatrix <- padj(drugResult)
    resultMatrix
  })
  
  # sets rownames to NA to avoid cluther of 1309 drugnames
  colnames(resultMatrix) <- it_codes
  rows <- rep("",1309)
  drugs <- colnames(drugTable)[3:ncol(drugTable)]
  
  # retrieving location of interesting drugs and set rownames for these locations
  int.drugs <- grep("trihexyphenidyl|propofol|spironolactone|clioquinol|guanfacine", drugs, perl = T)
  rows[int.drugs] <- drugs[int.drugs]
  rownames(resultMatrix) <- rows
  
  breakList = seq(0,(max(resultMatrix) + 1), 1)
  
  # splitting basename of file to determine file name and location
  str <- unlist(strsplit(wFile,"\\/"))
  baseFile <- tail(str,n =1)
  fileSnippets <- unlist(strsplit(baseFile, "\\."))
  tit <- paste("-log10 Padj value", fileSnippets[1], "-", fileSnippets[2], "/", fileSnippets[3], sep = " ")
  fileLoc <- paste("plots/gCMAP",fileSnippets[1], "", sep= "/" )
  pdfName <- paste(fileSnippets[1],fileSnippets[2],fileSnippets[3],"pdf", sep=".")
  
  print("start plot")
  # plotting heatmaps and saving them in location
  pdf(filename = paste(fileLoc, pdfName, sep= ""), onefile = F)
  pheatmap(resultMatrix[int.drugs,], main=tit,
           fontsize= 20, #,fontsize_col = 10, fontsize_row = 5,
           show_colnames = T, show_rownames = T,
           cellheight = 50, bg = "transparent",
           cluster_rows = T, cluster_cols = T,
           color = colorRampPalette(brewer.pal(n = 7, name = "YlOrBr"))(length(breakList)))
  dev.off()
}

# creates a comparison file displaying the results of the probability rank
# and the wilcox score probability score
createComparisonFile <- function(perc,geneDirection,setString) {
  probF <- paste("probability",geneDirection,setString, sep=".")
  probF <- paste("data/gcMAP/probability/",perc,"/",probF, sep="")
  load(probF)
  probFile <- exp
  probMatrix <- createResultMatrix(probFile)
  
  wilF <- paste("wilcox",geneDirection,setString, sep=".")
  wilF <- paste("data/gcMAP/wilcox/",perc,"/",wilF, sep="")
  load(wilF)
  wilFile <- exp
  wilMatrix <- createResultMatrix(wilFile)
  
  sapply(1:nrow(probMatrix), function(i){
    drug <- rownames(probMatrix)[i]
    probP <- probMatrix[i,]
    probRank <- -log(probP, base = 10)
    wilP <- wilMatrix[i,]
    wilRank <- -log(wilP, base = 10)
    compMatrix <- do.call(cbind, list(it_codes,probP,probRank,wilP,wilRank))
    colnames(compMatrix) <- c("IT_Names","Probability Padj","Probability Rank",
                              "Wilcox Padj","Wilcox Rank")
    #rownames(compMatrix) <- c("trihexyphenidyl","propofol",
    #                          "spironolactone","clioquinol","guanfacine")
    
    fileN <- paste("comparison",drug,geneDirection,setString, sep = ".")
    fileP <- paste("data/gcMAP/comparison/",perc,"/",sep="")
    checkDir(fileP)
    write.csv(compMatrix, paste(fileP, fileN, sep = ""),quote = F, row.names = F)
  })
  
}

# creates results matrix
# resultFile: input file 
# output: result
createResultMatrix <- function(resultFile){
  resultMatrix <- sapply(1:114,function(drugResult){
    drugResult = resultFile[[drugResult]]
    resultMatrix <- pval(drugResult)
    # Bonferri multiple testing correction
    resultMatrix <- resultMatrix * 114
  })
  # retrieving location of interesting drugs and set rownames for these locations
  drugs <- colnames(drugTable)[3:ncol(drugTable)]
  int.drugs <- grep("trihexyphenidyl|propofol|spironolactone|clioquinol|guanfacine", rownames(resultMatrix), perl = T)
  resultMatrix <- resultMatrix[int.drugs,]
  colnames(resultMatrix) <- it_codes
  resultMatrix
}


getItExprMatrix <- function(IT){
  percentage= ceiling(nrow(immuneTrait.subset) / 100 * 15)
  itValues <- immuneTrait.subset[,IT]
  immuneTrait.subset.order <- order(itValues,decreasing = T)
  
  lowValues <- rnaSeqSamples[rev(immuneTrait.subset.order)[1:percentage]]
  highValues <- rnaSeqSamples[immuneTrait.subset.order[1:percentage]]
  
  extremeValues <- c(lowValues,highValues)
  ItMatrix <- as.data.frame(geneLevelExpression[,extremeValues])
  varMeta <- matrix(c(rep("Control",15),rep("Case",15)),
                    byrow = T,
                    nrow = nrow(ItMatrix), 
                    ncol= ncol(ItMatrix))
  colnames(varMeta) <- c(paste("control",1:15,sep="_"),paste("case",1:15,sep="_"))
  varMeta <- as.data.frame(varMeta)
  
  annotatedIt <- AnnotatedDataFrame(as.data.frame(t(ItMatrix)), varMeta, dimLabels=c("gene id's","samples"))
  
  countSet <- newCountDataSet(ItMatrix, 
                             conditions = c(rep("low",15),rep("high",15)),
                             phenoData = annotatedIt
  )
  countSet
}


generateDrugGeneSetCollection <- function(drugColIndex){
  subsetLength <- nrow(drugTable) * 0.05
  drugCol <- drugTable[,drugColIndex]
  drugName <- colnames(drugTable)[drugColIndex]
  drug.sorted <- order(drugCol)
  
  drugGenesHi <- drugTable[,2][drug.sorted[0:subsetLength]]
  drugGenesLo <- drugTable[,2][rev(drug.sorted)[0:subsetLength]]
  
  drugGeneSetHi <- GeneSet(unique(drugGenesHi),setIdentifier = paste(drugName,"Hi",sep="_"),
                           setName =paste(drugName,"High",sep=""))
  drugGeneSetLo <- GeneSet(unique(drugGenesLo),setIdentifier = paste(drugName,"Lo",sep="_"),
                           setName =paste(drugName,"Low",sep=""))
  
  return(c(drugGeneSetLo,drugGeneSetHi))
}

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
  
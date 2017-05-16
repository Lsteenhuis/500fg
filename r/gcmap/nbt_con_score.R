load("data/gcMAP/nchannelSet")
source("data/nbt.3367-S3/CMapFxns.R")
drugTable <- read.csv("data/drugs.ens.csv", stringsAsFactors = F)
drugTable <- drugTable[isUnique(drugTable[,2]), ]
drugGenes <- drugTable[,2]


system.time(scores <- sapply(3:4,function(drugIndex){
  deg_set_len <- 250
  
  # create drug signature of metformin
  names(drugCol) <- drugGenes
  drugRank <- drugCol[which(shared %in% names(drugCol))]
  drug_signature <- matrix(nrow=length(shared),ncol =2)
  colnames(drug_signature) <- c("GeneID","value")
  drug_signature[,1] <- shared
  drug_signature[,2] <- drugRank
  drug_signature <- as.data.frame(drug_signature)
  rSc = rand_cmap_scores(drug_sig=drug_signature, m_genes=shared,
  deg_set_len=deg_set_len, nTrials=nTrials)
  
  scores <- sapply(1:114,function(i){
    # Read drug table, get unique gene ids
   drugCol <- drugTable[,drugIndex]
    # read the log_fc  changes from first IT
    profile <- assayDataElement(cde[,i], "z")
    profile[!complete.cases(profile)] <- 0
    profile[!is.finite(profile)] <- 0
    
    # get shared drug genes
    shared <- intersect(drugGenes,rownames(profile))
    
    # create it_signature of shared genes. 
    rowP <- rownames(profile)[which(shared %in% rownames(profile))]
    profile <- profile[which(shared %in% rownames(profile))]
    cell_signature <- data.frame(GeneID=rowP,value=profile,stringsAsFactors = F)
    cell_signature = cell_signature[order(cell_signature$value, 
                                          decreasing=TRUE), ]
    
    # get up and down regulated genes
    
    cell_up = cell_signature$GeneID[1:deg_set_len]
    gsel = (nrow(cell_signature) - deg_set_len + 1):nrow(cell_signature)
    cell_down = cell_signature$GeneID[gsel]
    

    
    # call score function from NBT
    score = cmap_score(cell_up, cell_down, drug_signature)
    if (score != 0) {
      score_pvalue = length(which(abs(rSc) >= abs(score)))  / nTrials
    }
    list(score=score, pvalue=score_pvalue)
  })
  cell_drug_scores = do.call(rbind, scores)
  ofile = paste("data/gcMAP/nbt_scores/D_", deg_set_len, "/", i, ".csv", sep="")
  write.csv(cell_drug_scores, file=ofile)
})
)
length(scores)

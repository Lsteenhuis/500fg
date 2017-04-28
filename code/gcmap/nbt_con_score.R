source("data/nbt.3367-S3/CMapFxns.R")
drugTable <- read.csv("data/drugs.ens.csv", stringsAsFactors = F)
drugTable <- drugTable[isUnique(drugTable[,2]), ]
drugGenes <- unique(drugTable[,2])


scores <- sapply(3:1309,2,function(drugCol){
  scores <- sapply(1:114,function(i){
    # Read drug table, get unique gene ids
   
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
    deg_set_len <- 250
    cell_up = cell_signature$GeneID[1:deg_set_len]
    gsel = (nrow(cell_signature) - deg_set_len + 1):nrow(cell_signature)
    cell_down = cell_signature$GeneID[gsel]
    
    # create drug signature of metformin
    names(drugCol) <- drugTable[,2]
    drugRank <- drugCol[which(shared %in% names(drugCol))]
    drug_signature <- matrix(nrow=length(shared),ncol =2)
    colnames(drug_signature) <- c("GeneID","value")
    drug_signature[,1] <- shared
    drug_signature[,2] <- drugRank
    drug_signature <- as.data.frame(drug_signature)
    
    # call score function from NBT
    score = cmap_score(cell_up, cell_down, drug_signature)
  })
})
length(scores)

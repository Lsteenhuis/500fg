#old
system.time({
  # for each PRS disease 
  as <- foreach (prsIndex = 1:100) %dopar% {
    # for each GWAS trait 
    as <- sapply(1:nrow(cytokineType),function(cytIndex){
      # creates a matrix with PRS and cytokine levels based on overlapping samples
      corMatrix <- cbind.data.frame(
        unlist(prsVal[prsIndex,]) ,
        unlist(cytokineType[cytIndex,])
      )
      # calculate correlation between PRS score and cytokine levels
      res <- rcorr(corMatrix[,1],corMatrix[,2], type = "spearman")
      #save P value, rho value, name of the cytokine, and prs name
      c(res$P[2,1],res$r[2,1],cytoNames[cytIndex],prsNames[prsIndex])
    })
    as <- t(as)
    colnames(as) <- c("pval","rho","celltype","gwas")
    as.data.frame(as,stringsAsFactors = F)
  }  
  #correlation.df <- rbind.fill(as)
})

#new 1
system.time({
  # for each PRS disease 
  as <- foreach (prsIndex = 1:100) %dopar% {
    # for each GWAS trait 
    as <- apply(monocytes,1,function(cytIndex){
      # creates a matrix with PRS and cytokine levels based on overlapping samples
      corMatrix <- cbind.data.frame(
        unlist(prsVal[prsIndex,]) ,
        unlist(cytIndex[2:ncol(monocytes)])
      )
      # calculate correlation between PRS score and cytokine levels
      res <- rcorr(corMatrix[,1],corMatrix[,2], type = "spearman")
      #save P value, rho value, name of the cytokine, and prs name
      c(res$P[2,1],res$r[2,1],cytIndex[1],prsNames[prsIndex])
    })
    as <- t(as)
    colnames(as) <- c("pval","rho","celltype","gwas")
    as.data.frame(as,stringsAsFactors = F)
  }  
  #correlation.df <- rbind.fill(as)
})

#new 2 Never use this wow so slow much ineffectiveness 
system.time({
  # for each PRS disease 
  as <- apply(prs,1,function(prsIndex){
    # for each GWAS trait 
    as <- apply(monocytes,1,function(cytIndex){
      # creates a matrix with PRS and cytokine levels based on overlapping samples
      corMatrix <- cbind.data.frame(
        unlist(prsIndex[2:ncol(prs)]) ,
        unlist(cytIndex[2:ncol(monocytes)])
      )
      # calculate correlation between PRS score and cytokine levels
      res <- rcorr(corMatrix[,1],corMatrix[,2], type = "spearman")
      #save P value, rho value, name of the cytokine, and prs name
      c(res$P[2,1],res$r[2,1],cytIndex[1],prsIndex[1])
    })
    as <- t(as)
    colnames(as) <- c("pval","rho","celltype","gwas")
    as.data.frame(as,stringsAsFactors = F)
  })
  #correlation.df <- rbind.fill(as)
})

cl <- makeCluster(mc <- getOption("cl.cores", 2))
cl <- makeCluster(4)
system.time({
  # for each PRS disease 
  as <- parRapply(cl,x = prs[1:100,] ,testF,monocytes,prs)
})
stopCluster(cl)

suppressPackageStartupMessages(library(Hmisc))

#correlation.df <- rbind.fill(as)
new1.df <- rbind.fill(as)

new.3df <- rbind.fill(as)
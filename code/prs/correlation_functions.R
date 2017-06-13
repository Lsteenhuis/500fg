# function for the correlation analysis
# returns list object with correlation scores and p values
# xmat <- data matrix X
# ymat <- data matrix Y
getCorMat <- function(xMat,yMat){
  
  if(length(grep("HV", rownames(xMat)) > 1)){
    xMat <- as.data.frame(t(xMat))
  }
  if(length(grep("HV", rownames(yMat)) > 1)){
    yMat <- as.data.frame(t(yMat))
  }
  xSamples <- colnames(xMat)
  ySamples <- colnames(yMat)
  
  intersectSample <- intersect(xSamples,ySamples)
  
  xMat <- xMat[,which(colnames(xMat) %in% intersectSample)]
  yMat <- yMat[,which(colnames(yMat) %in% intersectSample)]
  xMat <- xMat[,order(colnames(xMat))]
  yMat <- yMat[,order(colnames(yMat))]
  corRMat <- apply(xMat,1,function(x){
    corRes <- apply(yMat,1,function(y){
    
      res <- cor.test(as.numeric(x),
                      as.numeric(y),
                      method = "spearman")
      rRes <- res$estimate
    })
  })
  corPMat <- apply(xMat,1,function(x){
    corRes <- apply(yMat,1,function(y){
      res <- cor.test(as.matrix(as.numeric(x)),
                      as.matrix(as.numeric(y)),
                      method = "spearman")
      rRes <- res$p.value
    })
  })
  
  return(list(corPMat,corRMat))
  
}

# random permutation test
getRandomPvalues <- function(xMat,ySample, perm = 1000, observed, significanceThreshold=0.05){
  dist <- unlist(lapply(1:perm, FUN=function(x){
    xSample <- sample(xMat, size =430,replace = F)
    xSample <- xSample[sample(1:nrow(xMat), replace = F, size =1 ),]
    res <- cor.test(as.matrix(as.numeric(xSample)),
                    as.matrix(as.numeric(y)),
                    method = "spearman")$estimate
  }))
  pval<- sum(dist > rRes) / perm
  pval <- p.adjust(pval, method="BH")
}

# function for the creation of heatmaps with ggplot
createHeatMap <- function(cormod,xLab=NULL,yLab=NULL,lLab=NULL){
  minV <- min(cormod$value,na.rm = T)
  maxV <- max(cormod$value, na.rm = T)
  limitValue <- ifelse(abs(minV) > maxV,
                       floor(minV),
                       -ceiling(maxV)
  )
  myBreaks <- seq(limitValue,abs(limitValue), by = 1)
  
  p <- ggplot(na.omit(cormod), aes(Var2, Var1)) + 
    geom_tile(aes(fill = value), color="black") + 
    scale_fill_gradient2(low = "darkblue", mid = "white",high = "darkred",
                         limits = c(-5,5),
                         name=lLab,
                         breaks = myBreaks,
                         midpoint = 0, space = "Lab",
                         na.value = "grey50", guide = "colourbar")+
    xlab(xLab) +
    ylab(yLab) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  p
  
}


# creates abbreviations for diseases for the plots
getAbbr <- function(ibdDis){
  rownames(ibdDis) <- gsub(pattern = "Coronary_artery_disease_mi_additive_2015_26343387_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "CAD_mi_2015_add" )
  rownames(ibdDis) <- gsub(pattern = "Coronary_artery_disease_additive_2015_26343387_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "CAD_2015_add" )
  rownames(ibdDis) <- gsub(pattern = "Coronary_artery_disease_2011_21378988_hg18_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "CAD_2011_add" )
  rownames(ibdDis) <- gsub(pattern = "Coronary_artery_disease_recessive_2015_26343387_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "CAD_2015_rec" )
  rownames(ibdDis) <- gsub(pattern = "Coronary_artery_disease_2011_21378990_hg18_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "CAD_2011" )
  rownames(ibdDis) <- gsub(pattern = "Inflammatory_Bowel_Disease_EUR_2015_26192919_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "IBD" )
  rownames(ibdDis) <- gsub(pattern = "Crohns_disease_EUR_2015_26192919_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "Crohn" )
  
  rownames(ibdDis) <- gsub(pattern = "Ulcerative_Colitis_2011_21297633_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "UC_2011" )
  rownames(ibdDis) <- gsub(pattern = "Ulcerative_colitis_EUR_2015_26192919_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "UC_2015" )
  
  rownames(ibdDis) <- gsub(pattern = "Coronary_artery_disease_mi_additive_2015_26343387_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "CAD_mi_2015_add" )
  rownames(ibdDis) <- gsub(pattern = "Coronary_artery_disease_additive_2015_26343387_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "CAD_2015_add" )
  rownames(ibdDis) <- gsub(pattern = "Coronary_artery_disease_2011_21378988_hg18_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "CAD_2011_add" )
  rownames(ibdDis) <- gsub(pattern = "Coronary_artery_disease_recessive_2015_26343387_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "CAD_2015_rec" )
  rownames(ibdDis) <- gsub(pattern = "Coronary_artery_disease_2011_21378990_hg18_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "CAD_2011" )
  
  rownames(ibdDis) <- gsub(pattern = "T1D_meta_2015_25751624_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "T1D_meta" )
  rownames(ibdDis) <- gsub(pattern = "T1D_CC_2015_25751624_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "T1D_cc" )
  rownames(ibdDis) <- gsub(pattern = "Rheumatoid_Arthritis_2014_24390342_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "RA_2014" )
  rownames(ibdDis) <- gsub(pattern = "Rheumatoid_Arthritis_2010_20453842_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "RA_2010" )
  rownames(ibdDis) <- gsub(pattern = "Juvenile_Idiopathic_Arthritis_2013_23603761_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "JIA" )
  rownames(ibdDis) <- gsub(pattern = "Asthma_2010_860503_fixed_effects_hg18_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "Asthma_FE" )
  rownames(ibdDis) <- gsub(pattern = "Asthma_2010_860503_random_effects_hg18_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "Ashtma_RE" )
  rownames(ibdDis) <- gsub(pattern = "Primary_biliary_cirrhosis_2012_22961000_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "PBC_2012" )
  rownames(ibdDis) <- gsub(pattern = "Psoriasis_2012_23143594_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "Psoiasis_2012" )
  rownames(ibdDis) <- gsub(pattern = "Multiple_sclerosis_2013_24076602_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "MS_2013" )
  rownames(ibdDis) <- gsub(pattern = "Asthma_2010_860503_random_effects_hg18_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "Ashtma_RE" )
  rownames(ibdDis) <- gsub(pattern = "Primary_biliary_cirrhosis_2015_26394269_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "PBC_2015" )
  rownames(ibdDis) <- gsub(pattern = "Psoriasis_2012_23143594_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "Psoiasis_2012" )
  rownames(ibdDis) <- gsub(pattern = "Crohns_disease_EUR_2015_26192919_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "Crohns_2015" )
  rownames(ibdDis) <- gsub(pattern = "Celiac_disease_2011_22057235_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "Celiac_2011" )
  rownames(ibdDis) <- gsub(pattern = "Celiac_disease_2010_20190752_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "Celiac_2010" )
  rownames(ibdDis) <- gsub(pattern = "Multiple_sclerosis_2011_21833088_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "MS_2011" )
  rownames(ibdDis) <- gsub(pattern = "Inflammatory_Bowel_Disease_EUR_2015_26192919_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "IBD_2015" )
  rownames(ibdDis) <- gsub(pattern = "Systemic_lupus_erythematosus_2015_26502338_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "sle_2015" )
          
  rownames(ibdDis) <- gsub(pattern = "ALS_MLMA_2016_27455348_hg19.txt.gz_P5.0E-8",x = rownames(ibdDis), replacement = "ALS" )
  rownames(ibdDis) <- gsub(pattern = "Major_depression_2016_27479909_hg19.txt.gz_P5.0E.8",x = rownames(ibdDis), replacement = "depression" )
  rownames(ibdDis) <- gsub(pattern = "Alzheimers_disease_2013_24162737_hg19.txt.gz_P5.0E-8",x = rownames(ibdDis), replacement = "Alzheimer" )
  rownames(ibdDis) <- gsub(pattern = "Type_2_Diabetes_GWAS_metabo_2012_22885922_hg18_hg19.txt.gz_P5.0E-8",x = rownames(ibdDis), replacement = "T2D_GWAS" )
  rownames(ibdDis) <- gsub(pattern = "Type_2_Diabetes_2012_22885922_hg18_hg19.txt.gz_P5.0E-8",x = rownames(ibdDis), replacement = "T2D_2012" )
  rownames(ibdDis) <- gsub(pattern = "Type_2_Diabetes_metabo_2015_26551672_hg19.txt.gz_P5.0E-8",x = rownames(ibdDis), replacement = "T2D_2015" )
  rownames(ibdDis) <- gsub(pattern = "Type_2_Diabetes_2014_24509480_hg18_hg19.txt.gz_P5.0E-8",x = rownames(ibdDis), replacement = "T2D_2014" )
  
  ibdDis
}

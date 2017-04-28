setwd("~/git/500fg/")
geneFiles <- list.files("deseq_results/", pattern=".csv$",full.names = T)
drugData <- read.csv("matrices/hsa_ensembl/ensembl.drug.data.csv", stringsAsFactors = F)
deSeqR <-  read.csv("deseq_results/deseq.IT1.results.csv", stringsAsFactors = F,header=F)

###
# Names are derived by the example provided by ?phyper.
#
# deSeqGene -> genes that are found to be DE by DESeq2 (filtered by extremes of IT1)
# balls -> total number of human genes in CRCh37.p13
# ballsDrawn -> Number of genes that are found by DESeq2
# whiteBalls -> genes which are derived from the drug data base
# whiteBalls.drawn -> genes from deSeqGene which are also found in the drug data base
# blackballs -> number of total genes - number of genes found in the drug data base
##

#deseqGene <- deSeqR[,1]
balls = 64162
ballsDrawn = 2617
whiteBalls = 14281
whiteBalls.drawn = 1711
blackBalls = balls - whiteBalls
phyper( whiteBalls.drawn, whiteBalls, blackBalls, ballsDrawn )

#deseqGene <- deSeqR[,1]
# balls = 64162
# ballsDrawn = length(deseqGene)
# whiteBalls = length(unique(drugData[,2]))
# whiteBalls.drawn = length( intersect(deseqGene, unique(drugData[,2])) )
# blackBalls = balls - whiteBalls
# phyper( whiteBalls.drawn, whiteBalls, blackBalls, ballsDrawn )

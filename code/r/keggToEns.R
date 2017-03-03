#!/usr/bin/env Rscript
#####
# Date: 02-03-2017
# Author: Lars Steenhuis
# Script which adds a column with ensembl id's to supplementary table 3.
# also described how to retrieve the ensembl id list from supp table 3.
#####
setwd("/Volumes/MacOS/500fg/")
st = read.csv("data/SupplementaryTable3.csv", header=T)

# code to retrieve ensembl genes from kegg id's 
# keggId = st[,1]
# lapply(keggId, function(x) write.table( data.frame(x,stringsAsFactors = F)[1],
#                                        'data/test.csv'  , append= T, sep='\n',quote = F,row.names = F,col.names = F))
# convert output file via http://biodb.jp/
# grep -vwE "\d+\s+LRG_\d+" kegg.ensembl.tsv > keggtest.tsv
# ^ removes LRG lines from file
# cat kegg.ensembl.tsv | grep ENSG > jemoeke.txt
# removes lines which do not have an ensembl id

nt = read.csv("data/jemoeke.txt", sep="\t")

drugTable = apply(nt,1,function(x){
  return(st[st[,1]==as.numeric(x[1]),-1])
})
ensAdded = cbind(nt, matrix(unlist(drugTable), byrow=T, nrow=nrow(nt)))
colnames(ensAdded) = c(colnames(st)[1], "ensID", colnames(st)[-1])

write.csv("data/drugs.ens.csv", x = ensAdded, row.names = F,quote = F)

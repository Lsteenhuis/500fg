library(pheatmap)
library(RColorBrewer)
simTab <- read.table("/Volumes/MacOS/500fg/halla/prs.diseases.cyto.fdr.0.5/similarity_table.txt", sep="\t",header = T,
                     comment.char = "@",row.names = 1)
simTab[is.na(simTab)] <- 0

breakList = seq(-0.25,0.25, 0.01)

pheatmap(simTab,
         show_rownames = T,show_colnames = T,
         fontsize_row = 5,fontsize_col = 5,      
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breakList)),
         treeheight_row = 0,treeheight_col = 0,
         legend_labels = c("-0.2","-0.15","-0.1","-0.05","0","0.05", "-0.1","-0.15","-0.2")
         )


nonImmune <- c("Major_depression|Chronic|Alzheimers|Coronary")

simTab <- simTab[-grep(nonImmune,rownames(simTab)),] 

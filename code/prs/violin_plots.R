# script to create violin plots of the distribution of the pvalues from metabolites
library(RColorBrewer)
library(ggplot2)
library(grid)
################
# functions
################
# Function to plot multiple ggplots next to each other
# plotlist: list of plots
# filelist of files to plot
# cols: amount of columns in the plot
#layout: custom layout
# output: prints ggplots
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
# function to create violin plots from metabolite data
# grepString: string for the disease(s) which need to be included in the plot
# main: title of plot
# returns ggplot 
createViolin <- function(grepString,main){
  
  prsImmuno <- prs[grep(grepString,rownames(prs), ignore.case = T),]
  interSamples <- intersect(colnames(metaData), colnames(prsImmuno))
  prsImmuno <- prsImmuno[,interSamples]
  prsImmuno <- prsImmuno[,order(interSamples)]
  metaData <- metaData[,interSamples]
  metaData <- metaData[,order(interSamples)]
  
  corMat <- getCorMat(prsImmuno,metaData)
  # groups which will be included in the violin plot
  d1 ="VLDL"
  name1=rownames(corMat)[grep(d1, rownames(corMat))]
  d1="LDL"
  name2=rownames(corMat)[grep(d1, rownames(corMat))]
  name2 = name21[-grep("VLDL",name22)]
  d1="HDL"
  name3=rownames(corMat)[grep(d1, rownames(corMat))]
  

  # create on large vector of group values
  waardes <- c(as.vector(corMat[name1,]),
               as.vector(corMat[name2,]),
               as.vector(corMat[name3,]))
  
  # create vector based on group position and length
  # creates 2 dim df for violin plot to group based on variable 
  violin.df<-data.frame(
    variable=c(rep("VLDL",length(as.vector(corMat[name11,]))),
               rep("LDL",length(as.vector(corMat[name22,]))),
               rep("HDL",length(as.vector(corMat[name44,])))),
    value=waardes
  )
  
  violin.df <- as.data.frame(violin.df)
  color1 <-brewer.pal(7, "Paired")[1:nrow(violin.df)]
  colnames(violin.df) <- c("variable","value")
  p <- ggplot(violin.df, aes(x = variable, y = value, fill=variable),names=names(df.m)) + geom_violin(trim=FALSE)+
    labs(title=main,x="Metabolite", y = "Correlation Coefficients")+geom_boxplot(width=0.2)+ theme(legend.position="none")
  p
}
###################################################
# loading data
###################################################
metaData <- read.table("/Volumes/MacOS/500fg/500FG/data_for_lars/500FG_normalized_raw_metabolome.csv", 
                       row.names = 1, header = T, 
                       sep = ",")
metaData <- as.data.frame(t(metaData))
metaData <- metaData[grep("DL",rownames(metaData),ignore.case = T),]
metaData <- metaData[-grep("percent",rownames(metaData),ignore.case = T),]

prs <- as.data.frame(read.table(file = "/Volumes/MacOS/500fg/geneticRiskScore/output.test/TT/rawScoreMatrix.txt", 
                                header = T, row.names = 1,
                                sep = "\t", stringsAsFactors = F))
###################################################
par(mfrow=c(2,2))
rap <- createViolin("Rheumatoid_art", main="RA metabolites")
alp <- createViolin("Asthma", main="asthma metabolites")
multiplot(rap,alp,cols=2)




## ------------------------------------------------------------------------
plot.base <- function(p, ncol=1, xlab=NULL, ylab=NULL, main=NULL, horiz=T,
                      legend.title=NULL, legend.position="bottom",
                      legend.size=0.05, legend=T, legend.cex=1, label.cex=1,
                      axis.cex=1, main.cex=1, facet.cex=1, cols=NULL,
                      x.blank=F, y.blank=F, grid="xy",xlim=NULL, ylim=NULL) {
  # Flip x and y labels because ggplot flips them
  if (!horiz) {
    x.tmp   <- xlab
    xlab    <- ylab
    ylab    <- x.tmp
  }
  
  if (is.null(cols)) {
    cols <- colorRampPalette(c("#994d4d","#ffd9bf","#ff8800","#402200","#ffd940","#8c7723",
                               "#807d60","#74cc66","#00ff00","#005947","#40fff2","#00ccff",
                               "#739199","#0d2133","#b6c6f2","#23318c","#d9001d","#9173e6",
                               "#f200a2","#e6acd2","#59163a","#d9001d"))(ncol)}
  
  p <- p +  labs(x=xlab, y=ylab) +
    ggtitle(main) +
    #scale_fill_manual(values = cols, guide=legend) +
    theme(legend.position=legend.position,
          axis.title=element_text(size=rel(label.cex)),
          plot.title=element_text(size=rel(main.cex), hjust=0.5, vjust=0.5, face="bold"),
          strip.text=element_text(size=rel(facet.cex), face="bold"),
          axis.text=element_text(size=rel(axis.cex)),
          legend.text=element_text(size=rel(axis.cex)),
          legend.key.size=unit(x=legend.size, unit="npc"),
          legend.title=element_text(size=rel(label.cex)))
  
  if (grid=="x") {
    p <- p + theme(panel.grid.major.y=element_blank())
  } else if (grid=="y") {
    p <- p + theme(panel.grid.major.x=element_blank())
  }
  
  if (y.blank) {p <- p + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())}
  if (x.blank) {p <- p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())}
  if (horiz)   {p <- p + coord_flip(xlim=xlim, ylim=ylim)}
  legend=F
  if (legend)  {p <- p + guides(colour=guide_legend(title=legend.title), fill=guide_legend(title=legend.title))}
  
  return(p)
}

## ------------------------------------------------------------------------
bp.grouped <- function(data, group.x=NULL, group.x2=NULL, group.y=NULL, group.y2=NULL, id=NULL,
                       cat=NULL, grouped=T, scales="free_y", space="free", ...) {
  
  data.vec          <- as.vector(as.matrix(data))
  n.group           <- ncol(data)
  n.samp            <- nrow(data)
  if (is.null(group.x)){
    group.x         <- sapply(colnames(data), function(name) {return(rep(name, n.samp))})
    group.x         <- as.character(group.x)
  }
  if (is.null(id)){id <- rep(rownames(data), n.group)}
  df.plot           <- data.frame(row.names = 1:length(data.vec))
  df.plot$data      <- data.vec
  df.plot$id        <- factor(id, levels=unique(id))
  df.plot$cat       <- factor(id, levels=unique(id))
  df.plot$group.x   <- factor(group.x, levels=unique(group.x))
  df.plot$sign      <- ifelse(sign(df.plot$data) == 1, "positive", "negative")
  df.plot           <- df.plot[order(df.plot$group.x),]
  if (!is.null(group.y)){df.plot$group.y   <- factor(group.y, levels=unique(group.y))}
  if (!is.null(group.y2)){df.plot$group.y2 <- factor(group.y2, levels=unique(group.y2))}
  if (!is.null(group.x2)){df.plot$group.x2 <- factor(group.x2, levels=unique(group.x2))}
  
  if (!is.null(cat)){df.plot$cat           <- factor(cat, levels=unique(cat))}
  
  p <- ggplot(df.plot, aes(x=cat, y=data, fill=sign)) +  geom_bar(stat="identity")
  p <- plot.base(p, ncol=length(unique(df.plot$cat)), ...)
  
  if (grouped) {
    if (!is.null(group.x) && !is.null(group.y) && is.null(group.x2) && is.null(group.y2)){
      
      p <- p + facet_grid(group.y ~ group.x, scales=scales, space=space)
      
    } else if (!is.null(group.x) && !is.null(group.y) && !is.null(group.y2)){
      
      p <- p + facet_grid(group.y2 + group.y ~ group.x, scales=scales, space=space)
      
    } else if (!is.null(group.x) && !is.null(group.x2) && !is.null(group.y)&& !is.null(group.y2)){
      
      p <- p + facet_grid(group.y2 + group.y ~ group.x2 + group.x, scales=scales, space=space)
      
    } else if (!is.null(group.x) && !is.null(group.x2) && !is.null(group.y)){
      
      p <- p + facet_grid(group.y ~ group.x2 + group.x, scales=scales, space=space)
      
    } else if(!is.null(group.x)){
      
      p <- p + facet_grid(. ~ group.x, scales=scales, space=space)
      
    } else if (!is.null(group.y)) {
      
      p <- p + facet_grid(group.y ~ ., scales=scales, space=space)
    }
  }
  return(p)
}


## ------------------------------------------------------------------------
hm <- function(data, main=NULL, dendrogram="none", cols=NULL, symkey=F, symbreaks=F,  ...) {
  
  if (is.null(cols)) {
    cols <- colorRampPalette(c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603", "#7F2704"))(100)
  }
  
  if (dendrogram == "none") {
    lmat <- rbind(c(6,6,6), c(5, 4, 2), c(6, 1, 3), c(6, 6, 6))
    lhei <- c(1, 2.5, 5, 1)
    lwid <- c(1, 10, 3)
  } else if (dendrogram == "both") {
    lmat <- NULL
    lhei <- NULL
    lwid <- NULL
  } else if (dendrogram == "column") {
    lmat <- rbind(c(6,6,6), c(2, 3, 4), c(5, 1, 5), c(6, 6, 6))
    lhei <- c(1, 2.5, 5, 1)
    lwid <- c(1, 10, 3)
  } else if (dendrogram == "row") {
    lmat <- rbind(c(6,6,6), c(5, 4, 5), c(2, 1, 3), c(6, 6, 6))
    lhei <- c(1, 1.5, 5, 1)
    lwid <- c(3, 10, 3)
  }
  
  heatmap.2(data,
            scale="none",
            col=cols,
            trace='none',
            symkey=symkey,
            symbreaks=symbreaks,
            dendrogram=dendrogram,
            density.info='histogram',
            denscol="black",
            keysize=1,
            #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
            key.par=list(mar=c(3.5,3.5,3,4)),
            # lmat -- added 2 lattice sections (5 and 6) for padding
            lmat=lmat,
            lhei=lhei,
            lwid=lwid,
            cex.lab=1.5,
            srtCol=45,
            ...
  )
  
  title(main, line=1)
}



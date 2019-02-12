# get colors

GetColor <- function(method,score,color,score.values,n,q){
  switch(method,
         
         "confidence" = {
           color.value <- "black"
           switch (score,
                   "1" = {color.value = "grey"},
                   "2" = {color.value = "royalblue1"},
                   "3" = {color.value = "royalblue3"},
                   "4" = {color.value = "royalblue4"}
           )
         },
         
         "ploidy" = {
           # need redo
           if(score==0){
             color.value = color[1]
           }else{
             color.value = color[score]
           }
         },
         
         "ploidy2" = {
           # need redo
           if(missing(score)){
             color.value <- color[length(score.values)+1]
           }else{
             color.value <- color[score]
           }
         },
         
         "cohort" = {
           if(missing(color)){
             cohort.size <- length(unique(cohorts))
             color1 <- colorRampPalette(c("red2","indianred4","royalblue4","steelblue1","chartreuse3","darkgreen"))(cohort.size)
           }else{
             switch (color,
                     "1" = {rainbow(cohort.size)},
                     "2" = {color1 <- colorRampPalette(c("seashell2","seagreen2","turquoise2","palevioletred2"))(cohort.size)},
                     "3" = {color1 <- colorRampPalette(c("gray7","mediumblue","deeppink4","sienna3"))(cohort.size)},
                     "4" = {color1 <- colorRampPalette(c("darkslateblue","darkslategray4","deeppink4","tan4","gray10"))(cohort.size)},
                     "5" = {color1 <- colorRampPalette(c("navajowhite3","orange3","olivedrab3"))(cohort.size)},
                     "6" = {color1 <- colorRampPalette(c("royalblue2","yellow1"))(cohort.size)},
                     "7" = {color1 <- gray.colors(cohort.size)}
             )
           }
         
         color <- palette(color1)
         if(q==TRUE){return(color1)}else{return(color)}
         },
        
         "length" = {
           color.value <- color
         },
         
         "factor" = {
           color.value <- "black"
         }
  )
  return(color.value)
}

# sorting methods

# define class CNV and twinCNV

CNV_single = setClass("CNV_single",
                      slots = list(
                        name = "character",
                        matrix = "GenomicRanges",
                        gene_name = "character"
                      ))

CNV_twin = setClass("CNV_twin",
                    slots = list(
                      name = "character",
                      matrix_1 = "GenomicRanges",
                      matrix_2 = "GenomicRanges",
                      gene_name_1 = "character",
                      gene_name_2 = "character"
                    ))

# set generic plotCNV and plotCNVs(plot single CNV and CNV twin)

setGeneric('plotCNVs', function(object, ...) standardGeneric('plotCNVs'))
setGeneric('plotCNV', function(object, ...) standardGeneric('plotCNV'))

# CNV plot first step

CNV.by.method <- function(CNV_1,gene.name,pids,title,legend,legend.names,
                          out.dir,file.type,pixel.per.cnv,color,display,
                          gene.anno,start.gene,end.gene){
  # solid parameters
  chrom = as.vector(seqnames(CNV_1))
  start.CNV=start(CNV_1)
  end.CNV=end(CNV_1)
  legend=2
  
  # initail
  index <- (end.CNV - start.CNV) < 10000000 # only events shorter than 10 M
  m <- sum(index)
  startPos <- start.CNV[index]
  endPos <- end.CNV[index]
  
  if(method==""){
    rescore <- rep(100000000,m)
    score.values <- as.character(sort(unique(rescore)))
    n <- length(unique(rescore))
  }
  
  if(method=="pblcbp"){
    if(missing(score))  # if no argument is given --> score is 4 (diploid)
    {score <- rep(100000000,length(index))}
    rescore <- unlist(lapply(score,map.ploidy.classes))[index]
    score.values <- as.character(sort(unique(rescore)))
    n <- length(unique(rescore))
  }
  
  if(1>0){
  
  if(missing(gene.name)){gene.name <- "geneX"}
  cnv.type <- gene.name 
  if(missing(title)){
    if(missing(pids)){
      title <- cnv.type
      }else{
        title <- paste(cnv.type,":",m,"events from",length(unique(pids[index])),"samples") } # normal title genereated when pids given
    }else{
      if(missing(pids)){
        title <- title
        }else{
          title <- paste(title,":",m,"events from",length(unique(pids[index])),"samples") 
          } # normal title genereated when pids given
  }
  if(missing(legend)){legend <- "missing"}
  if(missing(legend.names)){
    legend.names <- c("unkown ploidy","0")
  }else{
    legend.names <- legend.names
  }
  ## default/optional parameter for file.type (default = pdf) ----------------------------------------------------------------------------------------------------------
  if(missing(file.type)){ 
    file.type1 <- pdf
    plot.type <- file.type1
  }else{
      plot.type <- file.type
      }
  ## default/optional parameter for out.dir (defaultdirectory = "/Users/CNV.by.cohort") --------------------------------------------------------------------------------
  if(missing(out.dir)){
    if(missing(file.type)){
      out.dir <- paste("CNV.by.length","_",gene.name,".","pdf",sep = "")
    }else{
        out.dir <- paste("CNV.by.length","_",gene.name,".",substitute(file.type),sep = "") 
        }
  }else{
    if(missing(file.type)){
      out.dir <- paste(out.dir,"_",gene.name,".","pdf",sep="")
    }else{
        out.dir <- paste(out.dir,"_",gene.name,".",substitute(file.type),sep="") 
        }
  }
  out.fp <- out.dir
  if(missing(pixel.per.cnv)){
    pixel.per.cnv <- 200/m
    }  ## better a equation dependened on the number of CNVs (index!)
  ## color ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(missing(color)){
    color <- "steelblue3"
  }else{
    color <- color
  }
  
  ## where to plot?----------------------------------------------------------------------------------------------------------------------------------------------------
  if(missing(display)){
    display <- "missing"
  }
  
  if(missing(gene.anno)){
    gene.anno <- "missing"
  }
  
  if(missing(start.gene) & gene.anno == TRUE){
    print("start.gene argument is missing")
  }
  
  if(missing(end.gene) & gene.anno == TRUE){
    print("end.gene argument is missing")
  }
  }
  
  ## sorting ----------------------------------------------------------------------------------------------------------------------------------------------------------
  
  if(method=="sort.by.length"){
    sorting <- order(endPos - startPos) # sort by length
  }else if(method=="sort.by.cohort"){
    if(missing(cohort)){
      print("use CNV.by.ploidy or CNV.by.length functions")
      }else{
        cohort <- cohort[index]
        sorting <- order(cohort,endPos - startPos)
      }
  }else if(method=="sort.by.ploidy"){
    sorting <- order(rescore,endPos - startPos) # sort by score and then by length, sort by ploidy
  }else if(method=="sblcbc"){
    sorting <- order(endPos - startPos) # sort by length, color by cohort
    cohort <- cohort[index]
  }else if(method=="sblcbp"){
    sorting <- order(endPos - startPos) # sort by length, color by plotdy
  }
  paralist <- list("gene.name"=gene.name,"cnv.type"=cnv.type,"title"=title,"pids"=pids,"legend"=lengend,
                   "legend.names"=legend.names,"file.type"=file.type,"out.dir"=out.dir,"pixel.per.cnv"=pixel.per.cnv,
                   "color"=color,"sorting"=sorting,"start.gene"=start.gene,"end.gene"=end.gene,"gene.anno"=gene.anno,
                   "chrom"=chrom,"start.CNV"=start.CNV,"end.CNV"=end.CNV,"rescore"=rescore,
                    "index"=index,"m"=m,"startPos"=startPos,"endPos"=endPos)
  return(paralist)
}


# definition of focallity score
focallity.score <- function(m,ends,starts)
{
  mean.length <-  (sum(ends - starts))/m
  range <- range(ends - starts)
  range.length <- range[2]-range[1]
  #sd.score <- round(sd(range.length/(ends - starts)),2)
  #sd.score <- round(sd(mean.length/(ends - starts)),2)
  mean.score <- round(mean(range.length/(ends - starts)),2)
  #mean.score <- round(mean(mean.length/(ends - starts)),2)
  #f.score <- paste(mean.score,"±",sd.score)
  f.score <- mean.score
}


# implement plotCNVs

plotCnvs.cohort <- function(){paralist,
  chrom = unlist(paralist["chrom"])
  startPos = unlist(paralist["startPos"])
  endPos = unlist(paralist["endPos"])
  rescore = unlist(paralist["rescore"])
  score = unlist(paralist["score"])
  sorting = unlist(paralist["sorting"])
  title = unlist(paralist["title"])
  pixel.per.cnv = unlist(paralist["pixel.per.cnv"])
  plot.type = unlist(paralist["plot.type"])
  # getcolor
  lengend = unlist(paralist["legend"])
  lengend.names = unlist(paralist["legend.names"])
  color = unlist(paralist["color"])
  score.values = unlist(paralist["score.values"])
  n = unlist(paralist["n"])
  display = unlist(paralist["display"])
  gene.anno = unlist(paralist["gene.anno"])
  cnv.type = unlist(paralist["cnv.type"])
  start.gene = unlist(paralist["start.gene"])
  end.gene = unlist(paralist["end.gene"])
  
  chroms <- chrom[sorting]
  starts <- startPos[sorting]
  ends <- endPos[sorting]
  cohorts <- rescore[sorting] #?
  cohorts <- droplevels.factor(cohorts, exclude = if(anyNA(levels(cohorts)))NULL else NA)  ## erase factor levels = 0 (turns out very important for color plotting)
  cnv.number <-  length(chroms) # number of lines in input
  chromWidth <- round((pixel.per.cnv * cnv.number) * 0.1)
  f.score <- focallity.score(m=length(starts),starts = starts,ends = ends)
  if (length(unique(chroms)) > 1){
    print(unique(chroms))
    print("More than one chromosome id - use other function")
    return()
  }
  y <- lengthChromosome(chroms[1],"bases") + 10000000
  plot.new()
  par(c(5,3,4,4))
  pixelPerChrom <- chromWidth + (pixel.per.cnv)*(cnv.number+1)+10 # determines space between chromsomes
  x.size <- pixelPerChrom
  y.size <- y+100
  plot(c(0,x.size),c(0,y.size),type="n",xaxt="n",yaxt="n",xlab="SVs",ylab="Chromosomal location",main=title)
  chrStr <- paste("chr",toString(chroms[1]))
  text(c((chromWidth/2)),c(0),labels=c(chrStr))
  if(gene.anno == TRUE)    ###added RT gene.anno arg
  {
    m <- mean(start.gene,end.gene)
    text(c(-1),c(y-m+(y*0.035)),labels=c(cnv.type),cex=0.5)
    paintCytobands(chroms[1],pos=c(chromWidth,y),units="bases",width=chromWidth-7,orientation="v",legend=FALSE)
    rect(7,y-start.gene,chromWidth-1,y-end.gene,col="gray50", border = "gray50")
    lines(c(0.6,2.25),c(y-m+(y*0.02),y-m),col="gray50")
    lines(c(2.25,5),c(y-m,y-m),col="gray50")
  }else{
    paintCytobands(chroms[1],pos=c(chromWidth,y),units="bases",width=chromWidth,orientation="v",legend=FALSE)
  }
  
  plotCnv.cohort(chroms,starts,ends,y,chromWidth=chromWidth,pixel.per.cnv=pixel.per.cnv,cohorts=cohorts,startPoint=chromWidth,getColor=getColor.cohort,color=color)
  
  
  # legend position decision (top or bottom)
  centro <- c(125,93.3,91,50.4,48.4,61,59.9,45.6,49,40.2,53.7,35.8,17.9,17.6,19,36.6,24,17.2,26.5,27.5,13.2,14.7,60.6,12.5)*1000000
  length <- lengthChromosome(c(1:22,"X","Y"),"bases")/2
  genome <- data.frame(chromosome=c(1:22,"X","Y"), centromere=centro, length=length) # dataframe containing chromosome and centromere position info
  
  mean.pos <- mean(c(starts,ends)) # mean position of all CNV´s
  #centroo <- genome[genome$chromosome %in% chroms,] # centromere position in the current chromosome
  half.length <- genome[genome$chromosome %in% chroms,] # half of the length of the current chromosome
  
  # mean CNV is over the centromere -> legend is plotted bottomright
  if(mean.pos < half.length$length){ 
    xtr <- "bottomright"
    xtf <- c(4,24,20.5,3)
    text(c(pixelPerChrom/2),c(y-10),labels = paste("score: ",f.score),cex=0.75) # score on opposite
  }
  
  if(mean.pos > half.length$length){
    xtr <- "topright"
    xtf <- c(21.5,24,4,3)
    text(c(pixelPerChrom/2),c(10),labels = paste("score: ",f.score),cex=0.75)
  }
  
  # legend type decision ----------------------------------------------------------------------------
  if(legend=="missing" || legend==1){
    legend(xtr,legend=unique(cohorts),col=getColor.cohort(color=color,cohorts=cohorts,q=F),cex=0.75,pch=16) # normal legend
  }
  
  if(legend==2){
    par(new=T,mar=xtf )
    pie(table(cohorts),col=getColor.cohort(color=color,cohorts=cohorts,q=F),cex=0.52) # piechart legend
  }else{} # no legend
  
  if(display == TRUE){}else{dev.off()}
  
}


# implement plotCNV


plotCnv.cohort <- function(chroms,starts,ends,y,chromWidth,pixel.per.cnv,cohorts,startPoint,color,method)
{
  indX <- chroms == 'X'
  indY <- chroms == 'Y'
  len <- length(starts)
  GetColor(method=method,color=color,cohorts=cohorts,q=TRUE)
  #getColor(color=color,cohorts=cohorts,q=TRUE) # source in the color palette (argument from package function)
  
  # Autosomes
  for(index in 1:len){
    x <- startPoint + pixel.per.cnv*index
    lines(c(x,x),c(y-starts[index],y-ends[index]),col=cohorts[index],lwd=pixel.per.cnv)
  }
  
  # X Chromosome
  for(index in 1:len){
    if(indX[index] == TRUE){
      x <- startPoint +pixel.per.cnv*index
      lines(c(x,x),c(y-starts[index],y-ends[index]),col=cohorts[index],lwd=pixel.per.cnv)
    }
  }
  
  # Y Chromosome
  for(index in 1:len){
    if(indY[index] == TRUE){
      x <- startPoint + pixel.per.cnv*index
      lines(c(x,x),c(y-starts[index],y-ends[index]),col=cohorts[index],lwd=pixel.per.cnv)  ## RT cohorts for scores
    }
  }
}

# implement pie plot

PlotPie <- function(method,feature,index,t="Overview"){
  switch (method,
    "standard" = {
     p <-  pie(sort(table(feature)),main=t)
    },
    "index" = {
     p <-  pie(sort(table(feature[index])),main=t)
    }
  )
  return(p)
}

PlotCohortLength <- function(data.feature,group.feature,t="Cohort_Length"){
  l <- length(group.feature)
  p <- boxplot(data.feature ~ group.feature,main=t)
  return(p)
}

# external color lengends

resolution.zscore <- 100

colors.zscore <- rev(terrain.colors(2*resolution.zscore + 1))

legend.col <- function(col, lev){
  
  opar <- par
  n <- length(col)
  bx <- par("usr")
  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
              bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n
  
  xx <- rep(box.cx, each = 2)
  
  par(xpd = TRUE)
  for(i in 1:n){
    #print(paste("Legend Box",i))
    
    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])
    
  }
  #print("After Box")
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(0, 100),
       yaxt = "n", ylab = "",
       xaxt = "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = .25)
  par <- opar
  #print("Finisched Legend")
}

get.background.rank <- function(gene,gene.symbol,score){
  ranking <- rank(-score)
  view <- data.frame(gene.symbol,score,ranking)
  background.rank = view$ranking[view$gene.symbol == gene]
  return(background.rank)
}




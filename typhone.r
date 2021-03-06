# dependencies
library(png)
library(grid)
library(GenomicRanges)
library(quantsmooth)
library(gridExtra)
# test data set generation

test = 0
if(test==0){

CNV <- data.frame(Chromosome=c(rep(12,2000)),
                  Start=c(round(runif(n=2000,min=15000000,max=35000000))),
                  End=c(round(runif(n=2000,min=1500000,max=35000000))),
                  Score=c(round(runif(n=2000,min=1,max=12))),Gene=c(rep("gene",2000)),
                  Cohort=sample(c("BRCA","AML","CML","CRC","GLIOMA"),size=2000,
                                prob=c(0.45,0.2,0.05,0.15,0.15),replace = T),
                  PID=c(stringi::stri_rand_strings(n=2000, length=9, pattern = "[A-Za-z0-9]")))

CNV <- cbind(CNV,length=CNV$Start-CNV$End)
CNV <- CNV[CNV$length <= 0,]
CNV <- makeGRangesFromDataFrame(CNV , keep.extra.columns = TRUE)
KRAS <- GRanges(seqnames =Rle(12) , ranges=IRanges(start=25204789,end=25250936))
CNV.KRAS <- subsetByOverlaps(CNV,KRAS)

CNV1 <- data.frame(Chromosome=c(rep(12,2000)),
                   Start=c(round(runif(n=2000,min=15000000,max=35000000))),
                   End=c(round(runif(n=2000,min=1500000,max=35000000))),
                   Score=sample(c(1:6),size=2000,
                                prob=c(0.01,0.1,0.04,0.15,0.65,0.15),replace = TRUE),
                   Gene=rep("KRAS",2000),
                   Cohort=sample(c("BRCA","AML","CML","CRC","GLIOMA"),size=2000,
                                 prob=c(0.45,0.2,0.05,0.15,0.15),replace = T),
                   PID=c(stringi::stri_rand_strings(n=2000, length=9, pattern = "[A-Za-z0-9]")))

CNV2 <- data.frame(Chromosome=c(rep(12,2000)),
                   Start=c(round(runif(n=2000,min=13000000,max=34000000))),
                   End=c(round(runif(n=2000,min=1300000,max=34000000))),
                   Score=sample(c(1:6),size=2000,
                                prob=c(0.15,0.5,0.1,0.05,0.15,0.05),replace = TRUE),
                   Gene=rep("STK38L",2000),
                   Cohort=sample(c("BRCA","AML","CML","CRC","GLIOMA"),size=2000,
                                 prob=c(0.55,0.1,0.15,0.05,0.15),replace = T),
                   PID=c(stringi::stri_rand_strings(n=2000, length=9, pattern = "[A-Za-z0-9]")))

# subsetting positive length
CNV1 <- cbind(CNV1,length=CNV1$Start-CNV1$End)
CNV1 <- CNV1[CNV1$length <= 0,]

CNV2 <- cbind(CNV2,length=CNV2$Start-CNV2$End)
CNV2 <- CNV2[CNV2$length <= 0,]

# creating GenomicRange data
library(GenomicRanges)
CNV1 <- makeGRangesFromDataFrame(CNV1 , keep.extra.columns = TRUE)
CNV2 <- makeGRangesFromDataFrame(CNV2 , keep.extra.columns = TRUE)
KRAS <- GRanges(seqnames =Rle(12) , ranges=IRanges(start=25204789,end=25250936))
STK38L <- GRanges(seqnames =Rle(12) , ranges=IRanges(start=27396901,end=27478892))

# subsetting overlays with KRAS region
CNV.KRAS <- subsetByOverlaps(CNV1,KRAS)
CNV.STK38L <- subsetByOverlaps(CNV2,STK38L)


#

test_cnv <- new("CNV_single",name="CNV_test",matrix=CNV.KRAS,gene_name="KRAS")
test_cnv_twin <- new("CNV_twin",name="Twin_Test",matrix_1=CNV.KRAS,matrix_2=CNV.STK38L,gene_name_1="KRAS",gene_name_2="STK38L")
#bb <- plotCnvs.cohort(paralist=para1,SaveAsObject = SaveAsObject)

}


# map classes

MapPloidyClasses <- function(v){
  class = 1
  if(v<5){
  switch(v,
         "1" = {class = 2},
         "2" = {class = 3},
         "3" = {class = 4},
         "4" = {class = 4}
         )
  }else if(v>=5 & v<=8){
    class = 5
  }else if(v>=9 & v<=99999999){
    class = 6
  }else{class = v}
  return(class)
}

# get colors

GetColor <- function(method,score,color,score.values,n,greyscale,cohorts){
  if(missing(score)){score=0}
  if(missing(score.values)){score.values=0}
  if(missing(n)){n=0}
  if(missing(greyscale)){greyscale=FALSE}
  if(missing(cohorts)){cohort=c("all_patients")}
  
  
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
           cohort.size <- length(unique(cohorts))
           if(missing(color)){
             color1 <- colorRampPalette(c("red2","indianred4","royalblue4","steelblue1","chartreuse3","darkgreen"))(cohort.size)
           }else if(color>0 & color<8 & is.integer(color)==TRUE){
             switch (color,
                     "1" = {color1 <- rainbow(cohort.size)},
                     "2" = {color1 <- colorRampPalette(c("seashell2","seagreen2","turquoise2","palevioletred2"))(cohort.size)},
                     "3" = {color1 <- colorRampPalette(c("gray7","mediumblue","deeppink4","sienna3"))(cohort.size)},
                     "4" = {color1 <- colorRampPalette(c("darkslateblue","darkslategray4","deeppink4","tan4","gray10"))(cohort.size)},
                     "5" = {color1 <- colorRampPalette(c("navajowhite3","orange3","olivedrab3"))(cohort.size)},
                     "6" = {color1 <- colorRampPalette(c("royalblue2","yellow1"))(cohort.size)},
                     "7" = {color1 <- gray.colors(cohort.size)}
             )
           }else if(greyscale==TRUE){
             color1 <- gray.colors(cohort.size)
             color1 <- palette(color1)
         }else{
           color1 <- colorRampPalette(c("red2","indianred4","royalblue4","steelblue1","chartreuse3","darkgreen"))(cohort.size)
           }
           color.value <- color1
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

setGeneric('TornadoPlots', function(object, ...) standardGeneric('TornadoPlots'))
setGeneric('plotCNVs', function(object, ...) standardGeneric('plotCNVs'))
setGeneric('plotCNV', function(object, ...) standardGeneric('plotCNV'))

# setMethod for TornadoPlots
setMethod("TornadoPlots",signature("CNV_single"),function(object,gene.name,pids,title,legend,legend.names,
                                  out.dir,file.type,pixel.per.cnv,color,display,
                                  gene.anno,start.gene,end.gene,color.method,sort.method,SaveAsObject){
  paralist0 <- CNV.by.method(object,gene.name,pids,title,legend,legend.names,
                             out.dir,file.type,pixel.per.cnv,color,display,
                             gene.anno,start.gene,end.gene,color.method,sort.method)
  if(SaveAsObject==TRUE){
    plot0 <- plotCnvs.cohort(paralist=paralist0,SaveAsObject = SaveAsObject)
  }else{
    print("Output image is saved!!")
  }
})

setMethod("TornadoPlots",signature("CNV_twin"),function(){
  paralist0 <- PlotTwinsInit()
  if(SaveAsObject==TRUE){
    plot0 <- PlotTwins(paralist=paralist0,SaveAsObject=SaveAsObject)
  }else{
    print("Output image is saved!!")
  }
})



# CNV plot first step

CNV.by.method <- function(CNV.input,gene.name,pids,title,legend,legend.names,
                          out.dir,file.type,pixel.per.cnv,color,display,
                          gene.anno,start.gene,end.gene,sort.method,color.method){
  
  CNV_1 <- CNV.input@matrix
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
  
  if(missing(sort.method) & missing(color.method)){sort.method = "length"}
  if(missing(sort.method)){sort.method = color.method}
  if(missing(color.method)){color.method = sort.method}
  
  score = CNV_1$Score
  cohort = CNV_1$Cohort
  pids = CNV_1$PID
  
  
  
  if(sort.method=="length" & color.method=="length"){
    rescore <- rep(100000000,m)
    score.values <- as.character(sort(unique(rescore)))
    n <- length(unique(rescore))
  }
  
  if(sort.method=="length" & color.method=="ploidy"){
    if(missing(score)){score <- rep(100000000,length(index))}  # if no argument is given --> score is 4 (diploid)
    rescore <- unlist(lapply(score,MapPloidyClasses))[index]
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
  
  if(sort.method=="length"){
    sorting <- order(endPos - startPos) # sort by length
  }else if(sort.method=="cohort"){
    if(missing(cohort)){
      print("use CNV.by.ploidy or CNV.by.length functions")
      }else{
        cohort <- cohort[index]
        sorting <- order(cohort,endPos - startPos)
      }
  }else if(sort.method=="ploidy"){
    sorting <- order(rescore,endPos - startPos) # sort by score and then by length, sort by ploidy
  }else if(sort.method=="length" & color.method=="cohort"){
    sorting <- order(endPos - startPos) # sort by length, color by cohort
    cohort <- cohort[index]
  }else if(sort.method=="length" & color.method=="ploidy"){
    sorting <- order(endPos - startPos) # sort by length, color by ploidy
  }
  if(missing(start.gene)){start.gene <- "geneX"}
  if(missing(end.gene)){end.gene <- "geneX"}
  
  file.type="default"
  out.dir="default"
  out.fp <- out.dir
  
  paralist <- list("gene.name"=gene.name,"cnv.type"=cnv.type,"title"=title,"legend"=legend,
                   "legend.names"=legend.names,"file.type"=file.type,"out.dir"=out.dir,"pixel.per.cnv"=pixel.per.cnv,
                   "color"=color,"sorting"=sorting,"start.gene"=start.gene,"end.gene"=end.gene,"gene.anno"=gene.anno,
                   "chrom"=chrom,"start.CNV"=start.CNV,"end.CNV"=end.CNV,"rescore"=rescore,
                    "index"=index,"m"=m,"startPos"=startPos,"endPos"=endPos,
                   "score"=score,"cohort"=cohort,"pids"=pids,
                   "sort.method"=sort.method,"color.method"=color.method)
  return(paralist)
}





PlotTwinsInit <- function(twin.cnv,
                          pids,file.type,out.dir,
                          out.fp,title,pixel.per.cnv,plot.type,
                          method,legend,legend.names,color,score.values_1,score.values_2,
                          n_1,n_2,display,gene.anno,cnv.type_1,cnv.type_2,start.gene_1,start.gene_2,
                          end.gene_1,end.gene_2){
  CNV_1 <- twin.cnv@matrix_1
  CNV_2 <- twin.cnv@matrix_2
  
  chrom_1 <- as.vector(seqnames(CNV_1))
  start.CNV_1 <- start(CNV_1)
  end.CNV_1 <- end(CNV_1)
  score_1 <- CNV_1$Score
  
  chrom_2 <- as.vector(seqnames(CNV_2))
  start.CNV_2 <- start(CNV_2)
  end.CNV_2 <- end(CNV_2)
  score_2 <- CNV_2$Score
  
  gene.name_1 <- as.character(twin.cnv@gene_name_1)
  gene.name_2 <- as.character(twin.cnv@gene_name_2)
  
  index_1 <- (end.CNV_1 - start.CNV_1) < 10000000 # only events shorter than 10 M
  m_1 <- sum(index_1)
  startPos_1 <- start.CNV_1[index_1]
  endPos_1 <- end.CNV[index_1]
  rescore_1 <- rep(gene.name_1,m_1)
  score.values_1 <- as.character(sort(unique(rescore_1)))
  n_1 <- length(unique(rescore_1))
  
  index_2 <- (end.CNV_2 - start.CNV_2) < 10000000 # only events shorter than 10 M
  m_2 <- sum(index_2)
  startPos_2 <- start.CNV_2[index_2]
  endPos_2 <- end.CNV_2[index_2]
  rescore_2 <- rep(gene.name_2,m_2)
  score.values_2 <- as.character(sort(unique(rescore_2)))
  n_2 <- length(unique(rescore_2))
  
  ## default/optional parameter for gene name (default = "geneX")-------------------------------------------------------------------------------------------------------
  if(missing(gene.name_1))  # if no argument is given --> gene.name is "geneX"
  {gene.name_1 <- "geneX"}
  cnv.type_1 <- gene.name_1   # assign gene.name to cnv.type --> can be replaced later on
  
  if(missing(gene.name_2))  # if no argument is given --> gene.name is "geneX"
  {gene.name_2 <- "geneY"}
  cnv.type_2 <- gene.name_2   # assign gene.name to cnv.type --> can be replaced later on
  
  ## default/optional parameter for "pids" and "title" (default:"gene.name: m events from unkown amount of unique samples")---------------------------------------------------------
  if(missing(title)){
    if(missing(pids)){
      title <- paste(cnv.type_1,"&",cnv.type_2) 
    }else{
      title <- paste(cnv.type_1,":",m_1,"events from",length(unique(pids[index_1])),"samples") } # normal title genereated when pids given
  }else{
    if(missing(pids)){
      title <- title
    }else{
      title <- paste(title,":",m_1,"events from",length(unique(pids[index_1])),"samples") } # normal title genereated when pids given
  }
  
  ## default/optional parameter for legend (default = "missing" processed to "normal legend")
  if(missing(legend)){legend <- "missing"}
  
  ## default/optional parameter for legend.names
  if(missing(legend.names)){
    legend.names <- c(score.values_2[1:n_2]) # ??
  }else{
    legend.names <- c(legend.names)[1:n_1]
  }
  ## default/optional parameter for file.type (default = pdf) ----------------------------------------------------------------------------------------------------------
  if(missing(file.type)){
    file.type <- pdf
    plot.type <- pdf
  }else{
    plot.type <- file.type
  }
  
  ## default/optional parameter for out.dir (defaultdirectory = "/package/TornadoCNV") --------------------------------------------------------------------------------
  if(missing(out.dir)){
    if(missing(file.type)){
      out.dir <- paste("twins.by.length","_",gene.name_1,"&",gene.name_2,".","pdf",sep = "") 
    }else{
      out.dir <- paste("twins.by.length","_",gene.name_1,"&",gene.name_2,".",sep = "") }
  }else{
    if(missing(file.type)){
      out.dir <- paste(out.dir,"_",gene.name_1,"&",gene.name_2,".","pdf",sep="") 
    }else{
      out.dir <- paste(out.dir,"_",gene.name_1,"&",gene.name_2,".",sep="") 
    }
  }
  out.fp <- out.dir
  
  ## default/optional parameter for pixel.per.cnv (default = 5)---------------------------------------------------------------------------------------------------------
  if(missing(pixel.per.cnv)){pixel.per.cnv <- 200/(m_1+m_2)}  ## better a equation dependened on the number of CNVs (index!)
  
  ## sorting ----------------------------------------------------------------------------------------------------------------------------------------------------------
  
  sorting_1 <- rev(order(endPos_1 - startPos_1)) # sort by length #### add ??
  sorting_2 <- order(endPos_2 - startPos_2)
  
  ## color ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(missing(color)){
    color <- "steelblue3"
  }else{
    color <- color
  }
  
  ## where to plot?----------------------------------------------------------------------------------------------------------------------------------------------------
  if(missing(display)){
    display <- ""
  }
  
  ## gene anno?----------------------------------------------------------------------------------------------------------------------------------------------------
  if(missing(gene.anno)){
    gene.anno <- ""
  }
  
  if(missing(start.gene_1)){start.gene_1 <- "geneX"}
  if(missing(end.gene_1)){end.gene_1 <- "geneX"}
  if(missing(start.gene_2)){start.gene_2 <- "geneY"}
  if(missing(end.gene_2)){end.gene_2 <- "geneY"}
  
  
  
  ## gene.anno TRUE but missing start/end.gene------------------------
  if(missing(start.gene_1) & gene.anno == TRUE){
    print("start.gene 1 argument is missing")
  }
  
  if(missing(end.gene_1) & gene.anno == TRUE){
    print("end.gene 1 argument is missing")
  }
  
  if(missing(end.gene_2) & gene.anno == TRUE){
    print("end.gene 2 argument is missing")
  }
  
  if(missing(end.gene_2) & gene.anno == TRUE){
    print("end.gene 2 argument is missing")
  }
  
  paralist <- list("gene.name_1"=gene.name_1,"gene.name_2"=gene.name_2,
                   "cnv.type_1"=cnv.type_1,"cnv.type_2"=cnv.type_2,
                   "title"=title,"legend"=legend,
                   "legend.names"=legend.names,"file.type"=file.type,"out.dir"=out.dir,"pixel.per.cnv"=pixel.per.cnv,
                   "color"=color,"sorting_1"=sorting_1,"sorting_2"=sorting_2,
                   "start.gene_1"=start.gene_1,"end.gene_1"=end.gene_1,"start.gene_2"=start.gene_2,"end.gene_2"=end.gene_2,
                   "gene.anno"=gene.anno,"n_1"=n_1,"n_2"=n_2,
                   "chrom_1"=chrom_1,"start.CNV_1"=start.CNV_1,"end.CNV_1"=end.CNV_1,"rescore_1"=rescore_1,
                   "chrom_2"=chrom_2,"start.CNV_2"=start.CNV_2,"end.CNV_2"=end.CNV_2,"rescore_2"=rescore_2,
                   "index_1"=index_1,"index_2"=index_2,"m_1"=m_1,"m_2"=m_2,
                   "startPos_1"=startPos_1,"endPos_1"=endPos_1,
                   "startPos_2"=startPos_2,"endPos_2"=endPos_2)
  
  
}


PlotTwins <- function(paralist,SaveAsObject = SaveAsObject){
  
  chroms_1 <- unlist(paralist["chrom_1"])
  starts_1 <- unlist(paralist["startPos_1"])
  ends_1 <- unlist(paralist["endPos_1"])
  scores_1 <- unlist(paralist["score_1"])
  f.score_1 <- focallity.score(m=length(starts_1),starts = starts_1,ends = ends_1)
  
  chroms_2 <- unlist(paralist["chrom_2"])
  starts_2 <- unlist(paralist["startPos_2"])
  ends_2 <- unlist(paralist["endPos_2"])
  scores_2 <- unlist(paralist["score_2"])
  f.score_2 <- focallity.score(m=length(starts_2),starts = starts_2,ends = ends_2)
  
  cnv.number <- (length(chroms_1)+length(chroms_2)) # number of lines in input
  chromWidth <- round((pixel.per.cnv * cnv.number) * 0.1)
  
  if (length(unique(chroms_1)) > 1){
    print(unique(chroms_1))
    print("More than one chromosome id - use other function")
    return()
  }
  
  
  if (length(unique(chroms_2)) > 1){
    print(unique(chroms_2))
    print("More than one chromosome id - use other function")
    return()
  }
  
  
  y <- lengthChromosome(chroms_1[1],"bases") + 10000000  ## are you sure??
  
  # plot parameters -----------------------------------------------------------------------------------------------------------------
  plot.new()
  png("t2.png",width = 1024,height=768,units = "px")
  
  par(c(5,3,4,4))
  pixelPerChrom_1 <-  (pixel.per.cnv)*(length(chroms_1)+1)
  pixelPerChrom_2 <-  (pixel.per.cnv)*(length(chroms_2)+1)
  pixelPerChrom <- chromWidth+pixelPerChrom_1+pixelPerChrom_2+10 # determines space between chromsomes
  
  x.size <- pixelPerChrom
  y.size <- y+100
  plot(c(0,x.size),c(0,y.size),type="n",xaxt="n",yaxt="n",xlab="CNVs",ylab="Chromosomal location",main=title)
  chrStr <- paste("chr",toString(chroms_1[1]))
  text(c(pixelPerChrom_1+(chromWidth/2)),c(0),labels=c(chrStr))
  
  if(gene.anno == TRUE){
    paintCytobands(chroms_1[1],pos=c(pixelPerChrom_1+chromWidth,y),units="bases",width=chromWidth,orientation="v",legend=FALSE)
    m_1 <- mean(start.gene_1,end.gene_1)
    m_2 <- mean(start.gene_2,end.gene_2)
    
    text(c(pixelPerChrom_1+chromWidth+15),c(y-m_1+(y*0.045)),labels=c(cnv.type_1),cex=0.7)
    rect(pixelPerChrom_1+1,y-m_1,pixelPerChrom_1+chromWidth-1,y-m_1,col="gray50", border = "gray50")
    lines(c(pixelPerChrom_1+chromWidth+7,pixelPerChrom_1+chromWidth+4),c(y-m_1+(y*0.03),y-m_1),col="gray50")
    lines(c(pixelPerChrom_1+chromWidth+4,pixelPerChrom_1+chromWidth+1),c(y-m_1,y-m_1),col="gray50")
    
    text(c(pixelPerChrom_1-15),c(y-m_2-(y*0.045)),labels=c(cnv.type_2),cex=0.7)
    rect(pixelPerChrom_1+1,y-m_2,pixelPerChrom_1+chromWidth-1,y-m_2,col="gray50", border = "gray50")
    lines(c(pixelPerChrom_1-7,pixelPerChrom_1-4),c(y-m_2-(y*0.03),y-m_2),col="gray50")
    lines(c(pixelPerChrom_1-4,pixelPerChrom_1-1),c(y-m_2,y-m_2),col="gray50")
    
  }else{
    paintCytobands(chroms_1[1],pos=c(pixelPerChrom_1+chromWidth,y),units="bases",width=chromWidth,orientation="v",legend=FALSE)
  }
  
  plotCnv(chroms_1,starts_1,ends_1,y,scores_1,pixel.per.cnv=pixel.per.cnv,method="by.length",color=color,score.values = score.values_1,n=n_1,startPoint=(pixelPerChrom_1),direction = "left")
  plotCnv(chroms_2,starts_2,ends_2,y,scores_2,pixel.per.cnv=pixel.per.cnv,method="by.length",color=color,score.values = score.values_2,n=n_2,startPoint=(pixelPerChrom_1+chromWidth),direction = "right")
  
  #plotCnv.cohort <- function(chroms,starts,ends,y,chromWidth,pixel.per.cnv,cohorts,startPoint,color,method)
  
  # legend parameters ------------------------------------------------------------------------------------------------------
  df.color.ploidy <- data.frame(color=color[1:n], # colors according to getColor.ploidy/2
                                score=score.values[1:n], # score according to
                                names=legend.names[1:n])
  
  data.score <- data.frame(score=c(sort(unique(scores)))) # unique and present scores of inout data
  dtt <- df.color.ploidy[df.color.ploidy$score %in% data.score$score,] # subset only present scores of input data
  color <- as.vector(dtt$color)
  labs <- as.vector(dtt$names)
  
  # legend position decision (top or bottom)
  centro <- c(125,93.3,91,50.4,48.4,61,59.9,45.6,49,40.2,53.7,35.8,17.9,17.6,19,36.6,24,17.2,26.5,27.5,13.2,14.7,60.6,12.5)*1000000
  length <- lengthChromosome(c(1:22,"X","Y"),"bases")/2
  genome <- data.frame(chromosome=c(1:22,"X","Y"), centromere=centro, length=length) # dataframe containing chromosome and centromere position info
  
  mean.pos <- mean(c(starts,ends)) # mean position of all CNV´s
  #centroo <- genome[genome$chromosome %in% chroms,] # centromere position in the current chromosome
  half.length <- genome[genome$chromosome %in% chroms,] # half of the length of the current chromosome
  
  if(mean.pos < half.length$length) {
    xtr <- "bottomleft"
    xtr2 <- "bottomright"
    xtf <- c(4,5,20.5,22)
    xtf2 <- c(4,24,20.5,3)
    text(c(pixelPerChrom_1/2),c(y-10),labels = paste("score: ",f.score_1),cex=0.75)
    text(c(pixelPerChrom_1+chromWidth+(pixelPerChrom_2/2)),c(y-10),labels = paste("score: ",f.score_2),cex=0.75)
    
  }    # mean start end smaller than subset chrom centromer
  if(mean.pos > half.length$length){
    xtr <- "topleft"
    xtr2 <- "topright"
    xtf <- c(21.5,5,4,22)
    xtf2 <- c(21.5,24,4,3)
    text(c(pixelPerChrom_1/2),c(10),labels = paste("score: ",f.score_1),cex=0.75)
    text(c(pixelPerChrom_1+chromWidth+(pixelPerChrom_2/2)),c(10),labels = paste("score: ",f.score_2),cex=0.75)
    
  }
  
  # legend type decision ----------------------------------------------------------------------------
  if(legend=="missing" || legend==1){
    legend(xtr,legend=labs,col=color,cex=0.75,pch=16) # normal legend
    legend(xtr2,legend=labs,col=color,cex=0.75,pch=16)
  }
  if(legend==2) {
    par(new=T,mar=xtf )
    pie(table(score_1),labels=labs,col=color,cex=0.52) # piechart legend
    par(new=T,mar=xtf2)
    pie(table(score_2),labels=labs,col=color,cex=0.52)
  } else{} # no legend
  
  dev.off()
  if(SaveAsObject==TRUE){
    img <- readPNG("t2.png")
    g <- rasterGrob(img, interpolate=TRUE)
    return(g)
  }
  
  
}



# definition of focallity score
focallity.score <- function(m,ends,starts){
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

plotCnvs.cohort <- function(paralist,SaveAsObject){
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
  legend = unlist(paralist["legend"])
  legend.names = unlist(paralist["legend.names"])
  color = unlist(paralist["color"])
  score.values = unlist(paralist["score.values"])
  n = unlist(paralist["n"])
  display = unlist(paralist["display"])
  gene.anno = unlist(paralist["gene.anno"])
  cnv.type = unlist(paralist["cnv.type"])
  start.gene = unlist(paralist["start.gene"])
  end.gene = unlist(paralist["end.gene"])
  sort.method = unlist(paralist["sort.method"])
  color.method = unlist(paralist["color.method"])
  score = unlist(paralist["score"])
  pids = unlist(paralist["pids"])
  cohort = unlist(paralist["cohort"])
  
  
  
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
  png("t1.png",width = 1024,height=768,units = "px")
  
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
  
  plotCnv.cohort(chroms,starts,ends,y,
                 chromWidth=chromWidth,pixel.per.cnv=pixel.per.cnv,
                 cohorts=cohort,startPoint=chromWidth,method=color.method,color=color)

  
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
  if(color.method=="cohort"){
    legend.color <- GetColor(method="cohort",color=color,cohorts=cohort)
  }
  
  
  if(legend=="missing" || legend==1){
    legend(xtr,legend=unique(cohorts),col=GetColor(color=color,cohorts=cohorts,q=F,method="by.cohort"),cex=0.75,pch=16) # normal legend
  }
  
  
  if(legend==2){
    par(new=T,mar=xtf )
    #par(new=T,mar=c(2,12,10,1))
    pie(table(cohorts),col=legend.color,cex=0.52) # piechart legend
  }
  
  dev.off()
  if(SaveAsObject==TRUE){
    img <- readPNG("t1.png")
    g <- rasterGrob(img, interpolate=TRUE)
    return(g)
  }
  
}


# implement plotCNV



plotCnv <- function(chroms,starts,ends,y,scores,pixel.per.cnv,method,color,score.values,n,startPoint,direction){
  
  indX <- chroms == 'X'
  indY <- chroms == 'Y'
  
  len <- length(starts)
  if(missing(direction)){direction="right"}
  if(direction=="right"){
  # Autosomes
    for(index in 1:len){
      x <- startPoint + pixel.per.cnv*index
      lines(c(x,x),c(y-starts[index],y-ends[index]),col="black",lwd=pixel.per.cnv)
    }
    
    # X Chromosome
    for(index in 1:len) {
      if(indX[index] == TRUE){
        x <- startPoint + pixel.per.cnv*index
        lines(c(x,x),c(y-starts[index],y-ends[index]),col="black",lwd=pixel.per.cnv)
      }
    }
    
    # Y Chromosome
    for(index in 1:len){
      if(indY[index] == TRUE){
        x <- startPoint + pixel.per.cnv*index
        lines(c(x,x),c(y-starts[index],y-ends[index]),col="black",lwd=pixel.per.cnv)
      }
    }
  }else{
    for(index in 1:len){
      x <- startPoint - pixel.per.cnv*index
      lines(c(x,x),c(y-starts[index],y-ends[index]),col="black",lwd=pixel.per.cnv)
    }
    
    # X Chromosome
    for(index in 1:len) {
      if(indX[index] == TRUE){
        x <- startPoint - pixel.per.cnv*index
        lines(c(x,x),c(y-starts[index],y-ends[index]),col="black",lwd=pixel.per.cnv)
      }
    }
    
    # Y Chromosome
    for(index in 1:len){
      if(indY[index] == TRUE){
        x <- startPoint - pixel.per.cnv*index
        lines(c(x,x),c(y-starts[index],y-ends[index]),col="black",lwd=pixel.per.cnv)
      }
    }
  }
}

plotCnv.cohort <- function(chroms,starts,ends,y,chromWidth,pixel.per.cnv,cohorts,startPoint,color,method){
  indX <- chroms == 'X'
  indY <- chroms == 'Y'
  len <- length(starts)
  
  color.value <- GetColor(method=method,color=color,cohorts=cohorts)
  cohort.list <- sort(unique(cohort))
  
  startPoint <- chromWidth
  
  # Autosomes
  for(index in 1:len){
    cohort.index <- match(cohort[index],cohort.list)
    x <- startPoint + pixel.per.cnv*index
    #lines(c(x,x),c(y-starts[index],y-ends[index]),col=cohort[index],lwd=pixel.per.cnv)
    lines(c(x,x),c(y-starts[index],y-ends[index]),
          col=color.value[cohort.index],lwd=pixel.per.cnv)
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




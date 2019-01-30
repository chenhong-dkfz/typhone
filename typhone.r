# get colors

GetColor <- function(method,score,color,score,values,n,q){
  switch(method,
         
         "confidence" = {
           color.value <- "black"
           switch (score,
                   "1" = {color.value = "grey"},
                   "2" = {color.value = "royalblue1"},
                   "3" = {color.value = "royalblue3"},
                   "4" = {color.value = "royalblue4"}
           )
           return(color.value)
         },
         
         "ploidy" = {
           # need redo
           if(score==0){
             color.value = color[1]
           }else{
             color.value = color[score]
           }
           return(color.value)
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
         },
         "cohort2" = {
           print("c2")
         },
         "length" = {
           print("l")
         },
         "factor" = {
           print("f")
         }
  ) 
}







getColor.confidence <- function(score)  #remnat not called in current package
{
  color.value <- "black"
  if (score == 1){
    color.value <- "grey"
  }
  if (score == 2){
    color.value <- "royalblue1"
  }
  if (score == 3){
    color.value <- "royalblue3"
  }
  if (score == 4){
    color.value <- "royalblue4"
  }
  return(color.value)
}

## for CNV.by.ploidy
getColor.ploidy <- function(score,color,score.values,n)
{
  if (score == 0){
    color.value <-  color[1]
  }
  for(i in 1:n)
  {if (score == i)
  {
    color.value <- color[i]
  }
  }
  return(color.value)
}

## for CNV.by.ploidy2
getColor.ploidy2 <- function(score,color,score.values,n)
{
  for(i in 1:n)
  {if (score == score.values[i])
  {
    color.value <- color[i] # variable
  }
  }
  if (missing(score)){
    color.value <- color[n+1] #unkown ploidy
  }
  return(color.value)
}

getColor.cohort <- function(cohorts,color,q)
{
  if(missing(color))
  {color1 <- colorRampPalette(c("red2","indianred4","royalblue4","steelblue1","chartreuse3","darkgreen"))(length(unique(cohorts)))
  } # default
  else
  {
    if(color==1)
    {color1 <- rainbow(length(unique(cohorts)))
    } # rainbow <- 1
    if(color==2)
    {color1 <- colorRampPalette(c("seashell2","seagreen2","turquoise2","palevioletred2"))(length(unique(cohorts)))
    } # light <- 2
    if(color==3)
    {color1 <- colorRampPalette(c("gray7","mediumblue","deeppink4","sienna3"))(length(unique(cohorts)))
    } # dark <- 3
    if(color==4)
    {color1 <- colorRampPalette(c("darkslateblue","darkslategray4","deeppink4","tan4","gray10"))(length(unique(cohorts)))
    } # darkmix <- 4
    if(color==5)
    {color1 <- colorRampPalette(c("navajowhite3","orange3","olivedrab3"))(length(unique(cohorts)))
    } # desert <- 5
    if(color==6)
    {color1 <- colorRampPalette(c("royalblue2","yellow1"))(length(unique(cohorts)))
    } # btoy <- 6
    if(color==7)
    {color1 <-gray.colors((length(unique(cohorts))))
    } # greyscale <- 7
  }
  color <- palette(color1)
  if(q==TRUE)
  {
    return(color1)
  }
  else
  {
    return(color)
  }
}

getColor.length <- function(score,color,score.values,n)
{
  color.value <- color
}

getColorFactor <- function(score) #remnat not called in current package
{
  color.value <- "black"
  
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



# set generic plotCNV(plot single CNV and CNV twin)

# implement of plotCNV

# set generic plotCNVs(plot single CNV and CNV twin)


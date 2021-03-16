library(data.table)
library(ggplot2)
library(stringr)
library(methods)
library(plotly)
library(gtools)

setwd('~/Documents/PING/PING_manuscript/manuscript_scripts/')

referenceDirectory <- 'reference_imputation'

alignedFilledReferencePath <- normalizePath(file.path(referenceDirectory,'KIR_gen_onelines_aligned_imputed.fasta'),mustWork=T)
baseReferencePath <- normalizePath(file.path(referenceDirectory,'KIR_gen_onelines.fasta'),mustWork=T)


kirLocusList <- c('KIR2DL1','KIR2DL2','KIR2DL3','KIR2DL4','KIR2DL5A','KIR2DL5B',
                  'KIR2DS1','KIR2DS2','KIR2DS3','KIR2DS4','KIR2DS5','KIR2DP1',
                  'KIR3DL1','KIR3DL2','KIR3DL3','KIR3DS1','KIR3DP1')

## This function returns a list of dataframes of allele sequences found in the reference fasta, KIR2DL5A and KIR2DL5B are different
read.kir_allele_sequence_from_reference_fasta <- function(fasta_path, kirLocusList){
  
  ## Make sure the fasta can be found
  fasta_path <- normalizePath(fasta_path, mustWork=T)
  
  ## Initialize a list of loci to store theh sequence strings from the file
  locusList <- as.list(kirLocusList)
  names(locusList) <- kirLocusList
  
  ## Initialize a list to store the sequence strings from the file
  for(locusName in names(locusList)){
    locusList[[locusName]] <- list()
  }
  
  ## Read in the fasta file and store the allele names and sequences
  for(currentLine in readLines(fasta_path)){
    alleleNameBool <- grepl('>',currentLine,fixed=T)
    
    if(alleleNameBool){
      alleleName <- strsplit(currentLine, '>',fixed=TRUE)[[1]][2]
      locusName <- strsplit(alleleName,'*',fixed=T)[[1]]
      
    }else{
      alleleSeq <- currentLine
      locusList[[locusName]][[alleleName]] <- alleleSeq
      #alleleSeqList[[alleleName]] <- alleleSeq
    }
  }
  
  return(locusList)
}

alignedFilledAlleleList <- read.kir_allele_sequence_from_reference_fasta(alignedFilledReferencePath, kirLocusList)
baseAlleleList <- read.kir_allele_sequence_from_reference_fasta(baseReferencePath, kirLocusList)


basePercVect <- unlist(lapply(baseAlleleList,function(x){
  lengthList <- lapply(x,nchar)
  maxLength <- max(unlist(lengthList))
  (unlist(lengthList)/maxLength)*100
  }),recursive=T,use.names = F)

p1 <- plot_ly(x = basePercVect, type = "histogram")

## Format the layout of the graph
p1 <- layout(p1,
             title='Filled Reference',
             xaxis = list(title='Fill Percentage',showgrid=F,range=c(0,100)),
             yaxis = list(title='Number of alleles',showgrid=F,range=c(0,600)))
print(p1)


alignFillPercVect <- unlist(lapply(alignedFilledAlleleList,function(x){
  lengthList <- lapply(x,function(y){
    #cat('\n',str_count(y,fixed('N')))
    #newY <- str_remove(y,fixed('N'))
    return(nchar(y)-str_count(y,fixed('N')))
    })
  maxLength <- max(unlist(lengthList))
  (unlist(lengthList)/maxLength)*100
}),recursive=T,use.names = F)

p2 <- plot_ly(x = alignFillPercVect, type = "histogram")

## Format the layout of the graph
p2 <- layout(p2,
             title='Filled Reference',
             xaxis = list(title='Fill Percentage',showgrid=F,range=c(0,100)),
             yaxis = list(title='Number of alleles',showgrid=F,range=c(0,600)))
print(p2)

p <- plot_ly(alpha = 0.6) %>%
  add_histogram(x = round(basePercVect,digits=0)) %>%
  add_histogram(x = round(alignFillPercVect,digits=0)) %>%
  layout(barmode = "overlay",
         xaxis = list(title='Allele fill percentage',showgrid=F,range=c(0,100)),
         yaxis = list(title='Number of alleles',showgrid=F,range=c(0,600)))

print(p)

hist(round(basePercVect,digits=0),breaks = 20,main = 'Base allele completion percentage',xlab = 'Percentage fill',ylab='Number of alleles')
hist(round(alignFillPercVect,digits=0),breaks=20, main='Filled allele completion percentage',xlab='Percentage fill',ylab='Number of alleles',xlim=c(0,100))




ggplot(round(basePercVect,digits=0)) + geom_histogram()

dumDF <- as.data.frame(round(alignFillPercVect, digits=0))
colnames(dumDF)  <- 'basePer'

ggplot(dumDF, aes(x=basePer)) + geom_histogram()

ggsave(
  "ggtest.png",
  ggplot_alternative(),
  width = 3.3,
  height = 3.3,
  dpi = 1200
)

library(extrafont)
font_import()
loadfonts(device="pdf")
windowsFonts(Times=windowsFont("TT Times New Roman"))

ggplot_alternative <- function()
{
  ggplot(dumDF, aes(x=basePer)) + geom_histogram(binwidth=5, color='black',fill='white',boundary=0,size=0.35)+
    xlab('Percentage characterized') +
    ylab('Number of alleles') +
    theme_classic() +
    theme(axis.text.x = element_text(colour='black'),
          axis.text.y = element_text(colour='black')) +
    xlim(0,100)
}
#theme(axis.text.y = element_text(family, face, colour, size))
ggsave(
  "figure_3B.png",
  ggplot_alternative(),
  width = 3.25,
  height = 2,
  dpi = 1200
)


for(id in unique(fileVect)){
  cat(paste0('\n',id))
}



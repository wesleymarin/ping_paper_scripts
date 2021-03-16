library(data.table)
library(stringr)
library(methods)
library(plotly)
library(circlize)

setwd('~/Documents/PING/PING_manuscript/manuscript_scripts/')

source('support/general_functions.R')
source('support/general_functions.R')

## Set up results directory
resultsDirectory <- '3_script_results/Figures/'
dir.create(resultsDirectory,recursive = T,showWarnings = F)
resultsDirectory <- normalizePath(resultsDirectory,mustWork=T)

## Define path to reference allele file
kirReferenceFasta <- 'reference_imputation/KIR_gen_onelines_aligned_imputed.fasta'

## Read in reference alleles
kirLocusAlleleSequenceList <- read.kir_allele_sequence_from_reference_fasta(kirReferenceFasta,kirLocusList)



# Shared k-mer bar plot ---------------------------------------------------
#Figure 2A

kmerLengthVect <- c(50,150,250)
kmerList <- list()
for(kmerLength in kmerLengthVect){
  cat('\nKmer size =',kmerLength)
  locusKmerList <- generate.kir_locus_kmer_list(kirLocusAlleleSequenceList, kmerLength)
  
  matchedOffLocusSumList <- list()
  for(checkLocus in names(locusKmerList)){
    cat('\n',checkLocus,'\tUnique Kmers = ',length(locusKmerList[[checkLocus]]))
    matchedOffLocusBoolList <- sapply(locusKmerList[[checkLocus]], function(x){analyze.find_off_locus_matches(checkLocus,x,locusKmerList)})
    
    matchedOffLocusSumList[[checkLocus]] <- list('totalOff' = sum(matchedOffLocusBoolList == TRUE),'totalKmers' = length(matchedOffLocusBoolList))
  }
  kmerList[[as.character(kmerLength)]] <- matchedOffLocusSumList
}


identityData <- data.frame(matrix(0,nrow=length(names(matchedOffLocusSumList)),ncol=7))
colnames(identityData) <- c('locus','sharedRatio.50','sharedRatio.100','sharedRatio.150','sharedRatio.200',
                            'sharedRatio.250','sharedRatio.300')
rownames(identityData) <- names(matchedOffLocusSumList)

identityData[,1] <- names(kmerList[[as.character(kmerLength)]])

for(kmerLength in names(kmerList)){
  matchedOffLocusSumList <- kmerList[[kmerLength]]
  
  for(locusName in names(matchedOffLocusSumList)){
    offInt <- as.integer(matchedOffLocusSumList[[locusName]]$totalOff)
    totalInt <- as.integer(matchedOffLocusSumList[[locusName]]$totalKmers)
    
    identityData[locusName,paste0('sharedRatio.',kmerLength)] <- round(offInt/totalInt, 2)
  }
}

t <- list(
  family = "calibri",
  size = 48,
  color = 'black')

p <- plot_ly(identityData, x=~locus, y = ~sharedRatio.50, type = 'bar', name = '50-mer') %>%
  add_trace(y=~sharedRatio.150, name = '150-mer') %>%
  add_trace(y=~sharedRatio.250, name = '250-mer') %>%
  layout(title='',
         yaxis = list(title='Shared ratio', barmode='group', showgrid = TRUE, showline = FALSE, showticklabels = TRUE,gridwidth=1,gridcolor='black'),
         xaxis = list(title='', zeroline = TRUE, showline = FALSE, showticklabels = TRUE, showgrid = FALSE, tickangle=-90),
         font=t,
         legend=list(x=0.87,y=0.95))

print(p)
htmlwidgets::saveWidget(p, file=file.path(resultsDirectory,paste0('shared_ratio_bars.html')))

# KIR k-mer Circlize plot -----------------------------------------------------------
# Figure 2B

kmerLength <- 150

locusKmerList <- generate.kir_locus_kmer_list(kirLocusAlleleSequenceList, kmerLength)

uniqueKmerList <- unique(unlist(locusKmerList))

matchedLocusList <- sapply(uniqueKmerList, function(x) analyze.return_all_locus_matches(x,locusKmerList))

concatMatchedLocusList <- sapply(matchedLocusList,paste0,collapse='-')
locusMatchTable <- table(concatMatchedLocusList)
locusMatchTable <- locusMatchTable[order(locusMatchTable, decreasing=T)]

adjFrame <- data.frame(matrix(0,nrow=length(kirLocusList),ncol=length(kirLocusList)),row.names=kirLocusList)
colnames(adjFrame) <- kirLocusList

for(connectionName in names(locusMatchTable)){
  connectionValue <- locusMatchTable[[connectionName]]
  
  locusNameVect <- unlist(strsplit(connectionName,'-',fixed=T))
  if(length(locusNameVect) == 1){
    adjFrame[locusNameVect,locusNameVect] <- as.integer(connectionValue)
  }else if(length(locusNameVect) == 2){
    adjFrame[locusNameVect[1],locusNameVect[2]] <- as.integer(adjFrame[locusNameVect[1],locusNameVect[2]]) + as.integer(connectionValue)
    adjFrame[locusNameVect[2],locusNameVect[1]] <- as.integer(adjFrame[locusNameVect[2],locusNameVect[1]]) + as.integer(connectionValue)
  }else{
    combMat <- combn(locusNameVect,2)
    for(i in ncol(combMat)){
      combLocusNameVect <- combMat[,i]
      adjFrame[combLocusNameVect[1],combLocusNameVect[2]] <- as.integer(adjFrame[combLocusNameVect[1],combLocusNameVect[2]]) + as.integer(connectionValue)
      adjFrame[combLocusNameVect[2],combLocusNameVect[1]] <- as.integer(adjFrame[combLocusNameVect[2],combLocusNameVect[1]]) + as.integer(connectionValue)
    }
  }
}

grid.col <- c(KIR3DP1='#FFB300',
              KIR2DS5='#803E75',
              KIR2DL3='#FF6800',
              KIR2DP1='#A6BDD7',
              KIR2DS3='#C10020',
              KIR2DS2='#CEA262',
              KIR2DL4='#817066',
              KIR3DL3='#007D34',
              KIR3DL1='#F6768E',
              KIR3DS1='#00538A',
              KIR2DL2='#FF7A5C',
              KIR3DL2='#53377A',
              KIR2DS4='#B32851',
              KIR2DL1='#7F180D',
              KIR2DS1='#593315',
              KIR2DL5A='#F13A13',
              KIR2DL5B='#232C16')

colFrame <- data.frame(matrix(0,nrow=nrow(adjFrame),ncol=ncol(adjFrame)))
rownames(colFrame) <- rownames(adjFrame)
colnames(colFrame) <- colnames(adjFrame)

for(locusName in rownames(colFrame)){
  colFrame[locusName,] <- paste0(grid.col[[locusName]],'80')
  colFrame[locusName,locusName] <- paste0(grid.col[[locusName]],'00')
  adjFrame[locusName,locusName] <- as.integer(adjFrame[locusName,locusName]/2)
}

png(file=file.path(resultsDirectory,paste0('circlize_',kmerLength,'-mer.png')),width=4000,height=4000)

chordDiagram(as.matrix(adjFrame),self.link=1,symmetric = T, keep.diagonal = T,grid.col=grid.col,col=as.matrix(colFrame),
             annotationTrack = c('grid'), 
             annotationTrackHeight = convert_height(c(56), "mm"),
             scale=F,
             preAllocateTracks = list(track.height = convert_height(120,'mm')))

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex=8)
}, bg.border = NA)

dev.off()

#### /Figure 2

#### Generate percent identity tables
kmerLengthVect <- c(50,150,250)
for(kmerLength in kmerLengthVect){
  cat('\nKmer size =',kmerLength)
  locusKmerList <- generate.kir_locus_kmer_list(kirLocusAlleleSequenceList, kmerLength)
    
  matchedLocusTableList <- list()
  for(checkLocus in names(locusKmerList)){
    cat('\n',checkLocus)
    matchedOffLocusList <- sapply(locusKmerList[[checkLocus]], function(x) analyze.return_locus_matches(checkLocus,x,locusKmerList))
    matchedLocusTableList[[checkLocus]] <- table(unlist(matchedOffLocusList))
  }
  
  kmerMatchFrame <- data.frame(matrix(0,nrow=(length(names(matchedLocusTableList))+1),ncol=(1+length(names(matchedLocusTableList)))))
  colnames(kmerMatchFrame) <- c('N',names(matchedLocusTableList))
  rownames(kmerMatchFrame) <- c('locus',names(matchedLocusTableList))
  
  kmerMatchFrame[1,] <- c('N',names(matchedLocusTableList))
  
  for(locus in names(matchedLocusTableList)){
    matchedTable <- matchedLocusTableList[[locus]]
    kmerMatchFrame[locus,names(matchedTable)] <- matchedTable#round(matchedTable/sum(matchedTable), 3)*100
    kmerMatchFrame[locus,'N'] <- length(locusKmerList[[locus]])
  }
  
  write.csv(kmerMatchFrame,file=file.path(resultsDirectory,paste0('unique_',kmerLength,'-mer_match_table.csv')))
}

library(data.table)
library(stringr)
library(methods)
library(plotly)
library(circlize)

resultsDirectory <- 'discordance_analysis/'

# ----- Pull misaligned read counts from alignment files -----
'This sections needs to be run in the results directory for the synthetic sequence dataset. Due to the size of the SAM files, this data is not included.
We have included the output of this section, syn.incoming.csv and syn.outgoing.csv, which are the recorded read counts.'

syn.mismatchID.list <- run.read_list('~/synMismatchIDs.txt')
incoming.count.list <- list()
outgoing.count.list <- list()

for( currentLocus in names(syn.mismatchID.list) ){
  cat('\n',currentLocus)
  locus.mismatchID.vect <- syn.mismatchID.list[[currentLocus]]
  
  for( currentID in locus.mismatchID.vect ){
    cat('\n\t',currentID)
    currentSample <- sampleList[[paste0(currentID,'.paired')]]
    #exonNameVect <- grep('E',kirLocusFeatureNameList[[currentLocus]],value=T)
    #exonNameVect <- grep('PE',exonNameVect,value=T,invert = T)
    
    sample.samPath <- file.path( paste0( "/home/wmarin/3_condensed_gc_align/alignmentFiles/", currentSample$name, "/iterAlign/iter_1/",currentSample$name,".sam"))
    uniqueSamDT <- alleleSetup.process_samDT( sample.samPath, delIndex.list, processSharedReads = F, 2 )
    #currentSnpDF <- read.csv(sample.snpPath,
    #                         check.names=F,stringsAsFactors = F,row.names=1,header = T,colClasses = c("character"))
    
    locus.readName.vect <- uniqueSamDT$read_name[ grepl( currentLocus, uniqueSamDT$read_name, fixed=T ) ]
    
    outgoingReads.tab <- table( uniqueSamDT[read_name %in% locus.readName.vect,]$locus )
    incomingReads.tab <- table( tstrsplit( uniqueSamDT[ uniqueSamDT$locus == currentLocus, ]$read_name, '*', fixed=T)[[1]] )
    incoming.count.list[[ paste0(currentLocus,'_',currentID) ]] <- incomingReads.tab
    outgoing.count.list[[ paste0(currentLocus,'_',currentID) ]] <- outgoingReads.tab
    
    #incoming.count.list[ names( incomingReads.tab ) ] <- as.integer( incomingReads.tab ) + unlist( incoming.count.list[ names( incomingReads.tab ) ] )
    #outgoing.count.list[ names( outgoingReads.tab ) ] <- as.integer( outgoingReads.tab ) + unlist( outgoing.count.list[ names( outgoingReads.tab ) ] )
  }
}


output.incomingCount.list <- vector(mode = "list", length = length(names(kirLocusFeatureNameList)) )
names(output.incomingCount.list) <- names(kirLocusFeatureNameList)
output.outgoingCount.list <- vector(mode = "list", length = length(names(kirLocusFeatureNameList)) )
names(output.outgoingCount.list) <- names(kirLocusFeatureNameList)

output.incoming.df <- data.frame(matrix(0,nrow=length(names(kirLocusFeatureNameList)),ncol=length(names(kirLocusFeatureNameList))))
colnames(output.incoming.df) <- names(kirLocusFeatureNameList)
rownames(output.incoming.df) <- names(kirLocusFeatureNameList)

for( locus.sampleID in names(incoming.count.list) ){
  elem.vect <- strsplit( locus.sampleID, '_', fixed=T )[[1]]
  currentLocus <- elem.vect[1]
  sampleID <- elem.vect[2]
  
  locusCount.tab <- incoming.count.list[[locus.sampleID]]
  locusID.vect <- names(locusCount.tab)
  
  locusCount.vect <- setNames( as.vector( locusCount.tab ), locusID.vect )
  
  if('KIR2DL5A' %in% locusID.vect | 'KIR2DL5B' %in% locusID.vect){
    
    if( 'KIR2DL5A' %in% locusID.vect & 'KIR2DL5B' %in% locusID.vect ){
      locusCount.vect <- c(locusCount.vect, setNames((locusCount.vect['KIR2DL5A'] + locusCount.vect['KIR2DL5B']),'KIR2DL5'))
    }else{
      locusCount.vect <- c(locusCount.vect, setNames( locusCount.vect[grepl('KIR2DL5',names(locusCount.vect))], 'KIR2DL5' ))
    }
    
    locusCount.vect <- locusCount.vect[ !grepl('KIR2DL5.',names(locusCount.vect)) ]
    locusID.vect <- names(locusCount.vect)
  }
  
  output.incoming.df[currentLocus,locusID.vect] <- output.incoming.df[currentLocus,locusID.vect] + as.integer(locusCount.vect)
}


output.outgoing.df <- data.frame(matrix(0,nrow=length(names(kirLocusFeatureNameList)),ncol=length(names(kirLocusFeatureNameList))))
colnames(output.outgoing.df) <- names(kirLocusFeatureNameList)
rownames(output.outgoing.df) <- names(kirLocusFeatureNameList)

for( locus.sampleID in names(outgoing.count.list) ){
  elem.vect <- strsplit( locus.sampleID, '_', fixed=T )[[1]]
  currentLocus <- elem.vect[1]
  sampleID <- elem.vect[2]
  
  locusCount.tab <- outgoing.count.list[[locus.sampleID]]
  locusID.vect <- names(locusCount.tab)
  
  locusCount.vect <- setNames( as.vector( locusCount.tab ), locusID.vect )
  
  if('KIR2DL5A' %in% locusID.vect | 'KIR2DL5B' %in% locusID.vect){
    
    if( 'KIR2DL5A' %in% locusID.vect & 'KIR2DL5B' %in% locusID.vect ){
      locusCount.vect <- c(locusCount.vect, setNames((locusCount.vect['KIR2DL5A'] + locusCount.vect['KIR2DL5B']),'KIR2DL5'))
    }else{
      locusCount.vect <- c(locusCount.vect, setNames( locusCount.vect[grepl('KIR2DL5',names(locusCount.vect))], 'KIR2DL5' ))
    }
    
    locusCount.vect <- locusCount.vect[ !grepl('KIR2DL5.',names(locusCount.vect)) ]
    locusID.vect <- names(locusCount.vect)
  }
  
  output.outgoing.df[currentLocus,locusID.vect] <- output.outgoing.df[currentLocus,locusID.vect] + as.integer(locusCount.vect)
}

write.csv(output.incoming.df, file = '~/syn.incoming.csv')
write.csv(output.outgoing.df, file = '~/syn.outgoing.csv')



# ----- Graph mismapped read counts -----
incoming.df <- read.csv('discordance_analysis/syn.incoming.csv',check.names=F,row.names=1)
outgoing.df <- read.csv('discordance_analysis/syn.outgoing.csv',check.names=F,row.names=1)

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
              KIR2DL5='#F13A13')

colFrame <- data.frame(matrix(0,nrow=nrow(incoming.df),ncol=ncol(incoming.df)))
rownames(colFrame) <- rownames(incoming.df)
colnames(colFrame) <- colnames(incoming.df)

for(locusName in rownames(colFrame)){
  colFrame[locusName,] <- paste0(grid.col[[locusName]],'80')
  colFrame[locusName,locusName] <- paste0(grid.col[[locusName]],'00')
  incoming.df[locusName,locusName] <- 0#as.integer(incoming.df[locusName,locusName]/2)
}

png(file=file.path(resultsDirectory,paste0('syn.incomingMisalignedReads.png')),width=4000,height=4000)
chordDiagram(t(as.matrix(incoming.df)),
             self.link=1,
             keep.diagonal = T,
             grid.col=grid.col,
             col=as.matrix(colFrame),
             annotationTrack = c('grid'),
             annotationTrackHeight = 0.07,
             convert_height(c(56, 54), "mm"),
             scale=F,
             preAllocateTracks = list(track.height = convert_height(120,'mm')))

for( locusName in rownames(incoming.df) ){
  #circos.axis(sector.index=locusName,h='bottom',direction = "inside",
  #            labels.facing = "reverse.clockwise",labels.cex = 0.6)
  circos.axis( sector.index=locusName, labels.cex=5, major.tick.length = mm_y(15), lwd=mm_x(0.15))
}


circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(-0.3, 0.5),cex=8,font = 3)
}, bg.border = NA)

dev.off()



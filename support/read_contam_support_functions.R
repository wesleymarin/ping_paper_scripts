kir.locus.vect <- c("KIR3DP1","KIR2DS5","KIR2DL3","KIR2DP1","KIR2DS3","KIR2DS2","KIR2DL4","KIR3DL3","KIR3DL1","KIR3DS1","KIR2DL2","KIR3DL2","KIR2DS4","KIR2DL1","KIR2DS1","KIR2DL5")
kirLocusList <- c('KIR2DL1','KIR2DL2','KIR2DL3','KIR2DL4','KIR2DL5A','KIR2DL5B',
                  'KIR2DS1','KIR2DS2','KIR2DS3','KIR2DS4','KIR2DS5','KIR2DP1',
                  'KIR3DL1','KIR3DL2','KIR3DL3','KIR3DS1','KIR3DP1')

kirLocusFeatureNameList <- list()
kirLocusFeatureNameList[['KIR2DL1']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DL2']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DL3']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DL4']] <- c("5UTR","E1","I1","E2","I2","E3","I3/4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DL5']] <- c("5UTR","E1","I1","E2","I2","E3","I3/4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DP1']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DS1']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DS2']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DS3']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DS4']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DS5']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR3DL1']] <- c("5UTR","E1","I1","E2","I2","E3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR3DL2']] <- c("5UTR","E1","I1","E2","I2","E3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR3DL3']] <- c("5UTR","E1","I1","E2","I2","E3","I3","E4","I4","E5","I5/6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR3DP1']] <- c("5UTR","E1","I1","E2","I2","E3","I3","E4","I4","E5","3UTR")
kirLocusFeatureNameList[['KIR3DS1']] <- c("5UTR","E1","I1","E2","I2","E3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")


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
    cat('\n',alleleNameBool)
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

## Format reads for snp phasing (fill in deletion positions with '.')
allele.read_formatter <- function(readSeq, startPos, endPos, refAlleleName){
  
  ## Format the start position
  startPos <- as.numeric(startPos)
  
  ## Format the end position
  endPos <- as.numeric(endPos)
  
  ## Testing if there deletions found in the reference for the current read
  #if((endPos-startPos+1) > nchar(readSeq)){
  #  
  #  delIndexList <- which(startPos:endPos %in% currentDelIndex[[refAlleleName]])
  #  
  #  ## If the index list is empty, something went wrong
  #  if(length(delIndexList) == 0){
  #    stop('Something weird happened with the deletion index.')
  #  }
  #  
  #  ## Split apart the read sequence and add in deletion symbols
  #  subStringList <- c()
  #  for(delIndex in delIndexList){
  #    preString <- substr(readSeq, 1, delIndex-1)
  #    postString <- substr(readSeq, delIndex, nchar(readSeq))
  #    
  #    readSeq <- paste0(preString, '.', postString)
  #  }
  #}
  
  ## Create a full index that cooresponds to the sequence nucleotides and where they go
  fullIndex <- startPos:endPos
  
  ## Turn the sequence string into a list
  seqList <- strsplit(readSeq,'')[[1]]
  
  if(length(seqList) != length(fullIndex)){
    cat('\n',startPos,endPos,refAlleleName)
  }
  
  names(seqList) <- fullIndex
  
  seqListList <- list(seqList)
  return(seqListList)
}

## Returns a vector of allele names that matched checkSeq that are not from checkLocus
analyze.find_locus_matches <- function(offendingNucStr, kirLocusAlleleSequenceList){
  
  locusMatchList <- list()
  for(currentLocus in names(kirLocusAlleleSequenceList)){
    
    checkSeq <- offendingNucStr
    
    alleleMatchVect <- names(kirLocusAlleleSequenceList[[ currentLocus ]])[ grepl(checkSeq, kirLocusAlleleSequenceList[[ currentLocus ]]) ]
    
    if( length(alleleMatchVect) > 0 ){
      formattedAlleleVect <- paste0( currentLocus, '*', alleleMatchVect )
      locusMatchList[[currentLocus]] <- formattedAlleleVect
    }
    
  }
  return(locusMatchList)
}

generate.kmer_list <- function( idStr, snpKmerStr, kSize ){
  
  cat('\nCreating kmer library for:',idStr)
  
  ## Determine how many kmers will be made
  maxIter <- nchar(snpKmerStr)-kSize+1
  
  ## Generate the kmer vector
  alleleKmerList <- sapply(1:maxIter, function(x){substr(snpKmerStr,x,x+kSize-1)})
  
  names(alleleKmerList) <- paste0(idStr,'_',1:maxIter)
  
  alleleKmerList <- alleleKmerList[ !grepl( '*', alleleKmerList, fixed=T ) ]
  
  return( alleleKmerList )
}

## This function reads in a SAM file with no header to a data.table
read.bowtie2_sam_nohd <- function(currentSample, allele_alignment=F, rows_to_skip=0){
  cat("\n\nReading in",currentSample$samPath)
  
  ## Make sure the SAM file can be read in
  sam_path <- normalizePath(currentSample$samPath, mustWork=T)
  
  ## SAM files can have a variable number of column names, the col.names=1:25 is a possible
  ## point of failure if the SAM file has more than 25 columns
  #output.samTable <- read.table(sam_path, sep='\t', col.names=1:25, stringsAsFactors=F, check.names=F, fill=T, skip=rows_to_skip, nrows = 30000000)
  output.samTable <- read.table(sam_path, sep='\t', stringsAsFactors=F, check.names=F, fill=T, skip=rows_to_skip,col.names=1:25)
  
  ## Convert the dataframe to a datatable for faster access
  output.samTable <- as.data.table(output.samTable)
  
  ## Name the columns that are used for downstream analysis
  colnames(output.samTable)[1] <- 'read_name'
  colnames(output.samTable)[3] <- 'reference_name'
  colnames(output.samTable)[4] <- 'ref_pos'
  colnames(output.samTable)[10] <- 'read_seq'
  
  ## Convert the alignment scores into integers
  alignmentScoreList <- as.integer(tstrsplit(output.samTable$`12`, ':', fixed=TRUE)[[3]])
  
  ## Save alignment scores to their own columns
  cat('\nSetting alignment scores.')
  output.samTable[,alignment_score := alignmentScoreList]
  
  ## Convert the allele names into locus names
  if(allele_alignment){
    locusList <- tstrsplit(output.samTable$reference_name, '_', fixed=T)[[1]]
  }else{
    locusList <- tstrsplit(output.samTable$reference_name, '*', fixed=TRUE)[[1]]
  }
  locusList <- tstrsplit(output.samTable$reference_name, '*', fixed=TRUE)[[1]]
  locusList[locusList %in% 'KIR2DL5A'] = 'KIR2DL5'
  locusList[locusList %in% 'KIR2DL5B'] = 'KIR2DL5'
  
  ## Save locus names to their own columns
  cat('\nSetting locus names.')
  output.samTable[,locus:=locusList]
  
  ## Create new column names to store universal coordinates
  output.samTable$startPos <- 0
  output.samTable$endPos <- 0
  
  ## Create a new column to store read lengths
  cat('\nSetting read lengths.')
  output.samTable[,readLen := nchar(read_seq)]
  
  cat('\nRemoving read alignments that do not include alignment scores.')
  
  ## Discard any rows of the samTable that do not have alignment scores (can happen for various reasons)
  output.samTable <- output.samTable[!is.na(alignmentScoreList),]
  
  ## Check if there were any formatting errors with the alignment score conversions
  if(any(is.na(output.samTable$alignment_score))){
    stop("NA found in SAM file alignment scores for ", currentSample$name)
  }
  
  return(output.samTable)
}

## This function counts the number of header lines in a SAM file (using '@SQ')
samfile.count_header_lines <- function(currentSample){
  if(!file.exists(currentSample$samPath)){
    stop('This sam file does not exist')
  }
  
  headerLines <- 0
  con = file(currentSample$samPath, "r")
  
  while(TRUE){
    
    line = readLines(con, n = 1)
    #if(strsplit(line,'\t')[[1]][1] != "@SQ"){
    #  break
    #}
    
    if(substr(line,1,1) != '@'){
      break
    }
    
    headerLines <- headerLines + 1
  }
  close(con)
  
  return(headerLines)
}





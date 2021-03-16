
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

## Returns a vector of allele names that matched checkSeq that are not from checkLocus
analyze.find_off_locus_matches <- function(checkLocus, checkSeq, locusKmerList){
  offLocusMatch <- FALSE
  for(locus in kirLocusList){
    
    if(locus == checkLocus){
      next
    }
    
    if(checkSeq %in% locusKmerList[[locus]]){
      offLocusMatch <- TRUE
    }
  }
  return(offLocusMatch)
}

analyze.return_locus_matches <- function(checkLocus, checkSeq, locusKmerList){
  offLocusMatchVect <- c()
  for(locus in kirLocusList){
    
    if(locus == checkLocus){
      next
    }
    
    if(checkSeq %in% locusKmerList[[locus]]){
      offLocusMatchVect <- c(offLocusMatchVect, locus)
    }
  }
  
  if(length(offLocusMatchVect) == 0){
    return(checkLocus)
  }else{
    return(offLocusMatchVect)
  }
}

analyze.return_all_locus_matches <- function(checkSeq, locusKmerList){
  locusMatchVect <- c()
  for(locus in kirLocusList){
    
    if(checkSeq %in% locusKmerList[[locus]]){
      locusMatchVect <- c(locusMatchVect, locus)
    }
  }
  
  return(locusMatchVect)
}

## Returns a list of unqiue kmers for each locus
generate.kir_locus_kmer_list <- function(kirLocusAlleleSequenceList, kSize){
  
  output.locusUniqueKmerList <- as.list(names(kirLocusAlleleSequenceList))
  names(output.locusUniqueKmerList) <- names(kirLocusAlleleSequenceList)
  
  cat('\nCreating kmer library for:')
  
  ## Iterate through each locus to generate a vector of unique kmers for each locus
  for(locusName in names(kirLocusAlleleSequenceList)){
    cat('',locusName)
    
    ## Initialize the list element to store a vector
    output.locusUniqueKmerList[[locusName]] <- c()
    
    locusKmerVect <- c()
    
    ## Iterate through each allele of a locus to itersect any unique kmers
    for(alleleName in names(kirLocusAlleleSequenceList[[locusName]])){
      
      ## Pull out the allele sequence
      alleleSeq <- kirLocusAlleleSequenceList[[locusName]][[alleleName]]
      
      ## Determine how many kmers will be made
      maxIter <- nchar(alleleSeq)-kSize+1
      
      ## Generate the kmer vector
      alleleKmerVect <- unlist(sapply(1:maxIter, function(x){substr(alleleSeq,x,x+kSize-1)}))
      
      ## Itersect the allele kmer vector with the locus kmer vector, only keeping unique strings
      locusKmerVect <- unique(c(locusKmerVect, alleleKmerVect))
      
      ## Remove all kmers with ambigious characters
      locusKmerVect <- locusKmerVect[!grepl('N',locusKmerVect)]
    }
    
    output.locusUniqueKmerList[[locusName]] <- locusKmerVect
  }
  
  return(output.locusUniqueKmerList)
}

## Returns a list of named kmers for each locus
generate.kir_named_kmer_DT <- function(kirLocusAlleleSequenceList, kSize, inverseDeletionIndexList){
  
  output.locusUniqueKmerList <- as.list(names(kirLocusAlleleSequenceList))
  names(output.locusUniqueKmerList) <- names(kirLocusAlleleSequenceList)
  
  ## Initialize data table for storing kmers w/ metadata
  output.kmerTable <- data.table(
    locus = character(),
    allele = character(),
    pos = numeric(),
    seq = character(),
    convPos = numeric(),
    stringsAsFactors=F
  )
  
  cat('\nCreating kmer library for:')
  
  ## Iterate through each locus to generate a vector of unique kmers for each locus
  for(locusName in names(kirLocusAlleleSequenceList)){
    cat('',locusName)
    
    ## Iterate through each allele of a locus to itersect any unique kmers
    for(alleleName in names(kirLocusAlleleSequenceList[[locusName]])){
    
      ## Pull out the allele sequence
      alleleSeq <- kirLocusAlleleSequenceList[[locusName]][[alleleName]]
      
      ## Pull out the deletion index
      currentDelIndex <- inverseDeletionIndexList[[names(alleleSeq)]]
      
      ## Determine how many kmers will be made
      maxIter <- nchar(alleleSeq)-kSize+1
      
      ## Initialize data table for storing kmers w/ metadata
      kmerTable <- setNames(data.table(matrix('',nrow=maxIter,ncol=5),stringsAsFactors=F),
                            c('locus','allele','pos','seq','convPos'))
      
      ## Generate the kmer vector
      alleleKmerVect <- unlist(sapply(1:maxIter, function(x){substr(alleleSeq,x,x+kSize-1)}))
      
      ## Fill out the kmerTable
      kmerTable$seq <- alleleKmerVect
      kmerTable$pos <- 1:maxIter
      kmerTable$locus <- locusName
      kmerTable$allele <- alleleName
      
      ## Convert the position coordinates using the del index
      kmerTable[,convPos := mapply(function(x){sum(x >= currentDelIndex)+x}, pos)]
      
      ## Merge the allele kmer data table with the output data table
      output.kmerTable <- merge(output.kmerTable,kmerTable,all=T)
    }
  }
  
  return(output.kmerTable)
}
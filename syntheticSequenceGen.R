library(data.table)
library(stringr)

source('support/read_contam_support_functions.R')

# ----- SUPPORTING FUNCTIONS -----

## This function checks to make sure the output of system2 is valid
check.system2_output <- function(system2_output, system2_error){
  
  ## Checking if the attributes of system2_command are NULL, if not the command was not found
  if(!is.null(attributes(system2_output))){
    cat('\n',system2_output)
    stop(system2_error, '. Stopping program.')
  }
}

### Check to make sure bowtie2-build is accessible <- only needed when building a new reference index
artIllumina <- system2('which', c('art_illumina'), stdout=T, stderr=T)
check.system2_output(artIllumina, 'art_illumina not found')

## This function returns the list of alleles found in the reference fasta
general.read_fasta <- function(fasta_path){
  fasta_path <- normalizePath(fasta_path, mustWork=T)
  
  output.alleleList <- list()
  
  for(currentLine in readLines(fasta_path)){
    alleleNameBool <- grepl('>',currentLine,fixed=T)
    
    if(alleleNameBool){
      currentAllele <- strsplit(currentLine, '>',fixed=TRUE)[[1]][2]
      output.alleleList[[ currentAllele ]] <- ''
    }else{
      output.alleleList[[ currentAllele ]] <- paste0( output.alleleList[[currentAllele]], currentLine )
    }
  }
  
  return(output.alleleList)
}

## Returns boolean
is_nuc.noDel <- function(chr){
  nuc_list <- c('A', 'T', 'C', 'G')
  return(as.character(chr) %in% nuc_list)
}

## Returns boolean
is_nuc <- function(chr){
  nuc_list <- c('A', 'T', 'C', 'G', '.')
  return(as.character(chr) %in% nuc_list)
}

synAllele.format_sequence <- function(currentAllele, currentLocus, knownSnpDFList){
  cat('\n',currentAllele)
  Nindex <- grep('*',knownSnpDFList[[currentLocus]][currentAllele,],fixed=T)
  
  # Pull out currentAllele seq DF
  alleleSeqDF <- knownSnpDFList[[currentLocus]][currentAllele,]
  
  # Mark '*' positions
  Nindex <- grep('*',alleleSeqDF,fixed=T)
  
  seqN.list <- list()
  i <- 1
  startIter <- T
  
  for( curElem in Nindex ){
    if(startIter == T){
      seqN.list[[i]] <- curElem
      startIter  <- F
    }else if(curElem == (1+prevElem)){
      seqN.list[[i]] <- c(seqN.list[[i]], curElem)
    }else{
      i <- i+1
      seqN.list[[i]] <- curElem
    }
    
    prevElem <- curElem
  }
  
  # Replace '*' positions with random nucleotide found across other alleles
  #
  # UPDATE: this is now programed to go by tracts of sequential '*' positions, 
  # each tract is replaced independently by fully defined tract from another allele
  #
  # UPDATE: replacement is now done using tracts from a random fully-defined (or maximully defined) sequence, backup is 'N' tract 
  
  #maxPos.int <- ncol(knownSnpDFList[[currentLocus]])
  #fullAllele.vect <- names( which( apply(knownSnpDFList[[currentLocus]], 1, function(row) sum(sapply(row, is_nuc)) == maxPos.int ) ) )
  #donorAllele <- sample(fullAllele.vect,1)
  
  donorAllele <- sample( fillSequence.list[[currentLocus]],1 )
  
  possAllele.vect <- fillSequence.list[[currentLocus]]
  
  for(nTract in seqN.list){
    n.colname.vect <- colnames(knownSnpDFList[[currentLocus]])[nTract]
    
    # donorAlleleVect <- names( which( apply( knownSnpDFList[[currentLocus]][,n.colname.vect,drop=F], 1, function(x){
    #   all(is_nuc(x))
    # }) ) )
    # 
    # donorAllele <- sample(donorAlleleVect,1)
    
    if( any( '*' %in% knownSnpDFList[[currentLocus]][donorAllele,n.colname.vect] ) ){

      donorAlleleVect <- names( which( apply( knownSnpDFList[[currentLocus]][,n.colname.vect,drop=F], 1, function(x){
        all(is_nuc(x))
      }) ) )
      
      donorAlleleVect <- intersect(possAllele.vect, donorAlleleVect)
      
      if(length(donorAlleleVect) == 0){
        stop(currentAllele,nTract)
      }
      
      #donorAllele <- sample(donorAlleleVect,1)
    #  alleleSeqDF[1,n.colname.vect] <- 'N'
    }s
    
    alleleSeqDF[1,n.colname.vect] <- knownSnpDFList[[currentLocus]][donorAllele,n.colname.vect]
  }
  
  # Condense allele seq to a string
  alleleSeq <- paste0(alleleSeqDF[1,],collapse='')
  
  # Remove '.' from seq
  alleleSeq <- gsub('.','',alleleSeq,fixed=T)
  
  return(alleleSeq)
}

synAllele.write_fasta <- function( alleleSeq.list, fastaPath){
  i <- 1
  for( alleleID in names(alleleSeq.list) ){
    if( i == 1){
      cat(paste0('>',alleleID,'\n',alleleSeq.list[[alleleID]],'\n'),file=fastaPath,append=F)
    }else{
      cat(paste0('>',alleleID,'\n',alleleSeq.list[[alleleID]],'\n'),file=fastaPath,append=T)
    }
    i <- i+1
  }
  #cat('\n',file=fastaPath,append=T)
}

synAllele.initialize_info_file <- function(resultsDirectory){
  path <- file.path( resultsDirectory, 'synSeq.info.txt')
  textStr <- ''
  cat(textStr, file=path)
  return(path)
}

synAllele.write_info <- function(synSeq.ID, selectedHaps.vect, allele.vect, infoPath){
  hapStr <- paste0(selectedHaps.vect,collapse='_')
  alleleStr <- paste0(allele.vect,collapse='_')
  cat(paste0(synSeq.ID,'\t',hapStr,'\t',alleleStr,'\n'),file=infoPath,append=T)
}

synAllele.run_art <- function(artIllumina, fastaPath, synSeq.ID, readDirectory, dp, rl){
  outFile <- file.path(readDirectory,paste0(synSeq.ID,'.paired_'))

  optionsCommand <- c('-ss HS25',
                      #'-ef',
                      '-na',
                      paste0('-i ',fastaPath),
                      paste0('-l ',rl),
                      paste0('-f ',dp),
                      '-m 200',
                      '-s 10',
                      paste0('-o ',outFile))
  
  cat('\n\n',artIllumina, optionsCommand)
  output.art <- system2(artIllumina, optionsCommand, stdout=T,stderr=T)
  check.system2_output(output.art, 'art_illumina command failed')
  return(NULL)
}

## Building up the run command
#optionsCommand <- c('view',paste0('-@', threads),
#                    currentSample$mhcSamPath, '-o', currentSample$mhcBamPath)

#cat('\n\n',samtools_command, optionsCommand)
#output.bamConv <- system2(samtools_command, optionsCommand, stdout=T, stderr=T)
#check.system2_output(output.bamConv, 'samtools sam to bam conversion failed')

# ----- RUN CODE ----

resultsDirectory <- 'synSeq_data_run5/'
dir.create(resultsDirectory)
readDirectory <- file.path(resultsDirectory,'reads')
dir.create(readDirectory)

synSeq.generate <- 50
synSeq.dp <- 50
synSeq.rl <- 150
synSeq.label <- 'hiseq'
#annotatedAlleleDirectory <- '/home/LAB_PROJECTS/PING2_PAPER/3_script_results/alleleFiles/'
annotatedAlleleDirectory <- 'synthetic_sequence/extended_SNP_files/'
kirGeneHapsFile <- 'synthetic_sequence/KIR_gene_haplotypes.csv'
kirHaps.DF <- read.csv(kirGeneHapsFile,stringsAsFactors = F, check.names=F,row.names=1)

## Load in the allele SNP dataframes
knownSnpDFList <- list()
for(locus in kir.locus.vect){
  alleleSnpDF <- read.csv( normalizePath(
    file.path(annotatedAlleleDirectory,
              paste0(locus, '_alleleSNPs.csv')
    )), check.names=F, row.names=1, stringsAsFactors = F,colClasses = 'character')
  
  knownSnpDFList[[locus]] <- alleleSnpDF
}

kirHapVect <- rownames(kirHaps.DF)

# Generate list of fully-defined (or max defined) sequence
maxPos.list <- sapply( names(knownSnpDFList), function(x) max( apply(knownSnpDFList[[x]], 1, function(row) sum(sapply(row, is_nuc)) ) ) )
fillSequence.list <- sapply( names(knownSnpDFList), function(x){
  names(which( apply(knownSnpDFList[[x]], 1, function(row) sum(sapply(row, is_nuc)) == maxPos.list[[x]] ) )) 
  namevect.2 <- rownames(knownSnpDFList[[x]])[ nchar( tstrsplit( rownames( knownSnpDFList[[x]] ), '*', fixed=T )[[2]] ) == 7 ]
  unique(c(namevect.1, namevect.2))
  })
sapply( fillSequence.list, length)
# Remove 3DL1*05901
fillSequence.list[['KIR3DL1']] <- fillSequence.list[['KIR3DL1']][fillSequence.list[['KIR3DL1']] != 'KIR3DL1*05901']

# Run code ----------------------------------------------------------------

# Initialize info file
infoPath <- synAllele.initialize_info_file(resultsDirectory)

set.seed(1)

for(synSeq.iter in 1:synSeq.generate){
#synSeq.iter <- 1
  synSeq.ID <- paste0('synSeq.',synSeq.label, '.dp',synSeq.dp,'.rl',synSeq.rl,'.',synSeq.iter)
  cat('\n',synSeq.ID)
  selectedHaps.vect <- sample(kirHapVect,2,replace=T)
  
  cat('\t',selectedHaps.vect,'\t')
  genoToBuild <- apply( kirHaps.DF[selectedHaps.vect,], 2, sum )
  
  master.alleleSeq.list <- list()
  for( currentLocus in names(genoToBuild) ){
    cat('',currentLocus)
    #currentLocus <- 'KIR3DP1'
    currentLocusCopy <- as.integer( genoToBuild[[currentLocus]] )
    
    if(currentLocusCopy == 0){
      next
    }
    
    currentLocusAlleleVect <- rownames(knownSnpDFList[[currentLocus]])
    currentLocusAlleleVect <- fillSequence.list[[currentLocus]]
    
    currentAlleleVect <- sample(currentLocusAlleleVect,currentLocusCopy,replace=F)
    
    alleleSeq.list <- sapply(currentAlleleVect, synAllele.format_sequence, currentLocus, knownSnpDFList)
    
    master.alleleSeq.list <- c(master.alleleSeq.list, alleleSeq.list)
  }
  
  
  fastaPath <- file.path(resultsDirectory,paste0(synSeq.ID,'.fa'))
  
  cat('\nWriting fasta and info.')
  synAllele.write_fasta(master.alleleSeq.list, fastaPath)
  synAllele.write_info(synSeq.ID, selectedHaps.vect, names(master.alleleSeq.list), infoPath)
  synAllele.run_art(artIllumina, fastaPath, synSeq.ID, readDirectory, synSeq.dp, synSeq.rl)
  cat('\n')
}
#synSeq.iter <- synSeq.iter + 1

# 
# 
# hapID <- names(hapSeq.list)[10]
# currentHapStr <- hapSeq.list[[hapID]]
# #currentLocus <- 'KIR3DL3'
# hapID
# 
# sapply(kir.locus.vect, function(currentLocus){
#   #cat('\n',currentLocus)
#   locusFeat.list <- sapply( kirLocusFeatureNameList[[ currentLocus ]], function(curFeat){
#     matchIndex <- tstrsplit( colnames( knownSnpDFList[[currentLocus]] ), '_', fixed=T )[[1]] == curFeat
#     
#     colnames( knownSnpDFList[[currentLocus]] )[matchIndex]
#   })
#   
#   exonFeatNameVect <- grep('E',kirLocusFeatureNameList[[ currentLocus ]],fixed=T,value=T)
#   exonFeatNameVect <- grep('PE', exonFeatNameVect, fixed=T, invert=T, value=T)
#   
#   featNameVect <- grep('/',kirLocusFeatureNameList[[ currentLocus ]],fixed=T,value=T,invert = T)
#   featNameVect <- grep('/',featNameVect,fixed=T,value=T,invert=T)
#   
#   locusExonMatch.list <- sapply( rownames(knownSnpDFList[[currentLocus]]), function(x){
#     #cat('\n\t',x)
#     exonMatch.list <- sapply(featNameVect, function(y){
#       #cat('\t',y)
#       nucStr <- paste0(knownSnpDFList[[currentLocus]][x,locusFeat.list[[y]]], collapse='')
#       
#       if( grepl('*', nucStr, fixed=T) ){
#         nucStr <- gsub('*','',nucStr,fixed=T)
#       }
#       
#       if( grepl('.', nucStr, fixed=T) ){
#         nucStr <- gsub('.','',nucStr,fixed=T)
#       }
#       
#       grepl(nucStr, currentHapStr,fixed=T)
#     })
#     
#     if(any(is.na(exonMatch.list))){
#       return(NA)
#     }else{
#       return(all(exonMatch.list))
#     }
#   })
#   
#   
#   names( which(locusExonMatch.list) )
# })
# 
# 
# 
# 
# 
# 

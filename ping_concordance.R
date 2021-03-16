
library(data.table)
library(stringr)
library(gtools)
library(ggplot2)
#install.packages(c('gtools','data.table','stringr'),dependencies = T)


kirLocusList <- c('KIR3DP1','KIR2DS5','KIR2DL3','KIR2DP1',
                  'KIR2DS3','KIR2DS2','KIR2DL4','KIR3DL3',
                  'KIR3DL1','KIR3DS1','KIR2DL2','KIR3DL2','KIR2DS4','KIR2DL1', 'KIR2DS1', 'KIR2DL5')

#synSeqAnswerKey <- '~/Documents/PING_projects/data/2020_11_16_synSeq_data/dp50rl150_seeded_alleleMatched/synSeq.info.txt'

resultsDirectory <- "concordance/"

## This function takes in a list and writes it to a text file
run.write_list <- function(listToWrite, filePath){
  firstLineBool = T
  
  for(listElemName in names(listToWrite)){
    listElem <- listToWrite[[listElemName]]
    
    write(listElemName, file=filePath, ncolumns=1, append=!firstLineBool, sep=' ')
    
    if(firstLineBool){
      firstLineBool = F
    }
    
    write(listElem, file=filePath, ncolumns=length(listElem), append=!firstLineBool, sep=' ')
  }
  
  return(filePath)
}

## This function takes in a text file reads it in as a list (its really only meant to work with run.write_list)
run.read_list <- function(filePath){
  #filePath <- file.path(assembledReferenceDirectory, paste0(currentSample$name,'assembled_deletion.index'))
  
  listToReturn <- list()
  
  conn <- file(filePath,open="r")
  linn <-readLines(conn)
  
  for (i in 1:length(linn)){
    if(i%%2 == 1){
      elemName <- linn[i]
    }else{
      elemValues <- linn[i]
      listToReturn[[elemName]] <- strsplit(elemValues, ' ', fixed=T)[[1]]
    }
  }
  close(conn)
  
  return(listToReturn)
}


kir.allele_resolution <- function(allele_name, res){
  
  ## Split the locus name from the allele number
  alleleLocusNumber <- str_split(allele_name, fixed('*'))[[1]]
  
  alleleLocus <- alleleLocusNumber[1]
  alleleNumber <- alleleLocusNumber[2]
  
  ## Only return the alleleLocus if the resolution is 0
  if(res == 0){
    return(alleleLocus)
  }
  
  ## Subset the allele_number based on res
  shortAlleleNumber <- substr(alleleNumber, 1, res)
  
  ## Return the shortened allele name
  return(paste0(alleleLocus,'*',shortAlleleNumber))
}

## Condenses allele fram into only positions that have multiple nucleotides
make_unique_pos_frame <- function(pos_frame){
  num_unique <- lapply(pos_frame, num_unique_nuc)
  unique_pos_frame <- pos_frame[,names(num_unique[num_unique > 1]), drop=F]
  return(unique_pos_frame)
}

## Returns boolean
is_nuc <- function(chr){
  nuc_list <- c('A', 'T', 'C', 'G', '.')
  return(as.character(chr) %in% nuc_list)
}

## Counts how many unique nucleotides are in a vector
num_unique_nuc <- function(vector_to_check){
  return(sum(is_nuc(names(table(vector_to_check)))))
}

# Reads in answer key (synSeq.info.txt)
synSeq.readAnswerKey <- function(keyFile){
  
  out.list <- list()
  con  <- file(keyFile, open = "r")
  
  while (length(oneLine <- readLines(con, n = 1)) > 0) {
    lineVect <- unlist((strsplit(oneLine, "\t",fixed=T)))
    sampleID <- lineVect[1]
    genoID <- lineVect[2]
    genoStr <- lineVect[3]
    
    genoVect <- unlist(strsplit(genoStr,'_',fixed=T))
    out.list[[sampleID]] <- list('genoID'=genoID,'genoVect'=genoVect)
  }
  
  close(con)
  
  return(out.list)
}

# Adds copy number to answer key
synSeq.keyCopyNumber <- function(rawKey.list, kirLocusList){
  
  for(i in 1:length(rawKey.list)){
    sampleID <- names(rawKey.list)[i]
    rawKey.list[[i]]$genoVect
    
    keyCopy.list <- sapply(kirLocusList, function(x){ sum(grepl(x,rawKey.list[[i]]$genoVect)) })
    
    rawKey.list[[i]] <- c(rawKey.list[[i]],list('copyList'=keyCopy.list))
  }
  
  return(rawKey.list)
}

# Adds formatted geno calls to answer key
synSeq.keyGeno <- function(rawKey.copy.list){
  
  for(i in 1:length(rawKey.copy.list)){
    sampleID <- names(rawKey.copy.list)[i]
    genoVect <- rawKey.copy.list[[i]]$genoVect
    copyList <- rawKey.copy.list[[i]]$copyList
    
    genoList <- list()
    
    for( currentLocus in names(copyList) ){
      alleleVect <- grep(currentLocus, genoVect, value = T)
      nullToAdd <- max(2-length(alleleVect), 0)
      nullAllele <- paste0(currentLocus,'*null')
      
      nullVect <- replicate(nullToAdd,nullAllele)
      if(length(nullVect) > 0){
        alleleVect <- c(alleleVect, nullVect)
      }
      
      names(alleleVect) <- paste0('allele',1:length(alleleVect))
      
      alleleList <- list(alleleVect)
      names(alleleList) <- currentLocus
      genoList <- c(genoList, alleleList)
    }
    
    #L1type <- grep('KIR3DL1',genoVect,value=T)
    #S1type <- grep('KIR3DS1',genoVect,value=T)
    
    # Combine KIR2DL23 typings
    L2type <- grep('KIR2DL2',genoVect,value=T)
    L3type <- grep('KIR2DL3',genoVect,value=T)
    
    L23Vect <- c(L2type, L3type)
    names(L23Vect) <- paste0('allele',1:length(L23Vect))
    L23List <- list(L23Vect)
    names(L23List) <- 'KIR2DL23'
    genoList <- c(genoList, L23List)
    
    rawKey.copy.list[[i]] <- c(rawKey.copy.list[[i]],list('genoList'=genoList))
    
  }
  return(rawKey.copy.list)  
}

#rawKey.list <- synSeq.readAnswerKey(synSeqAnswerKey)
#rawKey.copy.list <- synSeq.keyCopyNumber(rawKey.list, kirLocusList)
#genoCopyKey.list <- synSeq.keyGeno(rawKey.copy.list)


# ----- Read in results data -----
# Add NULL allele calls and format result genotype dataframe as data.table
synSeq.formatResultGeno <- function( genoCall.df, copyCall.df ){
  
  locusVect <- colnames(genoCall.df)
  sampleIDVect <- rownames(genoCall.df)
  
  genoCall.dt <- as.data.table(genoCall.df)
  genoCall.dt <- as.data.table( matrix( data=list(0),nrow=nrow(genoCall.df),ncol=ncol(genoCall.df)) )
  colnames(genoCall.dt) <- locusVect
  genoCall.dt$sampleID <- sampleIDVect
  
  for( currentLocus in locusVect ){ 
    #cat('\n\n',currentLocus)
    nullAllele <- paste0(currentLocus,'*null')
      
    for( currentID in sampleIDVect ){
      #cat('',currentID)
      currentGeno <- genoCall.df[currentID, currentLocus]
      
      # if(currentLocus == 'KIR2DL23'){
      #   currentCopy <- 2
      # }else if(currentLocus == 'KIR2DS35'){
      #   currentCopy <- as.integer(copyCall.df[currentID,'KIR2DS3']) + as.integer(copyCall.df[currentID,'KIR2DS5'])
      # }else if(currentLocus == 'KIR3DL1S1'){
      #   currentCopy <- as.integer(copyCall.df[currentID,'KIR3DL1']) + as.integer(copyCall.df[currentID,'KIR3DS1'])
      # }else{
      #   currentCopy <- as.integer( copyCall.df[currentID, currentLocus] )
      # }
      
      if( is.na(currentGeno) ){
        cat('\nEncountered NA',currentLocus, currentID)
        stop()
        # nullToAdd <- 2 
        # nullVect <- replicate(nullToAdd,nullAllele)
        # names(nullVect) <- paste0('allele',1:length(nullVect))
        # 
        # alleleMat <- t( as.matrix(nullVect) )
        # 
        # genoCall.dt[sampleID == currentID, currentLocus] <- list( alleleMat )
          
      }else{
        genoVect <- unlist( strsplit(currentGeno, ' ', fixed=T) )
        
        alleleMat <- t( sapply(1:length(genoVect), function(i){
          
          genoCall <- genoVect[i]
          
          alleleVect <- unlist( strsplit(genoCall, '+', fixed=T) )
          
          if( 'failed' %in% alleleVect ){
            alleleVect <- c('failed','failed')
          }
          
          # if( length(unique(alleleVect)) == 1 & currentCopy == 1 ){
          #   cat('\nAdding null',currentLocus,currentID)
          #   nullToAdd <- 1
          #   nullVect <- replicate(nullToAdd,nullAllele)
          #   alleleVect <- c(unique(alleleVect), nullVect)
          # }else if( length(alleleVect) == 1 & currentCopy > 1){
          #   alleleVect <- rep(alleleVect, currentCopy)
          # }
          
          names(alleleVect) <- paste0('allele',1:length(alleleVect))
          #alleleList <- list(alleleVect)
          return( alleleVect )
          
        }) )
        
        genoCall.dt[sampleID == currentID, currentLocus] <- list( alleleMat )
        
      }
      
    }
  }
 return(genoCall.dt) 
}


iterCallPath <- file.path(resultsDirectory,'380_genotypes/finalAlleleCalls_copyFixed.csv')
copyNumberPath <- file.path(resultsDirectory,'380_genotypes/manualCopyNumberFrame.csv')
kffPath <- file.path(resultsDirectory,'380_genotypes/kffPresenceFrame.csv')
validGenoCallPath <- file.path(resultsDirectory,'380_genotypes/INDIGO_controls_genoResults.csv')
validCopyCallPath <- file.path(resultsDirectory,'380_genotypes/INDIGO_controls_copyNumber.csv')

finalGeno.syn.path <- file.path(resultsDirectory,'syn_genotypes/finalAlleleCalls_1.csv')
copy.syn.path <- file.path(resultsDirectory,'syn_genotypes/manualCopyNumberFrame.csv')
valid.syn.path <- file.path(resultsDirectory,'syn_genotypes/synSeq.info.txt')


finalGeno.afr.path <- file.path(resultsDirectory,'afr_genotypes/finalAlleleCalls_feb16.csv')
finalCopy.afr.path <- file.path(resultsDirectory,'afr_genotypes/manualCopyNumberFrame.csv')
validGeno.afr.path <- file.path(resultsDirectory,'afr_genotypes/afr_validGenos.csv')
validCopy.afr.path <- file.path(resultsDirectory,'afr_genotypes/KhoeSan_paul_copy_number.csv')

# ----- Read in and format synthetic validation -----
validGeno.syn.list <- synSeq.readAnswerKey(valid.syn.path)

synSeq.genoListToDF <- function( validGeno.syn.list, kirLocusList ){
  
  output.validGeno.syn.DF <- as.data.frame( matrix('',nrow=length(validGeno.syn.list),ncol=length(kirLocusList) ), stringsAsFactors = F )
  rownames(output.validGeno.syn.DF) <- names(validGeno.syn.list)
  colnames(output.validGeno.syn.DF) <- kirLocusList
  
  output.validCopy.syn.DF <- as.data.frame( matrix('',nrow=length(validGeno.syn.list),ncol=length(kirLocusList) ), stringsAsFactors = F )
  rownames(output.validCopy.syn.DF) <- names(validGeno.syn.list)
  colnames(output.validCopy.syn.DF) <- kirLocusList
  
  for( sampleID in names(validGeno.syn.list) ){
    sample.genoVect <- validGeno.syn.list[[sampleID]]$genoVect
    
    for( currentLocus in kirLocusList){
      locus.genoVect <- grep(currentLocus, sample.genoVect, value=T)
      
      output.validCopy.syn.DF[sampleID, currentLocus] <- length(locus.genoVect)
      
      if(length(locus.genoVect) == 0){
        output.validGeno.syn.DF[sampleID,currentLocus] <- NA
      }else{
        output.validGeno.syn.DF[sampleID,currentLocus] <- paste0( locus.genoVect, collapse = '+' )
      }
    }
    
  }
  
  return(list('genoDF'=output.validGeno.syn.DF,'copyDF'=output.validCopy.syn.DF))
}

validGenoCopy.syn.list <- synSeq.genoListToDF( validGeno.syn.list, kirLocusList )

validGeno.syn.df <- validGenoCopy.syn.list$genoDF
validCopy.syn.df <- validGenoCopy.syn.list$copyDF


# Read in 380 validation
iterCall.df <- read.csv(iterCallPath,
                        stringsAsFactors=F,check.names=F,row.names=1)
copyCall.df <- read.csv(copyNumberPath,
                        stringsAsFactors=F,check.names=F,row.names=1)
kffCall.df <- read.csv(kffPath,
                       stringsAsFactors=F,check.names=F,row.names=1)
validGenoCall.df <- read.csv( validGenoCallPath,
                          stringsAsFactors = F,check.names=F, row.names=1)
validCopyCall.df <- read.csv( validCopyCallPath,
                              stringsAsFactors = F,check.names=F, row.names=1)

# modify sampleID's to match the answer key
#rownames(iterCall.df) <- unlist(strsplit(rownames(iterCall.df),'.paired'))
#rownames(copyCall.df) <- unlist(strsplit(rownames(copyCall.df),'.paired'))
#rownames(kffCall.df) <- unlist(strsplit(rownames(kffCall.df),'.paired'))
badValidID.vect <- c('IND00154')

goodSampleID.vect <- rownames(iterCall.df)[ !apply( iterCall.df, 1, function(x) all( is.na(x) ) ) ]

iterCall.df <- iterCall.df[goodSampleID.vect,]
copyCall.df <- copyCall.df[goodSampleID.vect,]
kffCall.df <- kffCall.df[goodSampleID.vect,]

goodValidID.vect <- intersect( rownames(validGenoCall.df), rownames(validCopyCall.df) )
goodValidID.vect <- goodValidID.vect[ !(goodValidID.vect %in% badID.vect) ]
validGenoCall.df <- validGenoCall.df[goodValidID.vect,]
validCopyCall.df <- validCopyCall.df[goodValidID.vect,]
validCopyCall.df$KIR2DL5 <- validCopyCall.df[,'KIR2DL5A'] + validCopyCall.df[,'KIR2DL5B']

pingFinalize.combineL23 <- function( iterCall.df ){
  
  L23.call.list <- sapply( rownames( iterCall.df ), function(x){
    L2.call <- iterCall.df[x,'KIR2DL2']
    L3.call <- iterCall.df[x,'KIR2DL3']
    
    L2.call.bool <- !is.na(L2.call) & L2.call != 'failed'
    L2.call.unresolved <- grepl('unresolved',L2.call)
    
    L3.call.bool <- !is.na(L3.call) & L3.call != 'failed'
    L3.call.unresolved <- grepl('unresolved',L3.call)
    
    if( L2.call.unresolved ){
      L2.call <- paste0('KIR2DL2*',L2.call)
    }
    if( L3.call.unresolved ){
      L3.call <- paste0('KIR2DL3*',L3.call)
    }
    
    if( L2.call.bool & L3.call.bool ){
      L2.call.vect <- strsplit(L2.call,' ',fixed=T)[[1]]
      L3.call.vect <- strsplit(L3.call,' ',fixed=T)[[1]]
      L23.call.mat <- expand.grid(L2.call.vect,L3.call.vect)
      L23.call.vect <- apply( L23.call.mat, 1, paste0, collapse='+' )
    }else if( L2.call.bool & !L3.call.bool ){
      
      L2.hetBool <- grepl('+',L2.call,fixed=T)
      
      if(!L2.hetBool){
        L2.callVect <- strsplit(L2.call,' ',fixed=T)[[1]]

        '
        Combinations could cause recursion limit problems for highly ambiguous calls
        '
        L23.call.mat <- combinations(length(L2.callVect),2,L2.callVect,repeats.allowed = T)
        L23.call.vect <- apply( L23.call.mat, 1, paste0, collapse='+')
      }else{
        L23.call.vect <- unlist( tstrsplit(L2.call,' ',fixed=T)[[1]] )
      }
      
    }else if( !L2.call.bool & L3.call.bool ){
      L3.hetBool <- grepl('+',L3.call,fixed=T)
      
      if(!L3.hetBool){
        L3.callVect <- strsplit(L3.call,' ',fixed=T)[[1]]
        L23.call.mat <- combinations(length(L3.callVect),2,L3.callVect,repeats.allowed = T)
        L23.call.vect <- apply( L23.call.mat, 1, paste0, collapse='+')
      }else{
        L23.call.vect <-  unlist( tstrsplit(L3.call,' ',fixed=T)[[1]] )
      }
    }else{
      L23.call.vect <- 'failed'
    }
    
    if( length(L23.call.vect) > 64 ){
      L23.call.str <- 'failed'
    }else{
      L23.call.str <- paste0(L23.call.vect, collapse=' ')
    }
    return(L23.call.str)
  })
  
  iterCall.df$KIR2DL23 <- NA
  iterCall.df[names(L23.call.list),'KIR2DL23'] <- unlist( L23.call.list )
  
  return(iterCall.df)
}

pingFinalize.combineS35 <- function( iterCall.df, copyCall.df ){
  S35.call.list <- sapply( rownames( iterCall.df ), function(x){
    cat('\n',x)
    S3.call <- iterCall.df[x,'KIR2DS3']
    S5.call <- iterCall.df[x,'KIR2DS5']
    
    S3.copy <- as.integer( copyCall.df[x,'KIR2DS3'] )
    S5.copy <- as.integer( copyCall.df[x,'KIR2DS5'] )
    
    S3.call.bool <- !is.na(S3.call) & S3.call != 'failed'
    S3.call.unresolved <- grepl('unresolved',S3.call)
    
    S5.call.bool <- !is.na(S5.call) & S5.call != 'failed'
    S5.call.unresolved <- grepl('unresolved',S5.call)
    
    if( S3.call.unresolved ){
      S3.call <- paste0('KIR2DS3*',S3.call)
    }
    if( S5.call.unresolved ){
      S5.call <- paste0('KIR2DS5*',S5.call)
    }
    
    if( S3.call.bool ){
      S3.het.bool <- grepl('+',S3.call,fixed=T)
    }
    
    if( S5.call.bool ){
      S5.het.bool <- grepl('+',S5.call,fixed=T)
    }
    
    if( S3.call.bool & S5.call.bool ){
      
      if( !S3.het.bool ){
        S3.allele1.vect <- strsplit(S3.call,' ',fixed=T)[[1]]
        
        if( S3.copy > 1 ){
          S3.call.mat <- combinations(length(S3.allele1.vect), S3.copy, S3.allele1.vect,repeats.allowed = T)
          S3.call.vect <- apply( S3.call.mat, 1, paste0, collapse='+')
        }else{
          S3.call.vect <- S3.allele1.vect
        }
      }else{
        S3.hetCall.vect <- tstrsplit(S3.call,' ',fixed=T)
        #S3.hetCall.list <- tstrsplit(S3.hetCall.vect,'+',fixed=T)
        #S3.allele1.vect <- S3.hetCall.list[[1]]
        #S3.allele2.vect <- S3.hetCall.list[[2]]
        
        #S3.call.vect <- apply( unique( expand.grid(S3.hetCall.list) ), 1, function(x) paste0(x,collapse='+'))
        S3.call.vect <- unlist( S3.hetCall.vect )
      }
      
      if( !S5.het.bool ){
        S5.allele1.vect <- strsplit(S5.call,' ',fixed=T)[[1]]
        if( S5.copy > 1 ){
          S5.call.mat <- combinations(length(S5.allele1.vect), S5.copy, S5.allele1.vect,repeats.allowed = T)
          S5.call.vect <- apply( S5.call.mat, 1, paste0, collapse='+')
        }else{
          S5.call.vect <- S5.allele1.vect
        }
      }else{
        S5.hetCall.vect <- tstrsplit(S5.call,' ',fixed=T)
        #S5.hetCall.list <- tstrsplit(S5.hetCall.vect,'+',fixed=T)
        #S5.allele1.vect <- S5.hetCall.list[[1]]
        #S5.allele2.vect <- S5.hetCall.list[[2]]
        
        #S5.call.vect <- apply( unique( expand.grid(S5.hetCall.list) ), 1, function(x) paste0(x,collapse='+'))
        S5.call.vect <- unlist( S5.hetCall.vect )
      }
      
      S35.call.mat <- expand.grid(S3.call.vect,S5.call.vect)
      S35.call.vect <- apply( S35.call.mat, 1, paste0, collapse='+')
    }else if( S3.call.bool & !S5.call.bool ){
      if( !S3.het.bool ){
        S3.allele1.vect <- strsplit(S3.call,' ',fixed=T)[[1]]
        
        if( S3.copy > 1 ){
          S3.call.mat <- combinations(length(S3.allele1.vect), S3.copy, S3.allele1.vect,repeats.allowed = T)
          S35.call.vect <- apply( S3.call.mat, 1, paste0, collapse='+')
        }else{
          S35.call.vect <- apply( expand.grid( S3.allele1.vect, 'KIR2DS35*null' ), 1, paste0, collapse='+' )
        }
      }else{
        S3.hetCall.vect <- tstrsplit(S3.call,' ',fixed=T)
        #S3.hetCall.list <- tstrsplit(S3.hetCall.vect,'+',fixed=T)
        
        #S35.call.vect <- apply( unique( expand.grid(S3.hetCall.list) ), 1, function(x) paste0(x,collapse='+') )
        S35.call.vect <- unlist( S3.hetCall.vect )
      }
    }else if(!S3.call.bool & S5.call.bool){
      if( !S5.het.bool ){
        S5.allele1.vect <- strsplit(S5.call,' ',fixed=T)[[1]]
        
        if( S5.copy > 1 ){
          S5.call.mat <- combinations(length(S5.allele1.vect), S5.copy, S5.allele1.vect,repeats.allowed = T)
          S35.call.vect <- apply( S5.call.mat, 1, paste0, collapse='+')
        }else{
          S35.call.vect <- apply( expand.grid( S5.allele1.vect, 'KIR2DS35*null' ), 1, paste0, collapse='+' )
        }
      }else{
        S5.hetCall.vect <- tstrsplit(S5.call,' ',fixed=T)
        #S5.hetCall.list <- tstrsplit(S5.hetCall.vect,'+',fixed=T)
        
        #S35.call.vect <- apply( unique( expand.grid(S5.hetCall.list) ), 1, function(x) paste0(x,collapse='+'))
        #S35.call.vect <- paste(S35.call.vect)
        S35.call.vect <- unlist( S5.hetCall.vect )
      }
    }else if(S3.copy == 0 & S5.copy == 0) {
      S35.call.vect <- 'KIR2DS35*null+KIR2DS35*null'
    }else{
      S35.call.vect <- 'failed'
    }
    
    if( length(S35.call.vect) > 64 ){
      S35.call.str <- 'failed'
    }else{
      S35.call.str <- paste0(S35.call.vect, collapse=' ')
    }
  })
  
  iterCall.df$KIR2DS35 <- NA
  iterCall.df[names(S35.call.list),'KIR2DS35'] <- unlist( S35.call.list )
  
  return(iterCall.df)
}

pingFinalize.combineL1S1 <- function( iterCall.df, copyCall.df ){
  L1S1.call.list <- sapply( rownames( iterCall.df ), function(x){
    S1.call <- iterCall.df[x,'KIR3DS1']
    L1.call <- iterCall.df[x,'KIR3DL1']
    
    S1.copy <- as.integer( copyCall.df[x,'KIR3DS1'] )
    L1.copy <- as.integer( copyCall.df[x,'KIR3DL1'] )
    
    S1.call.bool <- !is.na(S1.call) & S1.call != 'failed'
    S1.call.unresolved <- grepl('unresolved',S1.call)
    
    L1.call.bool <- !is.na(L1.call) & L1.call != 'failed'
    L1.call.unresolved <- grepl('unresolved',L1.call)
    
    if( S1.call.unresolved ){
      S1.call <- paste0('KIR3DS1*',S1.call)
    }
    if( L1.call.unresolved ){
      L1.call <- paste0('KIR3DL1*',L1.call)
    }
    
    if( S1.call.bool ){
      S1.het.bool <- grepl('+',S1.call,fixed=T)
    }
    
    if( L1.call.bool ){
      L1.het.bool <- grepl('+',L1.call,fixed=T)
    }
    
    if( S1.call.bool & L1.call.bool ){
      
      if( !S1.het.bool ){
        S1.allele1.vect <- strsplit(S1.call,' ',fixed=T)[[1]]
        
        if( S1.copy > 1 ){
          S1.call.mat <- combinations(length(S1.allele1.vect), S1.copy, S1.allele1.vect,repeats.allowed = T)
          S1.call.vect <- apply( S1.call.mat, 1, paste0, collapse='+')
        }else{
          S1.call.vect <- S1.allele1.vect
        }
      }else{
        S1.hetCall.vect <- tstrsplit(S1.call,' ',fixed=T)
        S1.call.vect <- unlist( S1.hetCall.vect )
      }
      
      if( !L1.het.bool ){
        L1.allele1.vect <- strsplit(L1.call,' ',fixed=T)[[1]]
        if( L1.copy > 1 ){
          L1.call.mat <- combinations(length(L1.allele1.vect), L1.copy, L1.allele1.vect,repeats.allowed = T)
          L1.call.vect <- apply( L1.call.mat, 1, paste0, collapse='+')
        }else{
          L1.call.vect <- L1.allele1.vect
        }
      }else{
        L1.hetCall.vect <- tstrsplit(L1.call,' ',fixed=T)
        L1.call.vect <- unlist( L1.hetCall.vect )
      }
      
      L1S1.call.mat <- expand.grid(S1.call.vect,L1.call.vect)
      L1S1.call.vect <- apply( L1S1.call.mat, 1, paste0, collapse='+')
    }else if( S1.call.bool & !L1.call.bool ){
      if( !S1.het.bool ){
        S1.allele1.vect <- strsplit(S1.call,' ',fixed=T)[[1]]
        
        if( S1.copy > 1 ){
          S1.call.mat <- combinations(length(S1.allele1.vect), S1.copy, S1.allele1.vect,repeats.allowed = T)
          L1S1.call.vect <- apply( S1.call.mat, 1, paste0, collapse='+')
        }else{
          L1S1.call.vect <- apply( expand.grid( S1.allele1.vect, 'KIR3DL1S1*null' ), 1, paste0, collapse='+' )
        }
      }else{
        S1.hetCall.vect <- tstrsplit(S1.call,' ',fixed=T)
        L1S1.call.vect <- unlist( S1.hetCall.vect )
      }
    }else if(!S1.call.bool & L1.call.bool){
      if( !L1.het.bool ){
        L1.allele1.vect <- strsplit(L1.call,' ',fixed=T)[[1]]
        
        if( L1.copy > 1 ){
          L1.call.mat <- combinations(length(L1.allele1.vect), L1.copy, L1.allele1.vect,repeats.allowed = T)
          L1S1.call.vect <- apply( L1.call.mat, 1, paste0, collapse='+')
        }else{
          L1S1.call.vect <- apply( expand.grid( L1.allele1.vect, 'KIR3DL1S1*null' ), 1, paste0, collapse='+' )
        }
      }else{
        L1.hetCall.vect <- tstrsplit(L1.call,' ',fixed=T)
        L1S1.call.vect <- unlist( L1.hetCall.vect )
      }
    }else if(S1.copy == 0 & L1.copy == 0) {
      L1S1.call.vect <- 'KIR3DL1S1*null+KIR3DL1S1*null'
    }else{
      L1S1.call.vect <- 'failed'
    }
    
    if( length(L1S1.call.vect) > 64 ){
      L1S1.call.str <- 'failed'
    }else{
      L1S1.call.str <- paste0(L1S1.call.vect, collapse=' ')
    }
    
    return(L1S1.call.str)
    
  })
  
  iterCall.df$KIR3DL1S1 <- NA
  iterCall.df[names(L1S1.call.list),'KIR3DL1S1'] <- unlist( L1S1.call.list )
  
  return(iterCall.df)
}

pingFinalize.otherLoci <- function( iterCall.df, copyCall.df ){
  
  lociVect <- colnames(iterCall.df)
  lociVect <- lociVect[ !lociVect %in% c('KIR3DL1S1','KIR2DL23','KIR2DS35','KIR3DL1','KIR3DS1','KIR2DS3','KIR2DS5','KIR2DL2','KIR2DL3')]
  
  for( currentLocus in lociVect ){
    locus.call.list <- sapply( rownames( iterCall.df ), function(x){
      locus.call <- iterCall.df[x,currentLocus]
      locus.copy <- as.integer( copyCall.df[x,currentLocus] )
      
      locus.call.bool <- !is.na(locus.call) & locus.call != 'failed' & locus.call != ''
      locus.call.unresolved <- grepl('unresolved',locus.call) | (grepl('R', locus.call) & !grepl('KIR', locus.call))
      
      if( locus.call.unresolved ){
        locus.call <- paste0(currentLocus,'*unresolved')
        #locus.call <- paste0(currentLocus,'*',locus.call)
      }
      
      if( locus.call.bool ){
        locus.het.bool <- grepl('+',locus.call,fixed=T)
      }
      
      if( locus.call.bool ){
        if( !locus.het.bool ){
          locus.allele1.vect <- strsplit(locus.call,' ',fixed=T)[[1]]
          
          if( locus.copy > 1 ){
            locus.call.mat <- combinations(length(locus.allele1.vect), locus.copy, locus.allele1.vect,repeats.allowed = T)
            locus.call.vect <- apply( locus.call.mat, 1, paste0, collapse='+')
          }else{
            locus.call.vect <- apply( expand.grid( locus.allele1.vect, paste0(currentLocus,'*null') ), 1, paste0, collapse='+' )
            #locus.call.vect <- locus.allele1.vect
          }
        }else{
          locus.hetCall.vect <- tstrsplit(locus.call,' ',fixed=T)
          #locus.hetCall.list <- tstrsplit(locus.hetCall.vect,'+',fixed=T)
          #locus.call.vect <- apply( unique( expand.grid(locus.hetCall.list) ), 1, function(x) paste0(x,collapse='+'))
          
          locus.call.vect <- unlist( locus.hetCall.vect )
        }
      }else if(locus.copy == 0){
        locus.call.vect <- paste0(currentLocus,'*null+',currentLocus,'*null')
      }else{
        locus.call.vect <- 'failed'
      }
      
      if( length(locus.call.vect) > 64 ){
        locus.call.str <- 'failed'
      }else{
        locus.call.str <- paste0(locus.call.vect, collapse=' ')
      }
      
      return(locus.call.str)
    })
    iterCall.df[names(locus.call.list),currentLocus] <- unlist( locus.call.list )
  }
  
  return(iterCall.df)
}

og.iterCall.df <- iterCall.df
# ----- Formatting Results Genotypes (carry directly over to PING_run) -----
mod.locusVect <- c("KIR3DP1","KIR2DP1","KIR2DS2","KIR2DL4","KIR3DL3","KIR3DL2","KIR2DS4",
                   "KIR2DL1","KIR2DS1","KIR2DL5","KIR2DL23","KIR2DS35","KIR3DL1S1")
#iterCall.df <- pingFinalize.combineL23( iterCall.df )
#iterCall.df <- pingFinalize.combineS35( iterCall.df, copyCall.df )
#iterCall.df <- pingFinalize.combineL1S1( iterCall.df, copyCall.df )
#iterCall.df <- pingFinalize.otherLoci( iterCall.df, copyCall.df )

iterCall.df <- iterCall.df[,mod.locusVect]


# ----- Formating validation genotypes -----
'Valid dataset had all 3DL3 recombinants removed (whole genotype removed). Any other recombinant gene was marked as failure'
pingValid.combineS35 <- function( iterCall.df, copyCall.df ){
  S35.call.list <- sapply( rownames( iterCall.df ), function(x){
    #cat('\n',x)
    S3.call <- iterCall.df[x,'KIR2DS3']
    S5.call <- iterCall.df[x,'KIR2DS5']
    
    S3.copy <- copyCall.df[x,'KIR2DS3']
    S5.copy <- copyCall.df[x,'KIR2DS5']
    
    if(S3.copy == 'R' | S5.copy == 'R'){
      return('failed')
    }else{
      S3.copy <- as.integer(S3.copy)
      S5.copy <- as.integer(S5.copy)
    }
    
    S3.call.bool <- !is.na(S3.call) & S3.call != 'failed' & S3.call != ''
    S3.call.unresolved <- grepl('new',S3.call)
    
    S5.call.bool <- !is.na(S5.call) & S5.call != 'failed' & S5.call != ''
    S5.call.unresolved <- grepl('new',S5.call)
    
    if( S3.call.unresolved ){
      S3.call <- paste0('KIR2DS3*',S3.call)
    }
    if( S5.call.unresolved ){
      S5.call <- paste0('KIR2DS5*',S5.call)
    }
    
    if( S3.call.bool ){
      S3.het.bool <- grepl('+',S3.call,fixed=T)
    }
    
    if( S5.call.bool ){
      S5.het.bool <- grepl('+',S5.call,fixed=T)
    }
    
    if( S3.call.bool & S5.call.bool ){
      
      if( !S3.het.bool ){
        S3.allele1.vect <- strsplit(S3.call,' ',fixed=T)[[1]]
        
        if( S3.copy > 1 ){
          S3.call.mat <- combinations(length(S3.allele1.vect), S3.copy, S3.allele1.vect,repeats.allowed = T)
          S3.call.vect <- apply( S3.call.mat, 1, paste0, collapse='+')
        }else{
          S3.call.vect <- S3.allele1.vect
        }
      }else{
        S3.hetCall.vect <- tstrsplit(S3.call,' ',fixed=T)
        #S3.hetCall.list <- tstrsplit(S3.hetCall.vect,'+',fixed=T)
        #S3.allele1.vect <- S3.hetCall.list[[1]]
        #S3.allele2.vect <- S3.hetCall.list[[2]]
        
        #S3.call.vect <- apply( unique( expand.grid(S3.hetCall.list) ), 1, function(x) paste0(x,collapse='+'))
        S3.call.vect <- unlist( S3.hetCall.vect )
      }
      
      if( !S5.het.bool ){
        S5.allele1.vect <- strsplit(S5.call,' ',fixed=T)[[1]]
        if( S5.copy > 1 ){
          S5.call.mat <- combinations(length(S5.allele1.vect), S5.copy, S5.allele1.vect,repeats.allowed = T)
          S5.call.vect <- apply( S5.call.mat, 1, paste0, collapse='+')
        }else{
          S5.call.vect <- S5.allele1.vect
        }
      }else{
        S5.hetCall.vect <- tstrsplit(S5.call,' ',fixed=T)
        #S5.hetCall.list <- tstrsplit(S5.hetCall.vect,'+',fixed=T)
        #S5.allele1.vect <- S5.hetCall.list[[1]]
        #S5.allele2.vect <- S5.hetCall.list[[2]]
        
        #S5.call.vect <- apply( unique( expand.grid(S5.hetCall.list) ), 1, function(x) paste0(x,collapse='+'))
        S5.call.vect <- unlist( S5.hetCall.vect )
      }
      
      S35.call.mat <- expand.grid(S3.call.vect,S5.call.vect)
      S35.call.vect <- apply( S35.call.mat, 1, paste0, collapse='+')
    }else if( S3.call.bool & !S5.call.bool ){
      if( !S3.het.bool ){
        S3.allele1.vect <- strsplit(S3.call,' ',fixed=T)[[1]]
        
        if( S3.copy > 1 ){
          S3.call.mat <- combinations(length(S3.allele1.vect), S3.copy, S3.allele1.vect,repeats.allowed = T)
          S35.call.vect <- apply( S3.call.mat, 1, paste0, collapse='+')
        }else{
          S35.call.vect <- apply( expand.grid( S3.allele1.vect, 'KIR2DS35*null' ), 1, paste0, collapse='+' )
        }
      }else{
        S3.hetCall.vect <- tstrsplit(S3.call,' ',fixed=T)
        #S3.hetCall.list <- tstrsplit(S3.hetCall.vect,'+',fixed=T)
        
        #S35.call.vect <- apply( unique( expand.grid(S3.hetCall.list) ), 1, function(x) paste0(x,collapse='+') )
        S35.call.vect <- unlist( S3.hetCall.vect )
      }
    }else if(!S3.call.bool & S5.call.bool){
      if( !S5.het.bool ){
        S5.allele1.vect <- strsplit(S5.call,' ',fixed=T)[[1]]
        
        if( S5.copy > 1 ){
          S5.call.mat <- combinations(length(S5.allele1.vect), S5.copy, S5.allele1.vect,repeats.allowed = T)
          S35.call.vect <- apply( S5.call.mat, 1, paste0, collapse='+')
        }else{
          S35.call.vect <- apply( expand.grid( S5.allele1.vect, 'KIR2DS35*null' ), 1, paste0, collapse='+' )
        }
      }else{
        S5.hetCall.vect <- tstrsplit(S5.call,' ',fixed=T)
        #S5.hetCall.list <- tstrsplit(S5.hetCall.vect,'+',fixed=T)
        
        #S35.call.vect <- apply( unique( expand.grid(S5.hetCall.list) ), 1, function(x) paste0(x,collapse='+'))
        #S35.call.vect <- paste(S35.call.vect)
        S35.call.vect <- unlist( S5.hetCall.vect )
      }
    }else if(S3.copy == 0 & S5.copy == 0) {
      S35.call.vect <- 'KIR2DS35*null+KIR2DS35*null'
    }else{
      S35.call.vect <- 'failed'
    }
    
    if( length(S35.call.vect) > 64 ){
      S35.call.str <- 'failed'
    }else{
      S35.call.str <- paste0(S35.call.vect, collapse=' ')
    }
  })
  
  iterCall.df$KIR2DS35 <- NA
  iterCall.df[names(S35.call.list),'KIR2DS35'] <- unlist( S35.call.list )
  
  return(iterCall.df)
}
pingValid.combineL1S1 <- function( iterCall.df, copyCall.df ){
  L1S1.call.list <- sapply( rownames( iterCall.df ), function(x){
    S1.call <- iterCall.df[x,'KIR3DS1']
    L1.call <- iterCall.df[x,'KIR3DL1']
    
    S1.copy <- copyCall.df[x,'KIR3DS1']
    L1.copy <- copyCall.df[x,'KIR3DL1']
    
    if(S1.copy == 'R' | L1.copy == 'R'){
      return('failed')
    }else{
      S1.copy <- as.integer(S1.copy)
      L1.copy <- as.integer(L1.copy)
    }
    
    S1.call.bool <- !is.na(S1.call) & S1.call != 'failed' & S1.call != ''
    S1.call.unresolved <- grepl('new',S1.call)
    
    L1.call.bool <- !is.na(L1.call) & L1.call != 'failed' & L1.call != ''
    L1.call.unresolved <- grepl('new',L1.call)
    
    if( S1.call.unresolved ){
      S1.call <- paste0('KIR3DS1*',S1.call)
    }
    if( L1.call.unresolved ){
      L1.call <- paste0('KIR3DL1*',L1.call)
    }
    
    if( S1.call.bool ){
      S1.het.bool <- grepl('+',S1.call,fixed=T)
    }
    
    if( L1.call.bool ){
      L1.het.bool <- grepl('+',L1.call,fixed=T)
    }
    
    if( S1.call.bool & L1.call.bool ){
      
      if( !S1.het.bool ){
        S1.allele1.vect <- strsplit(S1.call,' ',fixed=T)[[1]]
        
        if( S1.copy > 1 ){
          S1.call.mat <- combinations(length(S1.allele1.vect), S1.copy, S1.allele1.vect,repeats.allowed = T)
          S1.call.vect <- apply( S1.call.mat, 1, paste0, collapse='+')
        }else{
          S1.call.vect <- S1.allele1.vect
        }
      }else{
        S1.hetCall.vect <- tstrsplit(S1.call,' ',fixed=T)
        S1.call.vect <- unlist( S1.hetCall.vect )
      }
      
      if( !L1.het.bool ){
        L1.allele1.vect <- strsplit(L1.call,' ',fixed=T)[[1]]
        if( L1.copy > 1 ){
          L1.call.mat <- combinations(length(L1.allele1.vect), L1.copy, L1.allele1.vect,repeats.allowed = T)
          L1.call.vect <- apply( L1.call.mat, 1, paste0, collapse='+')
        }else{
          L1.call.vect <- L1.allele1.vect
        }
      }else{
        L1.hetCall.vect <- tstrsplit(L1.call,' ',fixed=T)
        L1.call.vect <- unlist( L1.hetCall.vect )
      }
      
      L1S1.call.mat <- expand.grid(S1.call.vect,L1.call.vect)
      L1S1.call.vect <- apply( L1S1.call.mat, 1, paste0, collapse='+')
    }else if( S1.call.bool & !L1.call.bool ){
      if( !S1.het.bool ){
        S1.allele1.vect <- strsplit(S1.call,' ',fixed=T)[[1]]
        
        if( S1.copy > 1 ){
          S1.call.mat <- combinations(length(S1.allele1.vect), S1.copy, S1.allele1.vect,repeats.allowed = T)
          L1S1.call.vect <- apply( S1.call.mat, 1, paste0, collapse='+')
        }else{
          L1S1.call.vect <- apply( expand.grid( S1.allele1.vect, 'KIR3DL1S1*null' ), 1, paste0, collapse='+' )
        }
      }else{
        S1.hetCall.vect <- tstrsplit(S1.call,' ',fixed=T)
        L1S1.call.vect <- unlist( S1.hetCall.vect )
      }
    }else if(!S1.call.bool & L1.call.bool){
      if( !L1.het.bool ){
        L1.allele1.vect <- strsplit(L1.call,' ',fixed=T)[[1]]
        
        if( L1.copy > 1 ){
          L1.call.mat <- combinations(length(L1.allele1.vect), L1.copy, L1.allele1.vect,repeats.allowed = T)
          L1S1.call.vect <- apply( L1.call.mat, 1, paste0, collapse='+')
        }else{
          L1S1.call.vect <- apply( expand.grid( L1.allele1.vect, 'KIR3DL1S1*null' ), 1, paste0, collapse='+' )
        }
      }else{
        L1.hetCall.vect <- tstrsplit(L1.call,' ',fixed=T)
        L1S1.call.vect <- unlist( L1.hetCall.vect )
      }
    }else if(S1.copy == 0 & L1.copy == 0) {
      L1S1.call.vect <- 'KIR3DL1S1*null+KIR3DL1S1*null'
    }else{
      L1S1.call.vect <- 'failed'
    }
    
    if( length(L1S1.call.vect) > 64 ){
      L1S1.call.str <- 'failed'
    }else{
      L1S1.call.str <- paste0(L1S1.call.vect, collapse=' ')
    }
    
    return(L1S1.call.str)
    
  })
  
  iterCall.df$KIR3DL1S1 <- NA
  iterCall.df[names(L1S1.call.list),'KIR3DL1S1'] <- unlist( L1S1.call.list )
  
  return(iterCall.df)
}
pingValid.otherLoci <- function( iterCall.df, copyCall.df ){
  
  lociVect <- colnames(iterCall.df)
  lociVect <- lociVect[ !lociVect %in% c('KIR3DL1S1','KIR2DS35','KIR3DL1S1','KIR2DS3','KIR2DS5')]
  
  for( currentLocus in lociVect ){
    #cat("\n",currentLocus)
    locus.call.list <- sapply( rownames( iterCall.df ), function(x){
      #cat('\n\t',x)
      locus.call <- iterCall.df[x,currentLocus]
      
      if(currentLocus == 'KIR2DL23'){
        L2.copy <- copyCall.df[x,'KIR2DL2']
        L3.copy <- copyCall.df[x,'KIR2DL3']
        
        if( L2.copy == 'R' | L3.copy == 'R' ){
          return('failed')
        }else{
          locus.copy <- as.integer(L2.copy) + as.integer(L3.copy)
        }
      }else{
        locus.copy <- copyCall.df[x,currentLocus]
        
        if( locus.copy == 'R' ){
          return('failed')
        }else{
          locus.copy <- as.integer(locus.copy)
        }
      }
      
      locus.call.bool <- !is.na(locus.call) & locus.call != 'failed' & locus.call != ''
      locus.call.unresolved <- grepl('new',locus.call)
      
      if( locus.call.unresolved ){
        locus.call <- paste0(currentLocus,'*',locus.call)
      }
      
      if( locus.call.bool ){
        locus.het.bool <- grepl('+',locus.call,fixed=T)
      }
      
      if( locus.call.bool ){
        if( !locus.het.bool ){
          locus.allele1.vect <- strsplit(locus.call,' ',fixed=T)[[1]]
          
          if( locus.copy > 1 ){
            locus.call.mat <- combinations(length(locus.allele1.vect), locus.copy, locus.allele1.vect,repeats.allowed = T)
            locus.call.vect <- apply( locus.call.mat, 1, paste0, collapse='+')
          }else{
            locus.call.vect <- apply( expand.grid( locus.allele1.vect, paste0(currentLocus,'*null') ), 1, paste0, collapse='+' )
            #locus.call.vect <- locus.allele1.vect
          }
        }else{
          locus.hetCall.vect <- tstrsplit(locus.call,' ',fixed=T)
          #locus.hetCall.list <- tstrsplit(locus.hetCall.vect,'+',fixed=T)
          #locus.call.vect <- apply( unique( expand.grid(locus.hetCall.list) ), 1, function(x) paste0(x,collapse='+'))
          
          locus.call.vect <- unlist( locus.hetCall.vect )
        }
      }else if(locus.copy == 0){
        locus.call.vect <- paste0(currentLocus,'*null+',currentLocus,'*null')
      }else{
        locus.call.vect <- 'failed'
      }
      
      if( length(locus.call.vect) > 64 ){
        locus.call.str <- 'failed'
      }else{
        locus.call.str <- paste0(locus.call.vect, collapse=' ')
      }
      
      return(locus.call.str)
    })
    iterCall.df[names(locus.call.list),currentLocus] <- unlist( locus.call.list )
  }
  
  return(iterCall.df)
}

og.validGenoCall.df <- validGenoCall.df
validGenoCall.df <- pingValid.combineS35( validGenoCall.df, validCopyCall.df )
validGenoCall.df <- pingValid.combineL1S1( validGenoCall.df, validCopyCall.df )
validGenoCall.df <- pingValid.otherLoci( validGenoCall.df, validCopyCall.df )

valid.mod.locusVect <- c("KIR2DP1","KIR2DS2","KIR2DL4","KIR3DL3","KIR3DL2","KIR2DS4",
                         "KIR2DL1","KIR2DS1","KIR2DL5","KIR2DL23","KIR2DS35","KIR3DL1S1")
validGenoCall.df <- validGenoCall.df[,valid.mod.locusVect]

# ----- Formatting for validation comparison -----
sharedSampleID.vect <- intersect(rownames(iterCall.df), rownames(validGenoCall.df))
iterCall.dt <- synSeq.formatResultGeno( iterCall.df[sharedSampleID.vect,], copyCall.df[sharedSampleID.vect,] )
validGenoCall.dt <- synSeq.formatResultGeno( validGenoCall.df[sharedSampleID.vect,], validCopyCall.df[sharedSampleID.vect,] )

# Copy concordance
pingValid.copyConcordance <- function( pingCopy.df, validCopy.df ){
  
  locusVect <- intersect( colnames(pingCopy.df), colnames(validCopy.df) )
  sampleVect <- intersect( rownames(pingCopy.df), rownames(validCopy.df) )
  
  out.copyConc.list <- list()
  
  for( currentLocus in locusVect ){
    sampleVect.noR <- sampleVect[ validCopy.df[sampleVect,currentLocus] != 'R' ]
    
    concVect <- as.integer( validCopy.df[sampleVect.noR,currentLocus] ) == as.integer( pingCopy.df[sampleVect.noR,currentLocus ])
    
    matchN <- sum(concVect)
    totalN <- length(concVect)
    
    if( any(!concVect) ){
      mismatchID.vect <- sampleVect.noR[ which( !concVect ) ]
    }else{
      mismatchID.vect <- NA
    }
    
    cat(paste0('\n',currentLocus,'\t',matchN,'\t',totalN,'\t',round(matchN/totalN,3)))
    out.copyConc.list[[currentLocus]] <- list('matchN'=matchN,'totalN'=totalN,mismatchIDVect=mismatchID.vect)
  }
  
  return(out.copyConc.list)
}

copyCon.list <- pingValid.copyConcordance( copyCall.df, validCopyCall.df )


# ----- Comparing genotypes to validation set -----
pingValid.genoConcordance <- function( genoCall.dt, validCall.dt, kirRes = 3, score.unresolved=T, score.copy4 = F){
  
  #genoKey.dt <- as.data.table( t( sapply( genoCopyKey.list, function(x) x$genoList ) ) )
  #genoKey.dt$sampleID <- names(genoCopyKey.list)
  
  sampleIDVect <- intersect( validCall.dt$sampleID, genoCall.dt$sampleID )
  
  orderedIDVect <- sampleIDVect[order(sampleIDVect)]
  locusVect <- intersect( colnames(validCall.dt), colnames(genoCall.dt) )
  locusVect <- setdiff( locusVect, 'sampleID' )
  
  out.list <- list()
  for(currentLocus in locusVect){
    
    locusConcordanceList <- lapply(orderedIDVect, function(currentID){
      keyList <- validCall.dt[sampleID == currentID,..currentLocus][[1]]
      callList <- genoCall.dt[sampleID == currentID,..currentLocus][[1]]
      
      if( 'failed' %in% keyList[[1]] ){
        return(list('matchN'=0,'totalN'=0,'mismatchID'=NA))
      }else if( any(grepl('new', keyList[[1]])) ){
        return(list('matchN'=0,'totalN'=0,'mismatchID'=NA))
      }else if( any( grepl('unresolved', keyList[[1]]) ) ){
        return(list('matchN'=0,'totalN'=0,'mismatchID'=NA))
      }
      
      if( !score.unresolved ){
        if( any( grepl('unresolved', callList[[1]]) ) ){
          return(list('matchN'=0,'totalN'=0,'mismatchID'=NA))
        }
      }
      
      if( !score.copy4 ){
        if( length(keyList[[1]]) > 3 ){
          return(list('matchN'=0,'totalN'=0,'mismatchID'=NA))
        }
      }
      
      concordanceList <- apply(callList[[1]], 1, function(alleleVect){
        
        alleleVect <- sapply(alleleVect, kir.allele_resolution, kirRes)
        keyAlleleVect <- sapply(keyList[[1]], kir.allele_resolution, kirRes)
        
        fw.matchN <- sum( alleleVect %in% keyAlleleVect )
        rv.matchN <- sum( keyAlleleVect %in% alleleVect )
        
        matchN <- min(fw.matchN, rv.matchN)
        
        totalN <- ncol(keyList[[1]])
        
        #matchN <- min( matchN, totalN )
        
        #callCopy <- ncol(callList[[1]])
        #copyDiff <- abs( totalN - callCopy )
        
        #matchN <- matchN - copyDiff
        #matchN <- max(0, matchN)
        
        #if(currentLocus == 'KIR2DL23'){
        #  copyNumber <- 2
        #}else if(currentLocus == 'KIR3DL1' | currentLocus == 'KIR3DS1'){
        #  copyNumber <- as.integer(genoCopyKey.list[[currentID]]$copyList[['KIR3DL1']])
        #  copyNumber <- copyNumber + as.integer(genoCopyKey.list[[currentID]]$copyList[['KIR3DS1']])
        #}else{
        #  copyNumber <- as.integer(genoCopyKey.list[[currentID]]$copyList[[currentLocus]])
        #}
        
        #if(copyNumber > 2){
        #  return(list('matchN'=0,'totalN'=0))
        #}else{
        if(matchN > totalN){
          cat('\n',currentID,currentLocus)
        }
        return(list('matchN'=matchN,'totalN'=totalN))
        #}
      })
      
      matchN <- max( sapply( concordanceList, function(x){x$matchN} ) )
      totalN <- max( sapply( concordanceList, function(x){x$totalN} ) )
      
      if( matchN != totalN ){
        mismatchID <- currentID
      }else{
        mismatchID <- NA
      }
      
      return(list('matchN'=matchN,'totalN'=totalN,'mismatchID'=mismatchID))
    })
    
    matchN <- sum( sapply(locusConcordanceList, function(x) x$matchN ) )
    totalN <- sum( sapply(locusConcordanceList, function(x) x$totalN ) )
    mismatchIDVect <- unique( sapply( locusConcordanceList, function(x) x$mismatchID ) )
    mismatchIDVect <- mismatchIDVect[ !is.na(mismatchIDVect) ]
    
    percMatch <- round(matchN/totalN,3)
    
    matchMessageStr <- paste(currentLocus, percMatch, totalN, sep='\t')
    cat('\n',matchMessageStr)
    
    out.list[[currentLocus]] <- list('N'=totalN,'match'=matchN,'percentage'=percMatch,'mismatchIDs'=mismatchIDVect)
  }
  
  return(out.list)
}

pingValid.genoConcordanceMakeup <- function( genoCall.dt, validCall.dt, kirRes = 3, score.copy4 = F){
  
  #genoKey.dt <- as.data.table( t( sapply( genoCopyKey.list, function(x) x$genoList ) ) )
  #genoKey.dt$sampleID <- names(genoCopyKey.list)
  
  sampleIDVect <- intersect( validCall.dt$sampleID, genoCall.dt$sampleID )
  
  orderedIDVect <- sampleIDVect[order(sampleIDVect)]
  locusVect <- intersect( colnames(validCall.dt), colnames(genoCall.dt) )
  locusVect <- setdiff( locusVect, 'sampleID' )
  
  out.list <- list()
  for(currentLocus in locusVect){
    
    locusConcordanceList <- lapply(orderedIDVect, function(currentID){
      keyList <- validCall.dt[sampleID == currentID,..currentLocus][[1]]
      callList <- genoCall.dt[sampleID == currentID,..currentLocus][[1]]
      
      if( 'failed' %in% keyList[[1]] ){
        return(list('matchN'=0,'totalN'=0,
                    'unresolvedN'=0,
                    'allele.mismatchN'=0,
                    'mismatchID'=NA))
      }else if( any(grepl('new', keyList[[1]])) ){
        return(list('matchN'=0,'totalN'=0,
                    'unresolvedN'=0,
                    'allele.mismatchN'=0,
                    'mismatchID'=NA))
      }else if( any( grepl('unresolved', keyList[[1]]) ) ){
        return(list('matchN'=0,'totalN'=0,
                    'unresolvedN'=0,
                    'allele.mismatchN'=0,
                    'mismatchID'=NA))
      }
      
      
      if( !score.copy4 ){
        if( length(keyList[[1]]) > 3 ){
          return(list('matchN'=0,'totalN'=0,
                      'unresolvedN'=0,
                      'allele.mismatchN'=0,
                      'mismatchID'=NA))
        }
      }
      
      concordanceList <- apply(callList[[1]], 1, function(alleleVect){
        
        alleleVect <- sapply(alleleVect, kir.allele_resolution, kirRes)
        keyAlleleVect <- sapply(keyList[[1]], kir.allele_resolution, kirRes)
        
        origin.alleleVect <- alleleVect
        
        matchVect <- c()
        for( compAllele in keyAlleleVect ){
          if(compAllele %in% alleleVect){
            compAllele.matchInd <- which( alleleVect == compAllele )[1]
            matchVect <- c(matchVect,compAllele)
            alleleVect <- alleleVect[-compAllele.matchInd]
          }
          
          if(length(alleleVect) == 0){
            break
          }
        }
        
        matchN <- length(matchVect)
        unresolvedN <- sum(grepl('*unr', origin.alleleVect,fixed=T))
        unresolvedN <- min( unresolvedN, length(keyAlleleVect) )
        allele.mismatchN <- length(keyAlleleVect) - length(matchVect) - unresolvedN
        
        if( allele.mismatchN < 0 ){
          cat('\nnegative mismatch',currentID,currentLocus)
        }
        

        #fw.matchN <- sum( alleleVect %in% keyAlleleVect )
        #rv.matchN <- sum( keyAlleleVect %in% alleleVect )
        
        #matchN <- min(fw.matchN, rv.matchN)
        
        totalN <- ncol(keyList[[1]])
        
        #fw.null.match <- sum( grepl('*null',alleleVect) )
        
        #fw.matchVect <- keyAlleleVect[ keyAlleleVect %in% alleleVect ]
        #rv.matchVect <- alleleVect[ alleleVect %in% keyAlleleVect ]
        
        #null.match <- sum(grepl('*null',matchVect))
        #allele.match <- length(matchVect) - null.match
        
        
        
        
        #matchN <- min( matchN, totalN )
        
        #callCopy <- ncol(callList[[1]])
        #copyDiff <- abs( totalN - callCopy )
        
        #matchN <- matchN - copyDiff
        #matchN <- max(0, matchN)
        
        #if(currentLocus == 'KIR2DL23'){
        #  copyNumber <- 2
        #}else if(currentLocus == 'KIR3DL1' | currentLocus == 'KIR3DS1'){
        #  copyNumber <- as.integer(genoCopyKey.list[[currentID]]$copyList[['KIR3DL1']])
        #  copyNumber <- copyNumber + as.integer(genoCopyKey.list[[currentID]]$copyList[['KIR3DS1']])
        #}else{
        #  copyNumber <- as.integer(genoCopyKey.list[[currentID]]$copyList[[currentLocus]])
        #}
        
        #if(copyNumber > 2){
        #  return(list('matchN'=0,'totalN'=0))
        #}else{
        if(matchN > totalN){
          cat('\noverN',currentID,currentLocus)
        }
        
        return(list('matchN'=matchN,'totalN'=totalN,
                    'unresolvedN'=unresolvedN,
                    'allele.mismatchN'=allele.mismatchN))
        #return(list('matchN'=matchN,'totalN'=totalN))
        #}
      })
      
      matchN <- max( sapply( concordanceList, function(x){x$matchN} ) )
      matchInd <- which( sapply( concordanceList, function(x) x$matchN == matchN ) )[1]
      
      matchResult <- concordanceList[[matchInd]]
      totalN <- matchResult$totalN
      unresolvedN <- matchResult$unresolvedN
      allele.mismatchN <- matchResult$allele.mismatchN
      #totalN <- max( sapply( concordanceList, function(x){x$totalN} ) )
      
      if( matchN != totalN ){
        mismatchID <- currentID
      }else{
        mismatchID <- NA
      }
      
      return(list('matchN'=matchN,'totalN'=totalN,
                  'unresolvedN'=unresolvedN,
                  'allele.mismatchN'=allele.mismatchN,
                  'mismatchID'=mismatchID))
      
      #return(list('matchN'=matchN,'totalN'=totalN,'mismatchID'=mismatchID))
    })
    
    matchN <- sum( sapply(locusConcordanceList, function(x) x$matchN ) )
    totalN <- sum( sapply(locusConcordanceList, function(x) x$totalN ) )
    unresolvedN <- sum( sapply(locusConcordanceList, function(x) x$unresolvedN ) )
    allele.mismatchN <- sum( sapply(locusConcordanceList, function(x) x$allele.mismatchN ) )

    mismatchIDVect <- unique( sapply( locusConcordanceList, function(x) x$mismatchID ) )
    mismatchIDVect <- mismatchIDVect[ !is.na(mismatchIDVect) ]
    
    percMatch <- round(matchN/totalN,3)
    
    matchMessageStr <- paste(currentLocus, percMatch, totalN, sep='\t')
    cat('\n',matchMessageStr)
    
    out.list[[currentLocus]] <- list('N'=totalN,'match'=matchN,'unresolved'=unresolvedN,
                                     'alleleMismatch'=allele.mismatchN,'mismatchIDs'=mismatchIDVect)
    #out.list[[currentLocus]] <- list('N'=totalN,'match'=matchN,'percentage'=percMatch,'mismatchIDs'=mismatchIDVect)
  }
  
  return(out.list)
}


#genoCon.list <- pingValid.genoConcordance( iterCall.dt, validGenoCall.dt, kirRes = 3 , score.unresolved = T, score.copy4=F)
genoConMakup.list <- pingValid.genoConcordanceMakeup( iterCall.dt, validGenoCall.dt, kirRes = 3, score.copy4 = F )

for( currentLocus in names( genoConMakup.list ) ){
  matchN <- genoConMakup.list[[currentLocus]]$match
  totalN <- genoConMakup.list[[currentLocus]]$N
  unresolvedN <- genoConMakup.list[[currentLocus]]$unresolved
  percMatch <- round(matchN / (totalN - unresolvedN),digits = 3)
  cat(paste0('\n',currentLocus,'\t',percMatch,'\t',(totalN - unresolvedN)))
}

# Graphing
gene <- sort( rep( names( genoConMakup.list ), 3) )
condition <- rep( c('match','unresolved','alleleMismatch'), length(names(genoConMakup.list)) )
value <- c()
for( currentLocus in unique(gene) ){
  for( currentCondition in unique(condition) ){
    value <- c(value, genoConMakup.list[[currentLocus]][[currentCondition]])
  }
}
plotData.380 <- data.frame(gene, condition, value)
plotData.380$condition <- factor( plotData.380$condition, levels=c('alleleMismatch','unresolved','match'))

# ----- European cohort concordance table -----
for( currentLocus in unique(plotData.380$gene) ){
  sub.df <- plotData.380[plotData.380$gene == currentLocus,c('condition','value'),drop=F]
  totalN <- sum(sub.df$value)
  valueVect <- sub.df$value / totalN
  cat(paste0('\n',currentLocus,'\t',paste0(round(valueVect,3),collapse='\t'),'\t',totalN))
}

g1 <- ggplot(plotData.380, aes(fill=condition, y=value, x=gene)) +
  geom_bar( position='fill',stat='identity') +
  theme_classic() + 
  scale_fill_manual(values=c('#b0b0b0','#595959','#000000')) +
  theme(axis.text.x = element_text(angle = 90),text = element_text(size=22))
g1
ggsave(file.path(resultsDirectory,"conc.380.percent.png"), width = 30, height = 15, units = "cm")

g2 <- ggplot(plotData.380, aes(fill=condition, y=value, x=gene)) +
  geom_bar( position='stack',stat='identity') +
  theme_classic() + 
  scale_fill_manual(values=c('#b0b0b0','#595959','#000000')) +
  theme(axis.text.x = element_text(angle = 90),text = element_text(size=22))
g2
ggsave(file.path(resultsDirectory,"conc.380.value.png"), width = 30, height = 15, units = "cm")

# ----- Discordance analysis -----

#View( iterCall.dt[sampleID %in% genoCon.list$KIR2DP1$mismatchIDs,] )
#View( validGenoCall.dt[sampleID %in% genoCon.list$KIR2DP1$mismatchIDs,] )

#View( iterCall.dt[sampleID %in% genoCon.list$KIR2DS2$mismatchIDs,] )


# ----- SYN SEQ CONCORDANCE -----
finalGeno.syn.df <- read.csv(finalGeno.syn.path,
                        stringsAsFactors=F,check.names=F,row.names=1)
rownames(finalGeno.syn.df) <- tstrsplit( rownames(finalGeno.syn.df), '.paired' )[[1]]

finalCopy.syn.df <- read.csv(copy.syn.path,
                        stringsAsFactors=F,check.names=F,row.names=1)
rownames(finalCopy.syn.df) <- tstrsplit( rownames(finalCopy.syn.df), '.paired' )[[1]]


validGeno.syn.df <- pingFinalize.combineL23( validGeno.syn.df )
validGeno.syn.df <- pingFinalize.combineS35( validGeno.syn.df, validCopy.syn.df )
validGeno.syn.df <- pingFinalize.combineL1S1( validGeno.syn.df, validCopy.syn.df )
validGeno.syn.df <- pingFinalize.otherLoci( validGeno.syn.df, validCopy.syn.df )
validGeno.syn.df <- validGeno.syn.df[,mod.locusVect]


sharedSampleID.syn.vect <- intersect(rownames(finalGeno.syn.df), rownames(validGeno.syn.df))

finalGeno.syn.dt <- synSeq.formatResultGeno( finalGeno.syn.df[sharedSampleID.syn.vect,], finalCopy.syn.df[sharedSampleID.syn.vect,] )
validGeno.syn.dt <- synSeq.formatResultGeno( validGeno.syn.df[sharedSampleID.syn.vect,], validCopy.syn.df[sharedSampleID.syn.vect,] )

## Concordance comparison
copyCon.syn.list <- pingValid.copyConcordance( finalCopy.syn.df, validCopy.syn.df )
#genoCon.syn.list <- pingValid.genoConcordance( finalGeno.syn.dt, validGeno.syn.dt, kirRes = 5, score.unresolved = F, score.copy4 = F )
genoConMakup.syn.list <- pingValid.genoConcordanceMakeup( finalGeno.syn.dt, validGeno.syn.dt, kirRes = 5, score.copy4 = F )


for( currentLocus in names( genoConMakup.syn.list ) ){
  matchN <- genoConMakup.syn.list[[currentLocus]]$match
  totalN <- genoConMakup.syn.list[[currentLocus]]$N
  unresolvedN <- genoConMakup.syn.list[[currentLocus]]$unresolved
  percMatch <- round(matchN / (totalN - unresolvedN),digits = 3)
  cat(paste0('\n',currentLocus,'\t',percMatch,'\t',(totalN - unresolvedN)))
}

# Graphing
gene <- sort( rep( names( genoConMakup.syn.list ), 3) )
condition <- rep( c('match','unresolved','alleleMismatch'), length(names(genoConMakup.syn.list)) )
value <- c()
for( currentLocus in unique(gene) ){
  for( currentCondition in unique(condition) ){
    value <- c(value, genoConMakup.syn.list[[currentLocus]][[currentCondition]])
  }
}
plotData.syn <- data.frame(gene, condition, value)
plotData.syn$condition <- factor( plotData.syn$condition, levels=c('alleleMismatch','unresolved','match'))

# ----- Synthetic dataset concordance table -----
for( currentLocus in unique(plotData.syn$gene) ){
  sub.df <- plotData.syn[plotData.syn$gene == currentLocus,c('condition','value'),drop=F]
  totalN <- sum(sub.df$value)
  valueVect <- sub.df$value / totalN
  cat(paste0('\n',currentLocus,'\t',paste0(round(valueVect,3),collapse='\t'),'\t',totalN))
}

g1 <- ggplot(plotData.syn, aes(fill=condition, y=value, x=gene)) +
  geom_bar( position='fill',stat='identity') +
  theme_classic() + 
  scale_fill_manual(values=c('#b0b0b0','#595959','#000000')) +
  theme(axis.text.x = element_text(angle = 90),text = element_text(size=22))
g1
ggsave(file.path(resultsDirectory,"conc.syn.percent.png"), width = 30, height = 15, units = "cm")


g2 <- ggplot(plotData.syn, aes(fill=condition, y=value, x=gene)) +
  geom_bar( position='stack',stat='identity') +
  theme_classic() + 
  scale_fill_manual(values=c('#b0b0b0','#595959','#000000')) +
  theme(axis.text.x = element_text(angle = 90),text = element_text(size=22))
g2
ggsave(file.path(resultsDirectory,"conc.syn.value.png"), width = 30, height = 15, units = "cm")




# ----- Read in and format AFR validation -----
mod.afr.locusVect <- c("KIR2DP1","KIR2DS2","KIR2DL4","KIR3DL3","KIR3DL2","KIR2DS4",
                       "KIR2DL1","KIR2DS1","KIR2DL5","KIR2DL23","KIR2DS35","KIR3DL1S1")
finalGeno.afr.df <- read.csv(finalGeno.afr.path,
                             stringsAsFactors=F,check.names=F,row.names=1)
finalCopy.afr.df <- read.csv(finalCopy.afr.path,
                             stringsAsFactors=F,check.names=F,row.names=1)
validGeno.afr.df <- read.csv(validGeno.afr.path,
                             stringsAsFactors=F,check.names=F,row.names=1)
validCopy.afr.df <- read.csv(validCopy.afr.path,
                             stringsAsFactors=F,check.names=F,row.names=1)

validGeno.afr.df <- pingFinalize.otherLoci( validGeno.afr.df, validCopy.afr.df )
validGeno.afr.df <- validGeno.afr.df[,mod.afr.locusVect]

sharedSampleID.afr.vect <- intersect(rownames(finalGeno.afr.df), rownames(validGeno.afr.df))

finalGeno.afr.dt <- synSeq.formatResultGeno( finalGeno.afr.df[sharedSampleID.afr.vect,], finalCopy.afr.df[sharedSampleID.afr.vect,] )
validGeno.afr.dt <- synSeq.formatResultGeno( validGeno.afr.df[sharedSampleID.afr.vect,], validCopy.afr.df[sharedSampleID.afr.vect,] )

## Concordance comparison
copyCon.afr.list <- pingValid.copyConcordance( finalCopy.afr.df, validCopy.afr.df )
#genoCon.afr.list <- pingValid.genoConcordance( finalGeno.afr.dt, validGeno.afr.dt, kirRes = 3, score.unresolved = F, score.copy4 = F )
genoConMakup.afr.list <- pingValid.genoConcordanceMakeup( finalGeno.afr.dt, validGeno.afr.dt, kirRes = 3, score.copy4 = F )

for( currentLocus in names( genoConMakup.afr.list ) ){
  matchN <- genoConMakup.afr.list[[currentLocus]]$match
  totalN <- genoConMakup.afr.list[[currentLocus]]$N
  unresolvedN <- genoConMakup.afr.list[[currentLocus]]$unresolved
  percMatch <- round(matchN / (totalN - unresolvedN),digits = 3)
  cat(paste0('\n',currentLocus,'\t',percMatch,'\t',(totalN - unresolvedN)))
}

# Graphing
gene <- sort( rep( names( genoConMakup.afr.list ), 3) )
condition <- rep( c('match','unresolved','alleleMismatch'), length(names(genoConMakup.afr.list)) )
value <- c()
for( currentLocus in unique(gene) ){
  for( currentCondition in unique(condition) ){
    value <- c(value, genoConMakup.afr.list[[currentLocus]][[currentCondition]])
  }
}
plotData.afr <- data.frame(gene, condition, value)

plotData.afr$condition <- factor( plotData.afr$condition, levels=c('alleleMismatch','unresolved','match'))

# ----- Khoisan dataset concordance table -----
for( currentLocus in unique(plotData.afr$gene) ){
  sub.df <- plotData.afr[plotData.afr$gene == currentLocus,c('condition','value'),drop=F]
  totalN <- sum(sub.df$value)
  valueVect <- sub.df$value / totalN
  cat(paste0('\n',currentLocus,'\t',paste0(round(valueVect,3),collapse='\t'),'\t',totalN))
}

g1 <- ggplot(plotData.afr, aes(fill=condition, y=value, x=gene)) +
  geom_bar( position='fill',stat='identity') +
  theme_classic() + 
  scale_fill_manual(values=c('#b0b0b0','#595959','#000000')) +
  theme(axis.text.x = element_text(angle = 90),text = element_text(size=22))
g1
ggsave(file.path(resultsDirectory,"conc.afr.percent.png"), width = 30, height = 15, units = "cm")


g2 <- ggplot(plotData.afr, aes(fill=condition, y=value, x=gene)) +
  geom_bar( position='stack',stat='identity') +
  theme_classic() + 
  scale_fill_manual(values=c('#b0b0b0','#595959','#000000')) +
  theme(axis.text.x = element_text(angle = 90),text = element_text(size=22))
g2
ggsave(file.path(resultsDirectory,"conc.afr.value.png"), width = 30, height = 15, units = "cm")



# ----- copy discordance correction -----
discord.idVect <- unique( unlist( lapply( copyCon.list, function(x) x$mismatchIDVect ), use.names = F ) )
discord.idVect <- discord.idVect[ !is.na(discord.idVect) ]

validLoci <- names(copyCon.list)

RvalueVect <- copyCall.df[discord.idVect,validLoci][validCopyCall.df[discord.idVect,validLoci] == 'R']
validCopyCall.df[discord.idVect,validLoci][validCopyCall.df[discord.idVect,validLoci] == 'R'] <- RvalueVect

copyCall.df[discord.idVect,validLoci] <- validCopyCall.df[discord.idVect,validLoci]
write.csv(copyCall.df,file='manualCopy_380_fixed.csv')




# ----- Unresolved/mismatch analysis -----
'2/25/2021 synthetic unresolved analysis'
currentLocus <- 'KIR3DP1'

syn.mismatchID.list <- lapply( genoConMakup.syn.list, function(x) x$mismatchIDs)
syn.mismatchID.list <- syn.mismatchID.list[ lapply( syn.mismatchID.list, length ) != 0 ]
run.write_list(syn.mismatchID.list, '~/Desktop/synMismatchIDs.txt')
run.read_list('~/Desktop/synMismatchIDs.txt')


cbind( finalGeno.syn.df[syn.mismatchID.list$KIR2DS35, 'KIR2DS35',drop=F], validGeno.syn.df[syn.mismatchID.list$KIR2DS35, 'KIR2DS35',drop=F] )



genoCon.syn.list$KIR2DL1$mismatchIDs # 3_copy mismatches 5_1101 mismatches
genoCon.afr.list$KIR3DL2$mismatchIDs # 9 / 14 *064 unresolved

S1.syn.mismatchVect <- genoCon.syn.list$KIR2DS1$mismatchIDs
L2.syn.mismatchVect <- genoCon.syn.list$KIR2DL2$mismatchIDs

intersect(L1.syn.mismatchVect, S1.syn.mismatchVect)
intersect(L1.syn.mismatchVect, L2.syn.mismatchVect)


#KIR2DL1 most connected with KIR2DS1, then KIR2DL2. Much lower connection to KIR2DS35, KIR2DP1 and KIR2DS4
KIR2DL1.syn.mismatchVect <- genoCon.syn.list$KIR2DL1$mismatchIDs # 23
KIR2DS1.syn.mismatchVect <- genoCon.syn.list$KIR2DS1$mismatchIDs # 20
KIR2DL23.syn.mismatchVect <- genoCon.syn.list$KIR2DL23$mismatchIDs # 24

KIR2DS35.syn.mismatchVect <- genoCon.syn.list$KIR2DS35$mismatchIDs # 23
KIR2DP1.syn.mismatchVect <- genoCon.syn.list$KIR2DP1$mismatchIDs # 8
KIR2DS4.syn.mismatchVect <- genoCon.syn.list$KIR2DS4$mismatchIDs # 8

length( intersect( KIR2DL1.syn.mismatchVect, KIR2DS1.syn.mismatchVect) ) # 13 / 20
length( intersect( KIR2DL1.syn.mismatchVect, KIR2DL23.syn.mismatchVect) ) # 16 / 24
length( intersect( KIR2DL1.syn.mismatchVect, KIR2DS35.syn.mismatchVect) ) # 12 / 23
length( intersect( KIR2DL1.syn.mismatchVect, KIR2DP1.syn.mismatchVect) ) # 4 / 8
length( intersect( KIR2DL1.syn.mismatchVect, KIR2DS4.syn.mismatchVect) ) # 3 / 8

KIR2DL1S1.leftoverVect <- setdiff( KIR2DL1.syn.mismatchVect, KIR2DS1.syn.mismatchVect )
setdiff( KIR2DL1S1.leftoverVect, KIR2DL23.syn.mismatchVect )

# 22 / 23 discordant KIR2DL1 samples are also discordant for either KIR2DL23 or KIR2DS1, the last discordant sample is also discordant for KIR2DP1


# KIR2DL4 misidentified delA
KIR2DL4.syn.mismatchVect <- genoCon.syn.list$KIR2DL4$mismatchIDs

KIR2DL4.snp.dir <- file.path(resultsDirectory,'syn_genotypes/syn_KIR2DL4_discordance')
KIR2DL4.snp.path.vect <- list.files(KIR2DL4.snp.dir,pattern='csv',full.names=T)

for( snpPath in KIR2DL4.snp.path.vect ){
  test.df <- read.csv(snpPath,header=T,check.names=F, row.names=1, stringsAsFactors = F,colClasses = 'character')
  cat('\n',snpPath,test.df[,c('E7_95')])
}

 
#test.df <- read.csv('/Users/wmarin/Documents/PING/PING_projects/data/2020_12_07_pingReborn_concordance//syn_genotypes/syn_KIR2DL4_discordance/final_KIR2DL4_synSeq.hiseq.dp50.rl150.40.paired_SNP.csv',header=T,check.names=F, row.names=1, stringsAsFactors = F,colClasses = 'character')
#test.df[,c('E7_95','E7_96','E7_97','E7_98','E7_99','E7_100','E7_101','E7_102','E7_103','E7_104','E7_105')]

# 15 / 17 discordant samples due to INS misplacement
# 2 / 17 discordant samples due to DEL misplacement

# KIR2DS35 discordance
KIR2DS35.discord.df <- cbind( validGeno.syn.df[KIR2DS35.syn.mismatchVect,'KIR2DS35'], finalGeno.syn.df[KIR2DS35.syn.mismatchVect,'KIR2DS35'] ) 
rownames(KIR2DS35.discord.df) <- KIR2DS35.syn.mismatchVect
colnames(KIR2DS35.discord.df) <- c('Valid','PING')
intersect(KIR2DS35.syn.mismatchVect, KIR2DL1.syn.mismatchVect)

refInfo.syn.dir <- file.path(resultsDirectory,'syn_genotypes/syn_refInfo')
refInfo.syn.pathVect <- list.files(refInfo.syn.dir,pattern='txt',full.names=F)

refInfo.syn.pathVect <- refInfo.syn.pathVect[ tstrsplit(refInfo.syn.pathVect, '.paired.refInfo.txt')[[1]] %in% KIR2DS35.syn.mismatchVect ]

for( refInfo.path in refInfo.syn.pathVect ){
  refAlleleVect.S5 <- paste0( grep('KIR2DS5',alleleSetup.readAnswerKey(file.path(refInfo.syn.dir,refInfo.path))[[1]][[1]],value = T), collapse=' ' )
  refAlleleVect.S3 <- paste0( grep('KIR2DS3',alleleSetup.readAnswerKey(file.path(refInfo.syn.dir,refInfo.path))[[1]][[1]],value = T), collapse=' ' )
  validAlleleVect <- validGeno.syn.df[strsplit( refInfo.path, '.paired.refInfo.txt' )[[1]],'KIR2DS35']
  pingAlleleVect <- finalGeno.syn.df[strsplit( refInfo.path, '.paired.refInfo.txt' )[[1]],'KIR2DS35']
  cat(paste0('\n',strsplit( refInfo.path, '.paired.refInfo.txt' )[[1]],'\t',refAlleleVect.S3,'\t',refAlleleVect.S5,'\t\t',validAlleleVect,'\t',pingAlleleVect))
}

alleleSetup.readAnswerKey <- function(keyFile){
  out.list <- list()
  con  <- file(keyFile, open = "r")
  
  while (length(oneLine <- readLines(con, n = 1)) > 0) {
    lineVect <- unlist((strsplit(oneLine, "\t",fixed=T)))
    sampleID <- lineVect[1]
    genoStr <- lineVect[2]
    
    genoVect <- unlist(strsplit(genoStr,'_',fixed=T))
    out.list[[sampleID]] <- list('genoVect'=genoVect)
  }
  
  close(con)
  
  return(out.list)
}

# KIR2DL23 discordance
finalGeno.syn.df[KIR2DL23.syn.mismatchVect,'KIR2DL23']
KIR2DS2.syn.mismatchVect <- genoCon.syn.list$KIR2DS2$mismatchIDs # 12
length(KIR2DS2.syn.mismatchVect)
setdiff( setdiff( setdiff( KIR2DL23.syn.mismatchVect, KIR2DL1.syn.mismatchVect ), KIR2DS1.syn.mismatchVect), KIR2DS2.syn.mismatchVect )

View( finalGeno.syn.df[setdiff( KIR2DL23.syn.mismatchVect, KIR2DL1.syn.mismatchVect ),] )


# KIR3DL1S1 discordance
KIR3DL1S1.syn.mismatchVect <- genoCon.syn.list$KIR3DL1S1$mismatchIDs # 12
finalGeno.syn.df[KIR3DL1S1.syn.mismatchVect,'KIR3DL1S1',drop=F]

KIR2DS4.syn.mismatchVect <- genoCon.syn.list$KIR2DS4$mismatchIDs # 8

finalGeno.syn.df[setdiff( KIR3DL1S1.syn.mismatchVect, KIR2DS4.syn.mismatchVect ),'KIR3DL1S1',drop=F] # 5 shared


# KIR3DP1 discordance
KIR3DP1.syn.mismatchVect <- genoCon.syn.list$KIR3DP1$mismatchIDs # 12


# KIR2DS1 discordance
KIR2DS1.syn.mismatchVect <- genoCon.syn.list$KIR2DS1$mismatchIDs # 20


# KIR3DL3 discordance
KIR3DL3.syn.mismatchVect <- genoCon.syn.list$KIR3DL3$mismatchIDs # 20




rerun.idVect <- c( "IND00190", "IND00213", "IND00223", "IND00365", "IND00374", "IND00016", "IND00310", "IND00063", "IND00127", "IND00177", "IND00229", "IND00293", "IND00166", "IND00178",
                   "IND00255", "IND00277", "IND00355", "IND00010", "IND00052", "IND00101", "IND00175", "IND00375", "IND00243", "IND00073", "IND00133", "IND00168", "IND00203",
                   "IND00161" ,"IND00284", "IND00114")



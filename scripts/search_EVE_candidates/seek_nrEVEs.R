### input: genome sequence (fasta file)
### output: a list of searched polyAs and TSDs (text file). Candidated nrEVEs (fasta file).

### read data (FASTA)
data = scan( "xxx.fasta", sep="	", what = "" )

### setting parameters
SUCCESSIVE = 7 # minimum number of times that the base A is successively appeared. 
TOTAL_COUNT = 8 # minimum total number of the base A
TSD_LENGTH = 7 # length of TSD
QUERY_START = -3 # start point of query for TSD search
UPPER_STREAM = 8000 # start point of upper stream for TSD search
UPPER_STREAM_END = 300 # end point of upper stream for TSD search
GAP = 0 # how many times can we tolerance other bases between the bases A and A?
TSD_mismatch = 0 # how many times can we tolerance the number of mismatches between the query and explored sequences on TSD search


RESULT = NULL
RESULT_A = NULL
RESULT_T = NULL
seq_genome = unlist( noquote( strsplit( data[-1], NULL ) ) )
seq_genome_rev = rev( seq_genome )
Chrom = data[1]
AllStart = 1
AllEnd = length(seq_genome)
StartPoint = NULL
EndPoint = NULL
StartPointT = NULL
EndPointT = NULL
count = 0
gap = 0
Successive = 0
successive = 0
countT = 0
gapT = 0
SuccessiveT = 0
successiveT = 0
for( i in 1:length(seq_genome) )
{
  #########################################################
  # search polyA
  #########################################################
  if( seq_genome[i] == "A" ){
    if( count == 0 ) startpoint = i
    count = count + 1
    gap = 0
  }
  
  if( i > 1 )
  {
    if( ( seq_genome[i] == "A" ) && ( seq_genome[i-1] == "A" ) ) successive = successive + 1
  }
  
  if( ( seq_genome[i] != "A" ) && ( count != 0 ) ){
    gap = gap + 1
    if( Successive <= successive ) Successive = successive
    successive = 0
    if( gap > GAP ){
      if( count >= TOTAL_COUNT ){
        if( ( Successive >= SUCCESSIVE ) || ( successive >= SUCCESSIVE ) ){
          StartPoint = c( StartPoint, startpoint )
          EndPoint = c( EndPoint, i )
          startpoint = count = gap = successive = Successive = 0
        }else{
          startpoint = count = gap = successive = Successive = 0
        }
        startpoint = count = Successive = 0
      }
      startpoint = count = Successive = 0
    }
  }
    
  #########################################################
  # search polyT
  #########################################################
  if( seq_genome_rev[i] == "T" ){
    if( countT == 0 ) startpointT = i
    countT = countT + 1
    gapT = 0
  }
  
  if( i > 1 )
  {
    if( ( seq_genome_rev[i] == "T" ) && ( seq_genome_rev[i-1] == "T" ) ) successiveT = successiveT + 1
  }
  
  if( ( seq_genome_rev[i] != "T" ) && ( countT != 0 ) ){
    gapT = gapT + 1
    if( SuccessiveT <= successiveT ) SuccessiveT = successiveT
    successiveT = 0
    if( gapT > GAP ){
      if( countT >= TOTAL_COUNT ){
        if( ( SuccessiveT >= SUCCESSIVE ) || ( successiveT >= SUCCESSIVE ) ){
          StartPointT = c( StartPointT, startpointT )
          EndPointT = c( EndPointT, i )
          startpointT = countT = gapT = successiveT = SuccessiveT = 0
        }else{
          startpointT = countT = gapT = successiveT = SuccessiveT = 0
        }
        startpointT = countT = SuccessiveT = 0
      }
      startpointT = countT = SuccessiveT = 0
    }
  }
}
  
### summarizing the results of polyA search
Result = matrix( 0, length(StartPoint), 2 )
if( length(Result) != 0 ){
  for( i in 1:length(StartPoint) )
  {
    NewPosition = paste( AllStart + ( StartPoint[i] - 1 ), AllStart + ( EndPoint[i] - 1 ), sep="-" )
    Result[ i, 1 ] = paste( Chrom, NewPosition, sep=":" )
    Result[ i, 2 ] = paste( seq_genome[StartPoint[i]: EndPoint[i]], collapse="" )
  }
  RESULT = rbind( RESULT, Result )
}
RESULT_A = rbind( RESULT_A, Result )
  
### summarizing the results of polyT search
ResultT = matrix( 0, length(StartPointT), 2 )
if( length(ResultT) != 0 ){
  for( i in 1:length(StartPointT) )
  {
    NewPosition = paste( AllEnd - ( EndPointT[i] - 1 ), AllEnd - ( StartPointT[i] - 1 ), sep="-" )
    ResultT[ i, 1 ] = paste( Chrom, NewPosition, sep=":" )
    ResultT[ i, 2 ] = paste( rev( seq_genome_rev[StartPointT[i]: EndPointT[i]] ), collapse="" )
  } 
  RESULT = rbind( RESULT, ResultT )
}
RESULT_T = rbind( RESULT_T, ResultT )

#########################################################
### TSD search
#########################################################

RESULT_TSD = NULL
RESULT_SEQ = NULL
All_suffix = seq( AllStart, AllEnd, by=1 )
seq_genome_length = length(seq_genome)

### TSD search from polyA results
for( itr in 1:nrow(RESULT_A) )
{
  seq_genome_num_A = RESULT_A[itr, ]
  Chrom_A = strsplit( seq_genome_num_A, ":" )[[1]][1]
  obj_A = strsplit( seq_genome_num_A, ":" )[[1]][2]
  AllStart_A = as.numeric( strsplit( obj_A, "-" )[[1]][1] )
  AllEnd_A = as.numeric( strsplit( obj_A, "-" )[[1]][2] )
  
  result_query = NULL
  if( ( AllStart <= ( AllStart_A - UPPER_STREAM - 25 ) ) && ( ( AllEnd_A  - ( GAP - QUERY_START ) + ( TSD_LENGTH - 1 ) ) <= seq_genome_length ) )
  {
    ### preparation of query
    query = seq_genome[ ( AllEnd_A  - ( GAP - QUERY_START ) ) : ( AllEnd_A  - ( GAP - QUERY_START ) + ( TSD_LENGTH - 1 ) ) ]
    query_posi = paste( AllEnd_A  - ( GAP - QUERY_START ), AllEnd_A  - ( GAP - QUERY_START ) + ( TSD_LENGTH - 1 ), sep="-" )
    query_posi = paste( Chrom_A, query_posi, sep=":" )
    query_posi = paste( query_posi, paste( query, collapse="" ), sep="," )
    
    if( ( is.na( sum( query == "N" ) ) == FALSE ) && ( ( AllEnd_A  - ( GAP - QUERY_START ) + ( TSD_LENGTH - 1 ) + 25 ) <= seq_genome_length) ){
      if( ( sum( query == "N" ) == 0 ) && ( sum( query == "a" ) == 0 ) && ( sum( query == "t" ) == 0 ) && ( sum( query == "g" ) == 0 ) && ( sum( query == "c" ) == 0 ) ){
        query_start = AllStart_A - UPPER_STREAM
        query_range = seq_genome[ ( AllEnd_A  - ( GAP - QUERY_START ) - 25 ) : ( AllEnd_A  - ( GAP - QUERY_START ) + ( TSD_LENGTH - 1 ) + 25 ) ]
        for( j in 1:( UPPER_STREAM - UPPER_STREAM_END ))
        {
          explore = seq_genome[ ( query_start + j - 1 )  : ( query_start + ( TSD_LENGTH - 1 ) + ( j - 1 ) )   ]
          explore_range = seq_genome[ ( ( query_start + j - 1 ) - 25 ) : ( ( query_start + ( TSD_LENGTH - 1 ) + ( j - 1 ) )  + 25)  ]
          temp = sum( query != explore )
          temp2 = mean( query_range == explore_range )
          temp3 = seq_genome[ (query_start+TSD_LENGTH-1+j) : (AllStart_A-1) ]
          if( ( temp <= TSD_mismatch ) && ( temp2 < 0.8 ) && ( mean(temp3!="N") >= 0.5 ) ){
            NewPosition = paste( query_start + j - 1, query_start + ( TSD_LENGTH - 1 ) + ( j - 1 ), sep="-" )
            Explore_posi = paste( Chrom_A, NewPosition, sep=":" )
            Explore = paste( explore, collapse="" )
            result_query = c( result_query, paste( Explore_posi, Explore, sep="," ) )
            NewPosition_seq = paste( query_start + TSD_LENGTH - 1 + j, AllStart_A - 1, sep="-" )
            Explore_posi_seq = paste( Chrom_A, NewPosition_seq, sep=":" )
            RESULT_SEQ = rbind( RESULT_SEQ, Explore_posi_seq )
            RESULT_SEQ = rbind( RESULT_SEQ, paste(seq_genome[ (query_start+TSD_LENGTH-1+j) : (AllStart_A-1) ], collapse="" ) )
          }
        }
        result_temp = cbind( rep( query_posi, length( result_query ) ), result_query )
        result_temp = cbind( rep( paste( seq_genome_num_A, collapse="," ), length( result_query ) ), result_temp )
        RESULT_TSD = rbind( RESULT_TSD, result_temp )
      }
    }
  }
  cat( " " )
  cat( itr )
}

### TSD search from polyT results
for( itr in 1:nrow(RESULT_T) )
{
  seq_genome_num_T = RESULT_T[itr, ]
  Chrom_T = strsplit( seq_genome_num_T, ":" )[[1]][1]
  obj_T = strsplit( seq_genome_num_T, ":" )[[1]][2]
  AllStart_T = as.numeric( strsplit( obj_T, "-" )[[1]][1] )
  AllEnd_T = as.numeric( strsplit( obj_T, "-" )[[1]][2] )
  
  result_query = NULL
  if( ( AllEnd >= ( AllEnd_T + UPPER_STREAM + 25 ) ) && ( AllStart_T  +  ( GAP - QUERY_START )  - ( TSD_LENGTH - 1 ) >= 1 ) )
  {
    ### preparation of query
    query = seq_genome[ ( AllStart_T  +  ( GAP - QUERY_START )  - ( TSD_LENGTH - 1 ) ) : ( AllStart_T  +  ( GAP - QUERY_START )  ) ]
    query_posi = paste( AllStart_T  +  ( GAP - QUERY_START )  - ( TSD_LENGTH - 1 ), AllStart_T  +  ( GAP - QUERY_START ), sep="-" )
    query_posi = paste( Chrom_T, query_posi, sep=":" )
    query_posi = paste( query_posi, paste( query, collapse="" ), sep="," )

    if( ( is.na( sum( query == "N" ) ) == FALSE ) && ( ( AllStart_T  +  ( GAP - QUERY_START )  - ( TSD_LENGTH - 1 ) - 25 ) >= 1  ) ){
      if( ( sum( query == "N" ) == 0 ) && ( sum( query == "a" ) == 0 ) && ( sum( query == "t" ) == 0 ) && ( sum( query == "g" ) == 0 ) && ( sum( query == "c" ) == 0 ) ){
        query_start = AllEnd_T + UPPER_STREAM
        query_range = seq_genome[ ( AllStart_T  +  ( GAP - QUERY_START )  - ( TSD_LENGTH - 1 ) - 25 ) : ( AllStart_T  +  ( GAP - QUERY_START )  + 25) ]
        for( j in 1:( UPPER_STREAM - UPPER_STREAM_END ))
        {
          explore = seq_genome[ ( query_start - ( TSD_LENGTH - 1 ) - j  ) : ( query_start - j ) ]
          explore_range = seq_genome[ ( ( query_start - ( TSD_LENGTH - 1 ) - j  ) - 25 ) : ( ( query_start - j ) + 25 ) ]
          temp = sum( query != explore )
          temp2 = mean( query_range == explore_range )
          temp3 = seq_genome[ (AllEnd_T + 1) : (query_start - ( TSD_LENGTH - 1 ) - j - 1) ]
          if( ( temp <= TSD_mismatch ) && ( temp2 < 0.8 ) && ( mean(temp3!="N") >= 0.5 ) ){
            NewPosition = paste( query_start - ( TSD_LENGTH - 1 ) - j , query_start - j, sep="-" )
            Explore_posi = paste( Chrom_T, NewPosition, sep=":" )
            Explore = paste( explore, collapse="" ) 
            result_query = c( result_query, paste( Explore_posi, Explore, sep="," ) )
            
            ### change complementary strand
            Complementary_Seq <- seq_genome[ (AllEnd_T + 1) : (query_start - ( TSD_LENGTH - 1 ) - j - 1) ]
            complementary_seq <- Complementary_Seq
            for(i in 1:length(Complementary_Seq))
            {
              if( Complementary_Seq[i] == "A" ) complementary_seq[i] = "T"
              if( Complementary_Seq[i] == "T" ) complementary_seq[i] = "A"
              if( Complementary_Seq[i] == "G" ) complementary_seq[i] = "C"
              if( Complementary_Seq[i] == "C" ) complementary_seq[i] = "G"
            }
            complementary_seq = rev(complementary_seq)

            NewPosition_seq = paste( query_start - ( TSD_LENGTH - 1 ) - j - 1, AllEnd_T + 1, sep="-" )
            Explore_posi_seq = paste( Chrom_T, NewPosition_seq, sep=":" )
            RESULT_SEQ = rbind( RESULT_SEQ, Explore_posi_seq )
            RESULT_SEQ = rbind( RESULT_SEQ, paste(complementary_seq, collapse="" ) )
          }
        }
        result_temp = cbind( rep( query_posi, length( result_query ) ), result_query )
        result_temp = cbind( rep( paste( seq_genome_num_T, collapse="," ), length( result_query ) ), result_temp )
        RESULT_TSD = rbind( RESULT_TSD, result_temp )
      }
    }
  }
  cat( " " )
  cat( itr )
}
colnames(RESULT_TSD) = c("polyA or polyT", "TSD_query", "TSD_search")

# change the bases a, t, g, c into N
extract_id = 0
for(i in 1:nrow(RESULT_TSD))
{
  temp = strsplit( RESULT_SEQ[2*i,1], NULL)[[1]]
  temp = gsub("a", "N", temp)
  temp = gsub("t", "N", temp)
  temp = gsub("g", "N", temp)
  temp = gsub("c", "N", temp)
  if( mean(temp=="N") >= 0.5 ) extract_id = c(extract_id, i)
}
if( sum(extract_id) != 0 )
{
  extract_id = extract_id[-1]
  RESULT_TSD = RESULT_TSD[ -extract_id, ]
  RESULT_SEQ = RESULT_SEQ[ -c(2*extract_id, 2*extract_id-1), ]
}

out_f <- "result.txt"
out_f_fasta <- "result.fasta"
### output of searched polyAs and TSDs
write.table(RESULT_TSD, out_f, sep="\t", quote=FALSE, row.names=FALSE)
### output of candidated nrEVEs
write.table(RESULT_SEQ, out_f_fasta, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

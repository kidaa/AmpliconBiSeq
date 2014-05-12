#concerts all C to T in a non CG context
bisConv=function(in.seq){
  
  if(class(in.seq)=="DNAStringSet"){
    C.pos=vmatchPattern('C', in.seq)
    CG.pos=vmatchPattern('CG', in.seq)
    x=do.call(c,
          lapply(seq(length(in.seq)),function(i){
            Cind=start(C.pos[[i]])
            CGind=start(CG.pos[[i]])
            convind=Cind[!Cind %in% CGind]
            DNAStringSet(replaceLetterAt(in.seq[[i]], convind, rep('T',length(convind))))
          })
    )
    names(x)=names(in.seq)
  }else if(class(in.seq)=="DNAString"){
    C.pos=matchPattern('C', in.seq)
    CG.pos=matchPattern('CG', in.seq)
    
    Cind=start(C.pos)
    CGind=start(CG.pos)
    convind=Cind[!Cind %in% CGind]
    x= replaceLetterAt(in.seq, convind, rep('T',length(convind)) )     
  }
  x
}

# internal function that calls primer3 and reads the output
# seq=sample.seq[[1]];name=names(sample.seq)[1]
.callP3Nread<-function(seq,size_range='151-500',Tm=c(55,57,58),name,
                       primer3="/work2/gschub/arnaud/softwares/primer3/primer3-2.3.4/bin/primer3_core",
                       thermo.param="/work2/gschub/arnaud/softwares/primer3/primer3-2.3.4/src/primer3_config/",
                       settings="/work2/gschub/arnaud/softwares/primer3/primer3-2.3.4/default_settings.txt"){
  
  # get C and CpGs in the sequence
  CG.pos=start(matchPattern('CG', seq) ) 
  C.pos=start(matchPattern('C', seq) )
  bis.seq=bisConv(seq) # convert sequence to bisulfite shit
  # exlcude CpGs
  excluded.regions=do.call(paste,(lapply(CG.pos, paste, '2', sep=',')))
  
  #print(excluded.regions)
  # make primer 3 input file
  p3.input=tempfile()
  p3.output=tempfile()
  write(
    paste( sprintf("SEQUENCE_ID=%s\n",name  ),
           sprintf("SEQUENCE_TEMPLATE=%s\n",as.character(bis.seq)),
           "PRIMER_TASK=pick_detection_primers\n",
           "PRIMER_PICK_LEFT_PRIMER=1\n" ,
           "PRIMER_PICK_INTERNAL_OLIGO=0\n",
           "PRIMER_PICK_RIGHT_PRIMER=1\n"  ,
           "PRIMER_EXPLAIN_FLAG=1\n"  ,
           "PRIMER_PAIR_MAX_DIFF_TM=3\n",
           sprintf("PRIMER_MIN_TM=%s\n" ,Tm[1]),
           sprintf("PRIMER_OPT_TM=%s\n" ,Tm[2]),
           sprintf("PRIMER_MAX_TM=%s\n" ,Tm[3]),
           sprintf("SEQUENCE_EXCLUDED_REGION=%s\n" ,excluded.regions),
           sprintf("PRIMER_PRODUCT_SIZE_RANGE=%s\n" ,size_range),
           
           sprintf("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s\n" ,thermo.param),
           "=",
           sep=''
    )
    ,
    p3.input
  )
  #call primer 3 and store the output in a temporary file
  
  try(system(
    paste(primer3 ,p3.input, "-p3_settings_file",settings, 
                ">", p3.output)
  ))
  
  #import and parse the output into a dataframe named designed.primers
  out=read.delim(p3.output, sep='=', header=FALSE)
  
  unlink(c(p3.input,p3.output) ) # delete temp files
  returned.primers=as.numeric(as.vector(out[out[,1]=='PRIMER_PAIR_NUM_RETURNED',][,2]))
  if (length(returned.primers)==0){ warning('primers not detected for ',name,call. = FALSE);return(NA)}
  if ((returned.primers)==0){ warning('primers not detected for ',name,call. = FALSE);return(NA)}
  if (returned.primers>0){
    designed.primers=data.frame()
    for (i in seq(0,returned.primers-1,1)){
      #IMPORT SEQUENCES
      id=sprintf(  'PRIMER_LEFT_%i_SEQUENCE',i)
      PRIMER_LEFT_SEQUENCE=as.character(out[out[,1]==id,][,2])
      id=sprintf(  'PRIMER_RIGHT_%i_SEQUENCE',i)
      PRIMER_RIGHT_SEQUENCE=as.character(out[out[,1]==id,][,2])
      
      #IMPORT PRIMING POSITIONS
      id=sprintf(  'PRIMER_LEFT_%i',i)
      PRIMER_LEFT=as.numeric(unlist(strsplit(as.vector((out[out[,1]==id,][,2])),',')))
      #PRIMER_LEFT_LEN=as.numeric(unlist(strsplit(as.vector((out[out[,1]==id,][,2])),',')))
      id=sprintf(  'PRIMER_RIGHT_%i',i)
      PRIMER_RIGHT=as.numeric(unlist(strsplit(as.vector((out[out[,1]==id,][,2])),',')))
      #IMPORT Tm
      id=sprintf(  'PRIMER_LEFT_%i_TM',i)
      PRIMER_LEFT_TM=as.numeric(as.vector((out[out[,1]==id,][,2])),',')
      id=sprintf(  'PRIMER_RIGHT_%i_TM',i)
      PRIMER_RIGHT_TM=as.numeric(as.vector((out[out[,1]==id,][,2])),',')
      
      res=out[grep(paste("_",i,"_",sep=""),out[,1]),]
      extra.inf=t(res)[2,,drop=FALSE]
      colnames(extra.inf)=sub( paste("_",i,sep=""),"",res[,1])
      extra.inf=extra.inf[,-c(4:9),drop=FALSE] # remove redundant columns
      extra.inf=apply(extra.inf,2,as.numeric)
      #Aggegate in a dataframe
      primer.info=data.frame(i,
                             PRIMER_LEFT_SEQUENCE,PRIMER_RIGHT_SEQUENCE,
                             PRIMER_LEFT_TM, PRIMER_RIGHT_TM,
                             PRIMER_LEFT_pos=PRIMER_LEFT[1],
                             PRIMER_LEFT_len=PRIMER_LEFT[2], 
                             PRIMER_RIGHT_pos=PRIMER_RIGHT[1], 
                             PRIMER_RIGHT_len=PRIMER_RIGHT[2],
                             t(data.frame(extra.inf))
                             
      )
      rownames(primer.info)=NULL
      designed.primers=rbind(designed.primers, primer.info)
      #print(primer.info)
    }
    
    
    
    #colnames(designed.primers)=c('PrimerID',
    #                             'Fwseq','Rvseq',
    #                             'FwTm','RvTm',
    #                             'FwPos','Fwlen',
    #                             'RvPos','Rvlen',
    #                             'fragLen' )
    
    #Rank acording to the number of CGs in the fragment
    total.CG=length(unlist(CG.pos))
    designed.primers$total.CG=total.CG
    covered.CG=unlist(
      lapply(seq(length(designed.primers$PRIMER_LEFT_pos)) ,function(i){
        sum(unlist(CG.pos) > designed.primers$PRIMER_LEFT_pos[i] & unlist(CG.pos) < designed.primers$PRIMER_RIGHT_pos[i]  )
      })
    )
    designed.primers$covered.CG=covered.CG
    #count the number of Cs covered
    FwC.covered=unlist(
      lapply(seq(length(designed.primers$PRIMER_LEFT_pos)) ,function(i){
        sum(unlist(C.pos) > designed.primers$PRIMER_LEFT_pos[i] & unlist(C.pos) <= designed.primers$PRIMER_LEFT_pos[i]+designed.primers$PRIMER_LEFT_len[i]  )
      })
    )
    
    RvC.covered=unlist(
      lapply(seq(length(designed.primers$PRIMER_RIGHT_pos)) ,function(i){
        sum( (unlist(C.pos) < designed.primers$PRIMER_RIGHT_pos[i]) & unlist(C.pos) >= designed.primers$PRIMER_RIGHT_pos[i]-designed.primers$PRIMER_RIGHT_len[i] )
      })
    )
    
    designed.primers$LeftC.covered=FwC.covered
    designed.primers$RightC.covered=RvC.covered
    
  }
  return(designed.primers)
}

#' Design primers for amplicons
#' 
#' Designs primers for Amplicon Bisulfite sequencing experiments
#' 
#'  @param target  \code{DNAStringSet} or \code{GRanges} object. If it is a GRanges
#'                object you need to provide BSgenome package with it. If it is a \code{DNAStringSet} object
#'                then names argument should contain chromomsome,start and end, by following
#'                "chr_start_end" naming convetion.
#'  @param genome default is NULL, it needs to be set only if target is set 
#'                to be a \code{GRanges} object
#'  @param tag extra information about amplicons, this will be concatanated to the
#'             names of BioStringSet object which is either created internally or provided
#'             by the user by supplying 'target' argument with a \code{BioStringSet} object
#'  @param primer3 the filesystem path to primer3 executable, default value NULL will use the version
#'         installed with the package
#'  @param settings text file for p3 settings. Default value NULL will use the default settings file
#'         installed with the package. 
#'  @param termoParam Location for thermodynamic parameters for primer3, this should designate a directory,
#'         not a single file. Default value NULL, will use the directory that is installed with
#'         the package.
#'  @param sizeRange a two element vector of integers for size range of the amplicons
#'  @param Tm a three element vector having minimum,optimum and maximum melting tempratures
#'            for primer desing (in that order).
#'  @param ncores number of cores to run on. default=1   
#'  
#'  @author Based on Arnaud Krebs' function, modified by Altuna Akalin   
#'   
#'   
#' @seealso \code{\link{primers2ranges}}, \code{\link{filterPrimers}}                     
#'    
#' @importFrom parallel mcmapply
#' @export
#' @docType methods         
designPrimers<-function(target,genome,tag=NULL,
                        primer3=NULL,
                        settings=NULL,
                        thermoParam=NULL,
                        sizeRange=c(151,500),
                        Tm=c(55,57,58),
                        ncores=4
                        ){
    if(is.null(primer3)){
      primer3=system.file("lib/x86_64/primer3_core", package="AmpliconBiSeq")
    }
    
    if(is.null(settings)){
      settings=system.file("lib/x86_64/default_settings.txt", package="AmpliconBiSeq")
      
    }
    
    if( is.null(thermoParam) ){
      thermoParam=system.file("lib/x86_64/primer3_config/", package="AmpliconBiSeq")
      
    }
    # get sequences if they are not provided by target argument
    if(class(target)=="GRanges" & !is.null(genome) ){
      cat("Extracting sequences from the genome...")
      sample.seq=getSeq(genome,target)
      names(sample.seq) = paste(seqnames(target),start(target),end(target),sep="_")
      names(sample.seq) = paste(names(sample.seq) ,tag,sep="|")
      cat("Sequence extraction completed.")
    }
    else if(class(target)=="GRanges" & ( !is.null(genome) | class(genome) != "BSgenome") ){
      stop("\nif 'target' is GRanges object, 'genome' can not be null and should be a BSgenome object\n")
    }
    else if(class(target)=="DNAStringSet" ){
      
        # check if the arguments have correct naming conventions if target is DNAStringSet
        if(!is.null(names(target)) & length(grep("_[[:digit:]]+_",names(target)[1]))>0  ) {
          sample.seq=target
        }else if( is.null(names(sample.seq)) ){
          stop("\nif argument 'target' is a DNAStringSet object it should have unique names assigned to each sequence\n",
               "you can set the names with names() function.\n",
               "Names must have 'chr_start_end|extraInfo|moreInfo|' convention\n",
               "where the first segment with underscores(_) denotes the location of the sequence in the genome\n",
               "and it should be seperated from the rest of the info with '|'\n")
        }else{
          stop("\nif argument 'target' is a DNAStringSet object it should have unique names assigned to each sequence\n",
               "you can set the names with names() function.\n",
               "Names must have 'chr_start_end|extraInfo|moreInfo|' convention\n",
               "where the first segment with underscores(_) denotes the location of the sequence in the genome\n",
               "and it should be seperated from the rest of the info with '|'\n")
        }
    }else{
      stop("\n'target' must be GRanges or DNAStringSet object")
    }
      
    #primer3="/work2/gschub/arnaud/softwares/primer3/primer3-2.3.4/bin/primer3_core"
    #thermo.param="/work2/gschub/arnaud/softwares/primer3/primer3-2.3.4/src/primer3_config/"
    #settings="/work2/gschub/arnaud/softwares/primer3/primer3-2.3.4/default_settings.txt"
    #.callP3Nread(seq,size_range='151-500',Tm=c(55,57,58),name, )
    size_range=paste(sizeRange[1],sizeRange[2],sep="-")
    cat("\nstarting primer3 jobs using",ncores,"core(s)\n")
    if(ncores==1){
      
      primers=mapply(.callP3Nread,seq=sample.seq,#[1:10],
                                  name=names(sample.seq),#[1:10]),
                     MoreArgs=list(primer3=primer3,size_range=size_range,Tm=Tm,
                                   thermo.param=thermoParam,settings=settings)
                     ,SIMPLIFY = FALSE)
    }else if(ncores>1){
      primers=parallel::mcmapply(.callP3Nread,
                                 seq=sample.seq,#[1:100],
                                 name=names(sample.seq),#[1:100]),
              MoreArgs=list(primer3=primer3,size_range=size_range,Tm=Tm,
                            thermo.param=thermoParam,settings=settings)
                                 ,SIMPLIFY = FALSE,mc.cores=ncores )
    }
    return(primers)
}


.filter.primers.list<-function(primers,minConPrimer,minCGonAmp){
  
  if(any(is.na(primers))){
    warning( "There are targets without primers\nfiltering those targets")
    primers=primers[ !is.na(primers) ]
  }
  
  lapply(primers,function(x) x[x$covered.CG>=minCGonAmp & x$LeftC.covered>= minConPrimer 
                               & x$RightC.covered>= minConPrimer ,]   )
}


.filter.primers.df<-function(primers,minConPrimer,minCGonAmp){
  
  x[x$covered.CG>=minCGonAmp & x$LeftC.covered>= minConPrimer 
                               & x$RightC.covered>= minConPrimer ,]   
}


#' filter primers based on thresholds
#' 
#' \code{filterPrimers} filters the primers based on covered CGs on amplicons
#' and number of non-CpG Cs on the primers
#' 
#' @param primers a list of primers returned by \code{designPrimers} or a data frame
#'                returned by \code{filterPrimers} function with as.data.frame=TRUE
#'                argument.
#' @param minConPrimer minimum number of non-CpG Cs on the primers
#' @param minCGonAmp   minimum number of CpGs covered on the targeted region
#' 
#' @author Based on Arnaud Krebs' function, modified by Altuna Akalin  
#' 
#' 
#' @examples
#'            data(bisPrimers)
#'           filt.primers=filterPrimers(bisPrimers,minConPrimer=1,minCGonAmp=1) 
#'           
#' @seealso \code{\link{primers2ranges}}, \code{\link{designPrimers}}
#'           
#' @export
#' @docType methods
filterPrimers<-function(primers,minConPrimer=1,minCGonAmp=2){
  
  if(class(primers)=="list"){
    .filter.primers.list(primers,minConPrimer,minCGonAmp)
      
  }else if(class(primers)=="data.frame"){
    .filter.primers.df(primers,minConPrimer,minCGonAmp)
  }
  
}



#' convert primer list to genomic intervals
#' 
#' function returns a GRanges or data frame object from a list of primers designed
#' by \code{designPrimers} function after calculating genomic location of the amplicon
#' targeted by the primers.
#' 
#' @param primers a list of primers returned by \code{designPrimers} function
#' @param as.data.frame logical indicating if a data frame should be returned
#'        instead of \code{GRanges} object.
#'        
#' @examples
#'  data(bisPrimers)
#'  gr.pr=primers2ranges(bisPrimers)
#'          
#' @seealso \code{\link{filterPrimers}}, \code{\link{designPrimers}}
#'        
#' @export
#' @docType methods
primers2ranges<-function(primers,as.data.frame=FALSE){
  
  if(any(is.na(primers))){
    warning( "There are targets without primers\nfiltering those before conversion")
    primers=primers[ !is.na(primers) ]
  }
  df=do.call("rbind",primers) # get primers to a df
  locs=gsub("\\|\\.+","",rownames(df)) # get the coordinates from list ids
  temp=do.call("rbind",strsplit(locs,"_")) #
  
  start=as.numeric(temp[,2])
  chr=as.character(temp[,1])
  
  
  
  amp.start= start + as.numeric(df$PRIMER_LEFT_pos)  
  amp.end  = start + as.numeric(df$PRIMER_RIGHT_pos)
  res=data.frame(chr=chr,start=amp.start,end=amp.end,df)
  #saveRDS(res,file="/work2/gschub/altuna/projects/DMR_alignments/all.designed.primers.to.amps.rds")
  if(as.data.frame)
  {
    return(res)
  }
  gr=GRanges(seqnames=res[,1],ranges=IRanges(res[,2],res[,3]) )
  values(gr)=res[,-c(1,2,3)]
  gr
}


#' check the quality of the spike-ins
#' 
#' This function provides methylation statistics for spike-ins. It plots a 
#' histogram or set of histograms from spike-in experiments which are helpful
#' to deduce conversion efficiency of the experiment.
#' 
#' @param proj \code{\link[QuasR]{qProject}} object from QuasR preferably 
#'    produced by ampBiSeqAlign
#' @param auxName string for which spike-ins should be checked, if equals "all"
#'        everything is checked
#' @param sampleName string for which sample the spike-ins should be checked
#' @param coverage minumum coverage before calculation of % methylation
#' @param targets named list of GRanges objects identifies targeted regions in
#'         given auxiliary genomes     
#' @param ... arguments to be passed to \code{\link{hist}} function
#' 
#' @examples # spikeCheck(proj,auxName="all",sampleName="ES.2")
#'           # spikeCheck(proj,auxName="T7",sampleName="ES.2")
#'           
#'           #la.gr=readRDS("/work2/gschub/Juliane/HTS/QuasRBiSeq/lamdaGR.rds")
#'           # par(mfrow=c(2,2))
#'           #spikeCheck(proj,auxName="all",sampleName="ES.2",targets=list(lambda=la.gr))
#'           
#' @return a list of percentages of methylated CpGs in all reads.
#'   
#' @importMethodsFrom Rsamtools scanBamHeader                  
#' @importFrom QuasR qMeth
#' 
#' @export
#' @docType methods
spikeCheck<-function(proj,auxName="all",sampleName=NULL,coverage=0,targets=NULL,...)
{  
  aux.files=alignments(proj)$aux
  if(nrow(proj@aux)==0){
    stop("The qProject object has no auxillary files\n",
         "This means probably no spike in experiments defined for the project!!")
  }
  if( is.null(sampleName) | (! sampleName %in% alignments(proj)$genome$SampleName ) ){
    stop("The 'sampleName' must be provided\n",
         "and it must be in the samples names of qProject object\n")
    
  }

  # filter QuasR object so it doesn't extract CpGs
  # from unnecessary samples
  proj2=proj
  proj2@alignments   =proj@alignments[proj@alignments$SampleName==sampleName,]
  proj2@reads        =proj@reads[proj@reads$SampleName==sampleName,]
  proj2@auxAlignments=proj@auxAlignments[,colnames(proj@auxAlignments)==sampleName,drop=FALSE]
  aux.files=alignments(proj2)$aux
  
  meths=list()
  for(i in 1:nrow(aux.files)){
    
    # if targets are provided
    if( is.list(targets) & (rownames(aux.files)[i] %in% names(targets)) ){
      
      #require(Rsamtools)
      chr.name=names( scanBamHeader(aux.files[i,])[[1]]$targets)# chrname of bam
      target=targets[[rownames(aux.files)[i]]] # get ranges for targeted regions
      
      # change the name in the GRanges object
      cur.name=as.character(unique(seqnames(target))[1])
      value=chr.name
      names(value)=cur.name
      target=renameSeqlevels(target,value)
      target=keepSeqlevels(target,chr.name)         
      meth <- qMeth(proj2,query=target, mode="CpG",reference=rownames(aux.files)[i] )
      
    }else{
      meth <- qMeth(proj2, mode="CpG",reference=rownames(aux.files)[i] )
    }
    
    Cs  =values(meth)[,paste(sampleName,"M",sep="_")]
    ToTs=values(meth)[,paste(sampleName,"T",sep="_")]
    Cs=Cs[ToTs>0]
    ToTs=ToTs[ToTs>0]
    
    id=rownames(aux.files)[i]
    meths[[id]][["Cs"]]=Cs
    meths[[id]][["ToTs"]]=ToTs
    
    # hist(Cs/ToTs,xlim=c(0,1))
  }
  
  if(auxName=="all"){
    
    #par(mfrow=c(2,ceiling(nrow(aux.files)/2) )  )
    
    result=list()
    for(i in 1:nrow(aux.files)){
      
      Cs=meths[[i]][["Cs"]]
      ToTs=meths[[i]][["ToTs"]]
      Cs  =Cs[ToTs>coverage]
      ToTs=ToTs[ToTs>coverage]
      vals=Cs/ToTs# meth values
      median1=median(vals)
      mean1  =mean(vals)
      hist(vals,xlim=c(0,1),xlab="Methylation Ratio",main=NULL,...)
      #hist(vals,xlim=c(0,1),xlab="Methylation Ratio",main=NULL)
      
      title(paste(rownames(aux.files)[i],":",sampleName,"Methylation Distribution"))
      mtext(side=3,paste("mean:",round(mean1,2),"median:",round(median1,2)) )
      
      result[[rownames(aux.files)[i]]]=sum(Cs)/sum(ToTs)
    }
     
  }else{
    Cs  =meths[[auxName]][["Cs"]]
    ToTs=meths[[auxName]][["ToTs"]]
    Cs  =Cs[ToTs>coverage]
    ToTs=ToTs[ToTs>coverage]
    vals=Cs/ToTs# meth values
    median1=median(vals)
    mean1  =mean(vals)
    hist(vals,xlim=c(0,1),xlab="Methylation Ratio",main=NULL,...)
    
    title(paste(auxName,":",sampleName,"Methylation Distribution"))
    mtext(side=3,paste("mean:",round(mean1,2),"median:",round(median1,2)) )
    
    result=list()
    result[[rownames(aux.files)[i]]]=sum(Cs)/sum(ToTs)
     
  }
  
  message("Ratio of methylated Cs to Total number of Cs\n\n")
  for(myname in names(result)){
    message(myname,": ",result[[myname]],"\n")
  };message("\n")
  return(result)
}
# ampliconView to get a matrix out of amplicons 



# calculate similarity of columns by various measures
# 
# Function takes input a matrix of 1s and 0s, and calculates similarity of the columns
# NA values are allowed but rows that have NA values will be ignored in pairwise comparisons.
# pearson for binary, phi coeff and Tanay's measure are the same for binary variables
# http://en.wikipedia.org/wiki/Phi_coefficient
csim<-function(x){
  
  ncy <- ncx <- ncol(x)
  if (ncx == 0) 
    stop("'x' is empty")
  rS <- matrix(0, nrow = ncx, ncol = ncy) # basic overlap similarity
  rT <- matrix(0, nrow = ncx, ncol = ncy) # Tanay correlation also Phi coeff
  rM <- matrix(0, nrow = ncx, ncol = ncy) # overlap similarity of bases at least one overlap
  m11<- matrix(0, nrow = ncx, ncol = ncy) # number of both CpGs are methylated
  m10<- matrix(0, nrow = ncx, ncol = ncy) # number of only 1st CpG methylated
  m00<- matrix(0, nrow = ncx, ncol = ncy) # number of reads both CpGs are unmethylated
  m01<- matrix(0, nrow = ncx, ncol = ncy) # number of only 2nd CpG methylated
  
  for (i in seq_len(ncx)) {
    for (j in seq_len(ncy)) {
      x2 <- x[, i]
      y2 <- x[, j]
      ok <- complete.cases(x2, y2)
      x2 <- x2[ok]
      y2 <- y2[ok]
      
      if (any(ok)){
        # my similarity
        rS[i, j] <- sum(x2==y2)/length(x2)
        
        # tanay coeff/phi
        n11=sum((y2+x2)==2)
        n10=sum( x2==1 & y2==0)
        n01=sum( x2==0 & y2==1)
        n00=sum((y2+x2)==0)
        N=length(x2)
        m1=mean(x2)
        m2=mean(y2)
        v1=m1*(1-m1)
        v2=m2*(1-m2)
        rT[i, j] <-  ( (n11/N) - (m1 * m2))/sqrt(v1*v2)
        
        # phi Coeff
        #row1=n11+n10
        #row2=n01+n00
        #col1=n11+n01
        #col2=n10+n00
        #rP[i,j] <- (n11*n00-n10*n01)/(sqrt(row1)*sqrt(row2)*sqrt(col1)*sqrt(col2))
        
        m11[i,j] <- n11 
        m01[i,j] <- n01 
        m10[i,j] <- n10 
        m00[i,j] <- n00 
        
        # methylation similarity
        # only look at methylated bases in at least one sample
        rM[i,j]<-n11/(n11+n10+n01)
      }else{
        rS[i,j]<-NA
        rT[i,j]<-NA
        rM[i,j]<-NA
        
        m11[i,j]<-NA 
        m01[i,j]<-NA 
        m10[i,j]<-NA 
        m00[i,j]<-NA 
      }
    }
  }
  
  rownames(rM) <- rownames(rS) <- rownames(rT) <- rownames(m11) <-rownames(m10) <- rownames(m01) <- rownames(m00) <- colnames(x)
  colnames(rM) <- colnames(rS) <- colnames(rT) <- colnames(m11) <-colnames(m10) <- colnames(m01) <- colnames(m00) <- colnames(x)
  return( list(similarity=rS,tanay=rT,msimilarity=rM,m11=m11,m10=m10,m01=m01,m00=m00) )
}


# function returns a matrix of locations and methylation statuses
# filter
# 
# 
# 
# a=readRDS("/work2/gschub/arnaud/amplicon_sequencing/primer_design/ampliconV2.gr.rds")
getCpGMatrix<-function(proj,range=a[2,],samp="amSeq",conv=NULL){
 # require(data.table)
  CpG=qMeth(proj, query=range,mode="CpG",reportLevel="alignment")
  #CpGm=tapply(CpG[[samp]]$meth ,list(CpG[[samp]]$aid ,CpG[[samp]]$Cid ),function(x) x )
  
  if(length(CpG[[samp]]$Cid)==0){return(NA)}
  
  # use data.table to get a 1,0 matrix of methylation profiles
  all.cids=unique(CpG[[samp]]$Cid) # get all possible CpG locations
  dt=data.table(meth=CpG[[samp]]$meth ,aid=CpG[[samp]]$aid ,cid=CpG[[samp]]$Cid) 
  
  # this function converts cids to columns
  myfun2<-function(x,all.cids){
    vec=rep(-1,length(all.cids))
    names(vec)=as.character(all.cids)
    b=as.list((vec))        
    b[ as.character(x$cid)]=as.double(x$meth)
    return(b)
  }
  dtm=dt[,myfun2(.SD,all.cids), by=aid]
  ronames=dtm$aid
  dtm[,aid:=NULL] # remove unwanted row
  CpGm=as.matrix(dtm) 
  CpGm[CpGm == -1]=NA # put NAs
  rownames(CpGm)=ronames

  
  #allCm=tapply(allC[[samp]]$meth,list(allC[[samp]]$aid,allC[[samp]]$Cid),function(x) x )

  #initialize values from conversion filtering
  f.conv.rate  =NULL
  conv.rate    =NULL
  max.numNonCG =NULL
  mean.numNonCG=NULL
  
  # remove reads that have fun!! (remove reads that are not converted)
  if( !is.null(conv) &  is.numeric(conv) ){
    
    conv=conv/100
    allC=qMeth(proj, query=range,mode="allC",reportLevel="alignment")
    
    
    CpGcor=c(CpG[[samp]]$Cid ,CpG[[samp]]$Cid+1) # get ids from CpG coordinates
    #rm(CpG);gc()
    
    # get non-CpG coordinates by removing stuff
    to.keep= (!allC[[samp]]$Cid %in% CpGcor)
    allC[[samp]]$meth  = allC[[samp]]$meth[to.keep ]
    allC[[samp]]$aid   = allC[[samp]]$aid[to.keep ]
    allC[[samp]]$strand= allC[[samp]]$strand[ to.keep ]
    allC[[samp]]$Cid   = allC[[samp]]$Cid[ to.keep ]
    
    
    # calculate conversion rate
    dt=data.table(meth=allC[[samp]]$meth,aln=allC[[samp]]$aid )
    rm(allC);gc()
    dt=dt[,list(conv=mean(meth),msum=sum(meth),cnum=length(meth)),by=aln] # get number of non-converted per read/amplicon
    c.rate=(1-dt$conv) # calculate conversion rate
    good.ids=dt$aln[c.rate>=conv] # get good alignment ids, that have conv rat
    
    # if no read has passed the conversion threshold
    if(length(good.ids)==0){
      CpGm=NA
    }else{
      CpGm=CpGm[rownames(CpGm) %in% good.ids,,drop=FALSE] # subset the CpG table
    }
    
    f.conv.rate=1-sum(dt[c.rate>=conv,]$msum)/sum(dt[c.rate>=conv,]$cnum)
    conv.rate  =1-sum(dt$msum)/sum(dt$cnum)
    max.numNonCG  =max(dt$cnum)
    mean.numNonCG  =mean(dt$cnum)
  }
  
  if(class(CpGm) != "matrix"){
    CpGm=as.matrix(rbind(CpGm))}
  
  # remove all NA columns from CpG matrix
  CpGm=CpGm[,!apply(CpGm,2, function(x) all(is.na(x)) ),drop=FALSE]
  if(dim(CpGm)[2]==0){CpGm=matrix(NA)}
  
  list(CpGm=CpGm, f.conv.rate=f.conv.rate, 
       conv.rate=conv.rate,
       mean.numNonCG=mean.numNonCG,
       max.numNonCG=max.numNonCG)
}
  
  


#' Makes an AmpliconViews object
#' 
#' Makes an AmpliconViews object from qProject object from QuasR.
#' 
#' @param proj qProject object from QuasR
#' @param range GRanges object for the amplicon locations
#' @param tag a character vector, containing a additional tags that describe the
#'        amplicons. The vector length should be same as \code{range} object.
#' @param sampleNames sample name as character
#' @param conv minimum conversion efficiency a numeric value between 0 and 100,
#'        reads below this conversion efficiency will be discarded
#' @param exp.var percentage of reads explained by meta-methylation patterns. This
#'        helps to infer top meta-methylation patterns. If set to 90, meta-patterns
#'        are ranked by their cluster sizes, and patterns that  attributes to 90 percent of the data
#'        are returned.
#' @param noise.tol percentage of variation to be considered as noise. This usually
#'        corresponds to PCA components that explain small percentage
#'        of the data to be removed. If equals to 'auto', the noise is automatically
#'        removed by removing components that explain low amount of variation.
#' 
#' @param example.size size of the example matrix returned. example matrix contains
#'        example reads from the experiment.
#' @param call.matrix \code{logical}, default FALSE, if TRUE methylation call matrix
#'        that contains methylation calls per fragment are returned. This might result
#'        in a large AmpliconViews object, should be set to FALSE normally.
#'
#'                     
#' @return returns an  \link{AmpliconViews-class}  object      
#'  
#' @examples
#'  # a=GRanges(seqnames=,range=IRanges() )
#'  # av=ampliconView(proj,range=a[2,],samp="amSeq",conv=80,exp.var=80,example.size=100)
#'  
#' @importFrom data.table data.table  
#' @importFrom QuasR qMeth
#' @importClassesFrom QuasR qProject
#' @export
#' @docType methods
AmpliconViews<-function(proj,range,tag=NULL,sampleNames,conv=NULL,exp.var=80,
                        noise.tol='auto',
                       example.size=100,call.matrix=FALSE,verbose=FALSE){
  if( (!is.null(tag))  & (length(tag) != length(range) )){
    stop("'tag' length should be equal to 'range' object ")
  }
  
  if(!is.null(tag)){
    tag=as.character(tag)
  }
  
  if(any(! sampleNames %in% proj@alignments$SampleName )){
    stop("\n",sampleNames[! sampleNames %in% proj@alignments$SampleName],
        "are not in aligned sample names in qProject object\n",
         "Please provide correct sample names")
  }
  
  av.obj=list() # output list with data
  for(samp in sampleNames){
  
    # filter QuasR object so it doesn't extract CpGs
    # from unnecessary sampleNames
    proj2=proj
    proj2@alignments   =proj@alignments[proj@alignments$SampleName==samp,]
    proj2@reads        =proj@reads[proj@reads$SampleName==samp,]
    proj2@auxAlignments=proj@auxAlignments[,colnames(proj@auxAlignments)==samp,drop=FALSE]
    
    
    for(i in 1:length(range)){
      
      # region info
      chr=as.character(seqnames(range[i,])[1])
      start=start(range[i,])[1]
      end=end(range[i,])[1]
      
      
      if(verbose){
        cat("accessing and analyzing region:",
            paste(chr,start ,end ,sep="."),"in sample",
            samp,"\n")
      }
      cpg=getCpGMatrix(proj2,range[i,],samp,conv)# get CpG meth matrix
      
      # if there is no overlap return NA
      if(all(is.na(cpg))){
        obj=list(meths=NA,covs=NA,smat=NA,mat=NA,
                 chr=chr,start=start,end=end,mf=NA,f.conv.rate=NA, 
                 conv.rate=NA,
                 mean.numNonCG=NA,
                 max.numNonCG=NA,tag=tag[i],col.sim=NA)
        av.obj[[samp]][[paste(chr,start,end,sep="_")]]=obj
        next
      }
    
      # if all the reads are lost due to conversion efficiency problems
      if( all(dim(cpg$CpGm)==c(1,1)) & is.na(cpg$CpGm[1,1]) ){
        obj=list(meths=NA,covs=NA,smat=NA,mat=NA,
                 chr=chr,start=start,end=end,mf=NA,f.conv.rate=cpg$f.conv.rate, 
                 conv.rate=cpg$conv.rate,
                 mean.numNonCG=cpg$mean.numNonCG,
                 max.numNonCG=cpg$max.numNonCG,tag=tag[i],col.sim=NA)
        av.obj[[samp]][[paste(chr,start,end,sep="_")]]=obj
        next
      }
      
      #calculate colum similarities
      colSim=csim(cpg$CpGm)
      
      #x=svd.bi(CpGm,exp.var=exp.var,nsmp=example.size) # get the SVD
      x=metaProfile(cpg$CpGm,pca.var='auto',noise.tol=noise.tol,meta.exp=exp.var,
                    example.size=example.size,call.matrix=call.matrix)
      
      if(class(cpg$CpGm) != "matrix"){
        cpg$CpGm=as.matrix(rbind(cpg$CpGm))}
      
      covs=colSums(!is.na(cpg$CpGm))
      meths=colMeans(cpg$CpGm,na.rm=TRUE)
      
      chr=as.character(seqnames(range[i,])[1])
      start=start(range[i,])[1]
      end=end(range[i,])[1]
      mf=x[names(x) != "o.mat"] # get the meta profiles from matrix factorization
      
      # if there is no example matrix to return
      if(is.null(x$o.mat)){
        smat=NA
      }else{
        smat=do.call("rbind",x$o.mat)
      }
      
      #
      if(call.matrix){cmat=cpg$CpGm[x$mat.order,]}else{cmat=NA}
      obj=list(meths=meths,covs=covs,smat=smat,mat=cmat,
              chr=chr,start=start,end=end,mf=mf,
              f.conv.rate=cpg$f.conv.rate, 
              conv.rate=cpg$conv.rate,
              mean.numNonCG=cpg$mean.numNonCG,
              max.numNonCG=cpg$max.numNonCG,tag=tag[i],
              col.sim=colSim)
      av.obj[[samp]][[paste(chr,start,end,sep="_")]]=obj
    }
  }
  
  new("AmpliconViews",sampleNames=sampleNames,amplicons=range ,data=av.obj)
}


plotEigens<-function(eigens){
  
pchs=rep(1,ncol(eigens))  
pchs[eigens[1,]==1]=19
plot(colnames(CpGm),rep(1,ncol(eigens)),ylim=c(nrow(eigens)+1,0),type="b",
     pch=pchs,yaxp=c(1, nrow(eigens), 1),ylab="Binarized Eigen Vectors",
     xlab="Chr coordinates")  
for(i in 2:nrow(eigens)){
  pchs=rep(1,ncol(eigens))  
  pchs[eigens[i,]==1]=19
  lines(colnames(CpGm),rep(i,ncol(eigens)),type="b",
       pch=pchs)   
  
}

}

# remove the noise from the matrix
# @param A matrix
# @param var is a numeric value or character string 'auto'
reduceMat <- function(A,var) {
  #Calculates the SVD
  #Approximate each result of SVD with the given dimension  
  #Create the new approximated matrix
  #return( u%*%d%*%t(v))
  pca=prcomp(A,center = FALSE, scale. = FALSE)
  
  #image(t(pca$x %*% t(pca$rotation)) )
  vars=pca$sdev**2/sum(pca$sdev**2)
  
  if(is.numeric(var)){
    dim.end=which(100*cumsum(vars)>=var)[1]
  }else if(var=="auto"){
    #dim.end=findElbow(vars)-1
    dim.end=findElbow(vars)
  }
  pca$rotation[,-c(1:dim.end)]=0
  return(pca$x %*% t(pca$rotation))
}


# paste the columns of a matrix to a string
# from plotrix package
pasteCols<-function(x,sep="") {
  pastestring<-paste("list(",paste("x","[",1:dim(x)[1],",]",
                                   sep="",collapse=","),")",sep="")
  return(do.call(paste,c(eval(parse(text = pastestring)),sep=sep)))
}

# get meta-profiles from a matrix of methylation statuses
# 
# 
# @param meta.exp numeric value between 0 and 100, threshold for what % of 
# @param pca.var numeric value between 0 and 100
# @param noise.tol numeric value between 0 and 100, character string 'auto' or NULL
# @param example.size numeric value of example set size 
# @param call.matrix if TRUE returns the PCA order of the input matrix
# 
# @return a list of with following slots 'metas':meta profiles, 'weight': numeric vector how many rows each
#         meta profile represent, 
#         'o.mat': list of matrices containing example rows for each meta-profile,
#         'imp.metas': most important meta-profiles, 
#          'meta.exp': numeric vector of how many rows each meta profile represent, 
#         'pca.dim': numeric, number of important PCA components 
metaProfile<-function(mat,meta.exp=80,pca.var="auto",noise.tol=NULL,example.size=100,call.matrix=FALSE){
  
  mat2=mat
  if(nrow(mat2) == 1){
    return(list(metas=NA,meta.weight=NA,o.mat=list(as.matrix(rbind(mat2))),
         imp.metas=NA,
         meta.exp=NA,pca.dim=NA))
  }
  
  if(ncol(mat2) == 1){
    if(call.matrix){
     return(list(metas=NA,meta.weight=NA,o.mat=list(as.matrix(rbind(mat2))),
                imp.metas=NA,
                meta.exp=NA,pca.dim=NA,pca.vars=NA,mat.order=mat2))
    }else{
      return(list(metas=NA,meta.weight=NA,o.mat=list(as.matrix(rbind(mat2))),
                  imp.metas=NA,
                  meta.exp=NA,pca.dim=NA,pca.vars=NA,mat.order=NA))
    }
  }
  
  # remove columns with only NAs
  mat2=mat2[,!apply(mat2,2, function(x) all(is.na(x)) ),drop=FALSE]
  
  # impute missing values
  means=round(colMeans(mat2,na.rm=TRUE))
  for(i in 1:ncol(mat2)){
    mat2[,i][is.na(mat2[,i])]=means[i]
  }
  
  
  mat2[mat2==0]=-1
  
  # reduce the noise
  if(is.numeric(noise.tol)){
    mat2=reduceMat(mat2,100-noise.tol)
  }else if (is.null(noise.tol)){
    mat2=mat2
  }else if(noise.tol=="auto"){
    mat2=reduceMat(mat2,"auto")  
  }else{
    stop("'noise.tol' can only be a numeric value, or string 'auto' or NULL")
  }
  
  # run SVD/PCA
  pca=prcomp(mat2,center = FALSE, scale. = FALSE)
  vars=pca$sdev**2/sum(pca$sdev**2)
  
  if(pca.var=="auto"){
    #dim.end=findElbow(vars)-1
    dim.end=findElbow(vars)
  }else if(is.numeric(pca.var)){
    dim.end=which(100*cumsum(vars)>=pca.var)[1]
  }else{
    stop("'pca.var' can only be a numeric value, or string 'auto'")}
  
  # discretize the PCA matrix
  subs=pca$x
  subs[subs <= 0.1 & subs >= -0.1]=0
  subs[subs > 0.1 ]=1
  subs[subs < -0.1 ]=-1
  if(call.matrix){mat.order=do.call("order",-as.data.frame(subs[,1:dim.end]))}
  else{
    mat.order=NA
  }
  # get most frequent pairings of discrete components
  # this is similar to clustering
  ids= pasteCols(t(subs[,1:dim.end]),sep=" ")
  #ids=paste(subs[,1],subs[,2],subs[,3])
  t.ids=sort(table(ids),decreasing = TRUE)
  t.ids=as.list(t.ids)
  mat.ind=lapply(names(t.ids),function(x) which(ids==x) )
  names(mat.ind)=names(t.ids)
  
  #get meta profiles by looking at clusters
  metas=c()
  for(i in 1:length(t.ids)){
    mprof=mat2[ids==names(t.ids)[i],]
    if(length(nrow(mprof)) ) {
      metas=rbind(metas,colMeans(mprof) )
    }else{
      metas=rbind(metas,mprof) 
      
    }
  }
  
  # combine patterns if same
  # similarities between meta profiles
  # subtract meta-profiles from the meta profile matrix
  s.metal=split(sign(metas),1:nrow(metas))
  sim=lapply(s.metal,function(x,y){which(rowSums(abs(sweep(y,2,x,FUN= "-")))==0)
  },sign(metas))

  #sim=apply(sign(metas),1,function(x,y){
  #  which(rowSums(abs(sweep(y,2,x,FUN= "-")))==0)
  #},sign(metas))
  
  ind=unique(sim[which(lapply(sim,length)>1)])
  
  b.metas=sign(metas)
  
  if(length(ind)){ #if there are duplicated meta-profiles merge them
    new.inds=list()
    new.t.ids=list()
    new.metas=c()
    for(i in 1:length(ind)){
      new.ind=unlist(lapply(names(t.ids)[ ind[[i]]],function(x) which(ids==x)))
      
      new.id=paste0(names(t.ids)[ ind[[i]]],collapse="|")
      new.t.ids[[new.id]]=length(new.ind)
      new.inds[[new.id]]=new.ind
      
      new.metas=rbind(new.metas,colMeans(b.metas[ind[[i]],])  )
    }
    
    # remove from old t.ids
    t.ids=as.list(t.ids[-unlist(ind)])
    b.metas=b.metas[-unlist(ind),]
    mat.ind=mat.ind[-unlist(ind)]
    
    t.ids=c(t.ids,new.t.ids)
    b.metas=rbind(b.metas,new.metas)
    mat.ind=c(mat.ind,new.inds)
  }
  
  # sort data structures based on cluster sizes
  b.metas=b.metas[order(-unlist(t.ids)),]
  mat.ind=mat.ind[order(-unlist(t.ids))]
  t.ids=sort(unlist(t.ids),decreasing=TRUE)
  
  mvars=100*(t.ids)/sum(t.ids)
  best.meta=1:which(cumsum(mvars) >=meta.exp)[1]
  
  if(nrow(mat)>example.size){
    o.mat=list()
    for(i in best.meta){
      o.mat[[i]]=mat[sample(mat.ind[[i]],example.size*mvars[[i]]/100),]
    }
  }else{
    o.mat=NULL
  }
  #o.mat[[i+5]]=mat[sample(mat.ind[[i]],example.size*mvars[[i]]/100),]
  if(class(b.metas) != "matrix"){
    b.metas=t(as.matrix(b.metas))
  }
  list(metas=b.metas,meta.weight=(t.ids),o.mat=o.mat,imp.metas=b.metas[best.meta,,drop = FALSE],
       meta.exp=mvars,pca.dim=dim.end,pca.vars=vars,mat.order=mat.order)
}


# finds elbow in explained variation by PCA  plots
findElbow<-function(vars){
  
  # if there are too few points return 1
  if( length(vars) <= 2){
    return(1)
  }
  
  nPoints = length(vars)
  allCoord <- cbind(1:nPoints,vars)              
  
  # pull out first point
  firstPoint = allCoord[1,];
  
  # get vector between first and last point - this is the line
  lineVec = allCoord[nrow(allCoord),] - firstPoint;
  
  # normalize the line vector
  lineVecN = lineVec / sqrt(sum(lineVec**2));
  
  # find the distance from each point to the line:
  # vector between all points and first point
  vecFromFirst =sweep(allCoord,2,firstPoint,FUN= "-")
  scalarProduct = vecFromFirst %*% lineVecN
  vecFromFirstParallel = do.call("rbind",lapply(scalarProduct,function(x) x*lineVecN))
  vecToLine = vecFromFirst - vecFromFirstParallel
  distToLine = sqrt(rowSums(vecToLine**2))
  distToLine = (distToLine-min(distToLine))/(max(distToLine)-min(distToLine))
  dists=abs(distToLine-max(distToLine))
  
  # this is to get the last component if there is no
  # elbow, where everything is on a straight line
  if( all(is.na(dists)) ){
    return(length(dists))
  }
  
  # this bit is to remove componets that are close
  # to elbow. If their distance to the line
  # is reasonably close to max dist and if they explain more
  cutoff=0.05
  if( any(dists<cutoff & dists != 0) ){
    
    my.order=order(dists)
    return(my.order[my.order %in%  which(dists<cutoff & dists != 0)][1] )
  }
  
  return(which.max(distToLine))
}


#' convert ampliconView objects to a table of methylation ratio
#' 
#' @param x an AmpliconViews object 
#' @param per.region logical, if FALSE base-pair methylation ratio (def:FALSE)
#' @param asGRanges logical, if TRUE object returned is a GRanges object (def:TRUE)
#' @param dup.resolve logical, if TRUE , and if per.region=FALSE, the duplicated
#'        CpGs will be resolved by taking the one with highest coverage. Since
#'        amplicon designs can overlap, there might be CpGs that are covered by 
#'        different amplicon designs at the same time, this will result in duplicated
#'        CpGs in the resulting object. (def:TRUE)
#' @param simple.tags logical. If TRUE, tags associated with amplicons will be output as a single 
#'        column. If FALSE, a tag column will be output for each experiment.
#' @param coverage.th default to 0. If this is set to a value larger than zero, the bases/regions that have
#'        have coverage below this value will have NA methylation and NA coverage.
#' 
#' @return data frame or GRanges object with locations of CpGs and methylation ratio and coverage
#' 
#' @usage methRatio(x,per.region=FALSE,asGRanges=TRUE,dup.resolve=TRUE,simple.tags=TRUE,coverage.th=0)
#' 
#' @examples 
#'  # methRatio(a.list,per.region=FALSE,asGRanges=TRUE,dup.resolve=TRUE)
#'  
#' @importFrom GenomicRanges GRanges
#' @importClassesFrom GenomicRanges GRanges
#'                   
#' @export
#' @docType methods
methRatio<-function(x,per.region=FALSE,asGRanges=TRUE,dup.resolve=TRUE,simple.tags=TRUE,coverage.th=0){
  #require(GenomicRanges)
  
  a.list=ampData(x)
  a.list=a.list[!is.na(a.list)] # remove NA amplicons
  
  # decide if a list contains multiple sampleNames
  if( names(a.list[[1]])[1] == "meths" ){
    
    a.list=list(sample=a.list)
  }
  # try to decide if the list is a list
  # from ampliconView function
  if(names(a.list[[1]][[1]])[1] != "meths"  )
  {stop("a.list provided doesn't look like a product of ampliconView function","\n")
  }
  
  
  df.list=list() # data.frame list, will be merged at the end for each sample
  
    
  
    for(i in 1:length(a.list))
    {   
      samp.name=names(a.list)[i]

      if(per.region){
    
        # if the tag is not null
        if(! is.null(a.list[[i]][[1]]$tag)){
          my.list=lapply(a.list[[i]], function(x)
                                                if(all(!is.na(x[["meths"]]))){
                                                data.frame(chr=x[["chr"]],start=x[["start"]],
                                                end=x[["end"]],
                                                #id=paste(x[["chr"]],x[["start"]],x[["end"]],sep="_"),
                                                coverage=mean(x[["covs"]],na.rm=TRUE), 
                                                methRatio=mean(x[["meths"]],na.rm=TRUE),
                                                tag=x$tag)} )
          names(my.list)=NULL
          res=do.call("rbind",my.list)
          colnames(res)[4:6]=paste(samp.name,colnames(res)[4:6],sep=".")
          df.list[[i]]=unique(res)
        }else{
          my.list=lapply(a.list[[i]], function(x)
            if(all(!is.na(x[["meths"]]))){
              data.frame(chr=x[["chr"]],start=x[["start"]],
                         end=x[["end"]],
                         #id=paste(x[["chr"]],x[["start"]],x[["end"]],sep="_"),
                         coverage=mean(x[["covs"]],na.rm=TRUE), 
                         methRatio=mean(x[["meths"]],na.rm=TRUE) )} )
          names(my.list)=NULL
          res=do.call("rbind",my.list)
          colnames(res)[4:5]=paste(samp.name,colnames(res)[4:5],sep=".")
          df.list[[i]]=unique(res)          
        }
            
      }else{
        
        if(!is.null(a.list[[i]][[1]]$tag) ){
          my.list=lapply(a.list[[i]], 
                       function(x){  # cat(paste(x[["chr"]],x[["start"]],x[["end"]],sep="_"),"\n")
                                     if(all(!is.na(x[["meths"]]))){  
                                      data.frame(chr=rep(x[["chr"]],length(x[["meths"]])),
                                                start=as.numeric(names(x[["meths"]])),
                                                end=as.numeric(names(x[["meths"]])),
                                                #id=paste(x[["chr"]],x[["start"]],x[["end"]],sep="_"),
                                                coverage=(x[["covs"]]), 
                                                methRatio=(x[["meths"]]),
                                                tag=x$tag  )}} )
          names(my.list)=NULL
          res=unique(do.call("rbind",my.list))
          if(dup.resolve){
            res=res[order(-res[,4]),]# order by coverage
            res=res[ !duplicated(res[,1:3]),] # remove the duplicated with lower cov
          }
          
          colnames(res)[4:6]=paste(samp.name,colnames(res)[4:6],sep=".")
          df.list[[i]]=unique(res)
        }else{
          my.list=lapply(a.list[[i]], 
                         function(x){  # cat(paste(x[["chr"]],x[["start"]],x[["end"]],sep="_"),"\n")
                           if(all(!is.na(x[["meths"]]))){  
                             data.frame(chr=rep(x[["chr"]],length(x[["meths"]])),
                                        start=as.numeric(names(x[["meths"]])),
                                        end=as.numeric(names(x[["meths"]])),
                                        #id=paste(x[["chr"]],x[["start"]],x[["end"]],sep="_"),
                                        coverage=(x[["covs"]]), 
                                        methRatio=(x[["meths"]]) )}} )
          names(my.list)=NULL
          res=unique(do.call("rbind",my.list))
          if(dup.resolve){
            res=res[order(-res[,4]),]# order by coverage
            res=res[ !duplicated(res[,1:3]),] # remove the duplicated with lower cov
          }
          colnames(res)[4:5]=paste(samp.name,colnames(res)[4:5],sep=".")
          df.list[[i]]=unique(res)
        } 
    }
     
    }
  

  
  if(length(df.list)>1){
    data <- Reduce(function(x, y) merge(x, y, all=T, 
          by=c("chr", "start", "end")), df.list, accumulate=F)
  }else{
    data <- df.list[[1]]
  }
  
  if(simple.tags & (!is.null(a.list[[i]][[1]]$tag) )  ){
    ti=grep(".tag" ,colnames(data))
    tvec=apply(data[ti],1, function(x) unique(x[!is.na(x)]) )
    if(class(tvec) != "list"   ){
      data=data[,-ti]
      data=cbind(data,tag=tvec)
    }
  }
  
  if(coverage.th>0){
    message("removing regions/bases with coverage below 'coverage.th'")
    cov.ind=grep("coverage",colnames(data))
    data[,cov.ind][data[,cov.ind]<coverage.th]=NA
    data[,cov.ind+1][is.na(data[,cov.ind])]=NA
  }
  # if wanted return GRanges
  if(asGRanges){
    res=GRanges(seqnames=as.character(data[,1]),
                ranges=IRanges(start=data[,2],end=data[,3]))
    
    values(res)=DataFrame(data[,-(1:3)]) 
    return(res)
  }
  return(data)
}




# _____________ UNUSED ____________________________
# _________________________________________________

# plot(dist(CpGm[1:100,]))

#  mat=CpGm[1:3000,]
#  my.dist=dist(mat)
#  clusters=kmeans(my.dist,centers=2)$cluster
#  colMeans(mat[clusters==1,],na.rm=TRUE)
#  colMeans(mat[clusters==2,],na.rm=TRUE)
#  plot(colMeans(mat[clusters==2,],na.rm=TRUE),type="b",ylim=c(0,1))
#  lines(colMeans(mat[clusters==1,],na.rm=TRUE),type="b" )
# smoothScatter( cmdscale(dist(CpGm[1:3000,])) )
# princomp(CpGm[1:3000,])


# do SVD on matrix of CpG methylation statuses
# 
# Takes a matrix from getCpGMatrix
# 
# @param CpGm CpG methylation matrix from amplicons, see getCpGMatrix()
# @param exp.var what portion of variation should be explained by components
# @param nsmp number of sample, to be plotted as example heatmap
# 
# @return a list including 
svd.bi<-function(CpGm,exp.var=0.8,nsmp=100)
{
  #exp.var=0.9 # variation to be explained
  #nsmp=200 # total number of examples to be plotted
  mat=CpGm
  
  # impute missing values
  means=round(colMeans(mat,na.rm=TRUE))
  for(i in 1:ncol(mat)){
    mat[,i][is.na(mat[,i])]=means[i]
  }
  mat[mat==0]=-1 # make 0s to -1 for svd
  
  pca=prcomp(mat,center = FALSE, scale. = FALSE) # compute SVD
  #screeplot(pca)
  #plot(pca$rotation[,1],type="b")
  #plot(pca$rotation[,2],type="b")
  #plot(pca$rotation[,3],type="b")
  cvars=cumsum(pca$sdev**2)/sum(pca$sdev**2) # get cumulative variance
  vars=(pca$sdev**2)/sum(pca$sdev**2) # get variance explained by each componnet
  n.comp= which(cvars>exp.var)[1] # get most explanatory components
  
  
  
  eigens=matrix(0,nrow=n.comp,ncol=ncol(mat))
  colnames(eigens)=colnames(mat)# create empty binary eigen values
  o.mat=list()# output matrix
  mat2=CpGm # initialize matrices
  #mat2[mat2==-1]=0 # convert -1 to 0s
  X=pca$x # component matrix
  for(comp in 1:n.comp){
    # create binary eigen vectors
    eigens[comp,]=sign(mean(pca$x[,comp])*pca$rotation[,comp])    
    
    # get the top examples for componets
    my.order=order(-1*sign(mean(pca$x[,comp]))*X[,comp])
    my.order=my.order[1:round(length(my.order)/4)] # get top quantile
    
    # sample from top quantile
    my.order=my.order[sample(1:length(my.order),round(nsmp*vars[comp]))]
    
    # put the subsample in the output matrix
    o.mat[[comp]]=rbind(mat2[my.order,],rep(NA,ncol(mat)) )
    
    #update mat2 and X (components) so the ones that are extracted won't be 
    # extracted again
    mat2=mat2[-my.order,]
    X   =X[-my.order,]
  }
  eigens[eigens==-1]=0
  
  #
  # get the top examples for rest of the componets
  #
  ##my.order=order(X[,(comp+1)])
  ##my.order=my.order[1:round(length(my.order)/2)] # get above median
  # sample from top quantile
  ##my.order=my.order[sample(1:length(my.order),round(nsmp*(1-cvars[comp])) )]
  # put the subsample in the output matrix
  ##o.mat[[comp+1]]=mat2[my.order,]
  colnames(o.mat)=colnames(mat)
  list(out.mat=o.mat,eigens=eigens,var.exp=vars)
}

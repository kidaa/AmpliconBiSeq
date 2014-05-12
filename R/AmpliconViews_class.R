# classes and accessors for the classes

#' An S4 class for amplicon bisulfite sequencing
#'
#' This class is designed to hold methylat ion statistics and locations for 
#' targeted amplicons from a amplicon bisulfite sequencing experiment .
#' \code{\link[AmpliconBiSeq]{AmpliconViews}} function returns an object of \code{AmpliconViews} class.
#'          
#' @section Slots:
#' \describe{
#'    \item{\code{sampleNames}}{Names of samples in a vector}
#'    \item{\code{amplicons}}{GRanges object for locations of the amplicons}
#'    \item{\code{data}}{list structure holding the information on amplicons}
#'
#' }
#' 
#' @section Details:
#' \code{AmpliconViews} class is desingned to contain amplicon information such
#' as methylation, coverage, meta-patterns and similarities between CpG profiles.
#' 
#' 
#' @section Constructor:
#' see \code{\link{AmpliconViews}}
#' 
#' @section Subsetting:
#'  an AmpliconViews object containing multiple amplicons and samples can be 
#'  subsetted using \code{\link[AmpliconBiSeq]{getAmplicon}} function.
#'  
#' 
#' 
#' @section Accessors: 
#' The following functions provides access to data slots of methylDiff:
#' \code{\link[AmpliconBiSeq]{getSamples}},\code{\link[AmpliconBiSeq]{getAmpliconRanges}},
#' \code{\link[AmpliconBiSeq]{getAmplicon}}, \code{\link[AmpliconBiSeq]{getMetaMethylation}}
#' 
#' @examples
#' 
#' library(GenomicRanges)
#' 
#' 
#' @name AmpliconViews-class
#' @rdname AmpliconViews-class
#' @exportClass AmpliconViews
#' @docType class
# x=new("AmpliconViews",sampleNames=names(sample.list),amplicons=unique(amps) ,data=sample.list )
#setClass("AmpliconViews",representation(
#  sampleNames = "character", amplicons = "GRanges"),contains="list")
setClass("AmpliconViews",representation(
  sampleNames = "character", amplicons = "GRanges",data="list") )

# new("AmpliconViews")

##############################################################################
## ACESSOR FUNCTIONS FOR AmpliconViews OBJECT
##############################################################################

#' get sample Names from AmpliconViews object
#' 
#' The function returns the sample Names stored in any of the \code{\link{AmpliconViews}}
#' 
#' @param x an \code{\link{AmpliconViews}}
#' @usage getSampleNames(x)
#' @examples data(ampViewEx);getSampleNames(ampViewEx)
#' 
#' 
#' 
#' @return character vector for sample Names
#' @export
#' @docType methods
#' @rdname getSampleNames
setGeneric("getSampleNames",function(x) standardGeneric("getSampleNames") )

#' @rdname getSampleNames
#' @aliases getSampleNames,AmpliconViews 
setMethod("getSampleNames", signature(x="AmpliconViews"),
          function(x) x@sampleNames
          )


#' get locations of amplicons from AmpliconViews object
#' 
#' The function returns amplicon locations of  stored in any of the 
#' \code{\link{AmpliconViews}}
#' 
#' @param x an \code{\link{AmpliconViews}}
#' 
#' @usage getAmpliconRanges(x)
#' 
#' @return character vector for sample Names
#' @examples data(ampViewEx);getAmpliconRanges(ampViewEx)
#' 
#' @export
#' @docType methods
#' @rdname getAmpliconRanges
setGeneric("getAmpliconRanges",function(x) standardGeneric("getAmpliconRanges") )

#' @rdname getAmpliconRanges
#' @aliases getAmpliconRanges,AmpliconViews 
setMethod("getAmpliconRanges",signature(x="AmpliconViews"),
          function(x) x@amplicons
          )


#' get data from AmpliconViews object
#' 
#' The function returns data stored in any of the 
#' \code{\link{AmpliconViews}}
#' 
#' @param x an \code{\link{AmpliconViews}}
#' @usage ampData(x)
#' @examples data(ampViewEx);ampData(ampViewEx)
#' 
#' @return list of amplicon data
#' @export
#' @docType methods
#' @rdname ampData
setGeneric("ampData",function(x) standardGeneric("ampData") )

#' @rdname ampData
#' @aliases ampData,AmpliconViews 
setMethod("ampData",signature(x="AmpliconViews"),
          function(x) x@data
)


#' subset AmpliconViews object
#' 
#' The function returns data stored in any of the 
#' \code{\link{AmpliconViews}}
#' 
#' @param x an \code{\link{AmpliconViews}} object
#' @param sampleNames character vector of sample names
#' @param ampliconNames character vector of amplicon names
#' @usage getAmplicon(x,sampleNames=NULL,ampliconNames=NULL)
#' @examples data(ampViewEx)
#'           getAmplicon(ampViewEx,"mock4","chr18_69674375_69674775")
#' 
#' @return \code{\link{AmpliconViews}} object
#' @export
#' @docType methods
#' @rdname getAmplicon
setGeneric("getAmplicon",function(x,sampleNames=NULL,ampliconNames=NULL) standardGeneric("getAmplicon") )


#' @rdname getAmplicon
#' @aliases getAmplicon,AmpliconViews,character,character 
setMethod("getAmplicon", signature(x="AmpliconViews",sampleNames="character",ampliconNames="character"),
          function(x,sampleNames,ampliconNames){
            if(is.null(sampleNames) & !is.null(ampliconNames)){
              stop("'ampliconNames' must be provided, if 'sampleNames' is provived")
            }
            #if(length(ampliconNames)>1){
            #  stop("'ampliconNames' must be character class of length 1, can not be a vector or a list")
            #}
            
            sampleNames=sampleNames[sampleNames %in% getSampleNames(x) ]
            if(is.null(ampliconNames)){
              
              new("AmpliconViews",sampleNames=sampleNames,amplicons=getAmpliconRanges(x)
                    ,data=ampData(x)[sampleNames] )
            }
            
            my.ranges=getAmpliconRanges(x)
            my.ranges=my.ranges[which(getAmpliconNames(x) %in% ampliconNames),]
            new("AmpliconViews",sampleNames=sampleNames,amplicons=my.ranges,
                data=lapply(ampData(x)[sampleNames],function(x,y) 
                           x[names(x) %in% y], y=ampliconNames )
                )
            
}
)

#' get amplicon Names from AmpliconView object
#' 
#' @param x an \code{\link{AmpliconViews}} object
#' @usage getAmpliconNames(x)
#' @return a character vector
#' 
#' @export
#' @docType methods
#' @rdname getAmpliconNames
setGeneric("getAmpliconNames",function(x) standardGeneric("getAmpliconNames") )

#' @rdname getAmpliconNames
#' @aliases getAmpliconNames,AmpliconViews 
setMethod("getAmpliconNames", signature(x="AmpliconViews"),
          function(x){
            names(x@data[[1]])
          }
)


#' get meta-methylation info from AmpliconView object
#' 
#' @param x an \code{\link{AmpliconViews}} object
#' @param sampleName name of the sample
#' @param ampliconName name of the amplicon
#' @usage getMetaMethylation(x,sampleName,ampliconName)
#' 
#' @return a list containing meta-methylation information: patterns, explained percentage of data
#'         SVD(PCA) results
#' @export
#' @docType methods
#' @rdname getMetaMethylation
setGeneric("getMetaMethylation",function(x,sampleName ,ampliconName ) standardGeneric("getMetaMethylation") )

#' @rdname getMetaMethylation
#' @aliases getMetaMethylation,AmpliconViews,character,missing 
setMethod( "getMetaMethylation", signature(x="AmpliconViews",sampleName="character",ampliconName="character"),
          function(x,sampleName,ampliconName){
            
            if( length(ampliconName) > 1 | length(sampleName) > 1){
              stop("'ampliconName' and 'sampleName' must be character class of length 1, can not be a vector or a list")
            }
            ampData(x)[[sampleName]][[ampliconName]]$mf
})

#' get example methylation matrix
#' the methylation matrix of an ampliconBiSeq experiment
#' 
#' @param x an \code{\link{AmpliconViews}} object
#' @param sampleName name of the sample
#' @param ampliconName name of the amplicon
#' @usage getExampleMethMat(x,sampleName,ampliconName)
#' 
#' @export
#' @docType methods
#' @rdname getExampleMethMat
setGeneric("getExampleMethMat",function(x,sampleName,ampliconName) standardGeneric("getExampleMethMat") )
#' @rdname getExampleMethMat
#' @aliases getExampleMethMat,AmpliconViews,character,character 
setMethod("getExampleMethMat", signature(x="AmpliconViews",sampleName="character",
                                         ampliconName="character"),
          function(x,sampleName,ampliconName){
            if(length(ampliconName)>1 | length(sampleName)>1){
              stop("\n'ampliconName' and 'sampleName' must be character class of", 
                   "length 1, can not be a vector or a list")
            }
            ampData(x)[[sampleName]][[ampliconName]]$smat
          }
)

#' get methylation matrix
#' 
#' If \code{AmpliconViews} object has the methylation call matrix this function
#' returns it.
#' 
#' @param x an \code{\link{AmpliconViews}} object
#' @param sampleName name of the sample
#' @param ampliconName name of the amplicon
#' @usage getMethMat(x,sampleName,ampliconName)
#' 
#' 
#' @export
#' @docType methods
#' @rdname getMethMat
setGeneric("getMethMat",function(x,sampleName,ampliconName) standardGeneric("getMethMat") )

#' @rdname getMethMat
#' @aliases getMethMat,AmpliconViews,character,character 
setMethod("getMethMat", signature(x="AmpliconViews",sampleName="character",
                                         ampliconName="character"),
          function(x,sampleName,ampliconName){
            if(length(ampliconName)>1 | length(sampleName)>1){
              stop("'ampliconName' and 'sampleName' must be character class of length 1, can not be a vector or a list")
            }
            ampData(x)[[sampleName]][[ampliconName]]$mat
          }
)


#' get average methylation for bases in an amplicon
#' 
#' The function returns average methylation for each base for
#' a given amplicon in a given sample.
#' 
#' @param x an \code{\link{AmpliconViews}} object
#' @param sampleName name of the sample
#' @param ampliconName name of the amplicon
#' @usage getAvMeth(x,sampleName,ampliconName)
#' 
#' @export
#' @docType methods
#' @rdname getAvMeth
setGeneric("getAvMeth",function(x,sampleName,ampliconName) standardGeneric("getAvMeth") )

#' @rdname getAvMeth
#' @aliases getAvMeth,AmpliconViews,character,character 
setMethod("getAvMeth", signature(x="AmpliconViews",sampleName="character",
                                 ampliconName="character"),
          function(x,sampleName,ampliconName){
            if(length(ampliconName)>1 | length(sampleName)>1){
              stop("'ampliconName' and 'sampleName' must be character class of length 1, can not be a vector or a list")
            }
            ampData(x)[[sampleName]][[ampliconName]]$meths
            
          }
)


#'  get coverage for bases in an amplicon 
#' 
#' @param x an \code{\link{AmpliconViews}} object
#' @param sampleName name of the sample
#' @param ampliconName name of the amplicon
#' @usage getCoverage(x,sampleName,ampliconName)
#' 
#' @export
#' @docType methods
#' @rdname getCoverage
setGeneric("getCoverage",function(x,sampleName,ampliconName) standardGeneric("getCoverage") )

#' @rdname getCoverage
#' @aliases getCoverage,AmpliconViews,character,character 
setMethod("getCoverage", signature(x="AmpliconViews",sampleName="character",
                                 ampliconName="character"),
          function(x,sampleName,ampliconName){
            if(length(ampliconName)>1 | length(sampleName)>1){
              stop("'ampliconName' and 'sampleName' must be character class of length 1, can not be a vector or a list")
            }
            ampData(x)[[sampleName]][[ampliconName]]$covs
            
          }
)




#'  get similarity for base methylation profiles in an amplicon 
#' 
#'  The function returns different metrics for methylation profile
#'  similarity for bases in an amplicon
#' 
#' @param x an \code{\link{AmpliconViews}} object
#' @param sampleName name of the sample
#' @param ampliconName name of the amplicon
#' @usage getSimilarity(x,sampleName,ampliconName,method="basic")
#' 
#' @export
#' @docType methods
#' @rdname getSimilarity
setGeneric("getSimilarity",function(x,sampleName,ampliconName,method="basic") standardGeneric("getSimilarity") )

#' @rdname getSimilarity
#' @aliases getSimilarity,AmpliconViews,character,character 
setMethod("getSimilarity", signature(x="AmpliconViews",sampleName="character",
                                   ampliconName="character"),
          function(x,sampleName,ampliconName,method){
            if(length(ampliconName)>1 | length(sampleName)>1){
              stop("'ampliconName' and 'sampleName' must be character class of length 1, can not be a vector or a list")
            }
            if(method=="basic"){
              ampData(x)[[sampleName]][[ampliconName]]$col.sim$similarity
            }else if(method=="asymetric"){
              ampData(x)[[sampleName]][[ampliconName]]$col.sim$msimilarity
            }else if(method=="tanay"){
              ampData(x)[[sampleName]][[ampliconName]]$col.sim$tanay
            }else if(method=="doubleJaccard"){
              m11=ampData(x)[[sampleName]][[ampliconName]]$col.sim$m11
              m01=ampData(x)[[sampleName]][[ampliconName]]$col.sim$m01
              m10=ampData(x)[[sampleName]][[ampliconName]]$col.sim$m10
              m00=ampData(x)[[sampleName]][[ampliconName]]$col.sim$m00
              s1=m11/(m11+m01+m10)    
              s2=m00/(m00+m01+m10) 
              (s1+s2)/2
            }
          }
)



#' show method for methylKit classes
#' 
#' The show method works for \code{methylRaw},\code{methylRawList},
#' \code{methylBase} and \code{methylDiff} objects
#' 
#' @examples
#' data(ampViewEx)
#' show(ampViewEx)
#' @rdname show-methods
#' @aliases show,AmpliconViews
setMethod("show",  "AmpliconViews", function(object) {
  
  cat("AmpliconViews object with",length(getSampleNames(object)),
    "samples and",length(object@data[[1]]),"amplicons\n--------------\n")
  cat("Sample Names:\n",getSampleNames(object),"\n--------------\n")
  cat("Amplicon locations\n")
  show(getAmpliconRanges(object) )
})

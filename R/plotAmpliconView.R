#### functions for plotAmpliconView


# trial for Gviz-Grid integration!!!!
# 

#```````````````````````````````````````````````````````````````````````````````
# helper internal functions to for the plotAmpliconView                       ||
#_______________________________________________________________________________

# converts given matrix mat to selected color palette based on it is intensity 
convertToColors <- function(mat,col.select=colorRampPalette(c("black", "yellow","red","purple"))(50) ) {
  # Produce 'normalized' version of matrix, with values ranging from 0 to 1
  rng <- range(mat, na.rm = TRUE)
  m <- (mat - rng[1])/diff(rng)
  # Convert to a matrix of sRGB color strings
  m2 <- m; class(m2) <- "character"
  m2[!is.na(m2)] <- rgb(colorRamp(col.select)(m[!is.na(m)]), max = 255)
  m2[is.na(m2)] <- "transparent"
  return(m2)
}



# modified version of hexypolygon from hexbin, where it plots squares
# instead of polygons
sqpolygon<-function(x, y, hexC = sqcoords(dx, dy, n = 1), dx, dy = NULL, 
                    fill = 1, border = 0, hUnit = "native", ...) 
{
  n <- length(x)
  stopifnot(length(y) == n)
  stopifnot(is.list(hexC) && is.numeric(hexC$x) && is.numeric(hexC$y))
  if (hexC$no.sep) {
    n6 <- rep.int(4:4, n)
    if (!is.null(hUnit)) {
      grid.polygon(x = unit(rep.int(hexC$x, n) + rep.int(x, n6), hUnit), 
                   y = unit(rep.int(hexC$y, n) + rep.int(y, n6), hUnit), id.lengths = n6, 
                   gp = gpar(col = border, fill = fill))
    }
    else {
      grid.polygon(x = rep.int(hexC$x, n) + rep.int(x, 
                                                    n6), y = rep.int(hexC$y, n) + rep.int(y, n6), 
                   id.lengths = n6, gp = gpar(col = border, fill = fill))
    }
  }
  else {
    n7 <- rep.int(7:7, n)
    polygon(x = rep.int(hexC$x, n) + rep.int(x, n7), y = rep.int(hexC$y, 
                                                                 n) + rep.int(y, n7), ...)
  }
}

# calculates coordinates for squares to be used in sqpolygon
#
sqcoords<-function (dx=0.5, dy = NULL, n = 1, sep = NULL) 
{
  stopifnot(length(dx) == 1)
  if (is.null(dy)) 
    dy <- dx 
  if (is.null(sep)) 
    list(x = rep.int(c(dx, 0,-dx, 0), n), y = rep.int(c(0, -dy, 0,  dy), n), no.sep = TRUE)
  else list(x = rep.int(c(dx, 0,-dx, 0, sep), n), 
            y = rep.int(c(0    ,-dy, 0,  dy, sep), 
                        n), no.sep = FALSE)
}

# function makes a new viewPort for meta-meth profiles
# and plots them on the new viewport
# @param x coordinate of bottom left corner of the viewport
# @param y y coordinate of the bottom left corner of the viewport
# @param height height of the viewport 
# @param width of the viewport
# @param metas meta-profiles
# @param importance importance scores of meta-profiles
# @param win.range  original range of the plotted window in the chromosome
# @param org.positions original positions of the bases
# @param meta.th control how many meta-profiles are drawn 
#tgv=plotTracks(trackList)
#rwidth =dev.size("px")[1]
#rheight=dev.size("px")[2]

#legend.width= coords(tgv$titles)[1,3]/rwidth
#flank= coords(tgv$titles)[1,1]/rwidth

#win.range=c(av.obj$start, av.obj$end)

#width=1-legend.width
#x=legend.width+flank+0.02
#y=0.04
#org.positions=as.numeric(colnames(av.obj$mf$metas))
#metas=av.obj$mf$metas
#importance=av.obj$mf$meta.exp
#win.range=c(av.obj$start,av.obj$end)

grid.meta<-function(x,y,width,height,metas,importance,win.range,org.positions,
                    meta.th=80){
  metaViewport <- viewport(x = unit(x,"npc"),y = unit(y, "npc"), 
                           width =unit(width,"npc"),height = unit(height, "npc"),
                           just=c("left","bottom"),
                           yscale = c(5,0))
  pushViewport(metaViewport)
  grid.rect()
  
  
  n <- length(org.positions)
  orgCord=org.positions
  v_x <- (orgCord-win.range[1]+1)/(win.range[2]-win.range[1]+1) # get coordinates
  if(any(v_x<0)){stop("win.range is out of bounds for the matrix")}
  
  imp.index=which(cumsum(importance)>meta.th)[1]
  no.meta=ifelse(imp.index>4,4,imp.index) # number of meta profiles to be plotted
  for(i in 1:no.meta){
    
    
    grid.lines(x=c(v_x[1],v_x[length(v_x)]),y=unit(i,"native") )
    
    my.col  =rep("black",length(v_x))
    my.col[metas[i,]<0]="white"
    grid.circle(x=v_x,y=unit(i,"native"),r=0.02,gp=gpar(fill=my.col))
    grid.rect(x = unit(0, "npc"), y = unit(i-0.4,"native"),
              width = unit(0.5*importance[i]/100, "npc"), height = unit(0.25, "native"),
              just=c("left"),gp=gpar(fill="blue",col=NA))
    grid.text(paste(round(importance[i],1),"%"),x=  unit(0.5*importance[i]/100, "npc"),
              y = unit(i-0.4,"native"),
              just=c("left"),gp=gpar(fill="blue",col=NULL,cex=0.8))
  }
  
  
  
  upViewport()
}

##########################################################################

# function plots similarity heatmap LD on the open Grid page
# 
#rwidth =dev.size("px")[1]
#rheight=dev.size("px")[2]

#legend.width= coords(tgv$titles)[1,3]/rwidth
#flank= coords(tgv$titles)[1,1]/rwidth

#win.range=c(av.obj$start, av.obj$end)
# for test:
# grid.newpage()
# pushViewport(plotViewport(margins=c(5.1, 4.1, 4.1, 2.1), xscale = c(0.5,ncol(sim.mat)+.5 ), yscale= c(-ncol(sim.mat),2) ))
simHeat<-function(av.obj,legend.width=0.2, spacing=0.02,flank=0.01,
                  win.range=NULL,sim.measure="similarity",
                  na.threshold=10,
                  col.select=colorRampPalette(c("black", "yellow","purple"))(50)
){
  
  # get similarity data
  sim.mat  <- av.obj$col.sim[[sim.measure[1]]]
  if(sim.measure[[1]]=="dissimilarity"){
    sim.mat  <- av.obj$col.sim[["similarity"]]
    sim.mat=1-sim.mat
  }else{
    sim.mat  <- av.obj$col.sim[[sim.measure[1]]]
  }
  #mMeth=colMeans(mat,na.rm=TRUE)
  
  # put NA values to sim.mat if paired coverage is low
  if(!is.null(na.threshold)){
    cmat=av.obj$col.sim$m11+ av.obj$col.sim$m00+av.obj$col.sim$m10+av.obj$col.sim$m01
    sim.mat[cmat<na.threshold]=NA
  }
  sim.mat[is.nan(sim.mat)]=NA # remove NaN
  
  # trick to scale colors from 0->1 add a 0 to the diagonal
  sim.mat[1,1]=0
  sim.mat[2,2]=1
  
  col.mat=convertToColors(sim.mat,col.select) # convert to colors
  
  
  # figure out locations of cells of the 45 degree
  # rotated heatmap
  # we will create the heatmap with hexagonal cells
  len=ncol(sim.mat)
  run.len=len-1
  Y<-X<-mX<-mY<-c()
  for(i in 1:(len-1)  ){
    
    X=c(X,seq(i+0.5,by=0.5,len=run.len))
    #Y=c(Y,-(1:run.len) )
    Y=c(Y,-(seq(1,by=0.5,length.out=run.len)))
    mX=c(mX,rep(i,run.len))
    mY=c(mY,seq(to=len,by=1,length.out=run.len))
    
    
    run.len=run.len-1
  }
  
  
  #figure out the relative locations of CpGs
  if( is.null(colnames(sim.mat)) ){
    n <- ncol(sim.mat)
    orgCord=1:n
    v_x <- orgCord/n
    X_x <- seq(0, 1, len=n)
    
  }else{
    n <- ncol(sim.mat)
    orgCord=as.numeric(colnames(sim.mat))
    
    if(! is.null(win.range)){
      v_x <- (orgCord-win.range[1]+1)/(win.range[2]-win.range[1]+1)
      
      if(any(v_x<0)){stop("win.range is out of bounds for the matrix")}
      
    }else{
      v_x <- (orgCord-min(orgCord)+1)/(max(orgCord)-min(orgCord)+1)
    }
    
    X_x <- seq(0, 1, len=n)
  }
  
  
  
  # PART 0: PLOT HEATMAP
  #_____________________________________________________________________________
  
  heatvp<-viewport(x = unit((legend.width+spacing), "npc"), y = unit(0.5, "npc"),
                   width = unit(1-(legend.width+spacing+flank), "npc"), height = unit(1, "npc"),
                   just="left",
                   xscale = c(0.5,ncol(sim.mat)+.5 ),
                   yscale = c(-(ncol(sim.mat)+1)*0.5,0.5) )
  
  pushViewport(heatvp)
  grid.rect() 
  #grid.yaxis()
  # get colors from color matrix
  my.cols=apply(cbind(mX,mY),1, function(x,col.mat) col.mat[x[2],x[1]],col.mat=col.mat  )
  
  #plot the cells
  # +0.3 on Y controls the relative height of the heatmap within the plot
  sqpolygon(X,Y+0.3, hexC=sqcoords(dx = 0.5,dy=0.5, sep=NULL), border = "black", fill=my.cols)
  
  
  # PLOT segments 
  # PART 1: DRAW CPG locations
  #_____________________________________________________________________________
  y_a=0.5;y_b=0.4;y_c=-0.2;y_d=-0.4 # controls the segment heights 
  seg.grb=grid.segments(x0 = unit(1:n,"native"), 
                        x1 = unit(v_x,"npc"), y0 = unit(y_c,"native"), y1 = unit(y_b,"native"))
  grid.segments(x0 = unit(v_x,"npc"), 
                x1 = unit(v_x,"npc"), y0 = unit(y_b,"native"), y1 = unit(y_a,"native"))
  grid.segments(x0 = unit(1:n,"native"), 
                x1 = unit(1:n,"native"), y0 = unit(y_c,"native"), y1 = unit(y_d,"native"))
  #grid.circle( x = unit(v_x,"npc"), y = unit(1,"npc")+unit(1,"lines"),r=0.01)
  #lab.grb=grid.text(label = orgCord, x = unit(v_x,"npc"), y = unit(1,"npc")+unit(2.5,"lines"), rot = 90,
  #                  just="left",
  #                  gp=gpar(cex=0.6))
  
  
  # PART 2: HEATMAP LEGEND
  #_____________________________________________________________________________
  upViewport()
  annotationViewport <- viewport(x = unit(flank+3*legend.width/4,"npc"),y = unit(0.5, "npc"), 
                                 width = legend.width/4,height = unit(0.9, "npc"),
                                 just=c("left"),
                                 xscale = c(0.5,ncol(sim.mat)+.5 ))
  pushViewport(annotationViewport)
  grid.rect()
  grid.yaxis(gp=gpar(cex=0.7)    )
  
  #my.cols=colorRampPalette(c("black", "yellow","red","purple"))(50)
  grid.raster(rev(convertToColors(seq(0,1,by=0.1),col.select)), interpolate = TRUE,
              height=unit(1,"npc"),width=1)
  #grid.raster(rev(convertToColors(seq(0,1,by=0.1),my.cols )), interpolate = TRUE)
  upViewport()
  grid.text(sim.measure[[1]],x = unit(flank,"npc"), rot=90,just="left")
  
  upViewport()
}

#```````````````````````````````````````````````````````````````````````````````
# main exported function: plotAmpliconView                                    ||
#_______________________________________________________________________________

#' plot ampliconView object with Gviz
#' 
#' @param obj \code{\link{AmpliconViews}} object with one amplicon. 
#'               Use \code{\link{getAmplicon}} function to get one amplicon
#'               from one sample if necessary.
#' @param example.reads if TRUE (default: FALSE), example reads sampled from actual 
#'                      reads will be displayed. Each \code{AmpliconViews} object will
#'                      have a set of example reads that correspond to meta-methylation
#'                      profiles.
#' @param h.panel height of the panels, this should be a numeric vector with length
#'          equaling to the number of panels plotted
#' @param sim.heat logical, if TRUE a similarity heatmap is drawn under the tracks,
#'                  showing similarity of base methylation profiles
#' @param h.heat relative height of the similarity heatmap, should be a numeric value between
#'               0 and 0.5. Smaller the value, smaller the heatmap on the plot.
#' @param sim.col colors for similarity heatmap   
#'                default: colorRampPalette(c("black", "yellow","purple"))(50) 
#' @param sim.na.threshold  number of pair-wise observation needed to calculate similarity
#'                          scores from CpGs with low number of overlap may be unreliable          
#' @param sim.measure similarity measure, should be one of "similarity","tanay","msimilarity",
#' @param meta.exp % of data set explained by meta-profiles, the meta-profiles explaining 
#'                  data set will be plotted
#' @param newpage if TRUE the plot will start a fresh graphical device.
#' @param ...    Other arguments to Gviz::plotTracks  
#' 
#'         
#' @examples
#'       data(ampViewEx) # load example data
#'       myAmp=getAmplicon(ampViewEx,"mock4","chr18_69674375_69674775")
#'       plotAmpliconView(myAmp)
#'
#' @import Gviz
#' @import grid
#'                   
#' @export
#' @docType methods
plotAmpliconView<-function(obj,example.reads=FALSE,h.panel=NULL,
                           sim.heat=TRUE,h.heat=0.5,
                           sim.col=colorRampPalette(c("blue", "yellow","red"))(50) ,
                           sim.na.threshold=10,
                           sim.measure="similarity",newpage=TRUE,meta.exp=80,
                           ...){
  
  #require(Gviz)
  
  # sjekk object class and number of Amplicons in the object
  stopifnot(class(obj) == "AmpliconViews")
  if( length(ampData(obj))> 1 | length(ampData(obj)[[1]]) >1 ){
    stop("\nplotAmpliconView can be used to plot one amplicon at a time.","\n",
         "Subset your AmpliconViews object with getAmplicon() function","\n")
  }
  
  av.obj=ampData(obj)[[1]][[1]]
  
  #cov.th=0
  #col.retain=which(av.obj$covs>cov.th) # columns to be retained
  #library(gplots)
  v=GRanges(seqnames=av.obj$chr,
            IRanges(start=as.numeric(colnames(av.obj$mf$metas)),width=1))
  
  trackList=list()
  values(v)=av.obj$meths
  Mtrack <- DataTrack(v, name = "methylation rat.",type=c("h","g"),ylim=c(0,1),
                      lwd=5,h=2,v=0,baseline=0,lwd.baseline=2,col.baseline="black",
                      col.grid="gray"
  )
  values(v)=av.obj$covs
  Ctrack <- DataTrack(v, name = "coverages",type=c("h","g"),col="red",lwd=5,
                      col.baseline="black",baseline=0,lwd.baseline=2,
                      col.grid="gray",ylim=c(0,max(av.obj$covs))
                      ,h=2,v=0)
  trackList=c(trackList,Mtrack,Ctrack)
  sizes=c(0.25,0.25)
  
  
  #values(v)=av.obj$mf$imp.metas[i,]
  
  trackList= c(trackList, 
               AnnotationTrack(GRanges(seqnames=av.obj$chr,ranges=IRanges(1,2)),
                               name="meta-methylation",
                               background.panel = "#FFFEDB",baseline= 0,lwd.baseline=2,
                               col.baseline="black")
  )
  sizes=c(sizes,0.25)
  
  
  # add heatmap
  if(example.reads){
    values(v)=t(av.obj$smat) #; image(t(smat))                                           
    Htrack <- DataTrack(v, name = "example profiles",type="heatmap",
                        gradient=colorRampPalette(c("blue", "red"))( 2 ))
    trackList=c(trackList,Htrack)
    sizes=c(sizes,2)
  }
  
  
  desc=paste( getSampleNames(obj),"|",
             av.obj$chr,av.obj$start,av.obj$end,"|",av.obj$tag,sep=" ")
  axisTrack <- GenomeAxisTrack(name= NULL,col.title="black",
                               showTitle=TRUE,cex.title=0.4,
                               littleTicks=TRUE)
  #pdf("/work2/gschub/altuna/Rdev/ampliconBiSeq/trial.pdf",width=8, height=10)

  if(is.null(h.panel)){
    sizes=c(0.15,sizes)
  }
  else if( !is.null(h.panel) & length(c(axisTrack,trackList)) == length(h.panel) ){
    sizes=h.panel
  }else if(length(c(axisTrack,trackList)) != length(h.panel)){
    stop("\n'h.panel' length is not equal to genome tracks number\n",
         "When plotting by default arguments there are 4 tracks:\naxis, methylation ratio,coverage and meta-methylation tracks\n",
         "If you have added additional tracks you must give a number higher than 4")
  }
  
  if(sim.heat){
    
    if(newpage){grid.newpage()}
    
    # see if we are plotting in a layout
    lo.col=1 # no of cols
    lo.row=1 # no of rows
    cu.col=1 # curret column
    cu.row=1 # current row
    vpaths=(current.vpPath() )
    if(length(vpaths)-2 >0){
      seekViewport(vpaths[[length(vpaths)-2]]) # get the second to last VP
      
      if(!is.null(current.viewport()$layout$ncol) ){
        lo.col=current.viewport()$layout$ncol
        lo.row=current.viewport()$layout$nrow
      }
      seekViewport(vpaths[[length(vpaths)-1]])
      cu.col=(current.viewport())$layout.pos.col[1]
      cu.row=(current.viewport())$layout.pos.row[1]
      
    }
    
    
    # get the layout
    lo=grid.layout(nrow=2, ncol=1, 
                   just="left",heights=c(1-h.heat,h.heat) )
    pushViewport( viewport(layout=lo) )
    
    pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
    tgv=plotTracks(c(axisTrack,trackList),from=av.obj$start,
                   to=av.obj$end,sizes=sizes,add=TRUE,
                   col.title="black",col.axis="black",
                   background.title=rgb(0,0,0,0)# transperant
                   , ... )
    coord=coords(tgv$titles) # get title coordinates
    rwidth =dev.size("px")[1]/lo.col # current device size
    rheight=dev.size("px")[2]*(1-h.heat) # this time have to divide device size by 2
    flank= (coord[1,1]-((cu.col-1)*rwidth))/rwidth # flank between the plot and device edge
    
    # height of the meta profile plot, relative to main graph
    my.height=(coord[grep("meta",rownames(coord) ),4]-coord[grep("meta",rownames(coord) ),2])/rheight
    left.flank=(coord[grep("meta",rownames(coord) ),3]-coord[grep("meta",rownames(coord) ),1])/rwidth+ 0.02 +flank
    bott.flank=1-coord[grep("meta",rownames(coord) ),4]/rheight
    grid.meta(x=left.flank,y=bott.flank,
              width=1-left.flank-flank,
              height=my.height,
              metas=av.obj$mf$metas,
              importance=av.obj$mf$meta.exp,
              win.range=c(av.obj$start,av.obj$end),
              org.positions=as.numeric(names(av.obj$meths)),
              meta.th=meta.exp)  
    grid.text(desc,y=unit(0.98,"npc"),gp=gpar(cex=0.8) ) # plot title
    upViewport()
    
    ## HEAT MAP
    pushViewport(viewport(layout.pos.col=1, layout.pos.row=2))
    
    # get coordinates so grid knows where to put what
    rwidth =dev.size("px")[1]/lo.col
    rheight=dev.size("px")[2]
    
    coord[,c(1,3)]=coord[,c(1,3)]-((cu.col-1)*rwidth)
    legend.width= coord[1,3]/rwidth
    flank= (coord[1,1])/rwidth
    #flank= coord[1,1]/rwidth
    
    win.range=c(av.obj$start, av.obj$end)
    simHeat(av.obj,legend.width=legend.width, spacing=0.02,flank=flank,
            win.range=win.range,sim.measure=sim.measure,
            na.threshold=sim.na.threshold,
            col.select=sim.col )
    upViewport()
    
  }else{
    tgv=plotTracks(c(axisTrack,trackList),from=av.obj$start,
                   to=av.obj$end,sizes=sizes,...)
    coord=coords(tgv$titles) # get title coordinates
    rwidth =dev.size("px")[1] # current device size
    rheight=dev.size("px")[2] # this time have to divide device size by
    flank= coords(tgv$titles)[1,1]/rwidth # flank between the plot and device edge
    
    # height of the meta profile plot, relative to main graph
    my.height=(coord[grep("meta",rownames(coord) ),4]-coord[grep("meta",rownames(coord) ),2])/rheight
    
    
    left.flank=(coord[grep("meta",rownames(coord) ),3]-coord[grep("meta",rownames(coord) ),1])/rwidth+ 0.02 +flank
    bott.flank=1-coord[grep("meta",rownames(coord) ),4]/rheight
    grid.meta(x=left.flank,y=bott.flank,
              width=1-left.flank-flank,
              height=my.height,
              metas=av.obj$mf$metas,
              importance=av.obj$mf$meta.exp,
              win.range=c(av.obj$start,av.obj$end),
              org.positions=as.numeric(names(av.obj$meths)),
              meta.th=80)  
    grid.text(desc,y=unit(0.97,"npc"),gp=gpar(cex=0.8) )
    
    
  }
  #dev.off()
}

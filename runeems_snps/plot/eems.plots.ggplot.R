mcmcpath <- paste0("../data/barrier-schemeX-nIndiv300-nSites3000-EEMS-nDemes153-simno", 1:3)
mcmcpath = "/data/eems-project/human_origins/results/africa/output/300"
mcmcpath = "/media/peterb/pps1/human_origin_2/p2/world/output/350_run1"


plotpath <- "tmp"


require(ggplot2)
require(fields)
require(mapdata)
require(SDMTools)


############################################################
# Plotting
############################################################

filled.polys <- function(dimns, z, bins=13, zvar=NULL, ...){
    z <- c(z)
    zvar <- c(zvar)
    minValues <- apply(dimns$marks, 2, min, na.rm=T)
    maxValues <- apply(dimns$marks, 2, max, na.rm=T)
    xstep <- dimns$xmrks[2] - dimns$xmrks[1]
    ystep <- dimns$ymrks[2] - dimns$ymrks[1]
    BINWIDTH <- (diff(range(z))/bins) # reference from ggplot2 code
    boundaryValue <- min(z) - BINWIDTH * 1.5
    boundaryValue <- min(z)

    xl <- minValues[1] - xstep
    xr <- maxValues[1] + xstep
    yl <- minValues[2] - ystep
    yr <- maxValues[2] + ystep
    mzv <- min(zvar)

    b1<-data.frame(x=xl, y=dimns$ymrks, z=boundaryValue, zvar=mzv)
    b3<-data.frame(x=xr, y=dimns$ymrks, z=boundaryValue, zvar=mzv)
    b2<-data.frame(x=dimns$xmrks, y=yl, z=boundaryValue, zvar=mzv)
    b4<-data.frame(x=dimns$xmrks, y=yr, z=boundaryValue, zvar=mzv)
    boundary<-rbind(b1, b2, b3, b4)

    data_set <- data.frame(dimns$marks, z=log(z), zvar=zvar)
    names(data_set) <- c("x", "y", "z", "zvar")
    bounded_data <- rbind(data_set,boundary)
    
    P <- ggplot()
#    P <- ggplot(bounded_data, aes(x=y, y=x, z = z))
    P <- P + stat_contour(geom="polygon", aes(x=x, y=y, z=z, fill=..level..), 
                          data=bounded_data,
                          bins=13, na.rm=T) +
    coord_fixed(xlim=c(min(data_set$x),max(data_set$x)),
                                  ylim=c(min(data_set$y),max(data_set$y)))
#    P <- P + stat_contour(geom="polygon", aes(x=x, y=y, z=zvar,
#                                              alpha=..level..), 
#                          bins=10,
#                          data=bounded_data, fill="grey60") 
    P <- P + scale_alpha_continuous(range=c(.0,.9))
    P <- P + theme(panel.background = element_rect(fill = 'grey20' )) 
    return(P)
}


average.eems.contours.ggplot <- function(mcmcpath,dimns, plot.params,
                                  is.mrates, alphaplot=T) {
    if (is.mrates) {
        message('Plotting effective migration rates m : posterior mean and variance')
        files <- c('/mcmcmtiles.txt','/mcmcmrates.txt', 
                   '/mcmcxcoord.txt','/mcmcycoord.txt')
    } else {
        message('Plotting effective diversity rates q : posterior mean and variance')
        files <- c('/mcmcqtiles.txt','/mcmcqrates.txt', 
                   '/mcmcwcoord.txt','/mcmczcoord.txt')
    }
    mcmcpath <- check.files.at.path(files, mcmcpath)
    n_runs <- length(mcmcpath)
    if (n_runs==0) { return(0) }

    Z <- get.z(mcmcpath, dimns, is.mrates)
    Zmean <- Z[[1]]
    Zvar <- Z[[2]]
    

    one.eems.contour.ggplot(mcmcpath, dimns, Zmean, Zvar,
                            plot.params, is.mrates, alphaplot)
}

one.eems.contour.ggplot <- function(mcmcpath, dimns, Zmean, Zvar,
                                    is.mrates, alphaplot=T){


    demes <- read.demes(mcmcpath)
    edges <- read.edges(mcmcpath)
    sample.sizes <- read.sample.sizes(mcmcpath)

    P <- filled.polys(dimns, Zmean, zvar=Zvar)
    P <- P + add.map() +
        add.grid(demes, edges) +
        add.samples(demes, sample.sizes)


    #breaks <- quantile(Zmean, c(0:12/12))

    #P <- P + get.eems.colors(is.mrates, breaks=breaks, labels=as.numeric(breaks))

    return(P)

}

add.map <- function(wrap=T){
    worldmap <- map_data("world")
    #worldmap <- map_data("world")
    if(wrap){
    worldmap$long[worldmap$long < -20] <- worldmap$long[worldmap$long < -20] + 360
    max.grp <- max(worldmap$group)
    prev <- 1
    prev.long <- worldmap$long[1]

    for(i in 2:nrow(worldmap))
        if(worldmap$group[i] != worldmap$group[i-1]){
            prev <- i
            prev.long <- worldmap$long[i]
        } else {
            if( abs(worldmap$long[i] - prev.long) > 300){
                worldmap$lat[i] = NA
            }
        }

    }
    ggworld <- geom_path(data=worldmap, 
                            aes(x=long,  y=lat, group=group),
                            fill=NA, colour="black")
    return(ggworld)
}


add.grid <- function(demes, edges){
    x.start <- demes[edges[,1] ,1]
    y.start <- demes[edges[,1] ,2]
    x.end <- demes[edges[,2] ,1]
    y.end <- demes[edges[,2] ,2]
    dt <- data.frame(x.start, y.start, x.end, y.end)
    gs <- geom_segment(aes(x=x.start, y=y.start, xend = x.end, yend= y.end),
                       colour="black", size=1, data=dt,
                       alpha=0.2
                       )
}

add.samples <- function(demes, sample.sizes){
    dx <- demes[sample.sizes[,1],1]
    dy <- demes[sample.sizes[,1],2]
    max.ss <- max(sample.sizes[,2])
    ss <- sample.sizes[,2] 
    dt <- data.frame(dx, dy, ss)
    q <- geom_point(aes(x=dx, y=dy, size=ss), data=dt, col="black") 
    q2<- scale_size_continuous(limits=c(1,max.ss), range=c(4,10), guide=F) 
    return(list(q,q2))
}

get.eems.colors <- function(par=F, ...) {
    if(!par){
	eems.colors <- c("#994000","#CC5800","#FF8F33","#FFAD66","#FFCA99","#FFE6CC", ## orange sequence
			 "#ffffff",                                                   ## white
			 "#BBFFBB","#86FF86","#50FF50","#00BB00","#008600","#005000") ## green sequence
    } else {
	eems.colors <- c("#994000","#CC5800","#FF8F33","#FFAD66","#FFCA99","#FFE6CC", ## orange sequence
			 "#ffffff",                                                   ## white
			 "#CCFDFF","#99F8FF","#66F0FF","#33E4FF","#00AACC","#007A99") ## blue sequence
    }
    scale_fill_gradientn(colours=eems.colors,  ...
                         )
}
############################################################
#  Density calculation
############################################################

get.z <- function(mcmcpath, dimns, is.mrates, standardize=T){
    niter <- 0
    Zmean <- matrix(0,dimns$nxmrks,dimns$nymrks)
    for (path in mcmcpath) {
        if(standardize){
            rslt <- standardize.rates(path,dimns, is.mrates)
        } else{
            print("here")
            rslt <- read.rates(path,dimns, is.mrates)
        }
        niter <- niter + rslt$niter
        Zmean <- Zmean + rslt$Zvals
    }
    Zmean <- Zmean/niter
    Zvar <- matrix(0,dimns$nxmrks,dimns$nymrks)
    for (path in mcmcpath) {
        rslt <- standardize.rates.var(path,dimns,Zmean, is.mrates)
        Zvar <- Zvar + rslt$Zvar
    }
    Zvar <- Zvar/(niter-1)

    return(list(mean=Zmean, var=Zvar))
}

read.rates <- function(mcmcpath, dimns, is.mrates){
    voronoi <- read.voronoi(mcmcpath,is.mrates, do.log=F)
    rates <- voronoi$rates
    tiles <- voronoi$tiles
    xseed <- voronoi$xseed
    yseed <- voronoi$yseed
    Zvals <- matrix(0,dimns$nxmrks,dimns$nymrks)
    niter <- 100#length(tiles)
    count <- 0
    for (i in 1:niter) {
        print(i)
        now.tiles <- tiles[i]
        now.rates <- rates[(count+1):(count+now.tiles)]
        now.xseed <- xseed[(count+1):(count+now.tiles)]
        now.yseed <- yseed[(count+1):(count+now.tiles)]
        now.seeds <- cbind(now.xseed,now.yseed)
        zvals <- compute.contour.vals(dimns,now.seeds,now.rates, 
                                      standardize=F)
        Zvals <- Zvals + zvals
        count <- count + now.tiles
    }
    return(list(Zvals=Zvals,niter=niter))
}

standardize.rates <- function(mcmcpath,dimns,is.mrates) {
    voronoi <- read.voronoi(mcmcpath,is.mrates)
    rates <- voronoi$rates
    tiles <- voronoi$tiles
    xseed <- voronoi$xseed
    yseed <- voronoi$yseed
    Zvals <- matrix(0,dimns$nxmrks,dimns$nymrks)
    Zvals[!dimns$filter] <- NA
    niter <- 100#length(tiles)
    count <- 0
    for (i in 1:niter) {
        print(i)
        now.tiles <- tiles[i]
        now.rates <- rates[(count+1):(count+now.tiles)]
        now.xseed <- xseed[(count+1):(count+now.tiles)]
        now.yseed <- yseed[(count+1):(count+now.tiles)]
        now.seeds <- cbind(now.xseed,now.yseed)
        zvals <- compute.contour.vals(dimns,now.seeds,now.rates)
        Zvals <- Zvals + zvals
        count <- count + now.tiles
    }
    return(list(Zvals=Zvals,niter=niter))
}
standardize.rates.var <- function(mcmcpath,dimns,Zmean,is.mrates) {
    voronoi <- read.voronoi(mcmcpath,is.mrates)
    rates <- voronoi$rates
    tiles <- voronoi$tiles
    xseed <- voronoi$xseed
    yseed <- voronoi$yseed
    Zvar <- matrix(0,dimns$nxmrks,dimns$nymrks)
    niter <- 100#length(tiles)
    count <- 0
    for (i in 1:niter) {
        print(i)
        now.tiles <- tiles[i]
        now.rates <- rates[(count+1):(count+now.tiles)]
        now.xseed <- xseed[(count+1):(count+now.tiles)]
        now.yseed <- yseed[(count+1):(count+now.tiles)]
        now.seeds <- cbind(now.xseed,now.yseed)
        zvals <- compute.contour.vals(dimns,now.seeds,now.rates)
        Zvar <- Zvar + (zvals - Zmean)^2
        count <- count + now.tiles
    }
    return(list(Zvar=Zvar,niter=niter))
}


compute.contour.vals <- function(dimns,seeds,rates,
                                 standardize=T,
                                 distf=rdist) {
    ## Here 'seeds' stores the generator seeds of a Voronoi tessellation
    ## and 'rates' stores the log10-transformed rates of the tiles.
    ## If there are C seeds in the partition, then 'seeds' is a matrix
    ## with C rows and 2 columns and 'rates' is a vector with C elements
    distances <- distf(dimns$marks[dimns$filter,],seeds)
    closest <- apply(distances,1,which.min)
        zvals <- matrix(NA,dimns$nxmrks,dimns$nymrks)
        zvals[dimns$filter] <- rates[closest]
        print(length(zvals[dimns$filter]))
        if(standardize)  zvals <- zvals - mean(zvals, na.rm=T)
    return(zvals)
}

############################################################
# Read files
############################################################
check.files.at.path <- function(files, paths=".", strict=F){
    # checks if ech subfolder `paths` contains all the files in `files`.
    file_paths <- c(outer(paths, files, paste, sep=.Platform$file.sep))
    file_exist <- file.exists(file_paths)
    if( all(file_exist) ) return( paths )
    
    print(paste0(file_paths[!file_exist], " not found", collate="\n") )

    if( strict ){
        stop( )
    }
    
    f <- matrix(file_exist, nrow=length(paths))
    return(paths[ rowSums(f) == ncol(f)])
}


read.edges <- function(mcmcpath) {
    edges <- read.table(paste(mcmcpath,'/edges.txt',sep=''), 
                        colClasses=numeric())
    edges <- as.matrix(edges)
    ## Previously EEMS output the edges with one vertex per line and
    ## the six neighbors of each vertex listed in order
    ## Currently EEMS outputs the edges with one edge per line (as a
    ## pair of vertices) which is more general
    if (ncol(edges)==6) {
        ## Convert the old format to the new format
        edges0 <- edges
        edges <- matrix(0,nrow=sum(edges0>0),ncol=2)
        nv <- nrow(edges0)
        nn <- ncol(edges0)
        nodes <- 1:nv
        e <- 0
        for (a in 1:nv) {
        for (i in 1:nn) {
            b <- edges0[a,i]
            if (b %in% nodes) {
                e <- e + 1
                edges[e,1] = a
                edges[e,2] = b
            }
        } }
    }
    return(edges)
}

read.demes <- function(mcmcpath){
    demes <- scan(paste(mcmcpath,'/demes.txt',sep=''), what=numeric(),
                  quiet=TRUE)
    demes <- matrix(demes, ncol=2, byrow=T)
    return(demes)
}


read.sample.sizes <- function(mcmcpath){
    ipmap <- scan(paste(mcmcpath,'/ipmap.txt',sep=''),what=numeric(),quiet=TRUE)
    sizes <- table(ipmap)
    sampled.demes <- as.numeric(names(sizes))
    sizes <- as.numeric(sizes)
    return(cbind(sampled.demes, sizes))
}

read.dimns <- function(mcmcpath,nxmrks=NULL,nymrks=NULL) {
    # reads a dimns object from mcmcpath. This object represents the 
    # plotting surface and the number of points at which densities
    # are calculated.
    #
    # Parameters:
    # mcmcpath    : path to a folder with eems output
    # nxmrks, nymrks : number of points to calcualte density on
    #
    # Returns:
    # Object of type dimns
    
    outer <- scan(paste(mcmcpath,'/outer.txt',sep=''), what=numeric(), 
                  quiet=TRUE)
    outer <- matrix(outer, ncol=2, byrow=TRUE)

    xmin <- min(outer[,1])
    xmax <- max(outer[,1])
    ymin <- min(outer[,2])
    ymax <- max(outer[,2])
    ## Choose the number of interpolation in each direction
    if (is.null(nxmrks)&&is.null(nymrks)) {
        xy.asp.ratio <- (xmax-xmin)/(ymax-ymin)
        if (xy.asp.ratio>1) {
            nxmrks <- 100
            nymrks <- round(nxmrks/xy.asp.ratio)
        } else {
            nymrks <- 100
            nxmrks <- round(nymrks*xy.asp.ratio)
        }
    }

    ## The interpolation points are equally spaced
    xmrks <- seq(xmin,xmax,length=nxmrks)
    ymrks <- seq(ymin,ymax,length=nymrks)
    marks <- cbind(rep(xmrks,times=nymrks),rep(ymrks,each=nxmrks))

    pip <- pnt.in.poly(marks, outer)[,3]
    filter <- pip == 1
    return(list(nxmrks=nxmrks, xmrks=xmrks, xrange=c(xmin,xmax), 
                xspan=(xmax-xmin),
                nymrks=nymrks, ymrks=ymrks, yrange=c(ymin,ymax),
                yspan=(ymax-ymin),
                marks=marks, filter=filter, outer=outer))
}
read.voronoi <- function(mcmcpath, is.mrates, do.log=T) {
    if (is.mrates) {
        rates <- scan(paste(mcmcpath,'/mcmcmrates.txt',sep=''),what=numeric(),quiet=TRUE)
        tiles <- scan(paste(mcmcpath,'/mcmcmtiles.txt',sep=''),what=numeric(),quiet=TRUE)
        xseed <- scan(paste(mcmcpath,'/mcmcxcoord.txt',sep=''),what=numeric(),quiet=TRUE)
        yseed <- scan(paste(mcmcpath,'/mcmcycoord.txt',sep=''),what=numeric(),quiet=TRUE)
    } else {
        rates <- scan(paste(mcmcpath,'/mcmcqrates.txt',sep=''),what=numeric(),quiet=TRUE)
        tiles <- scan(paste(mcmcpath,'/mcmcqtiles.txt',sep=''),what=numeric(),quiet=TRUE)
        xseed <- scan(paste(mcmcpath,'/mcmcwcoord.txt',sep=''),what=numeric(),quiet=TRUE)
        yseed <- scan(paste(mcmcpath,'/mcmczcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    }
    if(do.log) rates <- log10(rates)
    return(list(rates=rates,tiles=tiles,xseed=xseed,yseed=yseed))
}

############################################################
#  Main script
############################################################
run <- function(mcmcpath, plotpath='tmp', alphaplot=T, ...){
    infiles <- c("ipmap.txt", "demes.txt", "edges.txt")
    mcmcpath <- check.files.at.path(infiles, mcmcpath)

    n_runs <- length(mcmcpath)

    if(n_runs == 0) return( 0 )
    message("processing the following EEMS output")
    message(paste0(mcmcpath,collate="\n"))
    
    dimns <- read.dimns(mcmcpath[1])


    average.eems.contours.ggplot(mcmcpath, dimns, plot.params=list(),
                                 is.mrates=T, alphaplot=alphaplot, ...)
    ggsave(paste0(plotpath,"-mrates.pdf"))
    average.eems.contours.ggplot(mcmcpath, dimns, plot.params=list(),
                                 is.mrates=F, alphaplot=alphaplot, ...)
    ggsave(paste0(plotpath,"-qrates.pdf"))
    dev.off()
}

require(snpStats) #for read.plnk function

source("wishart.r")


project.folder = "/data/eems-project/human_origins/projects/europe/input/"
delta.folder = "/data/eems-project/human_origins/projects/europe/output/300_2/"



total.nnsnp = 600000


get.snp.mat <- function( snp.id, data){
    snp <- data$genotypes[,snp.id]
    snp <- as(snp, "numeric")[,1]
    snp[is.na(snp)] <- 0
    outer(snp, snp, function(x,y)(x-y)*(x-y))
}

subset.plink = function(data, indiv){
    to.keep = sapply(indiv[,2], function(x) which(x ==data$fam$member))
    data$fam = data$fam[to.keep,]
    data$genotypes = data$genotypes[to.keep,]
    data
}

subset.maf = function(data, maf = 0){
    maf.obs <- col.summary(data$genotypes)$MAF
    to.keep = sapply(maf, function(x) which(x < maf.obs))
    data$map = data$map[to.keep,]
    data$genotypes = data$genotypes[,to.keep]
    data
}

load.snp.data = function(bed.file="chr2", indiv=NULL, min.maf=0){
    chr2 = read.plink(bed.file)
    if(!is.null(indiv)){
        chr2 = subset.plink(chr2, indiv) 
    }

    europe = subset.maf(chr2, min.maf)
    return(europe)
}

get.diff.matrices <- function(diffs, data, tot.n.snp=600000, max.n.snp=100){
    n.snp.to.test <- min( max.n.snp, ncol(data$genotypes))

    f.diff <- function(x, tot.n.snp){
        as.matrix(diffs - get.snp.mat(x, data)/tot.n.snp)
    }
    diffs.one.removed <- lapply(1:n.snp.to.test, f.diff, tot.n.snp)
}


get.obs.matrices <- function(n.posterior.samples=10){
    deltas<- load.posterior(delta.folder, n.posterior.samples)
}

do.full.analysis <- function(max.n.snp=1, tot.n.snp=600000,
                             n.posterior.samples=10){
    indiv <- read.table(paste0(project.folder, "europe.order"))
    diffs <- read.table(paste0(project.folder, "europe.diffs"))
    diffs <- as.matrix(diffs)
    thetas <- read.table(paste0(delta.folder, "mcmcthetas.txt"))
    ll <- read.table(paste0(delta.folder, "mcmcpilogl.txt"))[,2]
    sigma2 <- thetas[,1]
    df <- thetas[,2]
    data <- load.snp.data(indiv=indiv)

    deltas <- get.obs.matrices(n.posterior.samples)


    diffs.one <-get.diff.matrices(diffs, data, tot.n.snp, max.n.snp)
    ll.full <- d.eems(diffs,deltas, df, sigma2, n.sims=n.posterior.samples)
    ll.one.1 <- sapply(diffs.one, d.eems, deltas, df, sigma2,
                     n.sims=n.posterior.samples)
    ll.one.2 <- sapply(diffs.one, d.eems, deltas, df-1, sigma2,
                     n.sims=n.posterior.samples)
    ll.one.3 <- sapply(diffs.one, d.eems, deltas, df*(1-1/tot.n.snp), sigma2,
                     n.sims=n.posterior.samples)

    lratio.1 <- exp(apply(ll.one.1, 2, function(x) x-ll.full)) 
    lratio.2 <- exp(apply(ll.one.2, 2, function(x) x-ll.full)) 
    lratio.3 <- exp(apply(ll.one.3, 2, function(x) x-ll.full)) 
    lratio.sum.1 <- colSums(lratio.1)
    lratio.sum.2 <- colSums(lratio.2)
    lratio.sum.3 <- colSums(lratio.3)
    p.1 <- n.posterior.samples/lratio.sum.1
    p.2 <- n.posterior.samples/lratio.sum.2
    p.3 <- n.posterior.samples/lratio.sum.3
    return( cbind(p.1,p.2,p.3) )
}

d.eems <- function(diffs, deltas, df, sigma2, n.sims){
    wishartll(deltas, diffs, sigma2, df, n.sims)
}


L <- function(x){
    if (!is.null(dim(x) ))
        x = ncol(x)
    cbind(-1, diag(rep(1, x-1)))
}


precompute <- function(diffs, deltas){
    mat <- list()
    LD <- L(diffs)
    LDt <- t(LD)
    deltas = lapply(deltas, function(Delta)-LD %*% Delta %*% LDt)
    diffs = - LD %*% diffs %*% LDt
    mat$ldet.diffs <- log(det(diffs))
    mat$inv.diffs <- solve(diffs)

    mat$det.deltas <- lapply(deltas, det)
    mat$ldet.deltas <- lapply(mat$det.deltas, log)

    mat$inv.deltas <- lapply(deltas, solve)

    mat$tr.invd.diffs <- lapply(mat$inv.deltas, trace.product, diffs)
    mat$ltr.invd.diffs <- lapply(mat$tr.invd.diffs, log)
    mat
}
LDL <- function(d2){
    d2 <- -L(d2) %*% ( d2) %*% t(L(d2))
}

lllr.wishart <- function(snp, diffs, df, df.rm, precomputes,
                         n.snps=600000){
    mat <- precomputes
    p <- dim(mat$inv.diffs)[1]
    delta <- df - df.rm

    snpn <- snp/n.snps
    d2 <- LDL(diffs-snpn)
    log.det.d2 <- log(det(d2))

    gamma.ratio <- function(df, df.rm){
        l <- length(df)
        res <- c()
        for(i in 1:l){
        res <- c(res, multivariate_gamma_function(df[i]/2, p) -
                multivariate_gamma_function(df.rm[i]/2, p))
        }
        res
    }
    gammas <- mapply(gamma.ratio, df, df.rm)

    log.c <- log(2) * delta / 2
    log.c <- log.c + log(df) * (-p) * df / 2
    log.c <- log.c - log(df.rm) * (-p) * df.rm / 2

    lllratio <- gammas + log.c 
    
    lllratio <- lllratio - mats$ldet.diffs * (p - df - 1) /2 
    lllratio <- lllratio + log.det.d2 * (p - df.rm - 1) /2 

    add.det.deltas <- function( lllr, dets){ lllr - delta/2 * dets}
    lllratio <- add.det.deltas( lllratio, unlist(mat$ldet.deltas))

    add.traces <- function( lllr, tr, df){
        lllr - 1/2/df * tr
    }
    lllratio <- add.traces( lllratio, unlist(mat$ltr.invd.diffs), df)


    traces <- sapply(mat$inv.deltas, trace.product, LDL(snpn))
    add.traces2 <- function( lllr, tr1, tr2, df){
        lllr - 1/2/df * (tr1 + tr2)
    }
    lllratio <- add.traces2( lllratio, unlist(mat$ltr.invd.diffs), 
                       log(traces), df.rm)


    return(lllratio)
}


multivariate_gamma_function<- function(a, p){                
        x = log(pi)*p*(p-1)/4 + sum(lgamma(a-(0:(p-1))/2));      
    return (x)                                               
}                                                            
                                                             
trace.product<- function(m1, m2) sum(m1 * t(m2))

harmonic.mean <- function(a){
    1/mean(1/a)
}

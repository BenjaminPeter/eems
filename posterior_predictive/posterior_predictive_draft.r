require(snpStats) #for read.plnk function

#source("wishart.r")


project.folder = "/data/eems-project/human_origins/projects/europe/input/europe"
delta.folder = "/data/eems-project/human_origins/projects/europe/output/300_1/"



tot.n.snp = 600000


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

get.freq <- function(data){
    maf.obs <- col.summary(data$genotypes)$MAF
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

    deltas <- load.posterior(delta.folder, n.posterior.samples)


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

log.det <- function(M){
    c.M <- chol(M)
    q <- 2 * sum(log(diag(c.M)))
    return(q)
}

get.df.rm <- function(df, mode,...){
    if( mode == 1){
        return(df)
    } else if(mode ==2) {
        return(df-1)
    } else if(mode ==3)  {
        l <- list(...)

        return(df*(1-1/tot.n.snp))
    }

}

load.files <- function(bed.file='chr2', project.folder,
                       delta.folder,
                       tot.n.snp = tot.n.snp, df.mode =1, ...){

    indiv <- read.table(paste0(project.folder, ".order"))
    diffs <- read.table(paste0(project.folder, ".diffs"))
    diffs <- as.matrix(diffs)

    thetas <- read.table(paste0(delta.folder, "mcmcthetas.txt"))
    ll <- read.table(paste0(delta.folder, "mcmcpilogl.txt"))[,2]
    sigma2 <- thetas[,1]
    df <- thetas[,2]
    df.rm <- get.df.rm(df, df.mode, tot.n.snp=tot.n.snp)
    

    data <- load.snp.data(bed.file=bed.file, indiv=indiv, ...)
    print("loaded data")

    return(list("data"=data,
                "diffs"=diffs,
                "sigma2"=sigma2,
                "df"=df,
                "df.rm"=df.rm,
                "data"=data,
                "ll"=ll))

}

do.full.analysis <- function(bed.file='chr2', project.folder, 
                             delta.folder,
                             n.posterior.samples=100,
                             max.n.snp=200, tot.n.snp=600000,
                             df.mode=1, ...){

    f <- load.files(bed.file, project.folder, delta.folder,
                   tot.n.snp, df.mode)

    deltas <- load.posterior(delta.folder, n.posterior.samples)
    print("loaded posterior")

    precomputes <- precompute(f$diffs, deltas)
    print("precomputed matrices")

    max.n.snp <- min(max.n.snp, dim(f$data$genotypes)[2])
    fn <- function(id,...){
        snp <- get.snp.mat(id, f$data)
        if(id %% 100==1) print(id)
        lllr.wishart(snp, ...)
    }

    lllr <- sapply(1:max.n.snp, fn, diffs=f$diffs, df=f$df, df.rm=f$df.rm,
                   sigma2=f$sigma2, precomputes=precomputes,
                   n.snps=tot.n.snp)
    return(lllr)

}


precompute <- function(diffs, deltas){
    mat <- list()
    LD <- L(diffs)
    LDt <- t(LD)
    deltasl = lapply(deltas, function(Delta)-LD %*% Delta %*% LDt)
    diffsl = - LD %*% diffs %*% LDt
    mat$ldet.diffs <- log.det(diffsl)
    mat$inv.diffs <- solve(diffsl)

    mat$ldet.deltas <- sapply(deltasl, log.det)

    mat$inv.deltas <- lapply(deltasl, solve)

    mat$tr.invd.diffs <- sapply(mat$inv.deltas, trace.product, diffsl)
    mat
}
LDL <- function(d2){
    d2 <- -L(d2) %*% ( d2) %*% t(L(d2))
}

lllr.wishart <- function(snp, diffs, df, df.rm, sigma2,
                         precomputes,
                         n.snps=600000){
    mat <- precomputes
    p <- dim(mat$inv.diffs)[1]
    delta.df <- df - df.rm

    snpn <- snp/n.snps
    d2 <- LDL(diffs-snpn)
    #d2 <- as.matrix(forceSymmetric(d2))
    log.det.d2 <- log.det(d2)

    gamma.ratio <- function(df, df.rm){
        l <- length(df)
        res <- c()
        for(i in 1:l){
        res <- c(res, multivariate_gamma_function(df[i]/2, p) -
                multivariate_gamma_function(df.rm[i]/2, p))
        }
        res
    }
    gammas <- mapply(gamma.ratio, df, df.rm) #gamma term

    log.c <- log(2) * delta.df / 2 * p #2^(.) term

    # constant term from detSigma
    log.c <- log.c + p/2 * ( delta.df*log(sigma2) - 
                            df*log(df) + df.rm*log(df.rm)) 

    lllratio <- gammas + log.c  
    

    #det data term
    det.data.1 <- log.det.d2 * (df.rm - p - 1) /2 
    det.data <- mat$ldet.diffs * (df - p - 1) /2 
    lllratio <- lllratio + det.data.1 - det.data

    # detSigma term
    det.sigma <- delta.df/2 * mat$ldet.deltas
    lllratio <- lllratio + det.sigma

    #traces term
    obs.traces <- mat$tr.invd.diffs * df.rm/sigma2 /2 
    snp.traces <- sapply(mat$inv.deltas, trace.product, LDL(diffs-snpn))
    snp.traces <- snp.traces  * df/ sigma2 /2 
    lllratio <- lllratio - snp.traces + obs.traces
    #snp.traces <- sapply(mat$inv.deltas, trace.product, LDL(snpn))
    #lllratio <- lllratio + df.rm/df * (1 - snp.traces/obs.traces)

    return(lllratio)
}

dWishart2 <- function(mat, df, sigma2, snp=NULL){
    p <- nrow(mat$inv.diffs)
    if(is.null(snp)){
        ldX <- mat$ldet.diffs
        tr.Sigma.inv.X <- mat$tr.invd.diffs[[1]] * df /sigma2
    } else{
        X <- LDL(diffs - snp/600000)
        ldX <- log.det(X)
        tr.Sigma.inv.X <- trace.product(mat$inv.deltas[[1]], X) * df /sigma2
    }
    ldS <- mat$ldet.deltas[1] + log(sigma2/df) * p

    ll <- -(p+1) * ldX -
        tr.Sigma.inv.X -
        df*(ldS - ldX) -
        df * p * log(2)  
    det.data <- (df-p-1) *ldX
    print(det.data/2)
    ll <- (df-p-1) * ldX -
        tr.Sigma.inv.X -
        df*ldS -
        df * p * log(2)  
    ll <- ll/2
    ll <- ll - multivariate_gamma_function(df/2, p)
    ll

}


multivariate_gamma_function<- function(a, p){                
        x = log(pi)*p*(p-1)/4 + sum(lgamma(a-(0:(p-1))/2));      
    return (x)                                               
}                                                            
                                                             
trace.product<- function(m1, m2) sum(m1 * t(m2))

harmonic.mean <- function(a){
    1/mean(1/a)
}

load.posterior = function(eems.folder, max.matrices=NULL){
    d <- read.table(paste0(eems.folder,"Delta.txt"), nrow=1)
    w <- ncol(d)

    if( is.null(max.matrices)){
	d <- read.table(paste0(eems.folder,"Delta.txt"))
	max.matrices = nrow(d) / w
    } else {
	d <- read.table(paste0(eems.folder,"Delta.txt"), 
			nrow=w * max.matrices)
    }

    f = function(i)as.matrix(d[((i-1)*w+1):(i*w),])
    Delta = lapply(1:max.matrices, f)
    return(Delta)
}

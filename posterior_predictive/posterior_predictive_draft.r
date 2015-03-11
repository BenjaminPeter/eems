require(snpStats) #for read.plnk function



tot.n.snp = 600000


get.snp.mat <- function( snp.id, data){
    snp <- data$genotypes[,snp.id]
    snp <- as(snp, "numeric")[,1]
    snp[is.na(snp)] <- 0
    outer(snp, snp, function(x,y)(x-y)*(x-y))
}

get.snp.meta <- function( snp.id, data){
    snp <- data$map[snp.id,]
    return(c(snp$chromosome, snp$position))
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


old.do.full.analysis <- function(max.n.snp=1, tot.n.snp=600000,
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
    mat <- cbind(-1, diag(rep(1, x-1)))
    Matrix(mat, sparse=TRUE)
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

    
    proj.files <- load.proj.files(project.folder)
    delta.files <- load.delta.files(delta.folder)

    df.rm <- get.df.rm(delta.files$df, df.mode, tot.n.snp=tot.n.snp)

    data <- load.snp.data(bed.file=bed.file, indiv=proj.files$indiv, ...)
    print("loaded data")

    return(list("data"=data,
                "diffs"=proj.files$diffs,
                "sigma2"=delta.files$sigma2,
                "df"=delta.files$df,
                "df.rm"=df.rm,
                "deltas"=delta.files$deltas,
                "ll"=delta.files$ll))
}

load.proj.files <- function(project.folder){
    indiv <- read.table(paste0(project.folder, ".order"))
    diffs <- read.table(paste0(project.folder, ".diffs"))
    diffs <- as.matrix(diffs)
    return(list("indiv"=indiv,
                "diffs"=diffs))
}

load.delta.files <- function(delta.folder, max.n.samples=10000){
    thetas <- read.table(paste0(delta.folder, "/mcmcthetas.txt"))
    ll <- read.table(paste0(delta.folder, "/mcmcpilogl.txt"))[,2]
    to.keep <- ll != 0
    ll <- ll[to.keep]
    sigma2 <- thetas[to.keep,1]
    df <- thetas[to.keep,2]

    n.samples.available <- sum(to.keep)

    if(n.samples.available > max.n.samples){
        n.samples <- max.n.samples
        ll <- ll[1:n.samples]
        sigma2 <- sigma2[1:n.samples]
        df <- df[1:n.samples]
    } else{
        n.samples <- n.samples.available
    }
    deltas <- load.posterior(delta.folder, n.samples)
    print(c("loaded posterior", delta.folder))

    return(list("ll"=ll,
                "sigma2"=sigma2,
                "df"=df,
                "deltas"=deltas,
                "n.samples"=n.samples))
}


load.files.multiple <-function(bed.file='chr2', folder,
                               tot.n.snp=600000, df.mode=3, 
                               n.posterior.samples=10000, ... ){
    id <- basename(folder)
    project.folder <- paste0(folder, "/input/", id)
    output.folders <- paste0(folder, "/output/") 
    delta.folders <- list.files(output.folders, pattern="???_run*", full.names=T)

    ll=c()
    sigma2=c()
    df=c()
    deltas=list()
    n.samples.left <- n.posterior.samples
    for(d.folder in delta.folders){
        delta.files <- load.delta.files(d.folder, n.samples.left)
        ll <- c(ll, delta.files$ll)
        sigma2 <- c(sigma2, delta.files$sigma2)
        df <- c(df, delta.files$df)
        deltas <- c(deltas, delta.files$deltas)
        print(c(n.samples.left, delta.files$n.samples))
        n.samples.left <- n.samples.left - delta.files$n.samples
        print(c(n.samples.left) )
        if(n.samples.left <= 0){
            print("got max posterior samples")
            break;
        }
    }
    proj.files <- load.proj.files(project.folder)
    print(c("loaded proj", project.folder))
    data <- load.snp.data(bed.file=bed.file, indiv=proj.files$indiv, ...)
    print("loaded data")

    df.rm <- get.df.rm(df, df.mode, tot.n.snp=tot.n.snp)

    return(list("data"=data,
                "diffs"=proj.files$diffs,
                "sigma2"=sigma2,
                "df"=df,
                "df.rm"=df.rm,
                "deltas"=deltas,
                "ll"=ll))
    
}

do.folder.analysis <- function(bed.file='chr2', folder, 
                             n.posterior.samples=500,
                             max.n.snp=200, tot.n.snp=600000,
                             df.mode=3, write.freq=100,
			     normalize=F,...){
    
    f <- load.files.multiple(bed.file, folder,
                   tot.n.snp, df.mode, n.posterior.samples, ...)


    freqs <- get.freq(f$data)

    precomputes <- precompute(f$diffs, f$deltas, mode=2)
    f$deltas = NULL
    print("precomputed matrices")

    max.n.snp <- min(max.n.snp, dim(f$data$genotypes)[2])
    print(c("max.n.snp is", max.n.snp))


    if(normalize){
	    out.name <- sprintf("%s/out_%s.alln.txt", folder, bed)
	    out.name.probs <- sprintf("%s/out_%s.probsn.txt", folder, bed)
    } else {
	    out.name <- sprintf("%s/out_%s.all1.txt", folder, bed)
	    out.name.probs <- sprintf("%s/out_%s.probs1.txt", folder, bed)
    }

    print(c("out names are", out.name, out.name.probs))

    opt <- c()
    opt.meta <- c()
    append <- F
    for(i in 1:max.n.snp){
	snp <- get.snp.mat(i, f$data)   
        lllr <- lllr.wishart2(snp, df=f$df, sigma2=f$sigma2,
		      precomputes=precomputes,
		      normalize=normalize)
        opt <- rbind(opt, lllr)
        opt.meta <- rbind(opt.meta, 
                          c(freqs[i], get.snp.meta(i, f$data), mean(snp)))
        if( i %% write.freq == 0 || i == max.n.snp){
            print(i)

            probs <- apply(exp(opt), 1, harmonic.mean)
            probs <- cbind(probs, rowMeans(exp(opt)))
            probs <- cbind(probs, rowSums(opt))
            probs <- cbind(opt.meta, probs)

            write.table(probs, out.name.probs, quote=F, col.names=F, 
                        row.names=F, append=append)
            write.table(opt, out.name, quote=F, col.names=F, 
                        row.names=F, append=append)


            if(i == write.freq){
                append <- T
            }
            opt <- c()
            opt.meta <- c()
        }
    }
    return(list(opt, opt.meta))
    u
}

do.full.analysis <- function(bed.file='chr2', project.folder, 
                             delta.folder,
                             n.posterior.samples=100,
                             max.n.snp=200, tot.n.snp=600000,
                             df.mode=1, ...){
    
    f <- load.files(bed.file, project.folder, delta.folder,
                   tot.n.snp, df.mode)


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


precompute <- function(diffs, deltas, mode=2){
    mat <- list()
    LD <- L(diffs)
    LDt <- t(LD)
    diffsl = - LD %*% diffs %*% LDt
    mat$ldet.diffs <- log.det(diffsl)
    mat$p <- nrow(diffsl)

    if(mode==1){
    deltasl = lapply(deltas, function(Delta)-LD %*% Delta %*% LDt)
    mat$ldet.deltas <- sapply(deltasl, log.det)
    }

    if(mode==2){
        mat$ldiffsl <- diffsl
        mat$diffs <- diffs
    }

    mat$inv.deltas <- lapply(deltas, function(x)solve(LDL(x)))

    mat$tr.invd.diffs <- sapply(mat$inv.deltas, trace.product, diffsl)
    mat
}
LDL <- function(d2){
    d2 <- -L(d2) %*% ( d2) %*% t(L(d2))
}

lllr.wishart <- function(snp, diffs, df, df.rm, sigma2,
                         precomputes,
                         n.snps=600000){
    if(is.null(dim(snp))){stop("not valid snp matrix")}

    mat <- precomputes
    p <- dim(mat$inv.diffs)[1]
    delta.df <- df - df.rm

    snpn <- snp/n.snps
    diffs <- diffs*(1-1/n.snps)
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
    #print(sprintf("det_data:\t%0.10f\t%0.10f", det.data.1, det.data))
    lllratio <- lllratio + det.data.1 - det.data

    # detSigma term
    det.sigma <- delta.df/2 * mat$ldet.deltas
    lllratio <- lllratio + det.sigma

    #traces term
    obs.traces <- mat$tr.invd.diffs * df.rm/sigma2 /2 
    snp.traces <- sapply(mat$inv.deltas, trace.product, LDL(diffs-snpn))
    snp.traces <- snp.traces  * df/ sigma2 /2 
    #print(sprintf("traces:\t%0.10f\t%0.10f", obs.traces, snp.traces))
    lllratio <- lllratio - snp.traces + obs.traces
    #snp.traces <- sapply(mat$inv.deltas, trace.product, LDL(snpn))
    #lllratio <- lllratio + df.rm/df * (1 - snp.traces/obs.traces)

    return(lllratio)
}
lllr.wishart2 <- function(snp, df, sigma2,
                         precomputes,
                         n.snps=600000, 
                         normalize=F){
    if(is.null(dim(snp))){stop("not valid snp matrix")}

    mat <- precomputes
    p <- mat$p

    if(normalize){
        if(normalize==1){
            snp.norm <- snp / mean(snp) 
            target <- mean(precomputes$diffs)
            m1.norm <- precomputes$diffs / mean(precomputes$diffs)
            d2 <- (n.snps / (n.snps -1) ) * m1.norm - (1 / (n.snps - 1)) *snp.norm
            d2 <- d2 * target

            d2 <- LDL(d2)
            d2 <- forceSymmetric(d2)
        } else if(normalize==2){
            snp.norm <- snp - mean(snp) 
            d2 <- precomputes$diffs - (1 / (n.snps )) *snp.norm
            d2 <- LDL(d2)
            d2 <- forceSymmetric(d2)
        }

    }else{
        snpn <- snp/n.snps
        lsnpn <- LDL(snpn)
        d2 <- n.snps/(n.snps-1) * mat$ldiffsl - lsnpn
    }
    #d2 <- as.matrix(forceSymmetric(d2))
    log.det.d2 <- log.det(d2)

    #det data term
    det.data.1 <- log.det.d2 * (df - p - 1) /2 
    det.data <- mat$ldet.diffs * (df - p - 1) /2 

    #traces term
    data.traces <- mat$tr.invd.diffs * df/sigma2 /2 
    data.1.traces <- sapply(mat$inv.deltas, trace.product, d2)
    data.1.traces <- data.1.traces  * df/ sigma2 /2 

    lllratio <- det.data.1 - det.data - data.1.traces + data.traces

    return(lllratio)
}

lllr3 <- function(id, data, ...){
    snp <- get.snp.mat(id, data)
    lllr.wishart2(snp=snp, ...)
}

dWishart2 <- function(mat, df, sigma2, snp=NULL){
    p <- mat$p
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
    d <- read.table(paste0(eems.folder,"/Delta.txt"), nrow=1)
    w <- ncol(d)

    if( is.null(max.matrices)){
	d <- read.table(paste0(eems.folder,"/Delta.txt"))
	max.matrices = nrow(d) / w
    } else {
	d <- read.table(paste0(eems.folder,"/Delta.txt"), 
			nrow=w * max.matrices)
    }

    f = function(i)as.matrix(d[((i-1)*w+1):(i*w),])
    Delta = lapply(1:max.matrices, f)
    return(Delta)
}


args <- commandArgs(T)

if(length(args)>1){
    bed <- args[1]
    folder <- args[2]
    n.posterior.samples <- as.numeric(args[3])
    normalize <- T

   fff <-do.folder.analysis(bed=args[1], folder=args[2],
                         max.n.snp=100000,
                         n.posterior.samples=n.posterior.samples, min.maf=0.1, normalize=normalize) 

   out.name <- sprintf("%s/out_%s.all.txt", folder, bed)
   write.table(f, out.name, quote=F, col.names=F, row.names=F)
   probs <- apply(exp(fff), 2, harmonic.mean)
   out.name <- sprintf("%s/out_%s.probs.txt", folder, bed)
   write.table(probs, out.name, quote=F, col.names=F, row.names=F)
}

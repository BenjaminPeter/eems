require(mixAK) #dWISHART
source("/home/peterb/test/wishart/test_wishart.r")

L <- function(x){
    if (!is.null(dim(x) ))
        x = ncol(x)
    cbind(-1, diag(rep(1, x-1)))
}

J <- function(coords){
    u <- unique(coords)
    outer((10000*coords[,1]+coords[,2]), (10000*u[,1]+u[,2]), function(x,y) {
           (x==y)})   
}


wishartll <- function(Delta, diffs, sigma2, df, n_sims=NULL){
    if(is.null(n_sims)) n_sims = length(sigma2)

    LD <- L(diffs)
    LDt <- t(LD)
    LDeL = lapply(Delta, function(Delta)-LD %*% Delta %*% LDt)
    LDL = - LD %*% diffs %*% LDt

    f2 = function(i) c( dWISHART(LDL, df[i], sigma2[i]/df[i]*LDeL[[i]], log=T) , ll[2]) 
    q2 <- t(sapply(1:n_sims, f2))

}

load_posterior = function(eems_folder, max_matrices=NULL){
    d <- read.table(paste0(eems_folder,"Delta.txt"), nrow=1)
    w <- ncol(d)

    if( is.null(max_matrices)){
	d <- read.table(paste0(eems_folder,"Delta.txt"))
	max_matrices = nrow(d) / w
    } else {
	d <- read.table(paste0(eems_folder,"Delta.txt"), 
			nrow=w * max_matrices)
    }

    f = function(i)as.matrix(d[((i-1)*w+1):(i*w),])
    Delta = lapply(1:max_matrices, f)
    return(Delta)
}


test<-function(){
    ll <- read.table("mcmcpilogl.txt")
    thetas <- read.table("mcmcthetas.txt")
    sigma2 <- thetas[,1]
    df <- thetas[,2]

    n_sims = sum(ll[,1] != 0)
    

    d <- read.table("Delta.txt")
    w <- ncol(d)
    Delta <- list()
    f <- function(i)as.matrix(d[((i-1)*w+1):(i*w),])
    Delta = lapply(1:n_sims, f)
    Delta <- lapply(Delta, function(x)as.matrix(forceSymmetric(x)))
    LDeL = lapply(Delta, function(Delta)-L(Delta) %*% Delta %*% t(L(Delta)))
    LDeL <- lapply(LDeL, function(x)as.matrix(forceSymmetric(x)))
    diffs = as.matrix(read.table(
                "../barrier-schemeX-nIndiv300-nSites3000.diffs"))
    LDL <- -L(diffs) %*% diffs %*% t(L(diffs))
    f1 = function(i) c( dWishart(LDL, sigma2[i]/df[i]*LDeL[[i]], df[i]) , ll[i,2])
    f2 = function(i) c( dWISHART(LDL, df[i], sigma2[i]/df[i]*LDeL[[i]], log=T) , ll[i,2]) 
    q1 <- t(sapply(1:n_sims, f1))
        
    q2 <- t(sapply(1:n_sims, f2))

    return(list(q1,q2))
}

trace <- function(m)sum(diag(m))                             
dWishart <- function( X, Sigma, df, pseudo=F){               
                                                             
    n = nrow(X)                                              
    r <- rankMatrix(X)[1]                                    
    if(pseudo) {                                             
        q <- assert_pseudo(r, n, df)                         
    } else {                                                 
        q <- n                                               
    }                                                        
                                                             
    cholesky_sigma = chol(Sigma)                             
                                                             
    if (pseudo){                                             
        log_det_X = pseudologdet(X)                          
    } else{                                                  
        log_det_X = log.det(X)                               
    }                                                        
    log_det_S = cholesky.log.det(cholesky_sigma)             
    Sigma.inv_X = solve(Sigma) %*% X                         
                                                             
    ll = -(n+1) * log_det_X -                                
        trace(Sigma.inv_X) -                                 
        df*(log_det_S - log_det_X) -                         
        df * n * log(2)                                      
                                                             
    if (pseudo) {                                            
        ll = ll - df * (n-q) * log(pi)                       
    }                                                        
                                                             
    ll = ll / 2                                              
    ll = ll - multivariate_gamma_function(df/2, q)           
    ll                                                       
}                                                            

multivariate_gamma_function<- function(a, p){                
    x = log(pi)*p*(p-1)/4 + sum(lgamma(a-(0:(p-1))/2));      
    return (x)                                               
}                                                            
                                                             


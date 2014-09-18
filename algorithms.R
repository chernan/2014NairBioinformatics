
source('./basic_probabilistic_partitioning.R')

# Complete algorithm for partitioning with random seeds
# Default K is 2 in article
basic_algorithm <- function(data, K) {
    N=dim(data)[1]
    p=matrix(nrow=N, ncol=K)
    for(i in 1:K) {p[,i] = rbeta(N,N**-0.5,1)}
    c = (t(p) %*% data)/colSums(p)
    q=rep(1/K,K)
    
    res <- list(c=c,q=q,p=p)
    for(i in 1:20) {res <- em_basic(res$c,res$q,data)}
    
    return(res)
}

# Iterative partitioning - standard version for two classes
# Default K is 2 in article
iterative_standard_algorithm <- function(data, K) {
    N=dim(data)[1]; L=dim(data)[2]
    c=matrix(data=colMeans(data), nrow=1, ncol=L)
    flat=matrix(data= mean(data), nrow=1, ncol=L)
    q = 1

    res <- list(c=c,q=q,p=data)
    
    for (m in 1:(K-1)) {
        res$c = rbind(flat,res$c)
        res$q = c(1/m,res$q); res$q = res$q/sum(res$q)
        
        for(i in 1:20) {res <- em_basic(res$c,res$q,data)}
    }
    
    return(res)
}

# Iterative EM partitioning with untrained flat class
# Default K is 2 in article
iterative_untrained_algorithm <- function(data, K) {
    N=dim(data)[1];
    L=dim(data)[2]
    c = matrix(data=colMeans(data), nrow=1, ncol=L)
    flat = matrix(data= mean(data), nrow=1, ncol=L)
    q = 1
    
    res <- list(c=c,q=q,p=data)
    
    for (m in 1:(K-1)) {
        res$c = rbind(flat,res$c)
        res$q = c(1/m,res$q); res$q = res$q/sum(res$q)
        for(i in 1:20) {
            res$c[1,]=flat;
            res <- em_basic(res$c,res$q,data)}
    }
    
    return(res)
}
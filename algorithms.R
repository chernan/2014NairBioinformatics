# Code extracted from the supplementary material of the following Bioinformatics publication:
#     
# Nishanth Ulhas Nair,
# Sunil Kumar,
# Bernard M.E. Moret,
# and Philipp Bucher
# 
# Probabilistic partitioning methods to find significant patterns in ChIP-Seq data Bioinformatics (2014)
# 30 (17): 2406-2413 first published online May 7, 2014 
# doi:10.1093/bioinformatics/btu318 
#
# Supplementary material:
# http://bioinformatics.oxfordjournals.org/content/suppl/2014/05/07/btu318.DC1/nair2014_suppl.pdf
#
# Adapted by Celine HERNANDEZ (celine.hernandez@ens.fr) Thieffry Lab - IBENS - France
#

source('./partitioning_strategies.R')

##############################################################################################
#
# Functions with each algorithm (random seeds, iterative standard, iterative EM),
#  and *basic partitioning* : Expectation-Maximization step for basic ChIP-partitioning algorithm
#
##############################################################################################

# Complete algorithm for partitioning with random seeds
#
# Default K (number of partitions) is 2 in article
# Iteration number is 20 in article
#
random_algorithm.basic_partitioning <- function(data, K=2, iterations=20) {
    N=dim(data)[1]
    p=matrix(nrow=N, ncol=K)
    for(i in 1:K) {p[,i] = rbeta(N,N**-0.5,1)}
    c = (t(p) %*% data)/colSums(p)
    q=rep(1/K,K)
    
    res <- list(c=c,q=q,p=p)
    for(i in 1:iterations) {res <- em_basic(res$c,res$q,data)}
    
    return(res)
}

# Iterative partitioning - standard version
#
# Default K (number of partitions) is 2 in article
# Iteration number is 20 in article
#
iterative_std_algorithm.basic_partitioning <- function(data, K=2, iterations=20) {
    N=dim(data)[1]; L=dim(data)[2]
    c=matrix(data=colMeans(data), nrow=1, ncol=L)
    flat=matrix(data= mean(data), nrow=1, ncol=L)
    q = 1

    res <- list(c=c,q=q,p=data)
    
    for (m in 1:(K-1)) {
        res$c = rbind(flat,res$c)
        res$q = c(1/m,res$q); res$q = res$q/sum(res$q)
        
        for(i in 1:iterations) {res <- em_basic(res$c,res$q,data)}
    }
    
    return(res)
}

# Iterative EM partitioning with untrained flat class
#
# Default K (number of partitions) is 2 in article
# Iteration number is 20 in article
#
iterative_untrained_algorithm.basic_partitioning <- function(data, K=2, iterations=20) {
    N=dim(data)[1];
    L=dim(data)[2]
    c = matrix(data=colMeans(data), nrow=1, ncol=L)
    flat = matrix(data= mean(data), nrow=1, ncol=L)
    q = 1
    
    res <- list(c=c,q=q,p=data)
    
    for (m in 1:(K-1)) {
        res$c = rbind(flat,res$c)
        res$q = c(1/m,res$q); res$q = res$q/sum(res$q)
        for(i in 1:iterations) {
            res$c[1,]=flat;
            res <- em_basic(res$c,res$q,data)
        }
    }
    
    return(res)
}


##############################################################################################
#
# Functions with each algorithm (random seeds, iterative standard, iterative EM),
#  and variations of the EM algorithm - shape-based partitioning
#
##############################################################################################

# Complete algorithm for partitioning with random seeds
#
# Default K (number of partitions) is 2 in article
# Iteration number is 20 in article
#
random_algorithm.shape_partitioning <- function(data, K=2, iterations=20) {
    N=dim(data)[1]
    p=matrix(nrow=N, ncol=K)
    for(i in 1:K) {p[,i] = rbeta(N,N**-0.5,1)}
    c = (t(p) %*% data)/colSums(p)
    q=rep(1/K,K)
    
    res <- list(c=c,q=q,p=p)
    for(i in 1:iterations) {res <- em_shape(res$c,res$q,data)}
    
    return(res)
}

# Iterative partitioning - standard version
#
# Default K (number of partitions) is 2 in article
# Iteration number is 20 in article
#
iterative_std_algorithm.shape_partitioning <- function(data, K=2, iterations=20) {
    N=dim(data)[1]; L=dim(data)[2]
    c=matrix(data=colMeans(data), nrow=1, ncol=L)
    flat=matrix(data= mean(data), nrow=1, ncol=L)
    q = 1
    
    res <- list(c=c,q=q,p=data)
    
    for (m in 1:(K-1)) {
        res$c = rbind(flat,res$c)
        res$q = c(1/m,res$q); res$q = res$q/sum(res$q)
        
        for(i in 1:iterations) {res <- em_shape(res$c,res$q,data)}
    }
    
    return(res)
}

# Iterative EM partitioning with untrained flat class
#
# Default K (number of partitions) is 2 in article
# Iteration number is 20 in article
#
iterative_untrained_algorithm.shape_partitioning <- function(data, K=2, iterations=20) {
    N=dim(data)[1];
    L=dim(data)[2]
    c = matrix(data=colMeans(data), nrow=1, ncol=L)
    flat = matrix(data= mean(data), nrow=1, ncol=L)
    q = 1
    
    res <- list(c=c,q=q,p=data)
    
    for (m in 1:(K-1)) {
        res$c = rbind(flat,res$c)
        res$q = c(1/m,res$q); res$q = res$q/sum(res$q)
        for(i in 1:iterations) {
            res$c[1,]=flat;
            res <- em_shape(res$c,res$q,data)
        }
    }
    
    return(res)
}

##############################################################################################
#
# Functions with each algorithm (random seeds, iterative standard, iterative EM),
#  and variations of the EM algorithm : Shape-based partitioning with flips
#
##############################################################################################

# Complete algorithm for partitioning with random seeds
#
# Default K (number of partitions) is 2 in article
# Iteration number is 20 in article
# Argument q is now a matrix of probabilities with rows corresponding to classes and 
#  columns to flip states (1 when there is no flip and 2 when there is one). Equiprobability 
#  is set by default.
#
random_algorithm.shape_flip_partitioning <- function(data, K=2, q=matrix(rep(1/(K*2),K*2), ncol=2), iterations=20) {
    N=dim(data)[1]
    p=matrix(nrow=N, ncol=K)
    for(i in 1:K) {p[,i] = rbeta(N,N**-0.5,1)}
    c = (t(p) %*% data)/colSums(p)
    
    res <- list(c=c,q=q,p=p)
    for(i in 1:iterations) {res <- em_shape_flip(res$c,res$q,data)}
    
    return(res)
}

# Iterative partitioning - standard version
#
# Default K (number of partitions) is 2 in article
# Iteration number is 20 in article
#
iterative_std_algorithm.shape_flip_partitioning <- function(data, K=2, iterations=20) {
    N=dim(data)[1]; L=dim(data)[2]
    c=matrix(data=colMeans(data), nrow=1, ncol=L)
    flat=matrix(data= mean(data), nrow=1, ncol=L)
    q=1
    
    res <- list(c=c,q=q,p=data)
    
    for (m in 1:(K-1)) {
        res$c = rbind(flat,res$c)
        res$q = c(1/m,res$q); res$q = res$q/sum(res$q)
        
        for(i in 1:iterations) {res <- em_shape_flip(res$c,res$q,data)}
    }
    
    return(res)
}

# Iterative EM partitioning with untrained flat class
#
# Default K (number of partitions) is 2 in article
# Iteration number is 20 in article
#
iterative_untrained_algorithm.shape_flip_partitioning <- function(data, K=2, iterations=20) {
    N=dim(data)[1];
    L=dim(data)[2]
    c = matrix(data=colMeans(data), nrow=1, ncol=L)
    flat = matrix(data= mean(data), nrow=1, ncol=L)
    q=1
    
    res <- list(c=c,q=q,p=data)
    
    for (m in 1:(K-1)) {
        res$c = rbind(flat,res$c)
        res$q = c(1/m,res$q); res$q = res$q/sum(res$q)
        for(i in 1:iterations) {
            res$c[1,]=flat;
            res <- em_shape_flip(res$c,res$q,data)
        }
    }
    
    return(res)
}

##############################################################################################
#
# Functions with each algorithm (random seeds, iterative standard, iterative EM),
#  and variations of the EM algorithm : Shape-based partitioning with shifts
#
##############################################################################################

# Complete algorithm for partitioning with random seeds
#
# Default K (number of partitions) is 2 in article
# Iteration number is 20 in article
# max_shift is the maximal numer of shifts allowed. Should be at least 1.
#
# Argument q is now a matrix with rows corresponding to classes and columns to shift indices
#
random_algorithm.shape_shift_partitioning <- function(data, K=2, max_shift=1, iterations=20) {
    #     N=dim(data)[1]
    #     p=matrix(nrow=N, ncol=K)
    #     for(i in 1:K) {p[,i] = rbeta(N,N**-0.5,1)}
    #     c = (t(p) %*% data)/colSums(p)
    
    N=dim(data)[1]; L=dim(data)[2]
    p=matrix(nrow=N, ncol=K)
    for(i in 1:K) {p[,i] = rbeta(N,N**-0.5,1)}
    range=(max_shift+1):(L-max_shift);
    c = ((t(p) %*% data)/colSums(p))[,range] 
    
    # argument q is now a matrix with rows corresponding to classes and columns to shift indices.
    q=matrix(rep(1/max_shift*K,max_shift*K), ncol=max_shift)
    
    res <- list(c=c, q=q, p=data)
    for(i in 1:iterations) {
        res <- em_shape_shift(res$c,res$q,data)
    }
    
    return(res)
}

# Iterative partitioning - standard version
# 
# Default K (number of partitions) is 2 in article
# Iteration number is 20 in article
# max_shift is the maximal numer of shifts allowed. Should be at least 1.
# 
iterative_std_algorithm.shape_shift_partitioning <- function(data, K=2, max_shift=1, iterations=20) {
    N=dim(data)[1]; L=dim(data)[2]
    c=matrix(data=colMeans(data[,(max_shift+1):(L-max_shift)]), 
             nrow=1, ncol=L-2*max_shift)
    flat=matrix(data= mean(data), nrow=1, ncol=L-2*max_shift)
    q = array(1,dim=c(1,max_shift))
    
    res <- list(c=c,q=q,p=data)
    
    for (m in 1:(K-1)) {
        res$c = rbind(flat,res$c)
        res$q = rbind(matrix(1, ncol=max_shift),res$q); res$q = res$q/sum(res$q)
        
        for(i in 1:iterations) {res <- em_shape_shift(res$c,res$q,data)}
    }
    
    return(res)
}

# Iterative EM partitioning with untrained flat class
# 
# Default K (number of partitions) is 2 in article
# Iteration number is 20 in article
# max_shift is the maximal numer of shifts allowed. Should be at least 1.
# 
iterative_untrained_algorithm.shape_shift_partitioning <- function(data, K=2, max_shift=1, iterations=20) {
    N=dim(data)[1];
    L=dim(data)[2]
    c=matrix(data=colMeans(data[,(max_shift+1):(L-max_shift)]), 
             nrow=1, ncol=L-2*max_shift)
    flat=matrix(data= mean(data), nrow=1, ncol=L-2*max_shift)
    q = array(1,dim=c(1,max_shift))
    
    res <- list(c=c,q=q,p=data)
    
    for (m in 1:(K-1)) {
        res$c = rbind(flat,res$c)
        res$q = rbind(matrix(1, ncol=max_shift),res$q); res$q = res$q/sum(res$q)
        for(i in 1:iterations) {
            res$c[1,]=flat;
            res <- em_shape_shift(res$c,res$q,data)
        }
    }
    
    return(res)
}

##############################################################################################
#
# Functions with each algorithm (random seeds, iterative standard, iterative EM),
#  and variations of the EM algorithm : Shape-based partitioning with shifts and flips
#
##############################################################################################

# Complete algorithm for partitioning with random seeds
#
# Default K (number of partitions) is 2 in article
# Iteration number is 20 in article
# max_shift is the maximal numer of shifts allowed. Should be at least 2.
#
# q is now a 3-dimensional array with dimensions 1, 2 and 3 corresponding classes, shift indices and flip states, respectively.
#
random_algorithm.shape_shift_flip_partitioning <- function(
    data, K=2, max_shift=2, 
    q=array(1/(K*max_shift*2), dim = c(K, max_shift, 2)),     #NB: 2 stands for the 2 flip states
    iterations=20) {
    
    N=dim(data)[1]; L=dim(data)[2]
    p=matrix(nrow=N, ncol=K)
    for(i in 1:K) {p[,i] = rbeta(N,N**-0.5,1)}
    range=(max_shift+1):(L-max_shift);
    c = ((t(p) %*% data)/colSums(p))[,range] 

    res <- list(c=c,q=q,p=p)
    for(i in 1:iterations) {res <- em_shape_shift_flip(res$c,res$q,data)}
    
    return(res)
}

# Iterative partitioning - standard version
#
# Default K (number of partitions) is 2 in article
# Iteration number is 20 in article
# max_shift is the maximal numer of shifts allowed. Should be at least 2.
#
iterative_std_algorithm.shape_shift_flip_partitioning <- function(data, K=2, max_shift=2, 
                                                                  iterations=20) {
    N=dim(data)[1]; L=dim(data)[2]
    c=matrix(data=colMeans(data[,(max_shift+1):(L-max_shift)]), 
             nrow=1, ncol=L-2*max_shift)
    flat=matrix(data= mean(data), nrow=1, ncol=L-2*max_shift)
    
    #NB: 2 stands for the 2 flip states
    q = array(1,dim=c(1,max_shift,2))
    
    res <- list(c=c,q=q,p=data)
    
    library(abind)
    for (m in 1:(K-1)) {
        res$c = rbind(flat,res$c)
        res$q = abind(array(1,dim=c(1,max_shift,2)), res$q, along=1); res$q = res$q/sum(res$q)
        
        for(i in 1:iterations) {res <- em_shape_shift_flip(res$c,res$q,data)}
    }
    
    return(res)
}




# Iterative EM partitioning with untrained flat class
#
# Default K (number of partitions) is 2 in article
# Iteration number is 20 in article
# max_shift is the maximal numer of shifts allowed. Should be at least 2.
#
iterative_untrained_algorithm.shape_shift_flip_partitioning <- function(data, K=2, max_shift=2, 
                                                                        iterations=20) {
    N=dim(data)[1];
    L=dim(data)[2]
    c=matrix(data=colMeans(data[,(max_shift+1):(L-max_shift)]), 
             nrow=1, ncol=L-2*max_shift)
    flat=matrix(data= mean(data), nrow=1, ncol=L-2*max_shift)

    #NB: 2 stands for the 2 flip states
    q = array(1,dim=c(1,max_shift,2))
    
    res <- list(c=c,q=q,p=data)
    
    for (m in 1:(K-1)) {
        res$c = rbind(flat,res$c)
        res$q = abind(array(1,dim=c(1,max_shift,2)), res$q, along=1); res$q = res$q/sum(res$q)
        
        for(i in 1:iterations) {
            res$c[1,]=flat;
            res <- em_shape_shift_flip(res$c,res$q,data)
        }
    }
    
    return(res)
}


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


# Expectation-Maximization step for basic ChIP-partitioning algorithm
#
# c: a matrix containing the classes to be optimized. c[i,j] is the expected bin count value of class i at position j.
# q: a vector defining the prior probabilities of each class.
# data: a matrix containing the samples. data[i,j] is the bin count of sample i at position j.
# 
em_basic = function(c,q,data) {
    K=dim(c)[1]; N=dim(data)[1]
    l=matrix(nrow=N, ncol=K); p=matrix(nrow=N, ncol=K)
    for(i in 1:N) { for (j in 1:K) {
        l[i,j] = sum(dpois(data[i,], c[j,], log=T)) }}
    for(i in 1:N) {
        p[i,] = q*exp(l[i,]-max(l[i,])); p[i,] = p[i,]/sum(p[i,])}
    q = colMeans(p)
    c = (t(p) %*% data)/colSums(p)
    
    #c <<-c; q <<-q; p <<-p;
    return(list(c=c,q=q,p=p))
}



# Variations of the EM algorithm - shape-based partitioning
#
# c: a matrix containing the classes to be optimized. c[i,j] is the expected bin count value of class i at position j.
# q: a vector defining the prior probabilities of each class.
# data: a matrix containing the samples. data[i,j] is the bin count of sample i at position j.
# 
em_shape = function(c,q,data) {
    K=dim(c)[1]; N=dim(data)[1]
    l=matrix(nrow=N, ncol=K); p=matrix(nrow=N, ncol=K)
    for(i in 1:K) {c[i,]=c[i,]/mean(c[i,])}
    rm=rowMeans(data)
    for(i in 1:N) { for (j in 1:K) {
        l[i,j] = sum(dpois(data[i,], c[j,]*rm[i], log=T)) }}
    for(i in 1:N) {
        p[i,] = q*exp(l[i,]-max(l[i,])); p[i,] = p[i,]/sum(p[i,])}
    q = colMeans(p)
    c = (t(p) %*% data)/colSums(p)
    
    #c <<-c; q <<-q; p <<-p;
    return(list(c=c,q=q,p=p))
}

# Shape-based partitioning with flips
#
# c: a matrix containing the classes to be optimized. c[i,j] is the expected bin count value of class i at position j.
# q: a matrix with rows corresponding to classes and columns to flip states (1: no flip, 2: flip)
# data: a matrix containing the samples. data[i,j] is the bin count of sample i at position j.
# 
em_shape_flip = function(c,q,data) {
    K=dim(c)[1]; N=dim(data)[1]
    l=array(dim=c(N,K,2)); p=array(dim=c(N,K,2))
    for(i in 1:K) {c[i,]=c[i,]/mean(c[i,])}
    rm=rowMeans(data)
    for(i in 1:N) {
        for(j in 1:K) {l[i,j,1] = sum(dpois(data[i,], c[j,] *rm[i],log=T))}
        for(j in 1:K) {l[i,j,2] = sum(dpois(data[i,],rev(c[j,])*rm[i],log=T))}
    }
    for(i in 1:N) {
        p[i,,] = q*exp(l[i,,]-max(l[i,,])); p[i,,] = p[i,,]/sum(p[i,,])}
    q = apply(p, c(2,3), mean)
    c = (t(p[,,1]) %*% data) + (t(p[,,2]) %*% t(apply(data,1,rev)))
    c = c / apply(p,2,sum)
    
    #c <<-c; q <<-q; p <<-p;
    return(list(c=c,q=q,p=p))
}

# Shape-based partitioning with shifts
#
# c, exp_bin_counts: a matrix containing the classes to be optimized. c[i,j] is the expected bin count value of class i at position j.
# q: a matrix with rows corresponding to classes and columns to shift indices
# data: a matrix containing the samples. data[i,j] is the bin count of sample i at position j.
# 
em_shape_shift <- function(c, q, data) {
    K = dim(c)[1]; # nb of classes
    L = dim(c)[2]; # length of samples
    N = dim(data)[1]; # nb of samples
    S = dim(q)[2] # ncol shifts
    l = array(dim=c(N,K,S)); p=array(dim=c(N,K,S))
    for(i in 1:K) {
        c[i,]=c[i,]/mean(c[i,])
    }
    rm=matrix(nrow=N, ncol=S)
    for(k in 1:S) {rm[,k] = rowMeans(data[,k:(k+L-1)])}
    for(i in 1:N) { for (j in 1:K) { for (k in 1:S) {
        l[i,j,k]=sum(dpois(data[i,k:(k+L-1)], c[j,] *rm[i,k],log=T)) }}}
    for(i in 1:N) {
        p[i,,] = q*exp(l[i,,] - max(l[i,,])); p[i,,] = p[i,,]/sum(p[i,,])}
    q = apply(p, c(2,3), mean)
    c = 0; for(k in 1:S) {
        c = c + (t(p[,,k]) %*% data[,k:(k+L-1)])/colSums(p[,,k])}
    c = c / apply(p,2,sum)
    m=sum((1:S)*colSums(q)); s=sum(((1:S)-m)**2*colSums(q))**0.5
    for (i in 1:K) {
        q[i,] = sum(q[i,]) * dnorm(1:S,floor(S/2)+1,s) /
            sum(dnorm(1:S,floor(S/2)+1,s))
    }
    
    #c <<-c; q <<-q; p <<-p;
    return(list(c=c,q=q,p=p))
}

# Shape-based partitioning with flips and shifts
#
# c: a matrix containing the classes to be optimized. c[i,j] is the expected bin count value of class i at position j.
# q: a 3-dimensional array with dimensions 1, 2 and 3 corresponding classes, shift indices and flip states, respectively.
# data: a matrix containing the samples. data[i,j] is the bin count of sample i at position j.
# 
em_shape_shift_flip = function(c,q,data) {
    K=dim(c)[1]; L=dim(c)[2]; N=dim(data)[1]; S=dim(q)[2]
    l=array(dim=c(N,K,S,2)); p=array(dim=c(N,K,S,2))
    for(i in 1:K) {c[i,]=c[i,]/mean(c[i,])}
    rm=matrix(nrow=N, ncol=S)
    for(k in 1:S) {rm[,k] = rowMeans(data[,k:(k+L-1)])}
    for(i in 1:N) { for (j in 1:K) { for (k in 1:S) {
        l[i,j,k,1]=sum(dpois(
            data[i,k:(k+L-1)],c[j,]*rm[i,k],log=T))
        l[i,j,k,2]=sum(dpois(
            rev(data[i,])[k:(k+L-1)],c[j,]*rm[i,S-k+1],log=T))
    }}}
    for(i in 1:N) {
        p[i,,,] = q*exp(l[i,,,]-max(l[i,,,]))
        p[i,,,] = p[i,,,]/sum(p[i,,,])
    }
    q = apply(p, c(2,3,4), mean)
    c = 0; for(k in 1:S) {
        c = c + 
            (t(p[,,k,1]) %*% data[,k:(k+L-1)]) +
            (t(p[,,k,2]) %*% t(apply(data[,(S-k+1):(S+L-k)],1,rev)))
    }
    c = c / apply(p,2,sum)
    m= sum((1:S)*apply(q,2,sum));
    s=sum(((1:S)-m)**2*apply(q,2,sum))**0.5
    for (i in 1:K) {
        q[i,,1] = q[i,,2] = sum(q[i,,]) * dnorm(1:S,floor(S/2)+1,s) /
            sum(dnorm(1:S,floor(S/2)+1,s)) / 2
    }
    
    #c <<-c; q <<-q; p <<-p;
    return(list(c=c,q=q,p=p))
}

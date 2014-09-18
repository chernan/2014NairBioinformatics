# http://bioinformatics.oxfordjournals.org/content/suppl/2014/05/07/btu318.DC1/nair2014_suppl.pdf

# Expectation-Maximization step for basic ChIP-partitioning algorithm
#
# c: a matrix containing the classes to be optimized. c[i,j] is the expected bin count value of class i at position j.
# q: a vector defining the prior probabilities of each class.
# data: a matrix containing the samples. data[i,j] is the bin count of sample i at position j.
em_basic = function(c,q,data) {
    K=dim(c)[1]; N=dim(data)[1]
    l=matrix(nrow=N, ncol=K); p=matrix(nrow=N, ncol=K)
    for(i in 1:N) { for (j in 1:K) {
        l[i,j] = sum(dpois(data[i,], c[j,], log=T)) }}
    for(i in 1:N) {
        p[i,] = q*exp(l[i,]-max(l[i,])); p[i,] = p[i,]/sum(p[i,])}
    q = colMeans(p)
    c = (t(p) %*% data)/colSums(p)
    c <<-c; q <<-q; p <<-p;
}


# Complete algorithm for partitioning with random seeds
K=2; N=dim(data)[1]
p=matrix(nrow=N, ncol=K)
for(i in 1:K) {p[,i] = rbeta(N,N**-0.5,1)}
c = (t(p) %*% data)/colSums(p)
q=rep(1/K,K)

for(i in 1:20) {em_basic(c,q,data)}

# Iterative partitioning - standard version for two classes
K=2; N=dim(data)[1]; L=dim(data)[2]
c=matrix(data=colMeans(data), nrow=1, ncol=L)
flat=matrix(data= mean(data), nrow=1, ncol=L)
q = 1
for (m in 1:(K-1)) {
    c = rbind(flat,c)
    q = c(1/m,q); q = q/sum(q)
    for(i in 1:20) {em_basic(c,q,data)}
}

# Iterative EM partitioning with untrained flat class
K=2; N=dim(data)[1];
L=dim(data)[2]
c = matrix(data=colMeans(data), nrow=1, ncol=L)
flat = matrix(data= mean(data), nrow=1, ncol=L)
q = 1
for (m in 1:(K-1)) {
    c = rbind(flat,c)
    q = c(1/m,q); q = q/sum(q)
    for(i in 1:20) {
        c[1,]=flat;
        em_basic(c,q,data)}
}






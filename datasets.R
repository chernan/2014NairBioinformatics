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

##############################################################################################
#
# Simulated datasets
#
##############################################################################################


# Two classes without shifts and flips
# 
# coverage parameter f varied from 0.5 to 50 to generate data sets of different coverage
create_dataset1 <- function(f) {
    n_samples = 1000
    class1_m = 60
    class1_s = 10
    class2_m = 75
    class2_s = 3

    data1 = matrix(NA,nrow=n_samples, ncol=100)
    for(sample in 1:n_samples) {
        lambda = dnorm(1:100, class1_m, class1_s)
        lambda = f*lambda/sum(lambda)
        data1[sample,] = rpois(100,lambda)
    }
    data2 = matrix(NA,nrow=n_samples, ncol=100)
    for(sample in 1:n_samples) {
        lambda = dnorm(1:100, class2_m, class2_s)
        lambda = f*
            lambda/sum(lambda)
        data2[sample,] = rpois(100,lambda)
    }
    data = rbind( data1, data2)
    
    return(data)
}

# Two classes with flips
#
# Code to generate the data for the experiments summarized by Figure 2 and Table 3 of the main text.
create_dataset2 <- function() {
    n_samples = 1000
    class1_m = 40
    class1_s = 1
    class2_m = 70
    class2_s = 5
    f = 50
    data1 = matrix(NA,nrow=n_samples, ncol=100)
    for(sample in 1:n_samples) {
        lambda = dnorm(1:100, class1_m, class1_s)
        lambda = f*lambda/sum(lambda)
        data1[sample,] =
            rpois(100,lambda)
    }
    data2 = matrix(NA,nrow=n_samples, ncol=100 )
    for(sample in 1:n_samples) {
        lambda = dnorm(1:100, class1_m, class1_s)
        lambda = f*lambda/sum(lambda)
        data2[sample,] = rev(rpois(100,lambda))
    }
    data3 = matrix(NA,nrow=n_samples, ncol=100)
    for(sample in 1:n_samples) {
        lambda = dnorm(1:100, class2_m, class2_s)
        lambda = f*lambda/sum(lambda)
        data3[sample,] = rpois(100,lambda)
    }
    data4 = matrix(NA,nrow=n_samples, ncol=100)
    for(sample in 1:n_samples) {
        lambda =
            dnorm(1:100, class2_m, class2_s)
        lambda = f*lambda/sum(lambda)
        data4[sample,] = rev(rpois(100,lambda))
    }
    data = rbind( data1, data2, data3, data4)
    
    return(data)
}

# Four classes
#
# Code to generate the data for the experiments summarized by Figure S1 in Supplementary material.
create_dataset3 <- function() {
    n_samples = 1000
    class1_m = 40
    class1_s = 3
    class2_m = 40
    class2_s = 10
    class3_m = 70
    class3_s = 2
    class4_m = 60
    class4_s = 5
    f = 5
    data1 = matrix(NA,nrow=n_samples, ncol=100)
    for(sample in 1:n_samples) {
        lambda = dnorm(1:100, class1_m, class1_s)
        lambda = f*lambda/sum(lambda)
        data1[sample,] = rpois(100,lambda)
    }
    data2 = matrix(NA,nrow=n_samples, ncol=100)
    for(sample in 1:n_samples) {
        lambda = dnorm(1:100, class2_m, class2_s)
        lambda = f*lambda/sum(lambda)
        data2[sample,] = rpois(100,lambda)
    }
    data3 = matrix(NA,nrow=n_samples, ncol=100)
    for(sample in 1:n_samples) {
        lambda = dnorm(1:100, class3_m, class3_s)
        lambda = f*lambda/sum(lambda)
        data3[sample,] = rpois(100,lambda)
    }
    data4 = matrix(NA,nrow=n_samples, ncol=100)
    for(sample in 1:n_samples) {
        lambda = dnorm(1:100, class4_m, class4_s)
        lambda = f*lambda/sum(lambda)
        data4[sample,] = rpois(100,lambda)
    }
    data = rbind( data1, data2, data3, data4)
    
    return(data)
}

# Two classes characterized by co-localizing peaks of different width
#
# Code to generate the data for the experiments summarized by Figure S2 in Supplementary material.
create_dataset4 <- function() {
    n_samples = 1000
    class1_m = 60
    class1_s = 10
    class2_m = 60
    class2_s = 3
    f = 50
    data1 = matrix(NA,nrow=n_samples, ncol=100)
    for(sample in 1:n_samples) {
        lambda = dnorm(1:100, class1_m, class1_s)
        lambda = f*lambda/sum(lambda)
        data1[sample,] = rpois(100,lambda)
    }
    data2 = matrix(NA,nrow=n_samples, ncol=100)
    for(sample in 1:n_samples)
    {
        lambda = dnorm(1:100, class2_m, class2_s)
        lambda = f*lambda/sum(lambda)
        data2[sample,] = rpois(100,lambda)
    }
    data = rbind( data1, data2)
    
    return(data)
}

##############################################################################################
#
# ChIP-seq dataset : Kundaje et al. 2012
#
##############################################################################################

# "Parallel application of probabilistic partitioning and CAGT to the same data set
# We tested probabilistic partitioning on an example that was used in (Kundaje et al. 2012, 
# Genome Res. 22:1735) for introduction and illustration of the CAGT method. The data consist 
# of H3K27ac signal profiles around CTCF binding site in the K562 cell line. For probabilistic 
# partitioning, ChIP- Seq tags were centered by 90 bp and counted in bins of 10 bp. We then 
# used the shape-based version with flips but without shifts. For CAGT, we normalized the 
# ChIP-Seq data as described in the afore- mentioned paper. Of a total of 67’986 samples, we a 
# priori eliminated about two third because of low coverage. Of the remaining samples, CAGT 
# automatically rejected another half, retaining 10226 samples for classification. Probabilistic 
# partitioning was applied to the same subset of samples that was retained by CAGT. Both methods 
# were forced to return exactly 10 classes."

# They don't provide the code to do the counting.


##############################################################################################
#
# ExoProfiler dataset : read a count file as generated by ExoProfiler
#
##############################################################################################

# file : file path
# header : boolean, if the file contains a header line
# row.names : index of column containing row names, starts at 1
read_count_dataset <- function(file, header=FALSE, row.names=1) {
    counts <- as.matrix(read.table(file, header=header, row.names=row.names))
    return(counts)    
}

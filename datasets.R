
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






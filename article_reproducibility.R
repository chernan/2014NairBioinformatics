
source('./datasets.R')
source('./algorithms.R')

color_lines <- c('red', 'blue', 'green', 'orange')

##############################################################################################
#
# Simulated datasets
#
##############################################################################################


##############################
# Simulated dataset 1 : 2 classes

data <- create_dataset1(50)

# Basic partitioning with random seeds
nb_partitions <- 2
res <- random_algorithm.basic_partitioning(data, 2)

# Display results : contour plots
contour(t(data))
for(partition_index in 1:nb_partitions) {
    contour(t(data[res$p[, partition_index]>1/nb_partitions, ]), col=color_lines[partition_index])
}

# Display results : aggregation plot/footprint
plot(colSums(res$c), type='l')
for(partition_index in 1:nb_partitions) {
    lines(res$c[partition_index, ], type = 'l', col=color_lines[partition_index])
}

# # Display results : aggregation plot/footprint using result probabilities
# plot(colSums(data), type='l')
# for(partition_index in 1:nb_partitions) {
#     lines(colSums(data[res$p[, partition_index]>1/nb_partitions,]), type = 'l', col=color_lines[partition_index])
# }


##############################
# Simulated dataset 1 : 2 classes - less coverage (smaller parameter for create_dataset1() )

data <- create_dataset1(0.5)

# Basic partitioning with random seeds
nb_partitions <- 2
res <- random_algorithm.basic_partitioning(data, 2)

# Display results : contour plots
contour(t(data))

# Display results : aggregation plot/footprint
plot(colSums(res$c), type='l')
for(partition_index in 1:nb_partitions) {
    lines(res$c[partition_index, ], type = 'l', col=color_lines[partition_index])
}



##############################
# Simulated dataset 2 : 2 classes - with flip

data <- create_dataset2()

# Display dataset : contour plots
contour(t(data))

nb_partitions <- 2

# Basic partitioning with random seeds
res <- random_algorithm.basic_partitioning(data, nb_partitions)
plot(colSums(res$c), type='l')
for(partition_index in 1:nb_partitions) {
    lines(res$c[partition_index, ], type = 'l', col=color_lines[partition_index])
}


# Flip partitioning, 
nb_flip_states <- 2

# 2 classes, each one also flipped
# q: a matrix with rows corresponding to classes and columns to flip states
res <- random_algorithm.shape_flip_partitioning(data, nb_partitions,
                                                #q: 1/4, 1/4, 1/4, 1/4
                                                #class 1: 50% of total in which 50% not-flipped, 50% flipped
                                                #class 2: 50% of total in which 50% not-flipped, 50% flipped
                                                q=matrix(c(1/(nb_partitions*nb_flip_states),1/(nb_partitions*nb_flip_states),
                                                           1/(nb_partitions*nb_flip_states),1/(nb_partitions*nb_flip_states)), ncol=2))
plot(colSums(data), type='l') 
for(partition_index in 1:nb_partitions) {
    lines(colSums(data[res$p[, partition_index, 1]>1/nb_partitions,]), type = 'l', col=color_lines[partition_index])
    lines(colSums(data[res$p[, partition_index, 2]>1/nb_partitions,]), type = 'l', col=color_lines[partition_index])
}

# Just a try... 3 classes, only one flipped
# q: a matrix with rows corresponding to classes and columns to flip states
nb_partitions <- 3
res <- random_algorithm.shape_flip_partitioning(data, nb_partitions,
                                                #q: 1/4, 1/4, 1/4, 1/4
                                                #class 1: 33% of total in which 50% not-flipped, 50% flipped
                                                #class 2: 33% of total in which 100% not-flipped, 0% flipped
                                                #class 2: 33% of total in which 100% not-flipped, 0% flipped
                                                q=matrix(c(.165, .165,
                                                           .33, 0,
                                                           .33, 0), ncol=2))
plot(colSums(data), type='l') 
for(partition_index in 1:nb_partitions) {
    lines(colSums(data[res$p[, partition_index, 1]>1/nb_partitions,]), type = 'l', col=color_lines[partition_index])
    lines(colSums(data[res$p[, partition_index, 2]>1/nb_partitions,]), type = 'l', col=color_lines[partition_index])
}



##############################
# Simulated dataset 3 : 4 classes

data <- create_dataset3()

# Iterative standard partitioning
nb_partitions <- 4
res <- iterative_std_algorithm.basic_partitioning(data, nb_partitions)

# Display dataset : contour plots
contour(t(data))

# Display dataset : aggregation plots
plot(colSums(data), type='l')
lines(colSums(data[1:1000,]), type = 'l', col='red')
lines(colSums(data[1001:2000,]), type = 'l', col='blue')
lines(colSums(data[2001:3000,]), type = 'l', col='green')
lines(colSums(data[3001:4000,]), type = 'l', col='orange')

# Display results : aggregation plot/footprint
plot(colSums(res$c), type='l')
for(partition_index in 1:nb_partitions) {
    lines(res$c[partition_index, ], type = 'l', col=color_lines[partition_index])
}


##############################
# Simulated dataset 4 : 2 classes characterized by co-localizing peaks of different width

data <- create_dataset4()

# Basic partitioning with random seeds
nb_partitions <- 2
res <- random_algorithm.basic_partitioning(data, nb_partitions)

# Display dataset : contour plots
contour(t(data))

# Display results : aggregation plot/footprint
plot(colSums(res$c), type='l')
for(partition_index in 1:nb_partitions) {
    lines(res$c[partition_index, ], type = 'l', col=color_lines[partition_index])
}


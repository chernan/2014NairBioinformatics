
source('./datasets.R')
source('./algorithms.R')

color_lines <- c('red', 'blue', 'green', 'orange')

# 2 classes

data <- create_dataset1(50)
res <- basic_algorithm(data, 2)

contour(t(data))
contour(t(data[res$p[,1]>0.5,]))
contour(t(data[res$p[,2]>0.5,]))

plot(colSums(data), type='l')
lines(colSums(data[res$p[,1]>0.5,]), type = 'l', col='red')
lines(colSums(data[res$p[,2]>0.5,]), type = 'l', col='blue')


# 2 classes - less coverage

data <- create_dataset1(0.5)
res <- basic_algorithm(data, 2)

contour(t(data))
contour(t(data[res$p[,1]>0.5,]))
contour(t(data[res$p[,2]>0.5,]))

plot(colSums(res$c), type='l')
lines(colSums(res$c[1,]), type='l', col='red')
lines(colSums(res$c[2,]), type='l', col='blue')

plot(colSums(data), type='l')
lines(colSums(data[res$p[,1]>0.5,]), type = 'l', col='red')
lines(colSums(data[res$p[,2]>0.5,]), type = 'l', col='blue')


# 4 classes

data <- create_dataset3()
res <- iterative_standard_algorithm(data, 4)

contour(t(data))

plot(colSums(data), type='l')
lines(colSums(data[1:1000,]), type = 'l', col='red')
lines(colSums(data[1001:2000,]), type = 'l', col='blue')
lines(colSums(data[2001:3000,]), type = 'l', col='green')
lines(colSums(data[3001:4000,]), type = 'l', col='orange')

plot(colSums(res$c), type='l')
for(i in 1:4) {
    lines(res$c[i,], type = 'l', col=color_lines[i])
}

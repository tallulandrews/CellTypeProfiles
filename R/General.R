factor_counts <- function(labels) {
        labels <- factor(labels)
        x <- split(seq(length(labels)), labels)
        result <- sapply(x, function(a) length(a))
        return(result);
}

my_row_mean_aggregate <- function(mat, groups) {
        # Much faster version
        MAT <- as.matrix(mat)
        x <- split(seq(ncol(MAT)), groups)
        result <- sapply(x, function(a) rowMeans(MAT[,a]))

        return(result);
}

cosine_dist <- function(x,y) {x %*% y / sqrt(x%*%x * y%*%y)}

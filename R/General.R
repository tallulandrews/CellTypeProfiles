factor_counts <- function(vec) {
        vec <- factor(vec)
        x <- split(seq(length(vec)), vec)
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

my_row_var_aggregate <- function(mat, groups) {
        # Much faster version
        MAT <- as.matrix(mat)
        x <- split(seq(ncol(MAT)), groups)
        result <- sapply(x, function(a) matrixStats::rowVars(MAT[,a]))

        return(result);
}

#cosine_dist <- function(x,y) {x %*% y / sqrt(x%*%x * y%*%y)}

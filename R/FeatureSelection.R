# Generalizable - Step 0 = feature selection functions
kw.features <- function(expr_mat, clusters) {
        p.value <- apply(expr_mat, 1, function(x){ kruskal.test(x,clusters)$p.value})
        p.value[is.na(p.value )] = 1;
        q.value <- p.adjust(p.value, "bon")
        return(as.numeric(q.value < 0.05))
}

no.feature.selection <- function(expr_mat, clusters) { return(rep(1, times=length(expr_mat[,1]))) }

m3drop.features <- function(expr_mat, clusters) {
        require("M3Drop")
        features <- M3DropFeatureSelection(expr_mat, suppress.plot=TRUE);
        f_names <- rownames(features[features$q.value < 0.05,])
        return(as.numeric(rownames(expr_mat) %in% f_names))
} 

marker.features <- function(expr_mat, clusters) {
        # To Do
}


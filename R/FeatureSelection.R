# Generalizable - Step 0 = feature selection functions
kw.features <- function(expr_mat, clusters) {
        p.value <- apply(expr_mat, 1, function(x){ stats::kruskal.test(x,clusters)$p.value})
        p.value[is.na(p.value )] = 1;
        q.value <- p.adjust(p.value, "bon")
        return(as.numeric(q.value < 0.05))
}

no.feature.selection <- function(expr_mat, clusters) { return(rep(1, times=length(expr_mat[,1]))) }

m3drop.features <- function(expr_mat, clusters) {
#        require("M3Drop")
        features <- M3Drop::M3DropFeatureSelection(expr_mat, suppress.plot=TRUE);
        f_names <- rownames(features[features$q.value < 0.05,])
        return(as.numeric(rownames(expr_mat) %in% f_names))
} 

danb.features <- function(expr_mat, clusters) {
#        require("M3Drop")
	counts <- M3Drop::NBumiConvertToInteger(expr_mat);
	fit <- M3Drop::NBumiFitModel(counts);
	topX <- ceiling(0.05*length(counts[,1]))
	f_names <- names(M3Drop::NBumiFeatureSelectionCombinedDrop(fit)[1:topX])
        return(as.numeric(rownames(expr_mat) %in% f_names))
}

hvg.features <- function (expr_mat, clusters) {
#        require("M3Drop")
        features <- M3Drop::BrenneckeGetVariableGenes(expr_mat, suppress.plot=TRUE, fdr=0.05);
        f_names <- rownames(features[features$q.value < 0.05,])
        return(as.numeric(rownames(expr_mat) %in% f_names))
}

marker.features <- function(expr_mat, clusters) {
	markers = complex_markers(expr_mat, clusters, n_max=length(unique(clusters))-1, strict_only=FALSE)
	auc_filter = markers[,1] > 0.7;
	sig_filter = markers$q.value < 0.05;
	return(as.numeric(auc_filter & sig_filter));
}

marker.strict.features <- function(expr_mat, clusters) {
	markers = complex_markers(expr_mat, clusters, n_max=length(unique(clusters))-1, strict_only=TRUE)
	auc_filter = markers[,1] > 0;
	return(as.numeric(auc_filter));
}

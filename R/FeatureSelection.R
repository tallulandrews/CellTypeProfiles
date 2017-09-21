# Generalizable - Step 0 = feature selection functions
get_assignment_matrix <- function(expr_mat, clusters) {
	meds <- apply(expr_mat, 1, mean)
	high <- expr_mat > meds
	freq <- rowSums(high)/ncol(high)
	cluster_assign <- my_row_mean_aggregate(high, clusters);
	cluster_assign <- apply(cluster_assign > 0.5, 2, as.numeric);
	cluster_assign[freq > 0.5,] <- apply(cluster_assign[freq > 0.5,] < 0.5, 2, as.numeric); # inverse for negative markers
	return(cluster_assign);
}

kw.features <- function(expr_mat, clusters) {
	clusters <- factor(clusters);
        p.value <- apply(expr_mat, 1, function(x){ stats::kruskal.test(x,clusters)$p.value})
        p.value[is.na(p.value )] = 1;
        q.value <- p.adjust(p.value, "bon")
	f_names <- rownames(expr_mat)[q.value < 0.05];


	assignment <- get_assignment_matrix(expr_mat, clusters);
	rownames(assignment) <- rownames(expr_mat);
	f_rows = rownames(assignment) %in% f_names;
	assignment[!f_rows,] <- matrix(0, nrow=sum(!f_rows), ncol=length(unique(clusters)))

	return(assignment)
#        return(as.numeric(q.value < 0.05))
}

no.feature.selection <- function(expr_mat, clusters) { return(rep(1, times=length(expr_mat[,1]))) }

# General FS methods

m3drop.features <- function(expr_mat, clusters) {
#        require("M3Drop")
        features <- M3Drop::M3DropFeatureSelection(expr_mat, suppress.plot=TRUE);
        f_names <- rownames(features[features$q.value < 0.05,])

#	assignment <- get_assignment_matrix(expr_mat, clusters);
#	rownames(assignment) <- rownames(expr_mat);
#	f_rows = rownames(assignment) %in% f_names;
#	assignment[!f_rows,] <- matrix(0, nrow=sum(!f_rows), ncol=length(unique(clusters)))

	nclust <- length(unique(clusters))
	is.f <- as.numeric(rownames(expr_mat) %in% f_names)
	assignement <- matrix(rep(is.f, times=nclust), ncol=nclust)

	return(assignment)
#        return(as.numeric(rownames(expr_mat) %in% f_names))
} 

danb.features <- function(expr_mat, clusters) {
#        require("M3Drop")
	counts <- M3Drop::NBumiConvertToInteger(expr_mat);
	fit <- M3Drop::NBumiFitModel(counts);
	topX <- ceiling(0.05*length(counts[,1]))
	f_names <- names(M3Drop::NBumiFeatureSelectionCombinedDrop(fit)[1:topX])

#	assignment <- get_assignment_matrix(expr_mat, clusters);
#	rownames(assignment) <- rownames(expr_mat);
#	f_rows = rownames(assignment) %in% f_names;
#	assignment[!f_rows,] <- matrix(0, nrow=sum(!f_rows), ncol=length(unique(clusters)))

	nclust <- length(unique(clusters))
	is.f <- as.numeric(rownames(expr_mat) %in% f_names)
	assignement <- matrix(rep(is.f, times=nclust), ncol=nclust)

	return(assignment)
}

hvg.features <- function (expr_mat, clusters) {
#        require("M3Drop")
        features <- M3Drop::BrenneckeGetVariableGenes(expr_mat, suppress.plot=TRUE, fdr=0.05);
        f_names <- rownames(features[features$q.value < 0.05,])

#	assignment <- get_assignment_matrix(expr_mat, clusters);
#	rownames(assignment) <- rownames(expr_mat);
#	f_rows = rownames(assignment) %in% f_names;
#	assignment[!f_rows,] <- matrix(0, nrow=sum(!f_rows), ncol=length(unique(clusters)))

	nclust <- length(unique(clusters))
	is.f <- as.numeric(rownames(expr_mat) %in% f_names)
	assignement <- matrix(rep(is.f, times=nclust), ncol=nclust)

	return(assignment)
}

marker.features <- function(expr_mat, clusters) {
	# returns matrix of gene->cluster
	markers = complex_markers(expr_mat, clusters, n_max=length(unique(clusters))-1, strict_only=FALSE)
	auc_filter = markers$AUC > 0.7;
	sig_filter = markers$q.value < 0.05;
	ngroups = length(unique(clusters))
	out_mat <- markers[,(1:ngroups)+1]
	exclude <- !(auc_filter & sig_filter)
	out_mat[exclude,] <- matrix(0, nrow=sum(exclude), ncol=ngroups)
	on_freq <- rowSums(out_mat)/ncol(out_mat);
	out_mat[on_freq > 0.5,] <- abs(out_mat[on_freq > 0.5,]-1) # invert negative markers
	#return(as.numeric(auc_filter & sig_filter));
	return(out_mat);
}

marker.strict.features <- function(expr_mat, clusters) {
	# returns matrix of gene->cluster
	markers = complex_markers(expr_mat, clusters, n_max=length(unique(clusters))-1, strict_only=TRUE)
	auc_filter = markers[,1] > 0;
	ngroups = length(unique(clusters))
	out_mat <- markers[,(1:ngroups)+1]
	on_freq <- rowSums(out_mat)/ncol(out_mat);
	out_mat[on_freq > 0.5,] <- abs(out_mat[on_freq > 0.5,]-1) # invert negative markers
	#return(as.numeric(auc_filter));
	return(out_mat);
}

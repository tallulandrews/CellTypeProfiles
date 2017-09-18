clade_o_score <- function(heatout, ann) {
	require("PhyloMeasures") #http://onlinelibrary.wiley.com/doi/10.1111/ecog.01814/full
	require("ape")

	#ape::is.monophyletic

	true_levels = levels(factor(ann));
	my_tree = heatout$heatmap_out$colDendrogram
	my_phylo = ape::as.phylo.hclust(stats::as.hclust(my_tree))
	my_truth_matrix <- matrix(0, nrow=length(true_levels), ncol=dim(heatout$dist_mat)[1])
	colnames(my_truth_matrix) = colnames(heatout$dist_mat)
	rownames(my_truth_matrix) = true_levels
	monophy = vector(length=length(true_levels));
	for (i in 1:length(true_levels)) {
		my_truth_matrix[i,ann == true_levels[i]] = 1
		# Monophyletic
		monophy[i]= ape::is.monophyletic(my_phylo, colnames(my_truth_matrix)[ann==true_levels[i]])
	}  
	# lower score = better, negative means better than random chance.
	pheno_div = pd.query(my_phylo, my_truth_matrix, standardize=FALSE, seed=103)
	mpd = mpd.query(my_phylo, my_truth_matrix, standardize=FALSE, seed=103)

	score = mpd; score[monophy] = 0; score = sum(score);	

	return(list(table=cbind(pheno_div, mpd, monophy, rowSums(my_truth_matrix)), score=score))
}

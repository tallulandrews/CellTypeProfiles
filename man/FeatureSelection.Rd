\name{FeatureSelection}
\alias{kw.features}
\alias{m3drop.features}
\alias{no.feature.selection}
\alias{marker.features}
\alias{danb.features}
\alias{hvg.features}
\title{Calculate biological features for a dataset}
\description{
	Determines whether or not each gene is a relevant feature of the dataset, i.e. is responding to biological factors not just technical noise.
}
\usage{
	kw.features(expr_mat, clusters)
	m3drop.features(expr_mat, clusters)
	danb.features(expr_mat, clusters)
	hvg.features(expr_mat, clusters)
	marker.features(expr_mat, clusters)
	no.feature.selection(expr_mat, clusters)
}
\arguments{
  \item{expr_mat}{a numeric matrix of expression values. Columns = cells, Rows = genes.}
  \item{clusters}{a vector of cluster assignments for each cell.}
}
\details{
	These functions are to be passed to the "get_cluster_profile" function. They return a vector of 1s and 0s indicating whether each gene was identified as a feature using the respective method. m3drop is recommended for very large datasets, being both fast and consistent across datasets. marker.features will select genes giving the clearest separation of clusters but is also the slowest.

	\code{kw.features} uses the kruskal-wallis test (kruskal.test) to find DE genes between the clusters then applies a Bonferroni multiple-testing correction (p.adjust). Genes with q-values < 0.05 are deemed features of the dataset. Requires multiple cell-types in the dataset.

	\code{m3drop.features} uses M3Drop feature selection with a 5\% FDR multiple testing correction. 

	\code{danb.features} uses DANB dropout-based feature selection (see: M3Drop package) and selects the top 5\% of genes.

	\code{hvg.features} uses highly variable gene feature selection as described in Brennecke et al. (2013) implementated in the M3Drop package. Uses default parameters with 5\% FDR multiple testing correction.

	\code{marker.features} uses the complex_markers function to find markers of any set of clusters (n_max = n-1). Features are selected with AUC > 0.7 and q.value < 0.05 (Bonferroni). Uses strict_only=FALSE. Requires multiple cell-types in the dataset.

	\code{no.feature.selection} keeps all genes. 

}
\value{
	vector of 1s and 0s for each gene in the expression matrix, indicating whether that gene is a feature or not.
}
\examples{
# TODO
}
\references{Brennecke et al. (2013) Accounting for technical noise in single-cell RNA-seq experiments. Nature Methods 10, 1093-1095. doi:10.1038/nmeth.2645}
\seealso{get_cluster_profile, complex_markers, M3DropFeatureSelection, NBumiFeatureSelectionCombinedDrop, kruskal.test}
\keyword{feature selection}


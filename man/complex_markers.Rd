\name{complex_markers}
\alias{complex_markers}
\title{Identify Marker Genes}
\description{
	Identifies marker genes of 1 or more groups of cells.
}
\usage{
	complex_markers(expr_mat, labels, n_max=1, strict_only = FALSE, strict_BON=FALSE)
}
\arguments{
  \item{expr_mat}{a numeric matrix of expression values. Columns = cells, Rows = genes.}
  \item{labels}{a vector of group assignments for each cell.}
  \item{n_max}{maximum number of groups a gene can be a marker of}
  \item{strict_only}{if TRUE then only genes significantly (p < 0.05) better markers of their set of groups than the next-best set of groups. see Details.}
  \item{strict_BON}{if TRUE then applies a strict Bonferroni multiple testing correction otherwise applies a weaker version. see Details}
}
\details{
	This function uses the ability of genes to discriminate between groups of cells to identify marker genes. Each set of 1 to n_max 
	groups is considered and the AUC of the classification based on the expression of each individual gene is calculated. Each gene is 
	assigned as markers to the set of groups with the highest AUC. 

	If \code{strict_only=TRUE} then only genes where the 95\% CI of the AUC of the best set of groups does not overlap the AUC of the next-
	best set of groups are assigned as markers. This is equivalent to testing whether the gene is a significantly better marker of 
	the best set of groups than all other sets of groups. 

	Overall significance of markers (useful for when \code{strict_only=FALSE}) are calculated using \code{wilcox.test} of the assigned 
	set of groups vs the remaining cells. A Bonferroni correction is applied accounting for the total number of genes and total number 
	of groups in the original data. If \code{strict_BON=TRUE} then p-values are corrected accounting for all possible combinations of 
	groups (i.e. n_clusters choose n_max), this is most likely excessively harsh since the Bonferroni correction assumes independent 
	tests which is not the case for testing different combinations of groups.

	AUCs and 95\% CIs are calculated using the \code{pROC} package, p-values are calculated with \code{wilcox.test}.
	Genes which are not considered markers recieve an AUC = -1.

}
\value{
	a dataframe with one row per gene and  columns: 
		if n_max == 1 : AUC, Group, p.value, q.value, 
		when n_max > 1 then the data.frame contains one column for each original group with a "1" indicating the gene was a marker of that group, and 0 otherwise.
}
\examples{
# TODO
}
\seealso{get_combo_names, marker_heatmap, get_cluster_profile}
\keyword{feature selection}


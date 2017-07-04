\name{get_cluster_profile}
\alias{get_cluster_profile}
\title{Calculate Cell-type Profiles}
\description{
	Calculates cell-type profiles and features for a dataset.
}
\usage{
	get_cluster_profile(expr_mat, clusters, norm_method="CPM", is.log=FALSE, out.log=2, feature_selection=kw.features)
}
\arguments{
  \item{expr_mat}{a numeric matrix of expression values. Columns = cells, Rows = genes.}
  \item{clusters}{a vector of group assignments for each cell.}
  \item{norm_method}{what method to use to normalize the cell-type profiles. Default is counts per million for other options see: edgeR::calcNormFactors}
  \item{is.log}{FALSE if input matrix is not log transformed, otherwise the base of the log-transformation of the input matrix.}
  \item{out.log}{base of the log-transformation to apply to the normalized profiles, or FALSE if no log transformation is desired.}
  \item{feature_selection}{function to use to identify relevant features of the dataset.}
}
\details{
	Calculates cell-type profile as the mean unlogged expression of each gene within each cluster. These profiles are then normalized 
	according to \code{norm_method}. Normalization factors are calculated using \code{edgeR}. 

	The normalized profiles are then scaled so the median total expression per profile is 1,000,000 thus making profiles more 
	comparable across scRNASeq experimental protocols. Finally these profiles are log-transformed as specified.

	Feature selection is performed on the full unlogged expression matrix by the provided function. The feature selection method may 
	be one of those provided (see: "See Also") or a custom function. Requirements for this function are (1) to take the expression 
	matrix as the first argument, (2) take the cluster assignments as the second argument, (3) return a vector of 1s & 0s indicating
	whether each gene is or is not a feature. Use \code{no.feature.selection} to label all genes as features.
}
\value{
	a matrix of profiles for each cluster. The last column of this matrix is a vector of 1s & 0s indicating which genes are the most important features of the dataset.
}
\examples{
# TODO
}
\seealso{kw.features, m3drop.features, danb.features, hvg.features, marker.features, no.feature.selection}
\keyword{feature selection, cell-type profiles}


\name{cluster_profile_heatmap}
\alias{cluster_profile_heatmap}
\title{Plot a heatmap of marker gene expression}
\description{
	Identifies marker genes of 1 or more groups of cells.
}
\usage{
	cluster_profile_heatmap(corrected_mat, matches, features_only=TRUE, npermute=0, distfun=function(x){as.dist(1-cor(t(x), method="spearman"))}, hclustfun=function(x){hclust(x,method="complete")}, ann=NULL)
}
\arguments{
  \item{matches}{output from \code{combine_and_match_clusters}.}
  \item{corrected_mat}{Batch corrected matrix of expression profiles. Columns = cells, Rows = genes. Calculated by \code{glm_of_matches}.}
  \item{features_only}{Whether to plot expression of all genes or just the union of features from each dataset.}
  \item{npermute}{Number of permutations used to assess significance.}
  \item{distfun}{Distance function to use.}
  \item{hclustfun}{Clustering function to use.}
  \item{ann}{Known annotations of the cell-types to be plotted as a colour bar across the top of the heatmap.}
}
\details{
	Plots a heatmap, using heatmap3, of the batch-corrected cell-type profiles. 
	Includes upto three colour bars across the top illustrating:
		(1) Which clusters were initially matched by reciprocal 
			best-hits (Matched)
		(2) Clusters not more similar to each other than expected by
			 chance (Signif) based on permutations.
		(3) Known equivalent cell-types provided by \code{ann}.

	For details of how (1) is calculated see: \code{combine_and_match_clusters}.
	Expected similarity for (2) is calculated as the 5th percentile of the
	distances between cell-type profiles obtained by permuting the 
	expression levels across profiles. Distances for the permutations and
	used in the heatmap clustering are calculated with the provided distance
	function. Cell-type profiles are clustered on the basis of this distance
	using any hclust-like function provided by the user and the resulting
	dendrogram is cut at the permutation-derived expected similarity. 

	Gene expression levels are scaled (converted to Z-scores) prior to plotting.
}
\value{
	Invisibly, output from heatmap call.
}
\examples{
# TODO
}
\seealso{glm_of_matches, combine_and_match_clusters, get_cluster_profile}
\keyword{cell-type profiles}

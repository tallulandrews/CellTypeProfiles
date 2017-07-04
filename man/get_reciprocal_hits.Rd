\name{get_reciprocal_hits}
\alias{get_reciprocal_hits}
\title{Match cluster profiles}
\description{
	Find reciprocal best-hits between clusters of two different datasets.
}
\usage{
	get_reciprocal_hits(profiles1, profiles2, features=NA)
}
\arguments{
  \item{profiles1}{output from \code{get_cluster_profile} for the first dataset}
  \item{profiles2}{output from \code{get_cluster_profile} for the second dataset}
  \item{features}{custom features to use for matching profiles.}
}
\details{
	Matches the profiles of clusters from two different datasets to each other. Similarity between clusters is calculated using 
	Spearman correlations. Each cluster from each dataset is matched to the most similar cluster in the other dataset. Those 
	clusters that were each other's match, i.e. reciprocal best-hits, are identified. 

	if \code{features} is undefined then the union of features identified in each dataset and stored in the profile matrix are used. 
	Otherwise, the provided features are used. Genes not present in both datasets are excluded from consideration.
}
\value{
	a named list including:
		matches1 - the clusters of the second dataset which are the best match to each cluster of the first dataset
		matches2 - the clusters of the first dataset which are the best match to each cluster of the second dataset
		recip - a table of reciprocal best hits, the first column being the cluster of dataset1 and the second column being the cluster of dataset2
}
\examples{
# TODO
}
\seealso{combine_and_match_clusters, get_cluster_profile}
\keyword{cell-type profiles}


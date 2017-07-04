\name{combine_and_match_clusters}
\alias{combine_and_match_clusters}
\title{Match cluster profiles}
\description{
	Find reciprocal best-hits between clusters of two different datasets.
}
\usage{
	combine_and_match_clusters(profile_list, features=NA)
}
\arguments{
  \item{profile_list}{a list of output from \code{get_cluster_profile} for all datasets to merge.}
  \item{features}{custom features to use for matching profiles.}
}
\details{
	Clusters are matched across pairs of clusters using reciprocal best-hits using Spearman correlation similarity.
	These pair-wise matches are used to create a graph linking all clusters across all datasets. The connected components of
	this graph are assigned to "Matched" cluster labels.

	Profiles for all datasets are sorted alphabetically by gene-name and genes not present in all datasets are excluded. 
	These profiles are combined into a single matrix. Feature columns were removed from the profiles and combined into a 
	separate matrix. 

	If \code{features} is undefined then the union of features identified in each dataset and stored in the profile matrix are used 
	for calculating reciprocal best-hits. Otherwise, the provided features are used. 
}
\value{
	a named list including:
		profiles - a matrix of gene expression profiles for all clusters in all datasets. 
		labels - a vector of cluster assignments for all profiles after merging reciprocal best-hit clusters across datasets.
		dataset - a vector of dataset of origin assignments for all profiles.
		feature - a matrix of feature assignments from each dataset (one column or 1s/0s per dataset).
}
\examples{
# TODO
}
\seealso{combine_and_match_clusters, get_cluster_profile}
\keyword{cell-type profiles}


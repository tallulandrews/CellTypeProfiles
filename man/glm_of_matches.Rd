\name{glm_of_matches}
\alias{glm_of_matches}
\title{Use GLM to remove batch effects}
\description{
	Uses general linear models to estimate batch-effects across datasets using matched clusters.
}
\usage{
	glm_of_matches(matches)
}
\arguments{
  \item{matches}{output from \code{combine_and_match_clusters}.}
}
\details{
	This function estimates batch effects using a glm based on clusters matched across all datasets, if possible. All
	sets of clusters matched across datasets are used to train a glm to distinguish biological and dataset effects. These
	effects are then subtracted from the profiles of all clusters in all datasets.

	If matched clusters are insufficient to estimate batch effects then the glm uses profiles from all clusters for all 
	datasets to estimate batch effects. 
}
\value{
	a named list including:
		corrected_mat - the batch corrected matrix of profiles.
		model_effects - a matrix of estimated model effects for each gene.
}
\examples{
# TODO
}
\seealso{glm_of_matches, combine_and_match_clusters, get_cluster_profile}
\keyword{cell-type profiles}


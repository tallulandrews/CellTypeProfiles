\name{get_combo_names}
\alias{get_combo_names}
\title{Get multi-group names}
\description{
	Converts the gene to group assignment matrix from complex_markers into a vector of compound group names.
}
\usage{
	get_combo_names(marker_matrix)
}
\arguments{
  \item{marker_matrix}{Output from \code{complex_markers}.}
}
\details{
	This function takes the output from \code{complex_markers} and creates a compound group name for the set of groups each marker 
	gene is assigned to. Names of groups are concatenated with a "+". If a gene is a marker of more than the majority of groups it is 
	instead described as a negative marker (prefix "NOT") of the remaining groups.

}
\value{
	vector of marker gene group assignment names.
}
\examples{
# TODO
}
\seealso{complex_markers, marker_heatmap}


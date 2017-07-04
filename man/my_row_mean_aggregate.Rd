\name{my_row_mean_aggregate}
\alias{my_row_mean_aggregate}
\title{Row Means by Factor}
\description{
	Efficiently calculates the row means of a matrix split by a factor.
}
\usage{
	my_row_mean_aggregate(mat, groups)
}
\arguments{
  \item{mat}{a matrix}
  \item{groups}{a factor of length ncol(mat) to split the matrix by.}
}
\details{
	Calculates the row means for columns belonging to each level of the provided groups. Uses an approach much more efficient than \code{aggregate}. Requires at least 2 samples from each group.
}
\value{
	a matrix of dimension nrow(mat) x length(levels(groups)) of means for each row for each group.
}
\examples{
	groups <- c("A", "A", "B", "B", "C","C","C")
	vals <- matrix(rnorm(70), ncol=7)
	group_means = my_row_mean_aggregate(vals, groups)
}

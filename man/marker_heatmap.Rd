\name{marker_heatmap}
\alias{marker_heatmap}
\title{Plot a heatmap of marker gene expression}
\description{
	Identifies marker genes of 1 or more groups of cells.
}
\usage{
	marker_heatmap(marker_matrix, expr_mat, cell_colour_bars=NULL, min_AUC=0.7, max_q_value=1, min_cat_size=10)
}
\arguments{
  \item{marker_matrix}{Output from \code{complex_markers}.}
  \item{expr_mat}{a numeric matrix of expression values. Columns = cells, Rows = genes.}
  \item{cell_colour_bars}{a matrixor data.frame with nrow = number of cells where is column is a vector of colours to be plotted as bars on the heatmap, corresponds to ColSideColors.}
  \item{min_AUC}{Threshold for AUC for markers to be plotted.}
  \item{max_q_value}{Threshold for significance of markers to be plotted.}
  \item{min_cat_size}{Minimum number of markers assigned to a particular set of clusters for those markers to be plotted.}
}
\details{
	Plots a heatmap, using heatmap3, of marker genes identified by \code{complex_markers}. Markers are filtered according to user 
	provided parameters, and the scaled expression across cells is plotted. Expression levels are coloured from blue to red by 
	z-scores with 10 bins between -2 and 2. Genes and cells are clustered using Ward's hierarchical clustering (see: \code{hclust}). 
	Zero, one or more bars are plotted across the top according to \code{cell_colour_bars}. And a colour bar is plotted down the 
	left-hand side according to the set of groups each marker was assigned to (see: \code{get_combo_names}).
}
\value{
	Invisibly, a named list including: legend = a named vector of gene-group colours used in the plot, and gene_assignments = a vector of compound group names for each marker gene.
}
\examples{
# TODO
}
\seealso{get_combo_names, complex_markers, heatmap3, hclust}
\keyword{heatmap}


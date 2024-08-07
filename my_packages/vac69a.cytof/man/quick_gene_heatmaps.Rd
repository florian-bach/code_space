% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quick_gene_heatmaps.R
\name{quick_gene_heatmaps}
\alias{quick_gene_heatmaps}
\title{This function searches the  DoD_Baseline and T6_Baseline pairwise comparison results for all the search terms supplied
in the character vector; the objects DoD_Baseline_sig AND T6_Baseline_sig need to be defined in the environment for this
to work}
\usage{
quick_gene_heatmaps(search_results, sort_by)
}
\arguments{
\item{search_results}{character string; a list of lists containing gene names for making heatmaps}

\item{sort_by}{character string; the name of the comparison used to order genes in the heatmap must be in unique(all_unique_genes$file_name)
("C14_Baseline" "DoD_Baseline" "C56_Baseline" "T6_Baseline"  "T6_DoD")}
}
\value{
a named list of lists containing the data frames used for plotting. The heatmaps will be save in ./figures relative to the
working directory
}
\description{
This function searches the  DoD_Baseline and T6_Baseline pairwise comparison results for all the search terms supplied
in the character vector; the objects DoD_Baseline_sig AND T6_Baseline_sig need to be defined in the environment for this
to work
}

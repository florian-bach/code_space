##                                             IMPORTANT !!!                                    ###
###    IF THIS DESCRIPTION IS ALTERED YOU HAVE TO RUN DOCUMENT() AT THE END OF THIS FILE        ###


#'
#' This function searches the  DoD_Baseline and T6_Baseline pairwise comparison results for all the search terms supplied
#' in the character vector; the objects DoD_Baseline_sig AND T6_Baseline_sig need to be defined in the environment for this
#' to work

#' @param search_terms character string; one or more search terms

#' @return a named list of lists containing the duplicated-cleaned gene symbols that were found using the search term
#' @importFrom magrittr %>%
#' @export
#'
#'



search_genes <- function(x){

  all_results <- list()

  for (i in x) {

    search_dod <- DoD_Baseline_sig %>%
      dplyr::filter(grepl(i, DoD_Baseline_sig$Description, ignore.case = T)) %>%
      dplyr::arrange(log2FoldChange) %>%
      dplyr::select(Symbol)

    search_t6 <- T6_Baseline_sig %>%
      dplyr::filter(grepl(i, T6_Baseline_sig$Description, ignore.case = T)) %>%
      dplyr::arrange(log2FoldChange) %>%
      dplyr::select(Symbol)

    search_results <- rbind(search_dod, search_t6)
    search_results <- list(search_results[!duplicated(search_results$Symbol),])
    names(search_results) <- i

    all_results <- c(all_results, search_results)

  }
  return(all_results)
}

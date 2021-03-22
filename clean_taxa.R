#' Function for cleaning leading characters from taxa names
#'
#' @import phyloseq
#' @return Phyloseq object with leading x___ removed
#' @param phylo Phyloseq object. Can be generated from meta_mat_to_phyloseq.
#' Required
#' @export
#' @details  Removes x___ from all samples
#' @examples
#' TBD

clean_taxa = function(phylo) {

  cleaned = apply(tax_table(phylo),2 , function(x) gsub("*.__", "", x))
  tax_table(phylo) = cleaned
  return(phylo)

}

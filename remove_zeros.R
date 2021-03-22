#' Function for remove samples that are all zero
#'
#' @import phyloseq
#' @import dplyr
#' @return Phyloseq object with all zero samples removed
#' @param phylo Phyloseq object. Can be generated from meta_mat_to_phyloseq.
#' Required
#' @export
#' @details  Removes all zero samples from a phyloseq object
#' @examples
#' TBD

remove_zeros = function(phylo) {
  phylo %>% phyloseq::prune_samples(sample_sums(.) != 0,. )
}

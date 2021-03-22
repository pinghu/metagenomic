#' Conversion of matrix and meta data to phyloseq object
#'
#' @import dplyr
#' @import magrittr
#' @import phyloseq
#' @importFrom tidyr separate
#' @return A phyloseq containing otu data from the matrix and sample data from
#' the meta matrix
#' @param mat Output from read_micro. Raw matrix where the leading columns
#' are taxonomy starting with kingdom and ending with species. Values are
#' case-sensitive
#' @param meta Matrix of meta data corresponding to mat. Can have any
#' number of columns but each row should correspond to a row in mat
#' @export
#' @details Function takes an otu matrix (From read_micro) and meta data and
#' automatically converts it into a phyloseq object.
#' @examples
#' #Read in data

meta_mat_to_phyloseq = function(mat, meta) {
  meta_tmp = meta
  #Creating samples mat
  row.names(meta_tmp) = paste0("Sample", seq_len(nrow(meta_tmp)))
  samples = phyloseq::sample_data(meta_tmp)

  #Creating taxa mat
  row.names(mat) = paste0("OTU", seq_len(nrow(mat)))
  tax = mat %>% dplyr::select(kingdom:species) %>% as.matrix() %>% phyloseq::tax_table()

  #Creating objects for phyloseq

  otu = mat %>% dplyr::select(-(kingdom:species)) %>% as.matrix() %>% round(.) %>%
    phyloseq::otu_table(., taxa_are_rows = TRUE)

  colnames(otu) = paste0("Sample", seq_len(ncol(otu)))

  #Placing everything in phyloseq

  phylo = phyloseq::phyloseq(otu, tax, samples)
}

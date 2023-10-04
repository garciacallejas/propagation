
#' remove a species from a meta-adjacency matrix
#'
#' in this version, rownames/colnames are of the form "spXX-C", 
#' with C being the cell id. if a name "spXX" is given, all instances
#' of the species will be removed, i.e. it will be removed in all cells.
#'
#' @param landscape meta-adjacency matrix, with proper rownames and colnames
#' @param sp.name string giving name to remove. 
#'
#' @return meta-adjacency matrix with the species removed in rows and columns
#' @export
#'
#' @examples
remove_sp_from_landscape <- function(landscape, sp.name){
  
  sp.pos <- grepl(sp.name,rownames(landscape))
  landscape.trimmed <- landscape[!sp.pos,!sp.pos]
  
  return(landscape.trimmed)
}

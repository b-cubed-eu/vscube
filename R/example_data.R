#' Read an example dataset
#' @return occurrence.txt file
#' @export
#' @examples
#' example_data()
example_data <- function() {
file <- system.file("inst", "occurrence.txt", package = "vscube") 
return(file)
}

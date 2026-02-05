#' Read an example dataset
#' @return path to occurrence.txt
#' @export
#' @examples
#' example_file()
example_file <- function() {
file <- system.file("extdata", "occurrence.txt", package = "vscube")
  return(file)
}

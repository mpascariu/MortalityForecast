
#' onAttach
#' @param lib lib
#' @param pkg pkg
#' @name onAttach
#' @keywords internal
".onAttach" <- function(lib, pkg){
  packageStartupMessage("\nMortalityForecast: ...",
                        "\nAuthors    : Pascariu M.D.",
                        "\nLast Update: August 21, 2018\n")
}

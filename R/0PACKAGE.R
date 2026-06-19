
#' @note
#' 
#' This R package author has retired from academic research. 
#' 
#' Accordingly, this package should not be considered a validated tool for use in peer-reviewed publications or as the basis for grant applications.  
#' 
#' Backward compatibility with user-code published in the following publications is not maintained in versions `>= 0.3.0` of this package. The authors of those publications are the appropriate contacts for reproducibility inquiries. 
#' 
#' \doi{10.1093/bioinformatics/btaf182}
#' 
#' \doi{10.1093/bioinformatics/btaf430}
#' 
#' 
#' @importFrom cli col_magenta col_blue
#' 
'_PACKAGE'



# [~~Highlighted strike through test.~~]{style="background-color: #FFFF00"}
# [Highlighted test.]{style="background-color: #FFFF00"}

#' @importFrom cli col_red col_green style_bold cli_text

.onAttach <- \(libname, pkgname) {
  
  'Author of this R package has retired from academia.' |>
    col_blue() |>
    style_bold() |>
    packageStartupMessage()
  
  'Backward compatibility with user-code published in' |>
    col_green() |>
    style_bold() |>
    packageStartupMessage()
  
  'https://doi.org/10.1093/bioinformatics/btaf182' |>
    col_red() |>
    style_bold() |>
    packageStartupMessage()

  'https://doi.org/10.1093/bioinformatics/btaf430' |>
    col_red() |>
    style_bold() |>
    packageStartupMessage()

  'is not maintained in versions >= 0.3.0 of this package.' |>
    col_blue() |>
    style_bold() |>
    packageStartupMessage()
  
  'Author(s) of those publications are the appropriate contacts for reproducibility inquiries.' |>
    col_green() |>
    style_bold() |>
    packageStartupMessage()
  
}


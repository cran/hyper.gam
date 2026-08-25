
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

#' @importFrom cli col_red make_ansi_style style_bold cli_text

.onAttach <- \(libname, pkgname) {
  
  'Backward compatibility with user-code published in' |>
    col_blue() |>
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

  'is not maintained in versions >= 0.3.0 (June 2026) of this package.' |>
    col_blue() |>
    style_bold() |>
    packageStartupMessage()

  'The author of this package has retired from academic life and bears no ongoing obligation to any research program, institution, or funding agency. This package is no longer actively maintained as a research tool. Users are strongly advised against relying on this package for academic research, manuscript preparation, or grant applications. Reproducibility inquiries regarding previously published work citing this package should be directed to the authors of those publications.' |>
    make_ansi_style('grey70')() |>
    packageStartupMessage()
  
}


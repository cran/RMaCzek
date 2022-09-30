## This file is part of RMaCzek

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


#' @title function to load data from an mdt file (MaCzek 3.3 - http://www.antropologia.uw.edu.pl/MaCzek/maczek.html)
#' @author Piotr Jaskulski
#' @description This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. Please understand that there may still be bugs and errors. Use it at your own risk. We take no responsibility for any errors or omissions in this package or for any misfortune that may befall you or others as a result of its use. Please send comments and report bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .
#' @param filepath path to file *.mdt.
#' @return data.frame.
#'

read_maczek_file <- function(filepath) {
  con <- file(filepath, "r")
  num_line <- 1
  start_rl <- FALSE
  start_vl <- FALSE
  start_data <- FALSE
  rlabel <- c()
  vlabel <- c()
  tmp <- tempfile()

  while ( TRUE ) {
    line <- readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    line <- trimws(line)
    if (num_line == 1 && line != 'MaCzek DATA 3.3') {
      close(con)
      stop("Error: This is not MaCzek DATA ver. 3.3 file!")
    }
    num_line <- num_line + 1
    if (line == '[RECORD_LABEL]') {
      start_rl <- TRUE
    }
    else if (line == '[VARIABLE_LABEL]') {
      start_vl <- TRUE
      start_rl <- FALSE
    }
    else if (line == '[DATA]') {
      start_data <- TRUE
      start_vl <- FALSE
    }
    else {
      if (start_rl && line != '') {
        rlabel <- append(rlabel, line)
      }
      else if (start_vl && line != '') {
        vlabel <- append(vlabel, line)
      }
      else if (start_data && line != '') {
        cat(line, file = tmp, sep = '\n', append = TRUE)
      }
    }
  }
  close(con)

  mdt <- utils::read.csv(tmp, header = FALSE, sep = ';')
  file.remove(tmp)
  rownames(mdt) <- rlabel
  colnames(mdt) <- vlabel

  return(mdt)
}


## This file is part of RMaCzek

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


.get_parameters <- function(parameter, defaults) {
  defaults <- as.list(defaults)
  parameter <- as.list(parameter)

  ## add verbose
  if(is.null(defaults$verbose)) defaults$verbose <- FALSE

  if(length(parameter) != 0) {
    o <- pmatch(names(parameter), names(defaults))

    ## unknown parameter
    if(any(is.na(o))){
      warning(sprintf(ngettext(length(is.na(o)),
                               "Unknown parameter: %s",
                               "Unknown parameters: %s"),
                      paste(names(parameter)[is.na(o)],
                            collapse = ", ")), call. = FALSE, immediate. = TRUE)

      message("Available parameter (with default values):\n")
      message(rbind(names(defaults)," = ", gsub("\n"," ",as.character(defaults))),
          sep=c("\t"," ","\n"))
    }

    defaults[o[!is.na(o)]] <- parameter[!is.na(o)]
  }

  if(defaults$verbose) {
    message("Used parameters:\n")
    message(rbind(names(defaults)," = ",
              strtrim(gsub("\n"," ",as.character(defaults)), 50)),
        sep=c("\t"," ","\n"))
  }

  defaults
}


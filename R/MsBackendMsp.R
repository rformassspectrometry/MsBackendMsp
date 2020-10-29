#' @include hidden_aliases.R
NULL

#' @title MS data backend for msp files
#'
#' @aliases MsBackendMsp-class
#'
#' @description
#'
#' The `MsBackendMsp` class supports import of MS/MS spectra data from
#' files in Mascot Generic Format
#' ([msp](http://www.matrixscience.com/help/data_file_help.html))
#' files. After initial import, the full MS data is kept in
#' memory. `MsBackendMsp` extends the [MsBackendDataFrame()] backend
#' directly and supports thus the [applyProcessing()] function to make
#' data manipulations persistent. The backend does however not
#' support export to msp files yet.
#'
#' New objects are created with the `MsBackendMsp` function. The
#' `backendInitialize` method has to be subsequently called to
#' initialize the object and import MS/MS data from (one or more) msp
#' files.  Optional parameter `nonStop` allows to specify whether the
#' import returns with an error if one of the xml files lacks required
#' data, such as `mz` and `intensity` values (default `nonStop =
#' FALSE`), or whether only affected file(s) is(are) skipped and a
#' warning is shown (`nonStop = TRUE`). Note that any other error
#' (such as xml import error) will abort import regardless of
#' parameter `nonStop`.
#'
#' @param object Instance of `MsBackendMsp` class.
#'
#' @param files `character` with the (full) file name(s) of the msp file(s)
#'     from which MS/MS data should be imported.
#'
#' @param nonStop `logical(1)` whether import should be stopped if an
#'     xml file does not contain all required fields. Defaults to
#'     `nonStop = FALSE`.
#'
#' @param BPPARAM Parameter object defining the parallel processing
#'     setup to import data in parallel. Defaults to `BPPARAM =
#'     bpparam()`. See [bpparam()] for more information.
#'
#' @param ... Currently ignored.
#'
#' @author Laurent Gatto and Johannes Rainer
#'
#' @importClassesFrom Spectra MsBackendDataFrame
#'
#' @exportClass MsBackendMsp
#'
#' @name MsBackendMsp
#'
#' @examples
#'
#' ## Create an MsBackendHmdbXml backend and import data from test xml files.
#' fls <- dir(system.file("extdata", package = "MsBackendMsp"),
#'     full.names = TRUE, pattern = "msp$")
#' be <- backendInitialize(MsBackendMsp(), fls)
#' be
#'
#' be$msLevel
#' be$intensity
#' be$mz
NULL

setClass("MsBackendMsp",
         contains = "MsBackendDataFrame",
         prototype = prototype(spectraData = DataFrame(),
                               readonly = FALSE,
                               version = "0.1"))

#' @importMethodsFrom Spectra backendInitialize spectraData<- $<- $
#'
#' @importFrom BiocParallel bpparam
#'
#' @importMethodsFrom BiocParallel bplapply
#'
#' @importFrom methods validObject
#'
#' @exportMethod backendInitialize
#'
#' @rdname MsBackendMsp
setMethod("backendInitialize", signature = "MsBackendMsp",
          function(object, files, nonStop = FALSE, ..., BPPARAM = bpparam()) {
              if (missing(files) || !length(files))
                  stop("Parameter 'files' is mandatory for ", class(object))
              if (!is.character(files))
                  stop("Parameter 'files' is expected to be a character vector",
                       " with the files names from where data should be",
                       " imported")
              files <- normalizePath(files)
              if (any(!file.exists(files)))
                  stop("file(s) ",
                       paste(files[!file.exists(files)], collapse = ", "),
                       " not found")
              ## Import data and rbind.
              message("Start data import from ", length(files), " files ... ",
                      appendLF = FALSE)
              res <- bplapply(files, FUN = .read_msp,
                              nonStop = nonStop, BPPARAM = BPPARAM)
              message("done")
              res <- do.call(rbind, res)
              if (nonStop && length(files) > nrow(res))
                      warning("Import failed for ", length(files) - nrow(res),
                              " files")
              spectraData(object) <- res
              object$dataStorage <- "<memory>"
              object$centroided <- TRUE
              validObject(object)
              object
          })

#' @rdname MsBackendMsp
#'
#' @importFrom methods new
#'
#' @export MsBackendMsp
MsBackendMsp <- function() {
    new("MsBackendMsp")
}

#' @export
#' 
#' @rdname MsBackendMsp
spectraVariableMapping <- function(format = c("msp")) {
  ## In future eventually define that in a text file and import upon package
  ## init.
  switch(match.arg(format),
         "msp" = c(
          # minimal information
           name = "NAME:",
           accession = "DB#:",
           formula = "FORMULA:",
           inchikey = "INCHIKEY:",
           adduct = "PRECURSORTYPE:",
           exactmass = "EXACTMASS:",
           rtime = "RETENTIONTIME:",
           precursorMz = "PRECURSORMZ:"
          
         )
  )
}

#' @importMethodsFrom Spectra export
#'
#' @exportMethod export
#'
#' @rdname MsBackendMassbank
setMethod("export", "MsBackendMsp", function(object, x, file = tempfile(),
                                                  mapping = spectraVariableMapping(),
                                                  ...) {
  .export_msp(x = x, con = file, mapping = mapping)
})
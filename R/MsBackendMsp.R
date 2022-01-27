#' @include hidden_aliases.R
NULL

#' @title MS data backend for msp files
#'
#' @aliases MsBackendMsp-class
#'
#' @description
#'
#' The `MsBackendMsp` class supports import of MS/MS spectra data from
#' files in NIST MSP file format. `MsBackendMsp` extends the
#' [MsBackendDataFrame()] backend directly and supports thus the
#' [applyProcessing()] function to make data manipulations persistent.
#'
#' New objects are created with the `MsBackendMsp` function. The
#' `backendInitialize` method has to be subsequently called to
#' initialize the object and import MS/MS data from (one or more) msp
#' files.
#'
#' The `spectraVariableMapping` function allows to provide the mapping between
#' spectra variable names (i.e. the names that will be used for the spectra
#' variables in the [Spectra()] object) and the data field names of the
#' MSP file. Parameter `format` allows to select pre-defined mapping (e.g. for
#' MSP files from MoNA).
#' 
#' @param object Instance of `MsBackendMsp` class.
#'
#' @param file `character` with the (full) file name(s) of the msp file(s)
#'     from which MS/MS data should be imported or exported.
#'
#' @param format For `spectraVariableMapping`: `character(1)` specifying for
#'     which MSP *flavour* the mapping should be returned. Currently supported
#'     are: `format = "msp"` (generic MSP format).
#' 
#' @param mapping named `character` vector to rename MSP fields to spectra
#'     variables (see [spectraVariableMapping()]). This allows to correctly
#'     import also custom fields or data from files with different MSP
#'     *flavors*.
#' 
#' @param BPPARAM Parameter object defining the parallel processing
#'     setup to import data in parallel. Defaults to `BPPARAM =
#'     bpparam()`. See [bpparam()] for more information.
#'
#' @param x For `export`: a [Spectra()] object that should be exported to the
#'     specified MSP file.
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
#' ## Import spectra from a MSP file from LipidBlast
#' f <- system.file("extdata", "small-export-LipidBlast.msp",
#'     package = "MsBackendMsp")
#' be <- backendInitialize(MsBackendMsp(), f)
#' be
#'
#' be$msLevel
#' be$intensity
#' be$mz
#'
#' ## Default spectra variable mapping
#' spectraVariableMapping(MsBackendMsp())
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
          function(object, file, 
                   mapping = spectraVariableMapping(object), ...,
                   BPPARAM = bpparam()) {
              if (missing(file) || !length(file))
                  stop("Parameter 'file' is mandatory for ", class(object))
              if (!is.character(file))
                  stop("Parameter 'file' is expected to be a character vector",
                       " with the file names from where data should be",
                       " imported")
              file <- normalizePath(file)
              if (any(!file.exists(file)))
                  stop("file(s) ",
                       paste(file[!file.exists(file)], collapse = ", "),
                       " not found")
              ## Import data and rbind.
              message("Start data import from ", length(file), " files ... ",
                      appendLF = FALSE)
              res <- bplapply(file, FUN = readMsp, mapping = mapping,
                              BPPARAM = BPPARAM)
              message("done")
              res <- do.call(rbindFill, res)
              spectraData(object) <- res
              object$dataStorage <- "<memory>"
              ## object$centroided <- TRUE
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

#' @importMethodsFrom Spectra spectraVariableMapping
#' 
#' @exportMethod spectraVariableMapping
#'
#' @rdname MsBackendMsp
setMethod("spectraVariableMapping", "MsBackendMsp",
          function(object, format = c("msp")) {
              switch(match.arg(format),
                     "msp" = c(
                         name = "NAME",
                         accession = "DB#",
                         formula = "FORMULA",
                         inchikey = "INCHIKEY",
                         adduct = "PRECURSORTYPE",
                         exactmass = "EXACTMASS",
                         rtime = "RETENTIONTIME",
                         precursorMz = "PRECURSORMZ"
                     )
                     )
          })

#' @importMethodsFrom Spectra export
#'
#' @exportMethod export
#'
#' @rdname MsBackendMsp
setMethod("export", "MsBackendMsp", function(object, x, file = tempfile(),
                                             mapping = spectraVariableMapping(),
                                             ...) {
    if (missing(x))
        stop("Required parameter 'x' is missing. 'x' should be a 'Spectra' ",
             "object with the full spectra data.")
    if (!inherits(x, "Spectra"))
        stop("Parameter 'x' is supposed to be a 'Spectra' object with the full",
             " spectra data to be exported.")
    .export_msp(x = x, con = file, mapping = mapping)
})

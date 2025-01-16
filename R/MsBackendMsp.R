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
#' [Spectra::MsBackendDataFrame()] backend directly and supports thus the
#' [Spectra::applyProcessing()] function to make data manipulations persistent.
#'
#' New objects are created with the `MsBackendMsp()` function. The
#' `backendInitialize()` method has to be subsequently called to
#' initialize the object and import MS/MS data from (one or more) msp
#' files.
#'
#' The `MsBackendMsp` backend provides an `export()` method that allows to
#' export the data from the `Spectra` object (parameter `x`) to a file in
#' MSP format.
#'
#' Parameters to this function are:
#'
#' - `x`: the `Spectra` object that should be exported.
#' - `file`: `character(1)` with the desired file name.
#' - `mapping`: named `character` providing the mapping between spectra
#'   variables and MSP data fields. Defaults to
#'   `mapping = spectraVariableMapping(MsBackendMsp())`.
#' - `allVariables`: `logical(1)` whether all spectra variables in `x` should be
#'   exported or only those defined with `mapping`.
#' - `exportName`: `logical(1)` whether a `NAME` field should always be exported
#'   even if not provided in `x`.
#' 
#' See the package vignette for details and examples.
#'
#' The `spectraVariableMapping()` function allows to provide the mapping
#' between spectra variable names (i.e. the names that will be used for the
#' spectra variables in the [Spectra::Spectra()] object) and the data field
#' names of the MSP file. Parameter `format` allows to select pre-defined
#' mappings. Currently supported mapping flavors are:
#' 
#' - `format = "msp"`: default MSP field names. Should work with standard NIST
#'   MSP files or MSP files exported from MS-DIAL.
#' - `format = "mona"`: MSP file format from MoNA including LipidBlast.
#'
#' @note
#'
#' Format requirements/assumptions of MSP files:
#'
#' - Comment lines are expected to start with a `#`.
#' - Multiple spectra within the same MSP file are separated by an empty line.
#' - The first n lines of a spectrum entry represent metadata.
#' - Metadata is provided as "name: value" pairs (i.e. name and value separated
#'   by a ":").
#' - One line per mass peak, with values separated by a whitespace or tabulator.
#' - Each line is expected to contain at least the m/z and intensity values (in
#'   that order) of a peak. Additional values are currently ignored.
#' 
#' @param object Instance of `MsBackendMsp` class.
#'
#' @param file `character` with the (full) file name(s) of the msp file(s)
#'     from which MS/MS data should be imported or exported.
#'
#' @param format For `spectraVariableMapping()`: `character(1)` specifying for
#'     which MSP *flavour* the mapping should be returned. Currently supported
#'     are: `format = "msp"` (generic MSP format, for example for MS-DIAL MSP
#'     files) and `format = "mona"` (MSP files in MoNA flavour).
#' 
#' @param mapping named `character` vector to rename MSP fields to spectra
#'     variables. This allows to correctly import also custom fields or data
#'     from files with different MSP *flavors*.
#' 
#' @param allVariables `logical(1)` whether all spectra variables in `x`
#'     should be exported or only those defined with `mapping`.
#'
#' @param exportName `logical(1)` whether a `NAME` field should always be
#'     exported even if not provided in `x`.
#' 
#' @param BPPARAM Parameter object defining the parallel processing
#'     setup to import data in parallel. Defaults to `BPPARAM =
#'     SerialParam()`. See [BiocParallel::bpparam()] for more information.
#'     Parallel processing would make most sense for import from a large
#'     set of individual MSP files, but could also improve performance for
#'     import from a (very large) single MSP file.
#'
#' @param x For `export()`: a [Spectra::Spectra()] object that should be
#'     exported to the specified MSP file.
#' 
#' @param ... Currently ignored.
#'
#' @return `MsBackendMsp()` and `backendInitialize()` return an instance of a
#'     `MsBackendMsp` class. `spectraVariableMapping()` a named `character`
#'     vector with the mapping between spectra variables and MSP data fields.
#' 
#' @author Steffen Neumann, Michael Witting, Laurent Gatto and Johannes Rainer
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
#' ## precursor m/z are however all missing
#' be$precursorMz
#'
#' ## Default spectra variable mapping
#' spectraVariableMapping(MsBackendMsp())
#'
#' ## In fact, to read MSP files in "LipidBlast flavour" (same as MoNA) we
#' ## should use a different spectra variable mapping
#' spectraVariableMapping(MsBackendMsp(), "mona")
#'
#' ## Importing the data with this will correctly retrieve data
#' be <- backendInitialize(MsBackendMsp(), f,
#'     mapping = spectraVariableMapping(MsBackendMsp(), "mona"))
#' be$precursorMz
#'
#' ## Other fields are also correctly mapped, but might need to be converted
#' ## to e.g. numeric, such as "exactmass"
#' be$exactmass
#'
#' be$exactmass <- as.numeric(be$exactmass)
#'
#' be$adduct
#' be$formula
#'
#' ## Exporting Spectra objects in MSP format.
#' 
#' sps <- Spectra(be)
#' export(MsBackendMsp(), sps, file = stdout())
NULL

setClass("MsBackendMsp",
         contains = "MsBackendDataFrame",
         prototype = prototype(spectraData = DataFrame(),
                               readonly = FALSE,
                               version = "0.1"))

#' @importMethodsFrom Spectra spectraData<- $<- $
#'
#' @importMethodsFrom ProtGenerics backendInitialize
#' 
#' @importFrom BiocParallel bpparam
#'
#' @importFrom BiocParallel SerialParam
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
                   BPPARAM = SerialParam()) {
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
              if (length(file) > 1) {
                  message("Start data import from ", length(file)," files ... ",
                          appendLF = FALSE)
                  res <- bplapply(file, FUN = readMsp, mapping = mapping,
                                  BPPARAM = BPPARAM)
                  message("done")
                  res <- do.call(rbindFill, res)
              } else
                  res <- readMsp(file, mapping = mapping, BPPARAM = BPPARAM)
              spectraData(object) <- res
              object$dataStorage <- "<memory>"
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
          function(object, format = c("msp", "mona")) {
              switch(match.arg(format),
                     "msp" = c(
                         name = "NAME",
                         accession = "DB#",
                         formula = "FORMULA",
                         inchikey = "INCHIKEY",
                         adduct = "PRECURSORTYPE",
                         exactmass = "EXACTMASS",
                         rtime = "RETENTIONTIME",
                         precursorMz = "PRECURSORMZ",
                         adduct = "PRECURSORTYPE",
                         smiles = "SMILES",
                         inchi = "INCHI",
                         polarity = "IONMODE",
                         instrument = "INSTRUMENT"
                     ),
                     "mona" = c(
                         name = "Name",
                         synonym = "Synon",
                         accession = "DB#",
                         inchikey = "InChIKey",
                         adduct = "Precursor_type",
                         precursorMz = "PrecursorMZ",
                         polarity = "Ion_mode",
                         formula = "Formula",
                         exactmass = "ExactMass",
                         collision_energy_text = "Collision_energy",
                         msLevel = "Spectrum_type",
                         instrument = "Instrument",
                         instrument_type = "Instrument_type"
                     )
                     )
          })

#' @importMethodsFrom Spectra export
#'
#' @exportMethod export
#'
#' @rdname MsBackendMsp
setMethod("export", "MsBackendMsp",
          function(object, x, file = tempfile(),
                   mapping = spectraVariableMapping(MsBackendMsp()),
                   allVariables = TRUE, exportName = TRUE, ...) {
    if (missing(x))
        stop("Required parameter 'x' is missing. 'x' should be a 'Spectra' ",
             "object with the full spectra data.")
    if (!inherits(x, "Spectra"))
        stop("Parameter 'x' is supposed to be a 'Spectra' object with the full",
             " spectra data to be exported.")
    .export_msp(x = x, con = file, mapping = mapping,
                allVariables = allVariables, exportName = exportName)
})

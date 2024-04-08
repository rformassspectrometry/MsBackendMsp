#' @title Reading MSP files
#'
#' @description
#' 
#' The `readMsp()` function imports the data from a file in MGF format reading
#' all specified fields and returning the data as a [DataFrame()].
#'
#' Format constraints for MSP files:
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
#' @param f `character(1)` with the path to an MSP file.
#'
#' @param msLevel `numeric(1)` with the MS level. Default is 2. This value will
#'     be reported as the spectra's MS level **unless** the source MSP file
#'     defines the MS level.
#'
#' @param mapping named `character` vector to rename MSP fields to spectra
#'     variables (see [spectraVariableMapping()]). This allows to correctly
#'     import also custom fields or data from files with different MSP
#'     *flavors*.
#'
#' @param BPPARAM parallel processing setup. See [bpparam()] for more details.
#' 
#' @param ... Additional parameters, currently ignored.
#'
#' @return
#'
#' A `DataFrame` with each row containing the data from one spectrum
#' in the MSP file. m/z and intensity values are available in columns `"mz"`
#' and `"intensity"` in a list representation.
#'
#' @export
#'
#' @importFrom Spectra coreSpectraVariables
#' 
#' @importFrom S4Vectors DataFrame
#'
#' @importFrom IRanges NumericList
#'
#' @importFrom MsCoreUtils rbindFill
#'
#' @importFrom BiocParallel bpmapply
#'
#' @importFrom methods as
#'
#' @author Laurent Gatto, Steffen Neumann, Johannes Rainer
#'
#' @examples
#'
#' f <- system.file("extdata", "minimona.msp", package = "MsBackendMsp")
#'
#' readMsp(f)
readMsp <- function(f, msLevel = 2L,
                    mapping = spectraVariableMapping(MsBackendMsp()),
                    BPPARAM = SerialParam(), ...) {
    if (length(f) != 1L)
        stop("Please provide a single msp file.")
    
    msp <- scan(file = f, what = "", sep = "\n", quote = "",
                allowEscapes = FALSE, quiet = TRUE,
                blank.lines.skip = FALSE, strip.white = TRUE)
    
    ## Ignore comments
    cmts <- grep("^[#]", msp)
    if (length(cmts))
        msp <- msp[-cmts]

    ## Find individual records. Instead of grepping by NAME: we use blank
    ## lines. These are expected to separate entries.
    wsp <- grep("^[[:space:]]|(^$)", msp)
    begin <- c(1, wsp +1L)
    end <- c(wsp -1L, length(msp))
    keep <- begin < end # drop consecutive blank lines.
    begin <- begin[keep]
    end <- end[keep]

    sp <- bpmapply(begin, end, FUN = function(a, b) {
         .extract_msp_spectrum(msp[a:b], mapping = mapping)
    }, SIMPLIFY = FALSE, USE.NAMES = FALSE, BPPARAM = BPPARAM)
    res <- DataFrame(rbindFill(sp))

    spv <- coreSpectraVariables()
    spv <- spv[!names(spv) %in% c("mz", "intensity")]
    for (i in seq_along(res)) {
        if (all(lengths(res[[i]]) == 1))
            res[[i]] <- unlist(res[[i]])
        if (any(col <- names(spv) == colnames(res)[i]))
            res[[i]] <- suppressWarnings(as(res[[i]], spv[col][1]))
    }

    res$mz <- NumericList(res$mz, compress = FALSE)
    res$intensity <- NumericList(res$intensity, compress = FALSE)
    res$dataOrigin <- f
    if (!any(colnames(res) == "msLevel"))
        res$msLevel <- as.integer(msLevel)
    res
}

#' @param msp `character()` of lines defining a spectrum in msp
#'     format.
#'
#' @param mapping spectra variable mapping that allows renaming data fields.
#' 
#' @author Laurent Gatto, Johannes Rainer
#' 
#' @importFrom stats setNames
#'
#' @note
#'
#' https://www.nist.gov/system/files/documents/srd/NIST1aVer22Man.pdf
#' 
#' @noRd
.extract_msp_spectrum <- function(msp, mapping) {
    ## grep description
    desc.idx <- grep(":", msp)
    desc <- msp[desc.idx]

    l <- length(desc.idx)
    if (desc.idx[1L] != 1L || desc.idx[l] != l)
        stop("MSP format error. Make sure that data for multiple spectra are ",
             "separated by an empty new line and that general spectrum ",
             "metadata/information is provided in 'name: value' format ",
             "(i.e. name and value of the metadata field separated by ",
             "a \":\").", call. = FALSE)
    
    spec <- trimws(msp[-desc.idx], "left")
    
    pks <- strsplit(sub("^(\\t|[[:space:]]+)", "", spec),
                    "[[:space:]]+")
    anns <- lengths(pks) > 2
    if (any(anns)) {
        warning("Unexpected number of values per peak found. Keeping only the ",
                "first two values assuming they correspond to m/z and ",
                "intensity", call. = FALSE)
        pks[anns] <- lapply(pks[anns], function(z) z[1:2])
    }
    ms <- do.call(rbind, pks)
    mode(ms) <- "double"
    
    if (!length(ms))
        ms <- matrix(numeric(), ncol = 2L)
    else {
        if (anyNA(ms[, 1L]))
            stop("MSP format error. Missing values (NA) for m/z are not ",
                 "supported. Please ensure that no missing or non-numeric ",
                 "values are reported as the peaks' m/z values.", call. = FALSE)
        if (is.unsorted(ms[, 1L]))
            ms <- ms[order(ms[, 1L]), ]
    }
    
    r <- regexpr(":", desc, fixed = TRUE)
    desc <- setNames(trimws(substring(desc, r + 1L, nchar(desc)), "left"),
                     substring(desc, 1L, r - 1L))

    ## map fields to spectra variables
    idx <- match(names(desc), mapping)
    not_na <- !is.na(idx)
    if (any(not_na))
        names(desc)[not_na] <- names(mapping)[idx][not_na]

    ## Handle eventually duplicated names -> list
    if (anyDuplicated(names(desc))) {
        res <- split(unname(desc), names(desc))
        dups <- lengths(res) > 1L
        dup_res <- res[dups]
        res <- as.data.frame(res[!dups])
        for (name in names(dup_res))
            res <- do.call("$<-", list(res, name, unname(dup_res[name])))
    } else res <- as.data.frame(as.list(desc))
    
    ## Ensure correct data type
    ## polarity
    if (any(have <- colnames(res) == "polarity"))
        res[, have] <- .process_polarity(res[, have])
    if (any(have <- colnames(res) == "msLevel"))
        res[, have] <- .process_mslevel(res[, have])

    res$mz <- list(ms[, 1L])
    res$intensity <- list(ms[, 2L])
    res
}

.process_polarity <- function(x, input = TRUE) {
    if (input) {
        if (grepl("^(p|\\+)", x, ignore.case = TRUE))
            return(1L)
        if (grepl("^(n|-)", x, ignore.case = TRUE))
            return(0L)
        -1L
    } else {
        xnew <- rep(NA_character_, length(x))
        xnew[x == 1L] <- "Positive"
        xnew[x == 0L] <- "Negative"
        xnew
    }
}

#' @param x value to be formatted
#'
#' @param input `logical(1)` whether the data is imported or exported.
#'
#' @noRd
.process_mslevel <- function(x, input = TRUE) {
    if (input)
        as.integer(sub("ms", "", x, ignore.case = TRUE))
    else paste0("MS", x)
}

#' @description
#'
#' Function to export a `Spectra` object in MSP format to `con`.
#'
#' @param x `Spectra`
#'
#' @param con output file.
#'
#' @param mapping named `character` vector that maps from `spectraVariables`
#'    (i.e. `names(mapping)`) to the variable name that should be used in the
#'    MSP file.
#'
#' @param allVariables `logical(1)` whether all spectra variables in `x` should
#'    be exported or only those that are listed in `mapping`. Note that if
#'    `exportName = TRUE` a field *NAME* will be exported regardless of
#'    `mapping`.
#'
#' @param exportName `logical(1)` whether a NAME field will always be exported
#'    even if no such spectra variable is available in `x`.
#' 
#' @author Michael Witting, Johannes Rainer
#'
#' @importMethodsFrom Spectra spectraVariables spectraNames spectraData
#'
#' @importMethodsFrom ProtGenerics peaksData
#'
#' @noRd
#'
#' @examples
#'
#' spd <- DataFrame(msLevel = c(2L, 2L, 2L), rtime = c(1, 2, 3))
#' spd$mz <- list(c(12, 14, 45, 56), c(14.1, 34, 56.1), c(12.1, 14.15, 34.1))
#' spd$intensity <- list(c(10, 20, 30, 40), c(11, 21, 31), c(12, 22, 32))
#'
#' sps <- Spectra(spd)
#'
#' .export_msp(sps)
#'
#' ## Handling of variables with multiple entries
#' sps$synonym <- list(c("a", "b"), "d", c("e", "f", "g"))
#' .export_msp(sps)
#'
.export_msp <- function(x, con = stdout(),
                        mapping = spectraVariableMapping(MsBackendMsp()),
                        allVariables = TRUE, exportName = TRUE) {
    spv <- spectraVariables(x)
    spv <- spv[!(spv %in% c("dataOrigin", "dataStorage"))]
    if (!allVariables)
        spv <- spv[spv %in% names(mapping)]
    spd <- spectraData(x, spv)
    ## Process any known required data conversions
    if (any(spv == "msLevel"))
        spd$msLevel <- .process_mslevel(spd$msLevel, input = FALSE)
    if (any(spv == "polarity"))
        spd$polarity <- .process_polarity(spd$polarity, input = FALSE)
    idx <- match(colnames(spd), names(mapping))
    colnames(spd)[!is.na(idx)] <- mapping[idx[!is.na(idx)]]
    ## Force variable NAME:
    if (!any(tolower(colnames(spd)) == "name") && exportName)
        spd$NAME <- seq_len(nrow(spd))
    idx <- which(tolower(colnames(spd)) == "name")
    if (length(idx)) {
        idx <- idx[1L]
        spd <- spd[, c(idx, seq_len(ncol(spd))[-idx])]
    }

    ## Determine which columns contain list-like data (i.e. multiple entries).
    mult <- colnames(spd)[!vapply(spd, function(z) is.vector(z) & !is.list(z),
                                  logical(1))]
    for (m in mult) {
        spd[, m] <- vapply(
            spd[, m], function(z) paste0(z, collapse = paste0("\n", m, ": ")),
            character(1))
    }

    tmp <- lapply(colnames(spd), function(z) {
        paste0(z, ": ", spd[, z], "\n")
    })

    pks <- vapply(peaksData(x), function(z)
        paste0("Num Peaks: ", nrow(z), "\n",
               paste0(paste0(z[, 1], " ", z[, 2], "\n"), collapse = ""),
               collapse = ""),
        character(1))
    tmp <- do.call(cbind, c(tmp, list(pks)))
    tmp[grep(": NA\n", tmp, fixed = TRUE)] <- ""
    writeLines(apply(tmp, 1, paste0, collapse = ""), con = con)
}

#' @title Parse the comment field from a MoNA MSP file
#'
#' @description
#'
#' Parse comment field from MoNA.
#' 
#' @author Johannes Rainer
#'
#' @noRd
parseMoNaComment <- function(x, names = c("InChI", "author", "SMILES",
                                          "date", "cas", "kegg",
                                          "pubchem cid")) {
    ## extract value between "<name>= and ".
    names(names) <- names
    as.data.frame(lapply(names, function(z) {
        tmp <- sub(paste0(".*?\"", z, "=(.*?)\".*"), "\\1", x, perl = TRUE)
        tmp[tmp == x] <- NA_character_
        tmp
    }))
}


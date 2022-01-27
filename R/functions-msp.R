#' @title Reading MSP files
#'
#' The `readMsp` function imports the data from a file in MGF format reading
#' all specified fields and returning the data as a [DataFrame()].
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
#' @importFrom S4Vectors DataFrame
#'
#' @importFrom IRanges NumericList
#'
#' @importFrom MsCoreUtils rbindFill
#'
#' @importFrom methods as
#'
#' @author Laurent Gatto, Steffen Neumann, Johannes Rainer
#'
#' @examples
#'
#' fls <- dir(system.file("extdata", package = "MsBackendMsp"),
#'     full.names = TRUE, pattern = "msp$")[1L]
#'
#' readMsp(fls)
readMsp <- function(f, msLevel = 2L,
                    mapping = spectraVariableMapping(MsBackendMsp()), ...) {
    if (length(f) != 1L)
        stop("Please provide a single msp file.")
    
    msp <- scan(file = f, what = "",
                sep = "\n", quote = "",
                allowEscapes = FALSE,
                quiet = TRUE)
    
    ## Ignore comments
    cmts <- grep("^[#]", msp)
    if (length(cmts))
        msp <- msp[-cmts]

    ## Find individual records
    begin <- grep("^NAME:", msp, ignore.case = TRUE)
    end <- c(begin[-1] -1L, length(msp))
    
    n <- length(begin)
    sp <- vector("list", length = n)

    for (i in seq_along(sp))
        sp[[i]] <- .extract_msp_spectrum(msp[begin[i]:end[i]],
                                         mapping = mapping)
    res <- DataFrame(rbindFill(sp))

    spv <- Spectra:::.SPECTRA_DATA_COLUMNS
    spv <- spv[!names(spv) %in% c("mz", "intensity")]
    for (i in seq_along(res)) {
        if (all(lengths(res[[i]]) == 1))
            res[[i]] <- unlist(res[[i]])
        if (any(col <- names(spv) == colnames(res)[i]))
            res[[i]] <- as(res[[i]], spv[col][1])
    }

    res$mz <- IRanges::NumericList(res$mz, compress = FALSE)
    res$intensity <- IRanges::NumericList(res$intensity, compress = FALSE)
    res$dataOrigin <- f
    if (!any(colnames(res) == msLevel))
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
#' @noRd
.extract_msp_spectrum <- function(msp, mapping) {
    ## grep description
    desc.idx <- grep(":", msp)
    desc <- msp[desc.idx]
    spec <- msp[-desc.idx]

    ms <- do.call(rbind, strsplit(sub("^(\\t|[[:space:]]+)", "", spec),
                                  "[[:space:]]+"))
    mode(ms) <- "double"

    if (!length(ms))
        ms <- matrix(numeric(), ncol = 2L)
    else if (is.unsorted(ms[, 1L]))
        ms <- ms[order(ms[, 1L]), ]
    
    r <- regexpr(":", desc, fixed = TRUE)
    desc <- setNames(substring(desc, r + 2L, nchar(desc)),
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

    res$mz = list(ms[, 1L])
    res$intensity = list(ms[, 2L])
    res
}

.process_polarity <- function(x) {
    if (grepl("^(p|\\+)", x, ignore.case = TRUE))
        return(1L)
    if (grepl("^(n|-)", x, ignore.case = TRUE))
        return(0L)
    -1L
}

#' @description
#'
#' Function to export a `Spectra` object in .msp format to `con`.
#'
#' @param x `Spectra`
#'
#' @param con output file.
#'
#' @param mapping named `character` vector that maps from `spectraVariables`
#'    (i.e. `names(mapping)`) to the variable name that should be used in the
#'    MGF file.
#'
#' @author Michael Witting
#'
#' @importMethodsFrom Spectra spectraVariables spectraNames
#'
#' @importMethodsFrom Spectra spectraData peaksData
#'
#' @noRd
.export_msp <- function(x, con = stdout(), mapping = spectraVariableMapping()) {
    
    if (class(con) == "character" && file.exists(con)) {
        
        message("Overwriting ", con, "!")
        unlink(con)
        
    }
    
    if (class(con)[1] == "character") {
        con <- file(description = con, open = "at")
        on.exit(close(con))
    }
    
    # custom cat function for writing of content
    .cat <- function(..., file = con, sep = " ", append = TRUE) {
        cat(..., file = file, sep = sep, append = append)
    }
    
    
    # iterate over all spectra
    for(i in 1:length(x)) {
        
        spv <- spectraVariables(x[i])
        spd <- spectraData(x[i], spv[!(spv %in% c("dataOrigin", "dataStorage"))])
        idx <- match(colnames(spd), names(mapping))
        colnames(spd)[!is.na(idx)] <- mapping[idx[!is.na(idx)]]
        
        spp <- peaksData(x[i])
        
        # here list with stuff in right order
        entries <- .getEntries()
        
        for(entry in entries) {
            
            #print(entry)
            
            if(entry %in% colnames(spd)) {
                
                value <- spd[entry][[1]]
                
                .cat(entry, value, "\n")
            
            }
        }
        
        .cat("Num Peaks:", length(peaksData(x[i])[[1]][,1]), "\n")
        
        .cat(paste0(peaksData(x[i])[[1]][,1],
                    " ",
                    peaksData(x[i])[[1]][,2],
                    collapse = "\n"))
        
        .cat("\n\n\n")
    }
}


#' This list defines the order and fields used for export
#'
#' @noRd
.getEntries <- function() {
    
    c(
        # record specific information
        "NAME:",
        "DB#:",
        "INCHIKEY:",
        "PRECURSORTYPE:",
        "PRECURSORMZ:",
        "RETENTIONTIME:",
        "EXACTMASS:",
        "FORMULA:"
        )
}



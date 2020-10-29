##' @param f `character(1)` with the path to an msp file.
##' 
##' @param msLevel `numeric(1)` with the MS level. Default is 2.
##' 
##' @param ... Additional parameters, currently ignored.
##'
##' @importFrom S4Vectors DataFrame
##'
##' @importFrom IRanges NumericList
##' 
##' @author Laurent Gatto
##' 
##' @noRd

if (FALSE) {
    library(MSnbase)
    library(microbenchmark)
    f <- "/vol/bioinvindex/Submissions/MIC044/MTBLS1582/20203111526_spectra_SMpos.msp"
    f <- "/vol/R/BioC/devel/Spectra/MassBank_MSMS_Pos_Rev173_vs1.msp"
    f <- "/vol/R/BioC/devel/2018-03-09_14_10_01_pos_387850_nist_msms_HR.MSP"
    f <- "/vol/R/BioC/devel/MsBackendMsp/inst/extdata/Spectrum.msp"
    f <- "/vol/R/BioC/devel/MsBackendMsp/inst/extdata/Spectrum2.msp"
    f <- "/vol/R/BioC/devel/MsBackendMsp/inst/extdata/msdial_pos.msp"
    f <- "/vol/R/BioC/devel/MsBackendMsp/inst/extdata/first-export-LipidBlast.msp"
    massbank <- .read_msp(f)
    orgmassbank <- ReadMspFile(f)
    
    microbenchmark(.read_msp(f), times = 5)
}


.read_msp <- function(f, msLevel = 2L, ...) {
    if (length(f) != 1L)
        stop("Please provide a single msp file.")
    
    msp <- scan(file = f, what = "",
                sep = "\n", quote = "",
                allowEscapes = FALSE,
                quiet = TRUE
                )
    
    ## Ignore comments
    cmts <- grep("^[#]", msp)
    if (length(cmts))
        msp <- msp[-cmts]

    ## Find individual records
    begin <- grep("^NAME:", msp, ignore.case = TRUE)
    stopifnot(!isEmpty(begin)) ## NO! DONT STOP ON ERROR! FIND OUT HOW EXCEPTIONS WORK!
    
    end <- c(begin[-1], length(msp))
    
    n <- length(begin)
    sp <- vector("list", length = n)

    for (i in seq(along = sp)) 
        sp[[i]] <- .extract_msp_spectrum(msp[begin[i]:end[i]])

    res <- DataFrame(do.call(rbind, sp))

    for (i in seq_along(res)) {
        if (all(lengths(res[[i]]) == 1))
            res[[i]] <- unlist(res[[i]])
    }

    res$mz <- IRanges::NumericList(res$mz)
    res$intensity <- IRanges::NumericList(res$intensity)
    res$dataOrigin <- f
    res$msLevel <- as.integer(msLevel)
    res
}

##' @param msp `character()` of lines defining a spectrum in msp
##'     format.
##' 
##' @author Laurent Gatto
##' 
##' @importFrom stats setNames
##'
##' @noRd
.extract_msp_spectrum <- function(msp) {
    ## grep description
    desc.idx <- grep(":", msp)
    desc <- msp[desc.idx]
    spec <- msp[-desc.idx]

    ms <- do.call(rbind, strsplit(spec, "[[:space:]]+"))
    mode(ms) <- "double"

    if (!length(ms))
        ms <- matrix(numeric(), ncol = 2L)
    
    if(length(ms) > 2)
        ms <- ms[order(ms[, 1L]),]
    
    r <- regexpr(":", desc, fixed = TRUE)
    desc <- setNames(substring(desc, r + 2L, nchar(desc)), tolower(substring(desc, 1L, r - 1L)))
    name <- unname(desc["name"])
    accession <- unname(desc["db#"])
    formula <- unname(desc["formula"])
    inchikey <- unname(desc["inchikey"])
    adduct <- unname(desc["precursor_type"])
    exactmass <- as.numeric(unname(desc["exactmass"]))

    ## select only values of interest and convert to numeric
    voi <- c("retentiontime", "ionmode", "precursormz")
    desc <- setNames(as.numeric(desc[voi]), voi)
    desc[is.na(desc[voi])] <- 0L
    
    list(rtime = unname(desc["retentiontime"]),
         precursorMz = unname(desc["precursormz"]),
         precursorCharge = unname(as.integer(desc["ionmode"])),
         mz = ms[, 1L],
         intensity = ms[, 2L],
         accession = accession,
         name = name,
         formula = formula,
         inchikey = inchikey,
         adduct = adduct,
         exactmass = exactmass)
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
#' @importMethodsFrom Spectra spectraVariables spectraNames spectraData
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


#'
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



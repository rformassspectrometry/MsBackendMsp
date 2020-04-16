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
    begin <- grep("^NAME:", msp)
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

    r <- regexpr(":", desc, fixed = TRUE)
    desc <- setNames(substring(desc, r + 1L, nchar(desc)), substring(desc, 1L, r - 1L))
    name <- unname(desc["NAME"])
    
    ## select only values of interest and convert to numeric
    voi <- c("RETENTIONTIME", "IONMODE", "PRECURSORMZ")
    desc <- setNames(as.numeric(desc[voi]), voi)
    desc[is.na(desc[voi])] <- 0L
    
    list(rtime = unname(desc["RETENTIONTIME"]),
         precursorMz = unname(desc["PRECURSORMZ"]),
         precursorCharge = unname(as.integer(desc["IONMODE"])),
         mz = ms[, 1L],
         intensity = ms[, 2L],
         name = name)
}

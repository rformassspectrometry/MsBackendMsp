test_that("readMsp works", {
    ## MS-DIAL
    f <- system.file("extdata", "msdial_pos.msp", package = "MsBackendMsp")

    res <- readMsp(f, mapping = spectraVariableMapping(MsBackendMsp()))
    expect_s4_class(res, "DataFrame")
    expect_true(all(c("name", "precursorMz", "adduct", "rtime") %in%
                    colnames(res)))
    expect_s4_class(res$mz, "NumericList")
    expect_s4_class(res$intensity, "NumericList")
    expect_true(length(res$mz[1L]) == 1)

    f <- system.file("extdata", "spectrum.msp", package = "MsBackendMsp")
    f2 <- system.file("extdata", "spectrum2.msp", package = "MsBackendMsp")

    res <- readMsp(f, mapping = spectraVariableMapping(MsBackendMsp()))
    res2 <- readMsp(f2, mapping = spectraVariableMapping(MsBackendMsp()))
    expect_true(all(colnames(res) %in% colnames(res2)))
    res <- res[, colnames(res2)]
    res$msLevel <- NA_integer_
    res$dataOrigin <- "a"
    res2$dataOrigin <- "a"
    expect_equal(res, res2[1L, ])
    
    expect_error(readMsp(c(f, f2)), "Please provide a single msp file.")
})

test_that("readMsp works for all test files", {
    ## wrong header: not name: value pairs
    f <- system.file("extdata", "fail_2_spectrum.msp", package = "MsBackendMsp")
    expect_error(readMsp(f), "MSP format error")
    
    f <- system.file("extdata", "spectrum.msp", package = "MsBackendMsp")
    res <- readMsp(f)
    expect_true(is.numeric(res$rtime))
    expect_true(is.integer(res$msLevel))
    expect_true(is(res$mz, "NumericList"))
    expect_true(is(res$intensity, "NumericList"))
    
    f <- system.file("extdata", "fail_spectrum.msp", package = "MsBackendMsp")
    expect_warning(res_2 <- readMsp(f), "Unexpected")
    expect_equal(res$intensity, res_2$intensity)
    expect_equal(res$mz, res_2$mz)

    ## spectra entries are not separated by an empty line
    f <- system.file("extdata", "fail_spectrum2.msp", package = "MsBackendMsp")
    expect_error(readMsp(f), "MSP format error")

    f <- system.file("extdata", "first-export-LipidBlast.msp",
                     package = "MsBackendMsp")
    res <- readMsp(f)
    expect_true(nrow(res) == 1)
    expect_equal(lengths(res$mz), 2)
    expect_equal(lengths(res$intensity), 2)
    expect_true(is.integer(res$msLevel))
    
    f <- system.file("extdata", "minimona.msp", package = "MsBackendMsp")
    res <- readMsp(f)
    expect_true(nrow(res) > 1)
    expect_equal(length(res$mz[[1L]]), 8)
    expect_equal(length(res$intensity[[1L]]), 8)

    f <- system.file("extdata", "msdial_pos.msp", package = "MsBackendMsp")
    res <- readMsp(f)
    expect_true(nrow(res) > 1)
    expect_equal(length(res$mz[[1L]]), 1)
    expect_equal(length(res$intensity[[1L]]), 1)
    expect_equal(length(res$mz[[4L]]), 8)
    expect_equal(length(res$intensity[[4L]]), 8)

    f <- system.file("extdata", "small-export-LipidBlast.msp",
                     package = "MsBackendMsp")
    res <- readMsp(f)
    expect_true(nrow(res) == 5)
    expect_equal(length(res$mz[[1L]]), 2)
    expect_equal(length(res$intensity[[1L]]), 2)
    expect_equal(length(res$mz[[5L]]), 12)
    expect_equal(length(res$intensity[[5L]]), 12)

    f <- system.file("extdata", "spectrum2.msp", package = "MsBackendMsp")
    res <- readMsp(f)
    expect_true(nrow(res) == 2)
    expect_equal(length(res$intensity[[1L]]), 86)
    expect_equal(length(res$mz[[1L]]), 86)
    expect_true(is.integer(res$msLevel))
    expect_equal(res$msLevel, c(NA_integer_, 1L))
})

test_that(".extract_msp_spectrum works", {
    f <- system.file("extdata", "msdial_pos.msp", package = "MsBackendMsp")
    x <- scan(file = f, what = "",
              sep = "\n", quote = "",
              allowEscapes = FALSE,
              quiet = TRUE)
    x2 <- x[1:37]
    x <- x[15:22]

    mapping <- spectraVariableMapping(MsBackendMsp())
    res <- .extract_msp_spectrum(x, mapping = mapping, fixupNISTRI=FALSE)
    expect_true(is.data.frame(res))
    expect_true(nrow(res) == 1)

    ## Duplicated values
    x <- c("synonym: a", "synonym: b", "other: 1", "other: 2", "other: 3", x)
    res <- .extract_msp_spectrum(x, mapping = mapping, fixupNISTRI=FALSE)
    expect_true(sum(colnames(res) == "other") == 1)
    expect_true(sum(colnames(res) == "synonym") == 1)
    expect_true(length(res$other[[1L]]) == 3)
    expect_true(length(res$synonym[[1L]]) == 2)

    ## Spectra not separated by blank lines.
    expect_error(.extract_msp_spectrum(x2, mapping), "MSP format")

    f <- system.file("extdata", "fail_spectrum.msp", package = "MsBackendMsp")
    x <- scan(file = f, what = "",
              sep = "\n", quote = "",
              allowEscapes = FALSE,
              quiet = TRUE)
    expect_warning(.extract_msp_spectrum(x, mapping = mapping, fixupNISTRI=FALSE),
                   "Unexpected number of values")
})

test_that(".process_polarity works", {
    expect_equal(.process_polarity("some"), -1L)
    expect_equal(.process_polarity("Pos"), 1L)
    expect_equal(.process_polarity("+"), 1L)
    expect_equal(.process_polarity("Neg"), 0L)
    expect_equal(.process_polarity("-"), 0L)

    expect_equal(.process_polarity(1, input = FALSE), "Positive")
    expect_equal(.process_polarity(0, input = FALSE), "Negative")
    expect_equal(.process_polarity(2, input = FALSE), NA_character_)
})

test_that(".process_mslevel works", {
    expect_equal(.process_mslevel("MS 2"), 2L)
    expect_equal(.process_mslevel("ms1"), 1L)
    expect_equal(.process_mslevel(1, input = FALSE), "MS1")
})

test_that(".export_msp works", {
    spd <- DataFrame(msLevel = c(2L, 1L, 2L), rtime = c(1, 2, 3.2))
    spd$mz <- list(c(12, 14, 45, 56), c(14.1, 34, 56.1), c(12.1, 14.15, 34.1))
    spd$intensity <- list(c(10, 20, 30, 40), c(11, 21, 31), c(12, 22, 32))
    sps <- Spectra(spd)

    expect_output(.export_msp(sps[1], exportName = FALSE), "^msLevel: MS2")
    expect_output(.export_msp(sps[1], exportName = FALSE), "Num Peaks: 4")

    expect_output(.export_msp(sps[1], exportName = TRUE), "^NAME: 1")

    sps$multi <- list("a", c("b", "c"), "e")
    expect_output(.export_msp(sps), "multi: b\\nmulti: c")

    map <- c(msLevel = "MSL", rtime = "TIME", multi = "OTHER")
    expect_output(.export_msp(sps, mapping = map), "OTHER: b\\nOTHER: c")
    expect_output(.export_msp(sps, mapping = map), "TIME: 1")
    expect_output(.export_msp(sps, mapping = map), "TIME: 2")
    expect_output(.export_msp(sps, mapping = map), "TIME: 3.2")    
})

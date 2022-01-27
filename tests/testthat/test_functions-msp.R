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
    res$dataOrigin <- "a"
    res2$dataOrigin <- "a"
    expect_equal(res, res2[1L, ])
    
    expect_error(readMsp(c(f, f2)), "Please provide a single msp file.")
})

test_that(".expect_msp_spectrum works", {
    f <- system.file("extdata", "msdial_pos.msp", package = "MsBackendMsp")
    x <- scan(file = f, what = "",
              sep = "\n", quote = "",
              allowEscapes = FALSE,
              quiet = TRUE)
    x <- x[15:22]

    mapping <- spectraVariableMapping(MsBackendMsp())
    res <- .extract_msp_spectrum(x, mapping = mapping)
    expect_true(is.data.frame(res))
    expect_true(nrow(res) == 1)

    ## Duplicated values
    x <- c(x, "synonym: a", "synonym: b", "other: 1", "other: 2", "other: 3")
    res <- .extract_msp_spectrum(x, mapping = mapping)
    expect_true(sum(colnames(res) == "other") == 1)
    expect_true(sum(colnames(res) == "synonym") == 1)
    expect_true(length(res$other[[1L]]) == 3)
    expect_true(length(res$synonym[[1L]]) == 2)
})

test_that(".process_polarity works", {
    expect_equal(.process_polarity("some"), -1L)
    expect_equal(.process_polarity("Pos"), 1L)
    expect_equal(.process_polarity("+"), 1L)
    expect_equal(.process_polarity("Neg"), 0L)
    expect_equal(.process_polarity("-"), 0L)
})


test_that(".read_msp works", {
    fls <- dir(system.file("extdata", package = "MsBackendMsp"),
               full.names = TRUE, pattern = "msp$")

    expect_error(.read_msp(fls), "Please provide a single msp file.")

    res1 <- .read_msp(fls[1])
    res2 <- .read_msp(fls[2])

    # cns <- c("rtime", "scanIndex", "precursorMz", "precursorIntensity",
    #          "precursorCharge", "mz", "intensity", "title", "dataOrigin",
    #          "msLevel")
    # expect_identical(sort(names(res1)), sort(cns))
    # expect_identical(sort(names(res1)), sort(names(res2)))
    # expect_true(is(res1$mz, "NumericList"))
    # expect_true(is(res1$intensity, "NumericList"))
    # expect_equal(length(res1$intensity[[1]]), length(res1$mz[[1]]))
    # expect_equal(length(res1$intensity[[2]]), length(res1$mz[[2]]))
    # expect_equal(length(res1$intensity[[3]]), length(res1$mz[[3]]))
    # 
    # expect_identical(res1$title,
    #                  c("File193 Spectrum1719 scans: 2162",
    #                    "File193 Spectrum1944 scans: 2406",
    #                    "File193 Spectrum1968 scans: 2432"))
})

test_that(".read_lipidblast_msp works", {
    f <- system.file("extdata/small-export-LipidBlast.msp", package = "MsBackendMsp")
    res1 <- .read_msp(f)
    expect(all(dim(res1)==c(5,15)), "Dimensionality should be 5,15")
})

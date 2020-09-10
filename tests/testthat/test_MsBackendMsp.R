test_that("backendInitialize,MsBackendMsp works", {
    fls <- dir(system.file("extdata", package = "MsBackendMsp"),
               full.names = TRUE, pattern = "msp$")
    be <- MsBackendMsp()

    ## Import a single file.
    res1 <- backendInitialize(be, fls[1])
    n1 <- length(res1) ## 3

    # expect_identical(length(res1), n1)
    # expect_identical(res1$dataStorage, rep("<memory>", n1))
    # expect_identical(res1$dataOrigin, rep(normalizePath(fls[1]), n1))
    # expect_identical(res1$msLevel, rep(2L, n1))
    # expect_identical(lengths(res1$mz), c(14L, 21L, 14L))

    res2 <- backendInitialize(be, fls[2])
    n2 <- length(res2) ## 4
    
    ## Import multiple files.
    res_all <- backendInitialize(be, fls)
    # expect_true(length(res_all) == n1 + n2)
    # expect_identical(res_all[1]$mz, res1[1]$mz)
    # expect_identical(res_all[n1 + 1]$mz, res2[1]$mz)
    # expect_true(all(res_all$msLevel == 2L))
    # expect_identical(res_all$dataOrigin,
    #                  c(rep(normalizePath(fls[1]), n1),
    #                    rep(normalizePath(fls[2]), n2)))
    # expect_true(is.integer(res_all@spectraData$msLevel))

    ## TODO: Import with failing file.
    ## TODO: Import with failing file and nonStop = TRUE
    
    ## errors
    expect_error(backendInitialize(be), "'files' is mandatory")
    expect_error(backendInitialize(be, 4), "expected to be a character")
    expect_error(suppressWarnings(backendInitialize(be, "a")), "a not found")
})

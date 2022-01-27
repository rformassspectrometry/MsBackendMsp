test_that("backendInitialize,MsBackendMsp works", {
    fls <- system.file("extdata", "small-export-LipidBlast.msp",
                       package = "MsBackendMsp")
    be <- MsBackendMsp()

    ## Import a LipidBlast file.
    res <- backendInitialize(be, fls)
    expect_s4_class(res, "MsBackendMsp")
    expect_true(length(res) == 5L)
    expect_true(all(res$msLevel == 2L))
    expect_true(all(is.na(res$precursorMz)))

    ## Import MoNa
    f <- system.file("extdata", "minimona.msp", package = "MsBackendMsp")
    res <- backendInitialize(be, f)
    expect_s4_class(res, "MsBackendMsp")
    expect_true(length(res) == 30L)
    expect_true(all(res$msLevel == 2L))
    expect_true(all(is.na(res$precursorMz)))
    expect_true(is.list(res$Synon))
    
    ## Import MoNa and LipidBlast file
    res <- backendInitialize(be, c(fls, f))
    expect_s4_class(res, "MsBackendMsp")
    expect_true(length(res) == 35L)
    expect_true(all(res$msLevel == 2L))
    expect_true(all(is.na(res$precursorMz)))
    expect_true(length(grep("MoNA", res$accession)) == 30L)
    expect_true(length(grep("LipidBlast", res$accession)) == 5L)
    
    ## errors
    expect_error(backendInitialize(be), "'file' is mandatory")
    expect_error(backendInitialize(be, 4), "expected to be a character")
    expect_error(suppressWarnings(backendInitialize(be, "a")), "a not found")
})

test_that("spectraVariableMapping works", {
    res <- spectraVariableMapping(MsBackendMsp())
    expect_true(is.character(res))
    expect_true(length(res) > 0)
    expect_error(spectraVariableMapping(MsBackendMsp(), "other"), "should be")
})

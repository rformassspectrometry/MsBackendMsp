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

    res <- backendInitialize(be, fls,
                             mapping = spectraVariableMapping(be, "mona"))
    expect_true(all(!is.na(res$precursorMz)))
    expect_equal(polarity(res), c(1L, 1L, 1L, 0L, 0L))
    
    ## Import MoNa
    f <- system.file("extdata", "minimona.msp", package = "MsBackendMsp")
    res <- backendInitialize(be, f)
    expect_s4_class(res, "MsBackendMsp")
    expect_true(length(res) == 30L)
    expect_true(all(res$msLevel == 2L))
    expect_true(all(is.na(res$precursorMz)))
    expect_true(is.list(res$Synon))
        
    ## Import MoNa and LipidBlast file
    res <- backendInitialize(be, c(fls, f), mapping = spectraVariableMapping(MsBackendMsp(), "mona"))
    expect_s4_class(res, "MsBackendMsp")
    expect_true(length(res) == 35L)
    expect_true(all(res$msLevel == 2L))
    expect_equal(length(which(is.na(res$precursorMz))), 1)
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

    ## Import MoNa 2023,
    ## https://github.com/rformassspectrometry/MsBackendMsp/issues/14

    f <- system.file("extdata", "mini-MoNA-export-LC-MS-MS_Positive_Mode_20231016.msp", package = "MsBackendMsp")
    mona <- Spectra(f, source = MsBackendMsp(), mapping = spectraVariableMapping(MsBackendMsp(), "mona"))
    expect_true(all(!is.na(precursorMz(mona)) &  is.numeric(precursorMz(mona))))

})


test_that("RI parsing works", {
    f <- system.file("extdata", "mini-lib2nist.msp", package = "MsBackendMsp")
    nist <- Spectra(f, source = MsBackendMsp())
    expect_equal(as.numeric(nist$RI), c(572,2674))
})

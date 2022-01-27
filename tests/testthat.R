library(testthat)
library(MsBackendMsp)

## Run specific (backend) package tests
test_check("MsBackendMsp")

## Run additional tests from Spectra:
test_suite <- system.file("test_backends", "test_MsBackend",
                          package = "Spectra")

fls <- system.file("extdata", "minimona.msp", package = "MsBackendMsp")
be <- MsBackendMsp()
be <- backendInitialize(
    be, fls, mapping = spectraVariableMapping(MsBackendMsp(), "mona"))

res <- test_file(paste0(test_suite, "/test_spectra_variables.R"),
                 reporter = check_reporter(), stop_on_failure = TRUE)




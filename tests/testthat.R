library(testthat)
library(MsBackendMsp)

## Run specific (backend) package tests
test_check("MsBackendMsp")

## Run generic tests defined in Spectra package 
## for the msmslibraries profile

# Initialize backend and test data

fls <- dir(system.file("extdata", package = "MsBackendMsp"),
           full.names = TRUE, pattern = "msp$")
be <- MsBackendMsp()





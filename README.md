# Mass Spectrometry Data Backend for NIST MSP Files

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/RforMassSpectrometry/MsBackendMsp/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RforMassSpectrometry/MsBackendMsp/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov](https://codecov.io/gh/rformassspectrometry/MsBackendMsp/branch/main/graph/badge.svg?token=RS145O5ER1)](https://codecov.io/gh/rformassspectrometry/MsBackendMsp)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)
[![years in bioc](http://bioconductor.org/shields/years-in-bioc/MsBackendMsp.svg)](https://bioconductor.org/packages/release/bioc/html/MsBackendMsp.html)
[![Ranking by downloads](http://bioconductor.org/shields/downloads/release/MsBackendMsp.svg)](https://bioconductor.org/packages/stats/bioc/MsBackendMsp/)
[![build release](http://bioconductor.org/shields/build/release/bioc/MsBackendMsp.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/MsBackendMsp/)
[![build devel](http://bioconductor.org/shields/build/devel/bioc/MsBackendMsp.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/MsBackendMsp/)

The `MsBackendMsp` package provides functionality to import and handle
MS/MS spectrum data from MSP files.
The package defines the `MsBackendMsp` backend which can be used to
import and use MS2 spectrum data from msp files with the
[Spectra](https://github.com/rformassspectrometry/Spectra) R package.

For more information see the package vignette or the package
[homepage](https://rformassspectrometry.github.io/MsBackendMsp).


# Installation

The package can be installed with

```r
install.packages("BiocManager")
BiocManager::install("MsBackendMsp")
```


# Contributions

Contributions are highly welcome and should follow the [contribution
guidelines](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#contributions).
Also, please check the coding style guidelines in the [RforMassSpectrometry
vignette](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html).


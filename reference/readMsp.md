# Reading MSP files

The `readMsp()` function imports the data from a file in MGF format
reading all specified fields and returning the data as a
[`S4Vectors::DataFrame()`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html).

Format constraints for MSP files:

- Comment lines are expected to start with a `#`.

- Multiple spectra within the same MSP file are separated by an empty
  line.

- The first n lines of a spectrum entry represent metadata.

- Metadata is provided as "name: value" pairs (i.e. name and value
  separated by a ":").

- One line per mass peak, with values separated by a whitespace or
  tabulator.

- Each line is expected to contain at least the m/z and intensity values
  (in that order) of a peak. Additional values are currently ignored.

## Usage

``` r
readMsp(
  f,
  msLevel = 2L,
  mapping = spectraVariableMapping(MsBackendMsp()),
  BPPARAM = SerialParam(),
  return.type = c("DataFrame", "data.frame"),
  ...
)
```

## Arguments

- f:

  `character(1)` with the path to an MSP file.

- msLevel:

  `numeric(1)` with the MS level. Default is 2. This value will be
  reported as the spectra's MS level **unless** the source MSP file
  defines the MS level.

- mapping:

  named `character` vector to rename MSP fields to spectra variables
  (see `spectraVariableMapping()` help). This allows to correctly import
  also custom fields or data from files with different MSP *flavors*.

- BPPARAM:

  parallel processing setup. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more details.

- return.type:

  `character(1)` defining whether the results should be returned as a
  `DataFrame` (`return.type = "DataFrame"`, the default) or a
  `data.frame` (`return.type = "data.frame"`).

- ...:

  Additional parameters, currently ignored.

## Value

A `DataFrame` with each row containing the data from one spectrum in the
MSP file. m/z and intensity values are available in columns `"mz"` and
`"intensity"` in a list representation.

## Author

Laurent Gatto, Steffen Neumann, Johannes Rainer

## Examples

``` r

f <- system.file("extdata", "minimona.msp", package = "MsBackendMsp")

readMsp(f)
#> DataFrame with 30 rows and 19 columns
#>                       Name                                            Synon
#>                <character>                                           <list>
#> 1                Ritonavir                                    $:00in-source
#> 2                  Unknown                                    $:00in-source
#> 3                  Unknown                                    $:00in-source
#> 4                  Unknown                                    $:00in-source
#> 5                  Unknown                                    $:00in-source
#> ...                    ...                                              ...
#> 26  7-MAHPDA [DMED-FAHFA.. $:00 ms2,$:05 30V CID,$:07 In-Silico-Spect..,...
#> 27  8-MAHPDA [DMED-FAHFA.. $:00 ms2,$:05 30V CID,$:07 In-Silico-Spect..,...
#> 28  9-MAHPDA [DMED-FAHFA.. $:00 ms2,$:05 30V CID,$:07 In-Silico-Spect..,...
#> 29  10-MAHPDA [DMED-FAHF.. $:00 ms2,$:05 30V CID,$:07 In-Silico-Spect..,...
#> 30  11-MAHPDA [DMED-FAHF.. $:00 ms2,$:05 30V CID,$:07 In-Silico-Spect..,...
#>       accession               InChIKey  Instrument_type      Formula
#>     <character>            <character>      <character>  <character>
#> 1    MoNA000010 NCDNCNXCDXHOMX-XGKFQ.. Waters Synapt G2 C37H48N6O5S2
#> 2    MoNA000012 MXNRLFUSFKVQSK-QMMMG..               NA          H2O
#> 3    MoNA000013 JFLIEFSWGNOPJJ-JTQLQ..               NA   C13H16N2O4
#> 4    MoNA000014 MXNRLFUSFKVQSK-QMMMG..               NA    C9H20N2O2
#> 5    MoNA000015 JFLIEFSWGNOPJJ-JTQLQ..               NA   C13H16N2O4
#> ...         ...                    ...              ...          ...
#> 26   MoNA011542 MQTXDUJVZFSMBT-UHFFF..               NA           NA
#> 27   MoNA011543 CYZGZMTYUCCRMO-UHFFF..               NA           NA
#> 28   MoNA011544 GREOXQBULJWSQV-UHFFF..               NA           NA
#> 29   MoNA011545 GPDSFBGRMDEPAY-UHFFF..               NA           NA
#> 30   MoNA011546 QGMGRNCZWXKENS-UHFFF..               NA           NA
#>              MW         ExactMass               Comments   Num.Peaks
#>     <character>       <character>            <character> <character>
#> 1           720 720.3127606360001 "computed SMILES=OC(..           8
#> 2            18      18.010564684 "computed SMILES=O" ..          29
#> 3           264     264.111006992 "computed SMILES=O=C..          21
#> 4           188      188.15247788 "computed SMILES=O=C..          10
#> 5           264     264.111006992 "computed SMILES=O=C..          18
#> ...         ...               ...                    ...         ...
#> 26           NA                NA "molecular weight=53..           4
#> 27           NA                NA "molecular weight=53..           4
#> 28           NA                NA "molecular weight=53..           4
#> 29           NA                NA "molecular weight=53..           4
#> 30           NA                NA "molecular weight=53..           4
#>                              mz                      intensity Precursor_type
#>                   <NumericList>                  <NumericList>    <character>
#> 1   140.054,171.096,197.075,...       18.018,18.018,31.031,...             NA
#> 2   41.0418,42.0375,43.0201,... 2.036603,1.781709,0.795269,...         [M+H]+
#> 3   41.0442,56.0522,65.0415,...    1.41478,2.52135,1.09140,...         [M+H]+
#> 4   41.0424,56.0532,65.0403,...    6.06061,7.02479,4.68320,...         [M+H]+
#> 5   42.0032,42.0105,67.0334,...    5.66221,3.25771,1.37677,...         [M-H]-
#> ...                         ...                            ...            ...
#> 26  266.248,311.306,494.457,...    60.0601,40.0400,30.0300,...         [M+H]+
#> 27  266.248,311.306,494.457,...    60.0601,40.0400,30.0300,...         [M+H]+
#> 28  266.248,311.306,494.457,...    60.0601,40.0400,30.0300,...         [M+H]+
#> 29  266.248,311.306,494.457,...    60.0601,40.0400,30.0300,...         [M+H]+
#> 30  266.248,311.306,494.457,...    60.0601,40.0400,30.0300,...         [M+H]+
#>     PrecursorMZ Collision_energy    Ion_mode Spectrum_type
#>     <character>      <character> <character>   <character>
#> 1            NA               NA          NA            NA
#> 2      189.1603               NA          NA            NA
#> 3      265.1188            35 eV          NA            NA
#> 4      265.1188            45 eV          NA            NA
#> 5      263.1031           -35 eV          NA            NA
#> ...         ...              ...         ...           ...
#> 26    539.51462               NA           P            NA
#> 27    539.51462               NA           P            NA
#> 28    539.51462               NA           P            NA
#> 29    539.51462               NA           P            NA
#> 30    539.51462               NA           P            NA
#>                 dataOrigin   msLevel
#>                <character> <integer>
#> 1   /__w/_temp/Library/M..         2
#> 2   /__w/_temp/Library/M..         2
#> 3   /__w/_temp/Library/M..         2
#> 4   /__w/_temp/Library/M..         2
#> 5   /__w/_temp/Library/M..         2
#> ...                    ...       ...
#> 26  /__w/_temp/Library/M..         2
#> 27  /__w/_temp/Library/M..         2
#> 28  /__w/_temp/Library/M..         2
#> 29  /__w/_temp/Library/M..         2
#> 30  /__w/_temp/Library/M..         2
```

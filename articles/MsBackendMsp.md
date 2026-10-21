# MsBackendMsp

**Package**:
*[MsBackendMsp](https://bioconductor.org/packages/3.24/MsBackendMsp)*\
**Authors**: Neumann Steffen \[aut\] (ORCID:
<https://orcid.org/0000-0002-7899-7192>), Johannes Rainer \[aut, cre\]
(ORCID: <https://orcid.org/0000-0002-6977-7147>), Michael Witting
\[ctb\] (ORCID: <https://orcid.org/0000-0002-1462-4426>)\
**Compiled**: Wed Oct 21 14:21:52 2026

## Introduction

The *[Spectra](https://bioconductor.org/packages/3.24/Spectra)* package
provides a central infrastructure for the handling of Mass Spectrometry
(MS) data. The package supports interchangeable use of different
*backends* to import MS data from a variety of sources (such as mzML
files). The
*[MsBackendMsp](https://bioconductor.org/packages/3.24/MsBackendMsp)*
package adds support for files in NIST MSP format which are frequently
used to share spectra libraries and hence enhances small compound
annotation workflows using the `Spectra` and
*[MetaboAnnotation](https://bioconductor.org/packages/3.24/MetaboAnnotation)*
packages (Rainer et al. 2022). This vignette illustrates the usage of
the *MsBackendMsp* package and how it can be used to import and export
data in MSP file format.

## Installation

To install this package, start `R` and enter:

\
`if`` ``(``!`[`requireNamespace`](https://rdrr.io/r/base/ns-load.html)`(``"BiocManager"``, quietly ``=`` ``TRUE``)``)`\
`    `[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"BiocManager"``)`\
\
`BiocManager``::`[`install`](https://bioconductor.github.io/BiocManager/reference/install.html)`(``"MsBackendMsp"``)`

This will install this package and all eventually missing dependencies.

## MSP file format

NIST MSP file format is supported. Some (eventually more stringent)
requirements of the format are:

- Comment lines are expected to start with a `#`.
- Multiple spectra within the same MSP file are separated by an empty
  line.
- The first n lines of a spectrum entry represent metadata.
- Metadata is provided as “name: value” pairs (i.e. name and value
  separated by a “:”).
- One line per mass peak, with values separated by a whitespace or
  tabulator.
- Each line is expected to contain at least the m/z and intensity values
  (in that order) of a peak. Additional values are currently ignored.

An MSP file can define/provide data for any number of spectra, with no
limit on the number of spectra, number of peaks per spectra or number of
metadata lines.

## Importing MS/MS data from MSP files

The MSP file format allows to store MS/MS spectra (m/z and intensity of
mass peaks) along with additional annotations for each spectrum. A
single MSP file can thus contain a single or multiple spectra. Below we
load the package and define the file name of an MSP file which is
distributed with this package.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`MsBackendMsp`](https://github.com/RforMassSpectrometry/MsBackendMsp)`)`\
\
`nist`` ``<-`` `[`system.file`](https://rdrr.io/r/base/system.file.html)`(``"extdata"``, ``"spectrum2.msp"``, package ``=`` ``"MsBackendMsp"``)`

We next import the data into a `Spectra` object by specifying in the
constructor function the *backend* class which can be used to read the
data (in our case a `MsBackendMsp`).

\
`sp`` ``<-`` `[`Spectra`](https://rdrr.io/pkg/Spectra/man/Spectra.html)`(``nist``, source ``=`` `[`MsBackendMsp`](https://rformassspectrometry.github.io/MsBackendMsp/reference/MsBackendMsp.md)`(``)``)`

With that we have now full access to all imported spectra variables that
can be listed with the
[`spectraVariables()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
function.

\
[`spectraVariables`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)`(``sp``)`

    ##  [1] "msLevel"                 "rtime"                  
    ##  [3] "acquisitionNum"          "scanIndex"              
    ##  [5] "dataStorage"             "dataOrigin"             
    ##  [7] "centroided"              "smoothed"               
    ##  [9] "polarity"                "precScanNum"            
    ## [11] "precursorMz"             "precursorIntensity"     
    ## [13] "precursorCharge"         "collisionEnergy"        
    ## [15] "isolationWindowLowerMz"  "isolationWindowTargetMz"
    ## [17] "isolationWindowUpperMz"  "name"                   
    ## [19] "adduct"                  "instrument_type"        
    ## [21] "instrument"              "smiles"                 
    ## [23] "inchikey"                "inchi"                  
    ## [25] "formula"                 "PUBCHEMID"              
    ## [27] "SOURCE"                  "COMMENT"                
    ## [29] "peak_count"

Besides default spectra variables, such as `msLevel`, `rtime`,
`precursorMz`, we also have additional spectra variables such as the
`name` or `adduct` that are additional data fields from the MSP file.

\
`sp``$``msLevel`

    ## [1] NA  1

\
`sp``$``name`

    ## [1] "3-Hydroxy-3-(2-(2-hydroxyphenyl)-2-oxoethyl)-1,3-dihydro-2H-indol-2-one"
    ## [2] "5-(4-Ethoxybenzylidene)-2-(4-morpholinyl)-1,3-thiazol-4(5H)-one"

\
`sp``$``adduct`

    ## [1] "[M+H]+" "[M+H]+"

The NIST file format is however only loosely defined and variety of
*flavors* (or *dialects*) exist which define their own data fields or
use different names for the fields. The `MsBackendMsp` supports data
import/export from all MSP format variations by defining and providing
different mappings between MSP data fields and spectra variables. Also
user-defined mappings can be used, which makes import from any MSP
flavor possible. Pre-defined mappings between MSP data fields and
spectra variables (i.e. variables within the `Spectra` object) are
returned by the
[`spectraVariableMapping()`](https://rdrr.io/pkg/Spectra/man/spectraVariableMapping.html)
function.

\
[`spectraVariableMapping`](https://rdrr.io/pkg/Spectra/man/spectraVariableMapping.html)`(`[`MsBackendMsp`](https://rformassspectrometry.github.io/MsBackendMsp/reference/MsBackendMsp.md)`(``)``)`

    ##             name        accession          formula         inchikey 
    ##           "NAME"            "DB#"        "FORMULA"       "INCHIKEY" 
    ##           adduct        exactmass            rtime      precursorMz 
    ##  "PRECURSORTYPE"      "EXACTMASS"  "RETENTIONTIME"    "PRECURSORMZ" 
    ##           adduct           smiles            inchi         polarity 
    ##  "PRECURSORTYPE"         "SMILES"          "INCHI"        "IONMODE" 
    ##       instrument         ontology  instrument_type            links 
    ##     "INSTRUMENT"       "ontology" "INSTRUMENTTYPE"          "LINKS" 
    ##          comment           splash       peak_count 
    ##        "Comment"         "Splash"      "Num Peaks"

The names of this `character` vector represent the spectra variable
names and the values of the vector the MSP data fields. Note that by
default, also all data fields for which no mapping is provided are
imported (with the field name being used as spectra variable name).

This default mapping works well for MSP files from NIST or from other
tools such as MS-DIAL. MassBank of North America
[MoNA](https://mona.fiehnlab.ucdavis.edu/) however, uses a slightly
different format. Below we read the first 6 lines of a MSP file from
MoNA.

\
`mona`` ``<-`` `[`system.file`](https://rdrr.io/r/base/system.file.html)`(``"extdata"``, ``"minimona.msp"``, package ``=`` ``"MsBackendMsp"``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(`[`readLines`](https://rdrr.io/r/base/readLines.html)`(``mona``)``)`

    ## [1] "Name: Ritonavir"                      
    ## [2] "Synon: $:00in-source"                 
    ## [3] "DB#: MoNA000010"                      
    ## [4] "InChIKey: NCDNCNXCDXHOMX-XGKFQTDJSA-N"
    ## [5] "Instrument_type: Waters Synapt G2"    
    ## [6] "Formula: C37H48N6O5S2"

The first 6 lines from a NIST MSP file:

\
[`head`](https://rdrr.io/r/utils/head.html)`(`[`readLines`](https://rdrr.io/r/base/readLines.html)`(``nist``)``)`

    ## [1] "NAME: 3-Hydroxy-3-(2-(2-hydroxyphenyl)-2-oxoethyl)-1,3-dihydro-2H-indol-2-one"
    ## [2] "PRECURSORMZ: 284.0917"                                                        
    ## [3] "PRECURSORTYPE: [M+H]+"                                                        
    ## [4] "INSTRUMENTTYPE: IT/ion trap"                                                  
    ## [5] "INSTRUMENT: Thermo Finnigan LCQ Deca"                                         
    ## [6] "SMILES: NA"

MSP files with MoNA flavor use slightly different field names, that are
also not all upper case, and also additional fields are defined. While
it is possible to import MoNA flavored MSP files using the default
variable mapping that was used above, most of the spectra variables
would however not mapped correctly to the respective spectra variable in
the resulting `Spectra` object (e.g. the precursor m/z would not be
available with the expected spectra variable `$precursorMz`).

The
[`spectraVariableMapping()`](https://rdrr.io/pkg/Spectra/man/spectraVariableMapping.html)
provides however also the mapping for MSP files with MoNA flavor.

\
[`spectraVariableMapping`](https://rdrr.io/pkg/Spectra/man/spectraVariableMapping.html)`(`[`MsBackendMsp`](https://rformassspectrometry.github.io/MsBackendMsp/reference/MsBackendMsp.md)`(``)``, ``"mona"``)`

    ##                  name               synonym             accession 
    ##                "Name"               "Synon"                 "DB#" 
    ##              inchikey                adduct           precursorMz 
    ##            "InChIKey"      "Precursor_type"         "PrecursorMZ" 
    ##              polarity               formula             exactmass 
    ##            "Ion_mode"             "Formula"           "ExactMass" 
    ## collision_energy_text               msLevel            instrument 
    ##    "Collision_energy"       "Spectrum_type"          "Instrument" 
    ##       instrument_type 
    ##     "Instrument_type"

Using this mapping in the data import will ensure that the fields get
correctly mapped.

\
`sp_mona`` ``<-`` `[`Spectra`](https://rdrr.io/pkg/Spectra/man/Spectra.html)`(``mona``, source ``=`` `[`MsBackendMsp`](https://rformassspectrometry.github.io/MsBackendMsp/reference/MsBackendMsp.md)`(``)``,`\
`                   mapping ``=`` `[`spectraVariableMapping`](https://rdrr.io/pkg/Spectra/man/spectraVariableMapping.html)`(`[`MsBackendMsp`](https://rformassspectrometry.github.io/MsBackendMsp/reference/MsBackendMsp.md)`(``)``, ``"mona"``)``)`\
`sp_mona``$``precursorMz`

    ##  [1]       NA 189.1603 265.1188 265.1188 263.1031 263.1031 229.1552 312.1302
    ##  [9] 525.4990 525.4990 525.4990 525.4990 525.4990 525.4990 525.4990 525.4990
    ## [17] 525.4990 525.4990 525.4990 525.4990 539.5146 539.5146 539.5146 539.5146
    ## [25] 539.5146 539.5146 539.5146 539.5146 539.5146 539.5146

Note that in addition to the predefined variable mappings, it is also
possible to provide any user-defined variable mapping with the `mapping`
parameter thus enabling to import from MSP files with a highly
customized format.

Multiple values for a certain spectrum are represented as duplicated
fields in an MSP file. The `MsBackendMsp` supports also import of such
data. MoNA MSP files use for example multiple `"Synon"` fields to list
all synonyms of a compound. Below we extract such values for two spectra
within our `Spectra` object from MoNA.

\
`sp_mona``[``29``:``30``]``$``synonym`

    ## [[1]]
    ## [1] "$:00 ms2"                "$:05 30V CID"           
    ## [3] "$:07 In-Silico-Spectrum" "$:00in-source"          
    ## 
    ## [[2]]
    ## [1] "$:00 ms2"                "$:05 30V CID"           
    ## [3] "$:07 In-Silico-Spectrum" "$:00in-source"

In addition to importing data from MSP files, `MsBackendMsp` allows also
to **export** `Spectra` to files in MSP format. Below we export for
example the `Spectra` with data from MoNA to a temporary file, using the
default NIST MSP format.

\
`tmpf`` ``<-`` `[`tempfile`](https://rdrr.io/r/base/tempfile.html)`(``)`\
\
[`export`](https://rdrr.io/pkg/Spectra/man/Spectra.html)`(``sp_mona``, backend ``=`` `[`MsBackendMsp`](https://rformassspectrometry.github.io/MsBackendMsp/reference/MsBackendMsp.md)`(``)``, file ``=`` ``tmpf``,`\
`       mapping ``=`` `[`spectraVariableMapping`](https://rdrr.io/pkg/Spectra/man/spectraVariableMapping.html)`(`[`MsBackendMsp`](https://rformassspectrometry.github.io/MsBackendMsp/reference/MsBackendMsp.md)`(``)``)``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(`[`readLines`](https://rdrr.io/r/base/readLines.html)`(``tmpf``)``)`

    ## [1] "NAME: Ritonavir"                      
    ## [2] "msLevel: MSNA"                        
    ## [3] "synonym: $:00in-source"               
    ## [4] "DB#: MoNA000010"                      
    ## [5] "INCHIKEY: NCDNCNXCDXHOMX-XGKFQTDJSA-N"
    ## [6] "INSTRUMENTTYPE: Waters Synapt G2"

Or export the `Spectra` with data in NIST MSP format to a MSP file with
MoNA flavor.

\
`tmpf`` ``<-`` `[`tempfile`](https://rdrr.io/r/base/tempfile.html)`(``)`\
\
[`export`](https://rdrr.io/pkg/Spectra/man/Spectra.html)`(``sp``, backend ``=`` `[`MsBackendMsp`](https://rformassspectrometry.github.io/MsBackendMsp/reference/MsBackendMsp.md)`(``)``, file ``=`` ``tmpf``,`\
`       mapping ``=`` `[`spectraVariableMapping`](https://rdrr.io/pkg/Spectra/man/spectraVariableMapping.html)`(`[`MsBackendMsp`](https://rformassspectrometry.github.io/MsBackendMsp/reference/MsBackendMsp.md)`(``)``, ``"mona"``)``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(`[`readLines`](https://rdrr.io/r/base/readLines.html)`(``tmpf``)``)`

    ## [1] "Name: 3-Hydroxy-3-(2-(2-hydroxyphenyl)-2-oxoethyl)-1,3-dihydro-2H-indol-2-one"
    ## [2] "Spectrum_type: MSNA"                                                          
    ## [3] "Ion_mode: Positive"                                                           
    ## [4] "PrecursorMZ: 284.0917"                                                        
    ## [5] "Precursor_type: [M+H]+"                                                       
    ## [6] "Instrument_type: IT/ion trap"

Thus, this could also be used to convert between MSP files with
different flavors.

## Session information

\
[`sessionInfo`](https://rdrr.io/r/utils/sessionInfo.html)`(``)`

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] MsBackendMsp_1.17.1  Spectra_1.23.5       BiocParallel_1.47.0 
    ## [4] S4Vectors_0.51.10    BiocGenerics_0.59.12 generics_0.1.4      
    ## [7] BiocStyle_2.41.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] jsonlite_2.0.0         compiler_4.6.1         BiocManager_1.30.27   
    ##  [4] parallel_4.6.1         cluster_2.1.8.3        jquerylib_0.1.4       
    ##  [7] systemfonts_1.3.2      IRanges_2.47.5         textshaping_1.0.5     
    ## [10] yaml_2.3.12            fastmap_1.2.0          R6_2.6.1              
    ## [13] ProtGenerics_1.45.0    knitr_1.52             htmlwidgets_1.6.4     
    ## [16] MASS_7.3-66            bookdown_0.48          desc_1.4.3            
    ## [19] bslib_0.12.0           rlang_1.3.0            cachem_1.1.0          
    ## [22] xfun_0.62              fs_2.1.0               MsCoreUtils_1.25.4    
    ## [25] sass_0.4.10            otel_0.2.0             cli_3.6.6             
    ## [28] pkgdown_2.2.1.9000     digest_0.6.39          MetaboCoreUtils_1.21.1
    ## [31] lifecycle_1.0.5        clue_0.3-68            data.table_1.18.6.1   
    ## [34] evaluate_1.0.5         codetools_0.2-20       ragg_1.5.2            
    ## [37] rmarkdown_2.32         tools_4.6.1            htmltools_0.5.9

## References

Rainer, Johannes, Andrea Vicini, Liesa Salzer, et al. 2022. “A Modular
and Expandable Ecosystem for Metabolomics Data Annotation in R.”
*Metabolites* 12 (2): 173. <https://doi.org/10.3390/metabo12020173>.

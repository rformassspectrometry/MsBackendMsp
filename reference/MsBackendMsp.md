# MS data backend for msp files

The `MsBackendMsp` class supports import of MS/MS spectra data from
files in NIST MSP file format. `MsBackendMsp` extends the
[`Spectra::MsBackendDataFrame()`](https://rdrr.io/pkg/Spectra/man/MsBackend.html)
backend directly and supports thus the
[`Spectra::applyProcessing()`](https://rdrr.io/pkg/ProtGenerics/man/processingQueue.html)
function to make data manipulations persistent.

New objects are created with the `MsBackendMsp()` function. The
`backendInitialize()` method has to be subsequently called to initialize
the object and import MS/MS data from (one or more) msp files.

The `MsBackendMsp` backend provides an `export()` method that allows to
export the data from the `Spectra` object (parameter `x`) to a file in
MSP format.

Parameters to this function are:

- `x`: the `Spectra` object that should be exported.

- `file`: `character(1)` with the desired file name.

- `mapping`: named `character` providing the mapping between spectra
  variables and MSP data fields. Defaults to
  `mapping = spectraVariableMapping(MsBackendMsp())`.

- `allVariables`: `logical(1)` whether all spectra variables in `x`
  should be exported or only those defined with `mapping`.

- `exportName`: `logical(1)` whether a `NAME` field should always be
  exported even if not provided in `x`.

See the package vignette for details and examples.

The `spectraVariableMapping()` function allows to provide the mapping
between spectra variable names (i.e. the names that will be used for the
spectra variables in the
[`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
object) and the data field names of the MSP file. Parameter `format`
allows to select pre-defined mappings. Currently supported mapping
flavors are:

- `format = "msp"`: default MSP field names. Should work with standard
  NIST MSP files or MSP files exported from MS-DIAL.

- `format = "mona"`: MSP file format from MoNA including LipidBlast.

## Usage

``` r
# S4 method for class 'MsBackendMsp'
backendInitialize(
  object,
  file,
  mapping = spectraVariableMapping(object),
  ...,
  BPPARAM = SerialParam()
)

MsBackendMsp()

# S4 method for class 'MsBackendMsp'
spectraVariableMapping(object, format = c("msp", "mona"))

# S4 method for class 'MsBackendMsp'
export(
  object,
  x,
  file = tempfile(),
  mapping = spectraVariableMapping(MsBackendMsp()),
  allVariables = TRUE,
  exportName = TRUE,
  ...
)
```

## Arguments

- object:

  Instance of `MsBackendMsp` class.

- file:

  `character` with the (full) file name(s) of the msp file(s) from which
  MS/MS data should be imported or exported.

- mapping:

  named `character` vector to rename MSP fields to spectra variables.
  This allows to correctly import also custom fields or data from files
  with different MSP *flavors*.

- ...:

  Currently ignored.

- BPPARAM:

  Parameter object defining the parallel processing setup to import data
  in parallel. Defaults to `BPPARAM = SerialParam()`. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more information. Parallel processing would make most sense for
  import from a large set of individual MSP files, but could also
  improve performance for import from a (very large) single MSP file.

- format:

  For `spectraVariableMapping()`: `character(1)` specifying for which
  MSP *flavour* the mapping should be returned. Currently supported are:
  `format = "msp"` (generic MSP format, for example for MS-DIAL MSP
  files) and `format = "mona"` (MSP files in MoNA flavour).

- x:

  For `export()`: a
  [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  object that should be exported to the specified MSP file.

- allVariables:

  `logical(1)` whether all spectra variables in `x` should be exported
  or only those defined with `mapping`.

- exportName:

  `logical(1)` whether a `NAME` field should always be exported even if
  not provided in `x`.

## Value

`MsBackendMsp()` and `backendInitialize()` return an instance of a
`MsBackendMsp` class. `spectraVariableMapping()` a named `character`
vector with the mapping between spectra variables and MSP data fields.

## Note

Format requirements/assumptions of MSP files:

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

## Author

Steffen Neumann, Michael Witting, Laurent Gatto and Johannes Rainer

## Examples

``` r

## Import spectra from a MSP file from LipidBlast
f <- system.file("extdata", "small-export-LipidBlast.msp",
    package = "MsBackendMsp")
be <- backendInitialize(MsBackendMsp(), f)
be
#> MsBackendMsp with 5 spectra
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         2        NA        NA
#> 2         2        NA        NA
#> 3         2        NA        NA
#> 4         2        NA        NA
#> 5         2        NA        NA
#>  ... 32 more variables/columns.

be$msLevel
#> [1] 2 2 2 2 2
be$intensity
#> NumericList of length 5
#> [[1]] 80.08008 100
#> [[2]] 80.08008 100
#> [[3]] 80.08008 100
#> [[4]] 30.03003 50.05005 50.05005 50.05005 ... 70.07007 50.05005 50.05005 100
#> [[5]] 30.03003 50.05005 50.05005 50.05005 ... 70.07007 50.05005 50.05005 100
be$mz
#> NumericList of length 5
#> [[1]] 85.02841 232.1543
#> [[2]] 85.02841 246.17
#> [[3]] 85.02841 260.1856
#> [[4]] 152.9958 227.2016 283.2643 327.2329 ... 691.4344 703.5283 1451.996
#> [[5]] 152.9958 253.2173 255.233 281.2486 ... 673.4814 721.4814 1451.996

## precursor m/z are however all missing
be$precursorMz
#> [1] NA NA NA NA NA

## Default spectra variable mapping
spectraVariableMapping(MsBackendMsp())
#>            name       accession         formula        inchikey          adduct 
#>          "NAME"           "DB#"       "FORMULA"      "INCHIKEY" "PRECURSORTYPE" 
#>       exactmass           rtime     precursorMz          adduct          smiles 
#>     "EXACTMASS" "RETENTIONTIME"   "PRECURSORMZ" "PRECURSORTYPE"        "SMILES" 
#>           inchi        polarity      instrument 
#>         "INCHI"       "IONMODE"    "INSTRUMENT" 

## In fact, to read MSP files in "LipidBlast flavour" (same as MoNA) we
## should use a different spectra variable mapping
spectraVariableMapping(MsBackendMsp(), "mona")
#>                  name               synonym             accession 
#>                "Name"               "Synon"                 "DB#" 
#>              inchikey                adduct           precursorMz 
#>            "InChIKey"      "Precursor_type"         "PrecursorMZ" 
#>              polarity               formula             exactmass 
#>            "Ion_mode"             "Formula"           "ExactMass" 
#> collision_energy_text               msLevel            instrument 
#>    "Collision_energy"       "Spectrum_type"          "Instrument" 
#>       instrument_type 
#>     "Instrument_type" 

## Importing the data with this will correctly retrieve data
be <- backendInitialize(MsBackendMsp(), f,
    mapping = spectraVariableMapping(MsBackendMsp(), "mona"))
be$precursorMz
#> [1]  232.1543  246.1700  260.1856 1451.9962 1451.9962

## Other fields are also correctly mapped, but might need to be converted
## to e.g. numeric, such as "exactmass"
be$exactmass
#> [1] "232.1543346040907"  "246.16998466809073" "260.18563473209076"
#> [4] "1453.003526472"     "1453.003526472"    

be$exactmass <- as.numeric(be$exactmass)

be$adduct
#> [1] "[M]+"   "[M]+"   "[M]+"   "[M-H]-" "[M-H]-"
be$formula
#> [1] "[C11H22NO4]+" "[C12H24NO4]+" "[C13H26NO4]+" "C81H146O17P2" "C81H146O17P2"

## Exporting Spectra objects in MSP format.

sps <- Spectra(be)
export(MsBackendMsp(), sps, file = stdout())
#> NAME: ACar 4:0
#> msLevel: MS2
#> IONMODE: Positive
#> PRECURSORMZ: 232.15433
#> Comments: "SMILES=CCCC(=O)OC(CC(O)=O)C[N+](C)(C)C" "compound class=ACar" "computed SMILES=O=C(O)CC(OC(=O)CCC)C[N+](C)(C)C" "computed InChI=InChI=1S/C11H21NO4/c1-5-6-11(15)16-9(7-10(13)14)8-12(2,3)4/h9H,5-8H2,1-4H3/p+1" "retention time=0.51" "collision energy spread=15 V" "author=Tobias Kind, Hiroshi Tsugawa" "computed mass accuracy=2.3431646471704717" "computed mass error=5.439758187435473E-4" "SPLASH=splash10-001r-7090000000-aa12589a2481560ea0d5" "submitter=Tobias Kind (University of California, Davis)" "MoNA Rating=3.6363636363636367"
#> MW: 232
#> Num.Peaks: 2
#> DB#: LipidBlast000001
#> PRECURSORTYPE: [M]+
#> collision_energy_text: 45 V
#> EXACTMASS: 232.154334604091
#> FORMULA: [C11H22NO4]+
#> INCHIKEY: QWYFHHGCZUCMBN-UHFFFAOYSA-O
#> INSTRUMENT: SCIEX 5600
#> instrument_type: in-silico QTOF
#> synonym: [M]+
#> synonym: $:00in-source
#> Num Peaks: 2
#> 85.02841 80.08008
#> 232.1543 100
#> 
#> NAME: ACar 5:0
#> msLevel: MS2
#> IONMODE: Positive
#> PRECURSORMZ: 246.16998
#> Comments: "SMILES=CCCCC(=O)OC(CC(O)=O)C[N+](C)(C)C" "compound class=ACar" "computed SMILES=O=C(O)CC(OC(=O)CCCC)C[N+](C)(C)C" "computed InChI=InChI=1S/C12H23NO4/c1-5-6-7-12(16)17-10(8-11(14)15)9-13(2,3)4/h10H,5-9H2,1-4H3/p+1" "retention time=0.68" "collision energy spread=15 V" "author=Tobias Kind, Hiroshi Tsugawa" "computed mass accuracy=2.209496944908765" "computed mass error=5.439118187382519E-4" "SPLASH=splash10-000b-7090000000-2b596f4df94dfefba50b" "submitter=Tobias Kind (University of California, Davis)" "MoNA Rating=3.6363636363636367"
#> MW: 246
#> Num.Peaks: 2
#> DB#: LipidBlast000002
#> PRECURSORTYPE: [M]+
#> collision_energy_text: 45 V
#> EXACTMASS: 246.169984668091
#> FORMULA: [C12H24NO4]+
#> INCHIKEY: VSNFQQXVMPSASB-UHFFFAOYSA-O
#> INSTRUMENT: SCIEX 5600
#> instrument_type: in-silico QTOF
#> synonym: [M]+
#> synonym: $:00in-source
#> Num Peaks: 2
#> 85.02841 80.08008
#> 246.17 100
#> 
#> NAME: ACar 6:0
#> msLevel: MS2
#> IONMODE: Positive
#> PRECURSORMZ: 260.18563
#> Comments: "SMILES=CCCCCC(=O)OC(CC(O)=O)C[N+](C)(C)C" "compound class=ACar" "computed SMILES=O=C(O)CC(OC(=O)CCCCC)C[N+](C)(C)C" "computed InChI=InChI=1S/C13H25NO4/c1-5-6-7-8-13(17)18-11(9-12(15)16)10-14(2,3)4/h11H,5-10H2,1-4H3/p+1" "retention time=0.86" "collision energy spread=15 V" "author=Tobias Kind, Hiroshi Tsugawa" "computed mass accuracy=2.0902300356715053" "computed mass error=5.438478186761131E-4" "SPLASH=splash10-03dr-7090000000-a9ec485bc43949b4b278" "submitter=Tobias Kind (University of California, Davis)" "MoNA Rating=3.6363636363636367"
#> MW: 260
#> Num.Peaks: 2
#> DB#: LipidBlast000003
#> PRECURSORTYPE: [M]+
#> collision_energy_text: 45 V
#> EXACTMASS: 260.185634732091
#> FORMULA: [C13H26NO4]+
#> INCHIKEY: VVPRQWTYSNDTEA-UHFFFAOYSA-O
#> INSTRUMENT: SCIEX 5600
#> instrument_type: in-silico QTOF
#> synonym: [M]+
#> synonym: $:00in-source
#> Num Peaks: 2
#> 85.02841 80.08008
#> 260.1856 100
#> 
#> NAME: CL 72:6
#> msLevel: MS2
#> IONMODE: Negative
#> PRECURSORMZ: 1451.99625
#> Comments: "SMILES=CCCCCCCCCCCCCCCCCC(=O)OCC(COP(O)(=O)OCC(O)COP(O)(=O)OCC(COC(=O)CC\C=C/C\C=C/C\C=C/C\C=C/C\C=C/C\C=C/CC)OC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCCCC" "compound class=CL" "computed SMILES=O=C(OCC(OC(=O)CCCCCCCCCCCCC)COP(=O)(O)OCC(O)COP(=O)(O)OCC(OC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCCCC)CCC=CCC=CCC=CCC=CCC=CCC=CCC" "computed InChI=InChI=1S/C81H146O17P2/c1-5-9-13-17-21-25-29-32-35-36-37-38-41-43-47-50-54-58-62-66-79(84)91-71-76(97-80(85)67-63-59-55-51-45-28-24-20-16-12-8-4)73-95-99(87,88)93-69-75(82)70-94-100(89,90)96-74-77(98-81(86)68-64-60-56-52-48-44-40-34-31-27-23-19-15-11-7-3)72-92-78(83)65-61-57-53-49-46-42-39-33-30-26-22-18-14-10-6-2/h9,13,21,25,32,35,37-38,43,47,54,58,75-77,82H,5-8,10-12,14-20,22-24,26-31,33-34,36,39-42,44-46,48-53,55-57,59-74H2,1-4H3,(H,87,88)(H,89,90)/b13-9-,25-21-,35-32-,38-37-,47-43-,58-54-" "retention time=11.91" "collision energy spread=15 V" "author=Tobias Kind, Hiroshi Tsugawa" "computed mass accuracy=3.2506971738511424E-4" "computed mass error=-4.720000106317457E-7" "SPLASH=splash10-0w4i-0134901100-d27b9577a1230ec65177" "submitter=Tobias Kind (University of California, Davis)" "MoNA Rating=4.5"
#> MW: 1453
#> Num.Peaks: 10
#> DB#: LipidBlast395475
#> PRECURSORTYPE: [M-H]-
#> collision_energy_text: 45 V
#> EXACTMASS: 1453.003526472
#> FORMULA: C81H146O17P2
#> INCHIKEY: KWABUIIFXPUAEH-OXNFOSRESA-N
#> INSTRUMENT: SCIEX 5600
#> instrument_type: in-silico QTOF
#> synonym: CL 14:0-22:6-18:0-18:0
#> synonym: $:00in-source
#> Num Peaks: 10
#> 152.9958 30.03003
#> 227.2016 50.05005
#> 283.2643 50.05005
#> 327.2329 50.05005
#> 363.1942 70.07007
#> 419.2568 70.07007
#> 463.2255 70.07007
#> 691.4344 50.05005
#> 703.5283 50.05005
#> 1451.996 100
#> 
#> NAME: CL 72:6
#> msLevel: MS2
#> IONMODE: Negative
#> PRECURSORMZ: 1451.99625
#> Comments: "SMILES=CCCCCCCCCCCCCCCC(=O)OC(COC(=O)CCCCCCC\C=C/CCCCCCCC)COP(O)(=O)OCC(O)COP(O)(=O)OCC(COC(=O)CCCCCCCC\C=C/C\C=C/C\C=C/C\C=C/CC)OC(=O)CCCCCCC\C=C/CCCCCC" "compound class=CL" "computed SMILES=O=C(OCC(OC(=O)CCCCCCCCCCCCCCC)COP(=O)(O)OCC(O)COP(=O)(O)OCC(OC(=O)CCCCCCCC=CCCCCCC)COC(=O)CCCCCCCCC=CCC=CCC=CCC=CCC)CCCCCCCC=CCCCCCCCC" "computed InChI=InChI=1S/C81H146O17P2/c1-5-9-13-17-21-25-29-33-35-36-37-38-40-44-46-50-54-58-62-66-79(84)92-72-77(98-81(86)68-64-60-56-52-48-42-32-28-24-20-16-12-8-4)74-96-100(89,90)94-70-75(82)69-93-99(87,88)95-73-76(97-80(85)67-63-59-55-51-47-41-31-27-23-19-15-11-7-3)71-91-78(83)65-61-57-53-49-45-43-39-34-30-26-22-18-14-10-6-2/h9,13,21,25,28,32-35,37-39,75-77,82H,5-8,10-12,14-20,22-24,26-27,29-31,36,40-74H2,1-4H3,(H,87,88)(H,89,90)/b13-9-,25-21-,32-28-,35-33-,38-37-,39-34-" "retention time=11.65" "collision energy spread=15 V" "author=Tobias Kind, Hiroshi Tsugawa" "computed mass accuracy=3.2506971738511424E-4" "computed mass error=-4.720000106317457E-7" "SPLASH=splash10-0uyi-0157901100-d3b7dfa9dd7050026484" "submitter=Tobias Kind (University of California, Davis)" "MoNA Rating=4.5"
#> MW: 1453
#> Num.Peaks: 12
#> DB#: LipidBlast398514
#> PRECURSORTYPE: [M-H]-
#> collision_energy_text: 45 V
#> EXACTMASS: 1453.003526472
#> FORMULA: C81H146O17P2
#> INCHIKEY: WINXBODRGKNJDS-DUULZMMUSA-N
#> INSTRUMENT: SCIEX 5600
#> instrument_type: in-silico QTOF
#> synonym: CL 16:0-18:1-16:1-22:4
#> synonym: $:00in-source
#> Num Peaks: 12
#> 152.9958 30.03003
#> 253.2173 50.05005
#> 255.233 50.05005
#> 281.2486 50.05005
#> 331.2643 50.05005
#> 389.2098 70.07007
#> 391.2255 70.07007
#> 417.2411 70.07007
#> 467.2568 70.07007
#> 673.4814 50.05005
#> 721.4814 50.05005
#> 1451.996 100
#> 
```

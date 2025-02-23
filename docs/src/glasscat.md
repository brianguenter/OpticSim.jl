# AGFFileReader

This submodule is used to download, parse, install and manage AGF glass specifications for use in OpticSim.

The central configuration file for AGFFileReader is located at `src/AGFFileReader/data/sources.txt`, which ships with the
following default entries.

```
NIKON a49714470fa875ad4fd8d11cbc0edbf1adfe877f42926d7612c1649dd9441e75 https://www.nikon.com/business/components/assets/pdf/nikon_zemax_data.zip
OHARA b8a0bea2a59ec643998040255925c7e051ea5352ff239f1aca50268ac27bb7f7 https://oharacorp.com/wp-content/uploads/2024/02/OHARA_240131_CATALOG.zip
HOYA b02c203e5a5b7a8918cc786badf0a9db1fe2572372c1c163dc306b0a8a908854 http://www.hoya-opticalworld.com/common/agf/HOYA20210105.agf
SCHOTT e9aabbb8ebff116ba0c106a71afd86e72f2a397ac9bc447469129e325e795f4e https://www.schott.com/d/advanced_optics/6959f9a4-0e4f-4ef2-a302-2347468a82f5/1.31/schott-optical-glass-overview-zemax-format.zip
SUMITA c1093e42a1d08acbe30698aba730161e3b43f8f0d50533f65de8b6b11100fdc8 https://www.sumita-opt.co.jp/en/wp/wp-content/themes/sumita-opt-en/download-files.php files%5B%5D=new_sumita.agf&category=data
```

Each line corresponds to one AGF source, which is described by 2 to 4 space-delimited columns. The first column provides
the installed module name for the catalog, e.g. `AGFFileReader.NIKON`. The second column is the expected SHA256 checksum for
the AGF file.

The final two columns are optional, specifying download instructions for acquiring the zipped AGF files
automatically from the web. The fourth column allows us to use POST requests to acquire files from interactive sites.

When `] build OpticSim` is run, the sources are verified and parsed into corresponding Julia files. These are then
included in OpticSim via `AGFGlassCat.jl`. These steps are run automatically when the package is first installed using
`] add OpticSim`, creating a sufficient working environment for our examples and tests.

## Adding glass catalogs
`sources.txt` can be edited freely to add more glass catalogs. However, this is a somewhat tedious process, so we have a
convenience function for adding a locally downloaded AGF file to the source list.

```@docs
AGFFileReader.add_agf
```

## Using installed glasses

Glass types are accessed like so: `AGFFileReader.CATALOG_NAME.GLASS_NAME`, e.g.

```julia
AGFFileReader.SUMITA.LAK7
AGFFileReader.SCHOTT.PK3
```

All glasses and catalogs are exported in their respective modules, so it is possible to invoke `using` calls for convenience, e.g.

```julia
using OpticSim
AGFFileReader.SUMITA.LAK7
using AGFFileReader
SCHOTT.PK3
using OpticsSim.AGFFileReader.SCHOTT
N_BK7
```

Autocompletion can be used to see available catalogs and glasses. All catalog glasses are of type [`AGFFileReader.Glass`](@ref).
Note that special characters in glass/catalog names are replaced with `_`.
There is a special type and constant value for air: [`AGFFileReader.Air`](@ref).

[Unitful.jl](https://github.com/PainterQubits/Unitful.jl) is used to manage units, meaning any valid unit can be used for all arguments, e.g., wavelength can be passed in as μm or nm (or cm, mm, m, etc.).
Non-unitful options are also available, in which case units are assumed to be μm, °C and Atm for length, temperature and pressure respectively.

`TEMP_REF` and `PRESSURE_REF` are constants:

```julia
const TEMP_REF = 20.0 # °C
const PRESSURE_REF = 1.0 # Atm
```

## Types

```@docs
AGFFileReader.AbstractGlass
AGFFileReader.Glass
AGFFileReader.Air
AGFFileReader.Glass
```

## Functions

```@docs
AGFFileReader.index
AGFFileReader.absairindex
AGFFileReader.absorption
```

---

```@docs
AGFFileReader.glassfromMIL
AGFFileReader.modelglass
```

---

```@docs
AGFFileReader.glasscatalogs
AGFFileReader.glassnames
AGFFileReader.info
AGFFileReader.findglass
AGFFileReader.isair
```

---

```@docs
AGFFileReader.glassname
AGFFileReader.glassid
AGFFileReader.glassforid
```

---

```@docs
AGFFileReader.polyfit_indices
AGFFileReader.plot_indices
AGFFileReader.drawglassmap
```

---

```@docs
AGFFileReader.verify_sources!
AGFFileReader.verify_source
AGFFileReader.download_source
```

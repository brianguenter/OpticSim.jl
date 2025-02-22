# OpticSim.jl

OpticSim.jl is a [Julia](https://julialang.org/) package for geometric optics (ray tracing) simulation and optimization of complex optical systems.

It is designed to allow optical engineers to create optical systems procedurally and then to simulate them. Unlike Zemax, Code V, or other interactive optical design systems OpticSim.jl has limited support for interactivity.

A large variety of surface types are supported, and these can be composed into complex 3D objects through the use of constructive solid geometry (CSG). A complete catalog of optical materials is provided through the complementary GlassCat submodule.

This software provides extensive control over the modelling, simulation, visualization and optimization of optical systems. It is especially suited for designs that have a procedural architecture.

## Installation

Install the Julia programming language from the [official download page](https://julialang.org/downloads/).

OpticSim.jl is optimized for use with Julia 1.7.x; using other versions will probably not work.

The system will automatically download glass catalog (.agf) files from some manufacturers when the package is built for the first time. These files are in an industry standard format and can be downloaded from many optical glass manufacturers.

Here are links to several publicly available glass files:

* [NIKON](https://www.nikon.com/business/components/assets/pdf/nikon_zemax_data.zip) (automatically downloaded)
* [NHG](http://hbnhg.com/down/data/nhgagp.zip) (you have to manually download)
* [OHARA](https://oharacorp.com/wp-content/uploads/2024/02/OHARA_240131_CATALOG.zip) (automatically downloaded)
* [HOYA](https://hoyaoptics.com/wp-content/uploads/2019/10/HOYA20170401.zip) (automatically downloaded)
* [SUMITA](https://www.sumita-opt.co.jp/en/download/) (automatically downloaded)
* [SCHOTT](https://www.schott.com/advanced_optics/english/download/index.html) (automatically downloaded)

OpticSim.jl will generate a glass database from the available files in `deps/downloads/glasscat/` and store it in the file `AGFClassCat.jl`. See [GlassCat](@ref) for a detailed description, including instructions on how to add more catalogs.

Run this example to check that everything installed properly:

```@example
using OpticSim
Vis.draw(SphericalLens(Examples.Examples_N_BK7, 0.0, 10.0, 10.0, 5.0, 5.0))
Vis.save("assets/test_install.png") # hide
nothing # hide
```

![install test image](assets/test_install.png)



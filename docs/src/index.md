# OpticSim.jl

OpticSim.jl is a [Julia](https://julialang.org/) package for geometric optics (ray tracing) simulation and optimization of complex optical systems.

It is designed to allow optical engineers to create optical systems procedurally and then to simulate them. Unlike Zemax, Code V, or other interactive optical design systems OpticSim.jl has limited support for interactivity.

A large variety of surface types are supported, and these can be composed into complex 3D objects through the use of constructive solid geometry (CSG). A complete catalog of optical materials is provided through the complementary AGFFileReader package.

## Installation
Run this example to check that everything installed properly:

```@example
using OpticSim, AGFFileReader
SphericalLens(AGFFileReader.Examples_N_BK7, 0.0, 10.0, 10.0, 5.0, 5.0)
```
`AGFFileReader` includes these 4 glass types to get you started
```julia
Examples_N_BK7
Examples_N_SF14
Examples_N_SF2
Examples_N_SK16
```

For serious work you will need more glass types. The `AGFFileReader` package has a function `AGFFileReader.initialize_AGFFileReader()` to make a local database of glass types generated from publicly available glass catalogs on the web. You must call this function **before** you execute your code:
```julia
module YourModule
using OpticSim
import AGFFileReader

AGFFileReader.initialize_AGFFileReader()

#your code goes here

end #YourModule
```

The first call to `AGFFileReader.initialize_AGFFileReader()` will search for publicly available glass files on the web and install them into a database locally on your machine. This may take a while. Subsequent calls will detect the local database and skip the download step. See the documentation of `AGFFileReader` for more information.

In previous versions of `OpticSim` the glass catalog download, visualization, and repeating structures code was included in the `OpticSim` package. This caused excessive load times. 

These have been broken out into separate packages: `AGFFileReader`, `OpticSimVisualization`, and `OpticSimRepeatingStructures`. The last two packages are not yet working.





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
If you want to use downloaded glass files then you will have to call the function `AGFFileReader.initialize_AGFFileReader()` before you 
Now download glass files. Open a REPL in the `OpticSim.jl` project directory. At the REPL type

```julia
julia> using AGFFileReader

julia> initialize_AGFFileReader()
12-element Vector{Module}:
 AGFFileReader.NIKON
 AGFFileReader.HOYA
 AGFFileReader.ARTON
 AGFFileReader.BIREFRINGENT
 AGFFileReader.CDGM
 AGFFileReader.HERAEUS
 AGFFileReader.HIKARI
 AGFFileReader.ISUZU
 AGFFileReader.LZOS
 AGFFileReader.MISC
 AGFFileReader.UMICORE
 AGFFileReader.ZEON
```
In this case 12 glass files were downloaded. Print the properties of a random glass to test that the files were downloaded correctly:
```julia
julia> info(NIKON.BAF3)
Dispersion formula:                                Schott (1)
Dispersion formula coefficients:
     a₀:                                           2.438682      
     a₁:                                           0.001609961   
     a₂:                                           0.02688829    
     a₃:                                           -0.002446579  
     a₄:                                           0.0003965661  
     a₅:                                           -1.973046e-5  
Valid wavelengths:                                 0.4μm to 0.7μm
Reference temperature:                             20.0°C        
TCE (÷1e-6):                                       9.1
Ignore thermal expansion:                          false
Density (p):                                       3.28g/m³      
ΔPgF:                                              0.0009        
RI at sodium D-Line (587nm):                       1.58267       
Abbe Number:                                       46.476929     
Cost relative to N_BK7:                            ?
Status:                                            Obsolete (2)  
Melt frequency:                                    0
Exclude substitution:                              false
Transmission data:
     Wavelength   Transmission      Thickness
         0.29μm            0.0         10.0mm
          0.3μm            0.0         10.0mm
         0.31μm            0.0         10.0mm
         0.32μm            0.1         10.0mm
         0.33μm           0.51         10.0mm
         0.34μm           0.81         10.0mm
         0.35μm          0.924         10.0mm
         0.36μm          0.962         10.0mm
         0.37μm           0.98         10.0mm
         0.38μm          0.983         10.0mm
         0.39μm          0.991         10.0mm
          0.4μm          0.994         10.0mm
         0.42μm          0.996         10.0mm
         0.44μm          0.997         10.0mm
         0.46μm          0.999         10.0mm
         0.48μm          0.999         10.0mm
          0.5μm          0.999         10.0mm
         0.55μm          0.999         10.0mm
          0.6μm          0.999         10.0mm
         0.65μm          0.999         10.0mm
          0.7μm          0.999         10.0mm
```

![install test image](assets/test_install.png)



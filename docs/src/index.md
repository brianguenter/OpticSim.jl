# OpticSim.jl

OpticSim.jl is a [Julia](https://julialang.org/) package for geometric optics (ray tracing) simulation and optimization of complex optical systems.

It is designed to allow optical engineers to create optical systems procedurally and then to simulate them. Unlike Zemax, Code V, or other interactive optical design systems OpticSim.jl has limited support for interactivity.

A large variety of surface types are supported, and these can be composed into complex 3D objects through the use of constructive solid geometry (CSG). A complete catalog of optical materials is provided through the complementary AGFFileReader package.

## Installation

Run this example to check that everything installed properly:

```@example
using OpticSim
Vis.draw(SphericalLens(Examples.Examples_N_BK7, 0.0, 10.0, 10.0, 5.0, 5.0))
Vis.save("assets/test_install.png") # hide
nothing # hide
```

![install test image](assets/test_install.png)



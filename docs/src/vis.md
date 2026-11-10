# Visualization

OpticSim draws systems, surfaces, rays and more in 3D using [Makie](https://docs.makie.org/stable/), and makes 2D plots such as spot diagrams and surface sag maps.

The drawing functions live in the `Vis` module. They are implemented in a package extension that loads automatically when you load a Makie backend, so OpticSim itself stays fast to load:

```julia
using OpticSim
using GLMakie       # interactive window, with buttons to switch between 3D and 2D views and take a screenshot
# or
using CairoMakie    # static, high quality images; works without a GPU
```

Makie and OpticSim export some of the same names, such as `Sphere` and `translation`. To use OpticSim's versions without a prefix, write `import GLMakie` rather than `using GLMakie`.

In [Pluto](https://plutojl.org/) notebooks call `Vis.set_current_mode(:pluto)` so the drawing functions return the figure instead of opening a window.

```@setup vis
using OpticSim, OpticSim.Geometry, OpticSim.Emitters, AGFFileReader
import CairoMakie
using DataFrames: DataFrame
```

## Tracing rays through a system

The examples on this page use a Cooke triplet:

```@example vis
using OpticSim, AGFFileReader
using DataFrames: DataFrame

cooketriplet(::Type{T}=Float64, detpix::Int=1000) where {T<:Real} = AxisymmetricOpticalSystem{T}(
    DataFrame(
        SurfaceType=["Object", "Standard", "Standard", "Standard", "Stop", "Standard", "Standard", "Image"],
        Radius=[Inf, 26.777, 66.604, -35.571, 35.571, 35.571, -26.777, Inf],
        Thickness=[Inf, 4.0, 2.0, 4.0, 2.0, 4.0, 44.748, missing],
        Material=[AGFFileReader.Air, AGFFileReader.Examples_N_SK16, AGFFileReader.Air, AGFFileReader.Examples_N_SF2, AGFFileReader.Air, AGFFileReader.Examples_N_SK16, AGFFileReader.Air, missing],
        SemiDiameter=[Inf, 8.580, 7.513, 7.054, 6.033, 7.003, 7.506, 15.0]
    ),
    detpix,
    detpix
)
nothing # hide
```

[`Vis.drawtracerays`](@ref OpticSim.Vis.drawtracerays) draws the system and the paths of rays from a ray generator, colored by wavelength:

```@example vis
Vis.drawtracerays(cooketriplet(), trackallrays=true, test=true, numdivisions=100)
```

`view=:side` draws it like a lens layout: from the side, optical axis horizontal, light travelling left to right:

```@example vis
Vis.drawtracerays(cooketriplet(), trackallrays=true, test=true, numdivisions=100, view=:side)
```

[`Vis.drawcomparison`](@ref OpticSim.Vis.drawcomparison) draws several systems one above the other with the same rays, for example two designs, and [`Vis.scene`](@ref OpticSim.Vis.scene) makes an empty figure to build your own drawing in.

The image formed on the detector:

```@example vis
image = Vis.drawtraceimage(cooketriplet(Float64, 400), test=true)
Vis.showimage(image)
```

## Basic drawing

`Vis.draw` draws an object in a new figure and `Vis.draw!` adds more objects to it. This is all you need to build up a visualization piece by piece:

```@example vis
obj = (Sphere(0.5) ∩ Plane(0.0, 1.0, 0.0, 0.0, 0.1, 0.0))()
ray1 = Ray([0.0, -0.1, 1.0], [0.0, 0.0, -1.0])
ray2 = Ray([0.8, 0.0, 0.0], [-1.0, 0.0, 0.0])
Vis.draw(obj)
Vis.draw!(ray1, rayscale=0.2)
Vis.draw!(ray2, rayscale=0.2, color=:blue)
Vis.draw!(surfaceintersection(obj, ray1), color=:red)
Vis.draw!(surfaceintersection(obj, ray2), color=:green)
```

```@docs
OpticSim.Vis
OpticSim.Vis.draw
OpticSim.Vis.draw!
OpticSim.Vis.scene
OpticSim.Vis.update_camera_orientation!
OpticSim.Vis.save
OpticSim.Vis.showimage
```

## Analysis plots

A spot diagram of where rays from a collimated source hit the detector:

```@example vis
source = Sources.Source(
    transform=translation(0.0, 0.0, 10.0),
    origins=Origins.Hexapolar(8, 8.0, 8.0),
    directions=Directions.Constant(0.0, 0.0, -1.0),
)
Vis.spotdiag(cooketriplet(), source)
```

The optical path difference across the entrance pupil:

```@example vis
Vis.plotOPD(cooketriplet(), samples=20)
```

The sag of a Zernike surface:

```@example vis
zernike = ZernikeSurface(1.5, radius=5.0, conic=1.0, zcoeff=[(7, 0.2), (9, 0.2)], aspherics=[(10, -0.01)])
Vis.surfacesag(AcceleratedParametricSurface(zernike), (256, 256), (1.55, 1.55))
```

```@docs
OpticSim.Vis.drawtracerays
OpticSim.Vis.drawtraceimage
OpticSim.Vis.spotdiag
OpticSim.Vis.spotdiaggrid
OpticSim.Vis.plotOPD
OpticSim.Vis.surfacesag
OpticSim.Vis.drawcurves
OpticSim.Vis.drawcomparison
OpticSim.Vis.sideview!
OpticSim.Vis.eyebox_eval_eye
OpticSim.Vis.eyebox_eval_planar
```

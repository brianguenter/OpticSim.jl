# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    Vis

Drawing and plotting functions for OpticSim objects.

The functions are declared here but implemented in the `OpticSimMakieExt` package extension, which loads
automatically once a Makie backend is loaded, e.g.

```julia
using OpticSim, GLMakie   # interactive window
# or
using OpticSim, CairoMakie # static images, works without a GPU
Vis.draw(Examples.cooketriplet())
```
"""
module Vis

export draw, draw!, drawtracerays, drawtracerays!, drawtraceimage, spotdiag, spotdiag!, spotdiaggrid, plotOPD, plotOPD!, surfacesag, drawcurves

# Settings shared between OpticSim and the extension.
# modes:  nothing/:default -> open a window and add camera buttons
#         :pluto           -> return the figure instead of displaying it (Pluto notebooks)
#         :docs            -> return the figure instead of displaying it (Documenter)
const current_mode = Ref{Union{Nothing,Symbol}}(nothing)
get_current_mode() = current_mode[]
set_current_mode(mode) = (current_mode[] = mode)

const _visfunctions = (
    :scene, :draw, :draw!, :display, :save, :showimage, :update_camera_orientation!,
    :drawtracerays, :drawtracerays!, :drawtraceimage,
    :spotdiag, :spotdiag!, :spotdiaggrid, :plotOPD, :plotOPD!, :surfacesag,
    :eyebox_eval_eye, :eyebox_eval_planar, :drawcurves, :λtoRGB,
)

for f in _visfunctions
    @eval function $f end
end

# When a Vis function is called before a Makie backend is loaded, explain how to fix it.
function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, exc, _, _
        f = exc.f
        if any(s -> getfield(Vis, s) === f, _visfunctions) && Base.get_extension(parentmodule(Vis), :OpticSimMakieExt) === nothing
            print(io, "\nOpticSim.Vis drawing functions need a Makie backend. Run `using GLMakie` (interactive) or `using CairoMakie` (static images) first.")
        end
    end
end

"""
    scene(resolution = (1000, 1000))

Create a new Makie `Figure` with a 3D `LScene` and make it the current figure. Returns `(figure, lscene)`.
"""
scene

"""
    draw(ob; resolution = (1000, 1000), kwargs...)

Draw an object in a new figure and display it. `kwargs` depend on the object type.
"""
draw

"""
    draw!([lscene], ob; kwargs...)

Draw an object in an existing scene (the current one if `lscene` is omitted). `kwargs` depend on the object type.
"""
draw!

"""
    save(path::String)

Save the current figure to an image file.
"""
save

"""
    drawtracerays(system; raygenerator, test = false, trackallrays = false, colorbysourcenum = false, colorbynhits = false, rayfilter = onlydetectorrays, verbose = false, resolution = (1000, 1000), kwargs...)

Draw an optical system along with the rays traced through it.
If `trackallrays` is true the whole ray path is drawn, otherwise only the final segment hitting the detector.
Rays are colored by wavelength unless `colorbysourcenum` or `colorbynhits` is set.
Only ray paths that reach the detector are drawn; set `rayfilter = nothing` to draw all of them.

`drawtracerays!([lscene], system; ...)` adds to an existing scene, with `drawsys` and `drawgen` controlling whether the system and ray generator are drawn.
"""
drawtracerays

"""
    drawtraceimage(system; raygenerator, test = false, verbose = false)

Trace rays from `raygenerator` through `system`, show the detector image and return it.
"""
drawtraceimage

"""
    spotdiag(sys::CSGOpticalSystem, raygenerator::OpticalRayGenerator; color = nothing, label = nothing, size = (500, 500))

Plot a spot diagram: trace all rays from `raygenerator` through `sys` and plot where they hit the detector, relative to the centroid.
`spotdiag!(ax, sys, raygenerator)` adds to an existing 2D `Axis`.
"""
spotdiag

"""
    spotdiaggrid(sys, raygenerators::Vector{<:OpticalRayGenerator}; colorbynum = false, kwargs...)

Plot a grid of spot diagrams, one per ray generator, all sharing the same axis limits.
"""
spotdiaggrid

"""
    plotOPD(sys::AxisymmetricOpticalSystem; label = nothing, color = nothing, collimated = true, axis = 1, samples = 5, wavelength = 0.55, sourcepos = SVector(0.0, 0.0, 10.0))

Plot optical path difference across the entrance pupil. `plotOPD!(ax, sys; ...)` adds to an existing 2D `Axis`.
"""
plotOPD

"""
    surfacesag(object::Union{CSGTree,Surface}, resolution::Tuple{Int,Int}, halfsizes::Tuple{T,T}; offset = 10, position = SVector(0.0, 0.0, 10.0), direction = SVector(0.0, 0.0, -1.0), rotationvec = SVector(0.0, 1.0, 0.0))

Calculate and plot the surface sag of a [`Surface`](@ref) or [`CSGTree`](@ref) as a heatmap.
"""
surfacesag

"""
    drawcurves(curves::Spline...; numpoints = 200)

Draw 2D spline curves along with their control polygons.
"""
drawcurves

end # module Vis
export Vis

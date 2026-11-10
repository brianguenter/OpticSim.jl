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

export draw, draw!, drawtracerays, drawtracerays!, drawtraceimage, spotdiag, spotdiag!, spotdiaggrid, plotOPD, plotOPD!, surfacesag, drawcurves, drawcomparison, sideview!

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
    :eyebox_eval_eye, :eyebox_eval_planar, :drawcurves, :λtoRGB, :sideview!, :drawcomparison,
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
Unrecognized `kwargs` are passed on to the underlying Makie plotting function.

Objects that can be drawn, and their main keyword arguments:

| Object | Keyword arguments |
|:-------|:------------------|
| `Surface` | `numdivisions = 30`, `normals = false`, `normalcolor = :blue`, plus those of `TriangleMesh` |
| `TriangleMesh` | `color = :orange`, `shaded = true`, `wireframe = false`, `linewidth = 3`, `transparency = false`, `normals = false` |
| `CSGTree`, `CSGGenerator` | `numdivisions = 30`, plus those of `TriangleMesh` |
| several `Surface`s or `TriangleMesh`es | `colors = false`: color each one differently |
| `LensAssembly`, `CSGOpticalSystem`, `AxisymmetricOpticalSystem` | each element is colored differently |
| `BoundingBox` | drawn as a wireframe |
| `Ray`, `OpticalRay` | `color = :yellow`, `rayscale = 1.0`; optical rays can use `colorbysourcenum`, `colorbynhits` |
| `LensTrace` or a vector of them | `colorbysourcenum = false`, `colorbynhits = false`, otherwise colored by wavelength |
| `Interval`, `DisjointUnion`, `Intersection` | `normal = false` for intersections |
| a point, a vector of points | `markersize = 20`, `color = :black` |
| `(point1, point2)` line, or a vector of them | `color = :yellow` |
| emitter `Source`, `CompositeSource` | `debug = false`: also draw the local axes, sample positions and rays |
| a mesh file name (`.obj`, `.stl`, `.ply`, `.off`, `.2dm`) | `color = :gray`, `transform`, `scale = 1.0`, `shaded`, `wireframe`. Requires `MeshIO` to be loaded. |
"""
draw!

"""
    showimage(image::AbstractArray)

Show a 2D image, such as a detector image, in grayscale with row 1 at the top. Images with more than 2 dimensions get a slider for each extra dimension.
"""
showimage

"""
    update_camera_orientation!(lscene, ϕ, θ)

Point the camera at the scene from angles `0 ≤ ϕ ≤ 2π` and `-π/2 ≤ θ ≤ π/2`. For example `(0, 0)` looks along the x axis and `(π/2, 0)` along the y axis.
"""
update_camera_orientation!

"""
    eyebox_eval_eye(assembly::LensAssembly, raygen::OpticalRayGenerator, eye_rotation_x, eye_rotation_y, sample_points_x::Int, sample_points_y::Int; pupil_radius = 2.0, resolution::Int = 512, eye_transform = identitytransform())

Show the images formed when tracing `assembly` with a human eye, for an evenly sampled `sample_points_x × sample_points_y` grid of eye rotations in the ranges `-eye_rotation_x:eye_rotation_x` and `-eye_rotation_y:eye_rotation_y`.
`resolution` is the size of the (square) detector image.

`eye_transform` should place the eye correctly relative to the system when at 0 rotation. By default the eye looks along the positive z axis with the vertex of the cornea at the origin.

The result is a 4D image, shown in 2D with sliders to vary the eye rotation in x and y. The whole image should be visible for all rotations in the range.
"""
eyebox_eval_eye

"""
    eyebox_eval_planar(assembly::LensAssembly, raygen::OpticalRayGenerator, eyebox::Rectangle, sample_points_x::Int, sample_points_y::Int, vsize; pupil_radius = 2.0, resolution::Int = 512)

Show the images formed when tracing `assembly` for multiple pupil positions within a planar `eyebox`.
The rays hitting a circle of radius `pupil_radius` around each point of an even `sample_points_x × sample_points_y` grid on the `eyebox` are added to the sub-image for that grid point.

A paraxial lens with focal length 1mm is placed at the eyebox and a `vsize × vsize` mm detector 1mm behind it.
The normal of `eyebox` should point toward the system. Rays which miss the detector are ignored.
The pupil is always fully inside the eyebox, so the extreme sample position in u is `eyebox.halfsizeu - pupil_radius`, for example.

The result is a 4D image, shown in 2D with sliders to vary the pupil position in x and y.
"""
eyebox_eval_planar

"""
    save(path::String)

Save the current figure to an image file.
"""
save

"""
    drawtracerays(system; raygenerator, view = :perspective, test = false, trackallrays = false, colorbysourcenum = false, colorbynhits = false, rayfilter = onlydetectorrays, verbose = false, resolution, kwargs...)

Draw an optical system along with the rays traced through it.
If `trackallrays` is true the whole ray path is drawn, otherwise only the final segment hitting the detector.
Rays are colored by wavelength unless `colorbysourcenum` or `colorbynhits` is set.
Only ray paths that reach the detector are drawn; set `rayfilter = nothing` to draw all of them.

`lensalpha` sets the opacity of the lens elements (default 1). The default `raygenerator` is a 9 × 9 grid of collimated
rays at 0.55μm coming down the optical axis from z = 10: few enough that the lenses stay visible.

`view = :side` shows the system from the side, like a lens layout drawing: the optical axis horizontal, light travelling
left to right, orthographic projection, fitted to the system (see [`sideview!`](@ref OpticSim.Vis.sideview!)).
`resolution` defaults to `(1000, 1000)`, or `(1400, 500)` for the side view.

`drawtracerays!([lscene], system; ...)` adds to an existing scene, with `drawsys` and `drawgen` controlling whether the system and ray generator are drawn.
"""
drawtracerays

"""
    sideview!(lscene)

Show the contents of `lscene` from the side, like a lens layout drawing: looking along the x axis at the y-z plane, so
the optical axis (z) is horizontal and light travelling toward -z goes left to right, with orthographic projection,
zoomed to fit. Returns `lscene`.
"""
sideview!

"""
    drawcomparison(systems; raygenerator, titles = nothing, view = :side, resolution = (1400, 380 * length(systems)), kwargs...)

Draw several optical systems one above the other, tracing the same rays from `raygenerator` through each: for example
a design before and after optimization. `titles` labels each drawing, `view` is `:side` or `:perspective` (see
[`drawtracerays`](@ref OpticSim.Vis.drawtracerays)), and other keyword arguments such as `trackallrays`, `test` and
`colorbysourcenum` are passed on to `drawtracerays!`.

```julia
Vis.drawcomparison([before, after]; raygenerator=rays, titles=["Before", "After"], trackallrays=true, test=true)
Vis.save("comparison.png")
```
"""
drawcomparison

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

Calculate and plot the surface sag of a `Surface` or `CSGTree` as a heatmap.

Rays are shot in a grid of size `resolution` across a rectangle with half sizes `halfsizes`, centered at `position` with normal along `direction` and rotation defined by `rotationvec`.
`offset` is subtracted from the sag measurements so values are relative to the appropriate zero level.
"""
surfacesag

"""
    drawcurves(curves::Spline...; numpoints = 200)

Draw 2D spline curves along with their control polygons.
"""
drawcurves

end # module Vis
export Vis

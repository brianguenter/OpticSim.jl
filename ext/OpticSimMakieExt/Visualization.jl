# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# All drawing functions follow the pattern draw(obj) and draw!(lscene, obj). The first draws the object in a new
# figure and displays it, the second draws into an existing LScene. draw!(obj) draws into the current LScene.

const current_main_scene = Ref{Any}(nothing)   # the current Makie.Figure
const current_3d_scene = Ref{Any}(nothing)     # the current Makie.LScene

isnoninteractive() = Vis.get_current_mode() in (:pluto, :docs)

#############################################################################
## figures, display and saving

Vis.showimage(image::AbstractArray) = Vis.display(imagefigure(image))

function Vis.display(fig=current_main_scene[])
    if isnoninteractive()
        return fig
    end
    Makie.display(fig)
end

function Vis.scene(resolution=(1000, 1000); show_axis::Bool=true)
    @assert resolution[1] > 0 && resolution[2] > 0

    fig = Makie.Figure(size=resolution)
    ax = Makie.LScene(fig[1, 1]; show_axis, scenekw=(; camera=Makie.cam3d_cad!))
    current_main_scene[] = fig
    current_3d_scene[] = ax

    # in these modes skip the buttons since the figure is not interactive
    if isnoninteractive()
        return fig, ax
    end

    buttons_layout = Makie.GridLayout(fig[2, 1]; tellwidth=false)
    buttoncolor = RGBf(0.8, 0.8, 0.8)
    threedbutton = Makie.Button(buttons_layout[1, 1], label="3D", buttoncolor=buttoncolor, height=40, width=80)
    twodxbutton = Makie.Button(buttons_layout[1, 2], label="2D-x", buttoncolor=buttoncolor, height=40, width=80)
    twodybutton = Makie.Button(buttons_layout[1, 3], label="2D-y", buttoncolor=buttoncolor, height=40, width=80)
    savebutton = Makie.Button(buttons_layout[1, 4], label="Screenshot", buttoncolor=buttoncolor, height=40, width=160)

    Makie.on(_ -> Vis.update_camera_orientation!(ax, π / 4, π / 4), threedbutton.clicks)
    Makie.on(_ -> Vis.update_camera_orientation!(ax, π / 2, 0.0), twodybutton.clicks)
    Makie.on(_ -> Vis.update_camera_orientation!(ax, 0.0, 0.0), twodxbutton.clicks)
    Makie.on(_ -> Vis.save("screenshot.png"), savebutton.clicks)

    return fig, ax
end

"""
Set the camera position based on two angles 0 ≤ ϕ ≤ 2π and -π/2 ≤ θ ≤ π/2.
"""
function Vis.update_camera_orientation!(ax::Makie.LScene, ϕ, θ)
    cam3d = Makie.cameracontrols(ax.scene)
    Makie.update_cam!(ax.scene, cam3d, ϕ, θ)
end

function Vis.draw(ob; resolution=(1000, 1000), kwargs...)
    fig, lscene = Vis.scene(resolution)
    Vis.draw!(lscene, ob; kwargs...)
    Vis.display(fig)
    return isnoninteractive() ? fig : nothing
end

function Vis.draw!(ob; kwargs...)
    if current_3d_scene[] === nothing
        Vis.scene()
    end
    fig = current_main_scene[]
    Vis.draw!(current_3d_scene[], ob; kwargs...)
    Vis.display(fig)
    return isnoninteractive() ? fig : nothing
end

function Vis.save(path::String)
    current_main_scene[] === nothing && error("There is no current figure to save.")
    Makie.save(path, current_main_scene[])
end
Vis.save(::Nothing) = nothing

#############################################################################
## colors

function Vis.λtoRGB(λ::T, gamma::T=T(0.8)) where {T<:Real}
    wavelength = λ * 1000 # λ is in um, need in nm
    if (wavelength >= 380 && wavelength <= 440)
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation)^gamma
        G = 0.0
        B = (1.0 * attenuation)^gamma
    elseif (wavelength >= 440 && wavelength <= 490)
        R = 0.0
        G = ((wavelength - 440) / (490 - 440))^gamma
        B = 1.0
    elseif (wavelength >= 490 && wavelength <= 510)
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490))^gamma
    elseif (wavelength >= 510 && wavelength <= 580)
        R = ((wavelength - 510) / (580 - 510))^gamma
        G = 1.0
        B = 0.0
    elseif (wavelength >= 580 && wavelength <= 645)
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580))^gamma
        B = 0.0
    elseif (wavelength >= 645 && wavelength <= 750)
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation)^gamma
        G = 0.0
        B = 0.0
    else
        R = 0.0
        G = 0.0
        B = 0.0
    end
    return RGB(R, G, B)
end

const hsvcolormap = Makie.to_colormap(:hsv)
# hues spaced by the golden ratio, so consecutive indices (sources, fields, hit counts) get clearly different colors
indexedcolor(i::Int) = Makie.interpolated_getindex(hsvcolormap, mod(i * 0.381966, 1.0))
function indexedcolor2(i::Int)
    c = Makie.interpolated_getindex(hsvcolormap, 1.0 - rem(i / (2.1 * π), 1.0))
    return RGBf(0.5 * c.r, 0.5 * c.g, 0.5 * c.b)
end

"""
Return a value that can be passed as a Makie `shading` keyword argument.
"""
default_shading(shaded::Bool=true) = shaded ? Makie.automatic : Makie.NoShading

#############################################################################
## helpers

# Convert an OpticSim TriangleMesh into a GeometryBasics mesh Makie can draw. Returns nothing for an empty mesh.
function tomesh(tmesh::TriangleMesh)
    points, indices = makiemesh(tmesh)
    if isempty(points) || isempty(indices)
        return nothing
    end
    pts = [Point3f(p) for p in points]
    faces = [GeometryBasics.GLTriangleFace(indices[i, 1], indices[i, 2], indices[i, 3]) for i in axes(indices, 1)]
    return GeometryBasics.Mesh(pts, faces)
end

# Draw 3D arrows. `tipsize` is the length of the arrow head in scene units.
function drawarrows!(ax, origins::AbstractVector, directions::AbstractVector; color=:yellow, tipsize=0.05, kwargs...)
    isempty(origins) && return nothing
    Makie.arrows3d!(ax, Point3f.(origins), Point3f.(directions); kwargs...,
        color=color,
        markerscale=1.0,
        minshaftlength=0.0,
        tiplength=tipsize,
        tipradius=tipsize / 3,
        shaftradius=tipsize / 10,
    )
end

#############################################################################
## imported meshes (.obj, .stl, ...)

function Vis.draw!(ax::Makie.AbstractAxis, ob::AbstractString;
    color=:gray,
    debug::Bool=false,  # make sure debug does not end up in kwargs (Makie would error)
    linewidth=3,
    shaded::Bool=true,
    wireframe::Bool=false,
    transform::Transform{Float64}=identitytransform(Float64),
    scale::Float64=1.0,
    kwargs...
)
    if !any(endswith(lowercase(ob), x) for x in (".obj", ".ply", ".2dm", ".off", ".stl"))
        @error "Unsupported file type"
        return nothing
    end
    meshdata = FileIO.load(ob)
    if transform != identitytransform(Float64) || scale != 1.0
        coords = [Point3f(transform * (scale * SVector{3,Float64}(p))) for p in GeometryBasics.coordinates(meshdata)]
        meshdata = GeometryBasics.Mesh(coords, GeometryBasics.faces(meshdata))
    end
    Makie.mesh!(ax, meshdata; kwargs..., color=color, shading=default_shading(shaded), visible=shaded)
    if wireframe
        Makie.wireframe!(ax, meshdata; color=shaded ? (:black, 0.1) : color, linewidth=linewidth)
    end
    return nothing
end

#############################################################################
## geometry

function Vis.draw!(ax::Makie.AbstractAxis, surf::Surface{T};
    numdivisions::Int=30,
    normals::Bool=false,
    normalcolor=:blue,
    kwargs...
) where {T<:Real}
    mesh = makemesh(surf, numdivisions)
    if nothing === mesh
        return nothing
    end
    Vis.draw!(ax, mesh; kwargs..., normals=false)
    if normals
        ndirs = vec(Point3f.(samplesurface(surf, normal, numdivisions ÷ 10)))
        norigins = vec(Point3f.(samplesurface(surf, point, numdivisions ÷ 10)))
        drawarrows!(ax, norigins, ndirs; color=normalcolor, tipsize=0.2)
    end
    return nothing
end

function Vis.draw!(ax::Makie.AbstractAxis, tmesh::TriangleMesh{T};
    debug::Bool=false,  # make sure debug does not end up in kwargs (Makie would error)
    linewidth=3,
    shaded::Bool=true,
    wireframe::Bool=false,
    color=:orange,
    normals::Bool=false,
    normalcolor=:blue,
    transparency::Bool=false,
    kwargs...
) where {T<:Real}
    mesh = tomesh(tmesh)
    if mesh !== nothing
        Makie.mesh!(ax, mesh; kwargs..., color=color, shading=default_shading(shaded), transparency=transparency, visible=shaded)
        if wireframe
            Makie.wireframe!(ax, mesh; color=shaded ? (:black, 0.1) : color, linewidth=linewidth)
        end
    end
    if normals
        @warn "Normals being drawn from triangulated mesh, precision may be low"
        norigins = [Point3f(centroid(t)) for t in tmesh.triangles[1:10:end]]
        ndirs = [Point3f(normal(t)) for t in tmesh.triangles[1:10:end]]
        drawarrows!(ax, norigins, ndirs; color=normalcolor, tipsize=0.2)
    end
    return nothing
end

function Vis.draw!(ax::Makie.AbstractAxis, meshes::Vararg{S};
    colors::Bool=false,
    kwargs...
) where {T<:Real,S<:Union{TriangleMesh{T},Surface{T}}}
    for i in 1:length(meshes)
        col = colors ? indexedcolor2(i) : :orange
        Vis.draw!(ax, meshes[i]; kwargs..., color=col)
    end
end

function Vis.draw(meshes::Vararg{S}; resolution=(1000, 1000), kwargs...) where {T<:Real,S<:Union{TriangleMesh{T},Surface{T}}}
    fig, lscene = Vis.scene(resolution)
    Vis.draw!(lscene, meshes...; kwargs...)
    Vis.display(fig)
    return isnoninteractive() ? fig : nothing
end

Vis.draw!(ax::Makie.AbstractAxis, csg::CSGTree{T}; numdivisions::Int=30, kwargs...) where {T<:Real} = Vis.draw!(ax, makemesh(csg, numdivisions); kwargs...)
Vis.draw!(ax::Makie.AbstractAxis, csg::CSGGenerator{T}; kwargs...) where {T<:Real} = Vis.draw!(ax, csg(); kwargs...)

function Vis.draw!(ax::Makie.AbstractAxis, bbox::BoundingBox{T}; kwargs...) where {T<:Real}
    p1 = Point3f(bbox.xmin, bbox.ymin, bbox.zmin)
    p2 = Point3f(bbox.xmin, bbox.ymax, bbox.zmin)
    p3 = Point3f(bbox.xmin, bbox.ymax, bbox.zmax)
    p4 = Point3f(bbox.xmin, bbox.ymin, bbox.zmax)
    p5 = Point3f(bbox.xmax, bbox.ymin, bbox.zmin)
    p6 = Point3f(bbox.xmax, bbox.ymax, bbox.zmin)
    p7 = Point3f(bbox.xmax, bbox.ymax, bbox.zmax)
    p8 = Point3f(bbox.xmax, bbox.ymin, bbox.zmax)
    Makie.linesegments!(ax, [p1, p2, p2, p3, p3, p4, p4, p1, p1, p5, p2, p6, p3, p7, p4, p8, p5, p6, p6, p7, p7, p8, p8, p5]; kwargs...)
end

#############################################################################
## optics

# `alpha` < 1 draws the elements semi-transparent, so rays inside them stay visible
function Vis.draw!(ax::Makie.AbstractAxis, ass::LensAssembly{T}; alpha=1.0, kwargs...) where {T<:Real}
    for (i, e) in enumerate(elements(ass))
        c = indexedcolor2(i)
        Vis.draw!(ax, e; kwargs..., color=RGBA(c.r, c.g, c.b, alpha), transparency=alpha < 1)
    end
end

function Vis.draw!(ax::Makie.AbstractAxis, sys::CSGOpticalSystem{T}; alpha=1.0, kwargs...) where {T<:Real}
    Vis.draw!(ax, sys.assembly; alpha=alpha, kwargs...)
    Vis.draw!(ax, sys.detector; kwargs...)
end

Vis.draw!(ax::Makie.AbstractAxis, sys::AxisymmetricOpticalSystem{T}; kwargs...) where {T<:Real} = Vis.draw!(ax, sys.system; kwargs...)

onlydetectorrays(system::Q, tracevalue::LensTrace{T,3}) where {T<:Real,Q<:AbstractOpticalSystem{T}} = onsurface(detector(system), point(tracevalue))

defaultraygenerator() = Sources.Source(
    transform=translation(0.0, 0.0, 10.0),
    origins=Origins.RectGrid(10.0, 10.0, 25, 25),
    directions=Directions.Constant(0.0, 0.0, -1.0),
)

# fewer rays, one wavelength: enough to see the ray paths without hiding the lenses
defaultdrawingrays() = Sources.Source(
    transform=translation(0.0, 0.0, 10.0),
    origins=Origins.RectGrid(10.0, 10.0, 9, 9),
    directions=Directions.Constant(0.0, 0.0, -1.0),
    spectrum=Spectrum.DeltaFunction(0.55),
)

## ray generators

# Ray generators without a specific drawing method are not drawn.
Vis.draw!(::Makie.AbstractAxis, ::AbstractRayGenerator; kwargs...) = nothing

function Vis.drawtracerays(system::AbstractOpticalSystem{T};
    raygenerator::AbstractRayGenerator{T}=defaultdrawingrays(),
    view::Symbol=:perspective,
    resolution::Tuple{Int,Int}=view === :side ? (1400, 500) : (1000, 1000),
    verbose::Bool=false,
    kwargs...
) where {T<:Real}
    verbose && println("Drawing System...")
    # the 3D axis box is no use in a side view, and it upsets the side view's camera fit, so leave it out
    fig, lscene = Vis.scene(resolution; show_axis=view !== :side)
    Vis.drawtracerays!(lscene, system; raygenerator=raygenerator, verbose=verbose, drawsys=true, drawgen=true, kwargs...)
    view === :side && Vis.sideview!(lscene)
    Vis.display(fig)
    return isnoninteractive() ? fig : nothing
end

function Vis.sideview!(ls::Makie.LScene; distance=nothing)
    cam = Makie.cameracontrols(ls.scene)
    cam.settings.projectiontype[] = Makie.Orthographic
    cam.settings.center[] = false # keep this view instead of re-centering when the figure is rendered
    bb = Makie.data_limits(ls.scene)
    center = Makie.Vec3d(bb.origin .+ bb.widths ./ 2)
    distance = something(distance, sidefitdistance(ls))
    # look along +x at the y-z plane, y up: the optical axis is horizontal and light travelling toward -z goes left to right
    Makie.update_cam!(ls.scene, cam, center .+ Makie.Vec3d(distance, 0, 0), center, Makie.Vec3d(0, 1, 0))
    return ls
end

# camera distance that fits the contents of `ls` in a side view. Orthographic: the visible half-height is the camera
# distance, the half-width that times the panel aspect ratio. Fits the system length (z) and height (y) with a margin.
function sidefitdistance(ls::Makie.LScene)
    bb = Makie.data_limits(ls.scene)
    viewport = Makie.widths(ls.scene.viewport[])
    aspect = viewport[2] > 0 ? viewport[1] / viewport[2] : 1.0
    return max(0.55 * bb.widths[3] / aspect, 0.6 * bb.widths[2])
end

function Vis.drawcomparison(systems::AbstractVector{<:AbstractOpticalSystem};
    raygenerator::AbstractRayGenerator,
    titles=nothing,
    view::Symbol=:side,
    resolution::Tuple{Int,Int}=(1400, 380 * length(systems)),
    kwargs...
)
    fig = Makie.Figure(size=resolution)
    lscenes = Makie.LScene[]
    for (i, system) in enumerate(systems)
        titles === nothing || Makie.Label(fig[2i-1, 1], titles[i], fontsize=18, tellwidth=false)
        lscene = Makie.LScene(fig[2i, 1]; show_axis=false, scenekw=(; camera=Makie.cam3d_cad!))
        Vis.drawtracerays!(lscene, system; raygenerator=raygenerator, drawsys=true, kwargs...)
        push!(lscenes, lscene)
    end
    if view === :side
        # one scale for all the drawings, so they can be compared: the one that fits the largest
        distance = maximum(sidefitdistance, lscenes)
        foreach(ls -> Vis.sideview!(ls; distance), lscenes)
    end
    lscene = last(lscenes)
    current_main_scene[] = fig
    current_3d_scene[] = lscene
    Vis.display(fig)
    return isnoninteractive() ? fig : nothing
end

function Vis.drawtracerays!(system::AbstractOpticalSystem; kwargs...)
    if current_3d_scene[] === nothing
        Vis.scene()
    end
    Vis.drawtracerays!(current_3d_scene[], system; kwargs...)
end

function Vis.drawtracerays!(ax::Makie.AbstractAxis, system::AbstractOpticalSystem{T};
    raygenerator::AbstractRayGenerator{T}=defaultdrawingrays(),
    test::Bool=false,
    trackallrays::Bool=false,
    colorbysourcenum::Bool=false,
    colorbynhits::Bool=false,
    rayfilter::Union{Nothing,Function}=onlydetectorrays,
    verbose::Bool=false,
    drawsys::Bool=false,
    drawgen::Bool=false,
    numdivisions::Int=30, # mesh resolution of the system's surfaces
    lensalpha=1.0, # opacity of the lens elements
    kwargs...
) where {T<:Real}
    raylines = Vector{LensTrace{T,3}}(undef, 0)

    # system first, rays on top: CairoMakie paints plots in order without depth testing, so the rays are drawn across
    # the lenses, as in a lens layout drawing. GLMakie depth-tests, so there rays inside a lens are hidden by it.
    drawgen && Vis.draw!(ax, raygenerator; kwargs...)
    drawsys && Vis.draw!(ax, system; numdivisions=numdivisions, alpha=lensalpha, kwargs...)

    verbose && println("Tracing...")
    for (i, r) in enumerate(raygenerator)
        if i % 1000 == 0 && verbose
            print("\r $i / $(length(raygenerator))")
        end
        allrays = Vector{LensTrace{T,3}}(undef, 0)
        if trackallrays
            res = trace(system, r, trackrays=allrays, test=test)
        else
            res = trace(system, r, test=test)
        end

        if trackallrays && !isempty(allrays)
            # filter on the trace that is hitting the detector
            if rayfilter === nothing || rayfilter(system, allrays[end])
                append!(raylines, allrays)
            end
        elseif res !== nothing
            if rayfilter === nothing || rayfilter(system, res)
                push!(raylines, res)
            end
        end
    end
    verbose && print("\r")

    verbose && println("Drawing Rays...")
    Vis.draw!(ax, raylines; colorbysourcenum=colorbysourcenum, colorbynhits=colorbynhits, kwargs...)
end

function Vis.drawtraceimage(system::AbstractOpticalSystem{T};
    raygenerator::AbstractRayGenerator{T}=defaultraygenerator(),
    test::Bool=false,
    verbose::Bool=false
) where {T<:Real}
    resetdetector!(system)
    start_time = time()
    for (i, r) in enumerate(raygenerator)
        if i % 1000 == 0 && verbose
            dif = round(time() - start_time, digits=1)
            left = round((time() - start_time) / i * (length(raygenerator) - i), digits=1)
            print("\rTraced: $i / $(length(raygenerator))     Elapsed: $(dif)s     Left: $(left)s")
        end
        trace(system, r, test=test)
    end
    verbose && print("\r")
    tot = round(time() - start_time, digits=1)
    verbose && println("Completed in $(tot)s")
    Vis.showimage(detectorimage(system))
    return detectorimage(system)
end

# Figure showing a 2D image with row 1 at the top. Extra dimensions (e.g. the 4D eyebox images) get a slider each.
function imagefigure(image::AbstractArray; resolution=(800, 800))
    image = Array(image)
    fig = Makie.Figure(size=resolution)
    ax = Makie.Axis(fig[1, 1], aspect=Makie.DataAspect(), yreversed=true)
    if ndims(image) <= 2
        slice = Makie.Observable(Float32.(image))
    else
        extradims = size(image)[3:end]
        sg = Makie.SliderGrid(fig[2, 1], [(label="dim $(d + 2)", range=1:n, startvalue=1) for (d, n) in enumerate(extradims)]...)
        indices = [s.value for s in sg.sliders]
        slice = Makie.lift((idx...) -> Float32.(image[:, :, idx...]), indices...)
    end
    hm = Makie.heatmap!(ax, Makie.lift(permutedims, slice); colormap=:grays)
    Makie.Colorbar(fig[1, 2], hm)
    return fig
end

#############################################################################
## rays, lines and points

function Vis.draw!(ax::Makie.AbstractAxis, rays::AbstractVector{<:AbstractRay{T,N}}; kwargs...) where {T<:Real,N}
    for r in rays
        Vis.draw!(ax, r; kwargs...)
    end
end

function tracecolor(trace::LensTrace, colorbysourcenum::Bool, colorbynhits::Bool)
    if colorbysourcenum
        color = indexedcolor(sourcenum(trace))
    elseif colorbynhits
        color = indexedcolor(nhits(trace))
    else
        color = Vis.λtoRGB(wavelength(trace))
    end
    return RGBA(color.r, color.g, color.b, sqrt(power(trace)))
end

# All traces are drawn with a single linesegments plot, which is much faster than one plot per trace.
function Vis.draw!(ax::Makie.AbstractAxis, traces::AbstractVector{LensTrace{T,N}};
    colorbysourcenum::Bool=false,
    colorbynhits::Bool=false,
    debug::Bool=false,
    kwargs...
) where {T<:Real,N}
    isempty(traces) && return nothing
    points = Vector{Point3f}(undef, 2 * length(traces))
    colors = Vector{RGBA{Float32}}(undef, 2 * length(traces))
    for (i, t) in enumerate(traces)
        points[2i-1] = Point3f(origin(ray(t)))
        points[2i] = Point3f(point(intersection(t)))
        colors[2i-1] = colors[2i] = tracecolor(t, colorbysourcenum, colorbynhits)
    end
    Makie.linesegments!(ax, points; kwargs..., color=colors, transparency=true)
end

function Vis.draw!(ax::Makie.AbstractAxis, trace::LensTrace{T,N}; kwargs...) where {T<:Real,N}
    Vis.draw!(ax, [trace]; kwargs...)
end

function Vis.draw!(ax::Makie.AbstractAxis, r::OpticalRay{T,N}; colorbysourcenum::Bool=false, colorbynhits::Bool=false, kwargs...) where {T<:Real,N}
    if colorbysourcenum
        color = indexedcolor(sourcenum(r))
    elseif colorbynhits
        color = indexedcolor(nhits(r))
    else
        color = Vis.λtoRGB(wavelength(r))
    end
    Vis.draw!(ax, ray(r); kwargs..., color=RGBA(color.r, color.g, color.b, sqrt(power(r))), transparency=true, rayscale=power(r))
end

function Vis.draw!(ax::Makie.AbstractAxis, ray::AbstractRay{T,N};
    color=:yellow,
    debug::Bool=false,  # make sure debug does not end up in kwargs (Makie would error)
    rayscale=1.0,
    kwargs...
) where {T<:Real,N}
    tipsize = min(0.05, rayscale * 0.05)
    drawarrows!(ax, [origin(ray)], [rayscale * direction(ray)]; kwargs..., color=color, tipsize=tipsize)
end

function Vis.draw!(ax::Makie.AbstractAxis, du::DisjointUnion{T}; kwargs...) where {T<:Real}
    Vis.draw!(ax, intervals(du); kwargs...)
end

function Vis.draw!(ax::Makie.AbstractAxis, intervals::AbstractVector{Interval{T}}; kwargs...) where {T<:Real}
    for i in intervals
        Vis.draw!(ax, i; kwargs...)
    end
end

function Vis.draw!(ax::Makie.AbstractAxis, interval::Interval{T}; kwargs...) where {T<:Real}
    if !(interval isa EmptyInterval)
        l = lower(interval)
        u = upper(interval)
        if !(l isa RayOrigin)
            Vis.draw!(ax, l)
        else
            @warn "Negative half space, can't draw ray origin"
        end
        if !(u isa Infinity)
            Vis.draw!(ax, u)
        else
            @warn "Positive half space, can't draw end point"
        end
        if !(l isa RayOrigin) && !(u isa Infinity)
            Vis.draw!(ax, (point(l), point(u)); kwargs...)
        end
    end
end

function Vis.draw!(ax::Makie.AbstractAxis, intersection::Intersection; normal::Bool=false, kwargs...)
    Vis.draw!(ax, point(intersection); kwargs...)
    if normal
        Vis.draw!(ax, Ray(point(intersection), OpticSim.normal(intersection)); kwargs...)
    end
end

function Vis.draw!(ax::Makie.AbstractAxis, lines::AbstractVector{Tuple{P,P}}; kwargs...) where {T<:Real,P<:AbstractVector{T}}
    for l in lines
        Vis.draw!(ax, l; kwargs...)
    end
end

function Vis.draw!(ax::Makie.AbstractAxis, line::Tuple{P,P};
    debug::Bool=false,  # make sure debug does not end up in kwargs (Makie would error)
    color=:yellow,
    kwargs...
) where {T<:Real,P<:AbstractVector{T}}
    Makie.linesegments!(ax, [Point3f(line[1]), Point3f(line[2])]; kwargs..., color=color)
end

function Vis.draw!(ax::Makie.AbstractAxis, point::AbstractVector{T}; kwargs...) where {T<:Real}
    Vis.draw!(ax, [point]; kwargs...)
end

function Vis.draw!(ax::Makie.AbstractAxis, points::AbstractVector{P};
    debug::Bool=false,  # make sure debug does not end up in kwargs (Makie would error)
    markersize=20,
    color=:black,
    kwargs...
) where {T<:Real,P<:AbstractVector{T}}
    Makie.scatter!(ax, Point3f.(points); markersize=markersize, color=color, strokewidth=0, kwargs...)
end

#############################################################################
## 2D plots: spline curves, OPD, spot diagrams, surface sag

function Vis.drawcurves(curves::Vararg{Spline}; numpoints::Int=200, resolution=(800, 800))
    fig = Makie.Figure(size=resolution)
    ax = Makie.Axis(fig[1, 1], aspect=Makie.DataAspect())
    us = range(0.0, 1.0, length=numpoints + 1)
    for (i, curve) in enumerate(curves)
        first = i == 1
        pts = [Point2f(point(curve, u)[1:2]...) for u in us]
        Makie.lines!(ax, pts, color=first ? RGBf(0, 1, 1) : RGBf(1, 0, 1), linewidth=first ? 5 : 1)
        controlpoints = [Point2f(p[1:2]...) for p in euclideancontrolpoints(curve)]
        Makie.lines!(ax, controlpoints, color=first ? :black : RGBf(0, 1, 0), linewidth=1)
        Makie.scatter!(ax, controlpoints, color=first ? :black : :red, markersize=first ? 10 : 6)
    end
    Vis.display(fig)
    return fig
end

function Vis.plotOPD(sys::AxisymmetricOpticalSystem{T}; resolution=(800, 600), kwargs...) where {T<:Real}
    fig = Makie.Figure(size=resolution)
    ax = Makie.Axis(fig[1, 1], xlabel="Relative position on entrance pupil", ylabel="OPD (waves)")
    Vis.plotOPD!(ax, sys; kwargs...)
    return fig
end

Vis.plotOPD!(sys::AxisymmetricOpticalSystem; kwargs...) = Vis.plotOPD!(Makie.current_axis(), sys; kwargs...)

function Vis.plotOPD!(ax::Makie.Axis, sys::AxisymmetricOpticalSystem{T}; label=nothing, color=nothing, collimated::Bool=true, axis::Int=1, samples::Int=5, wavelength::T=0.55, sourcepos::SVector{3,T}=SVector{3,T}(0.0, 0.0, 10.0), kwargs...) where {T<:Real}
    sysdiam = semidiameter(sys)
    rays = Vector{OpticalRay{T,3}}(undef, 0)
    positions = Vector{T}(undef, 0)
    for v in 0:samples
        s = (2 * (v / samples) - 1.0) * sysdiam
        p = axis == 2 ? SVector{3,T}(0, s, 0) : SVector{3,T}(s, 0, 0)
        push!(positions, s)
        if collimated
            dir = -sourcepos
            push!(rays, OpticalRay(p - dir, dir, 1.0, wavelength))
        else
            push!(rays, OpticalRay(sourcepos, p - sourcepos, 1.0, wavelength))
        end
    end
    chiefray = OpticalRay(sourcepos, -sourcepos, 1.0, wavelength)
    chiefres = trace(sys, chiefray, test=true)
    xs = Vector{T}()
    opls = Vector{T}()
    for (i, r) in enumerate(rays)
        lt = trace(sys, r, test=true)
        if !(nothing === lt)
            push!(xs, positions[i])
            push!(opls, (pathlength(lt) - pathlength(chiefres)) / (wavelength / 1000)) # wavelength is um while OPD is mm
        end
    end
    return Makie.lines!(ax, xs, opls; color=nothing === color ? Vis.λtoRGB(wavelength) : color, label=label, kwargs...)
end

function Vis.spotdiag(sys::CSGOpticalSystem{T}, raygenerator::OpticalRayGenerator{T}; size=(500, 500), kwargs...) where {T<:Real}
    fig = Makie.Figure(size=size)
    ax = Makie.Axis(fig[1, 1], xlabel="x (mm)", ylabel="y (mm)", aspect=Makie.DataAspect())
    Vis.spotdiag!(ax, sys, raygenerator; kwargs...)
    return fig
end
Vis.spotdiag(sys::AxisymmetricOpticalSystem{T}, raygenerator::OpticalRayGenerator{T}; kwargs...) where {T<:Real} = Vis.spotdiag(sys.system, raygenerator; kwargs...)

Vis.spotdiag!(sys::AbstractOpticalSystem, raygenerator::OpticalRayGenerator; kwargs...) = Vis.spotdiag!(Makie.current_axis(), sys, raygenerator; kwargs...)
Vis.spotdiag!(ax::Makie.Axis, sys::AxisymmetricOpticalSystem{T}, raygenerator::OpticalRayGenerator{T}; kwargs...) where {T<:Real} = Vis.spotdiag!(ax, sys.system, raygenerator; kwargs...)

function Vis.spotdiag!(ax::Makie.Axis, sys::CSGOpticalSystem{T}, raygenerator::OpticalRayGenerator{T}; color=nothing, label=nothing, kwargs...) where {T<:Real}
    xs = Vector{T}()
    ys = Vector{T}()
    det = detector(sys)
    for r in raygenerator
        lt = trace(sys, r, test=true)
        if !(nothing === lt)
            push!(xs, dot(point(lt) - centroid(det), det.uvec))
            push!(ys, dot(point(lt) - centroid(det), det.vvec))
        end
    end
    # relative to centroid
    if !isempty(xs)
        xs .-= sum(xs) / length(xs)
        ys .-= sum(ys) / length(ys)
    end
    if color === nothing
        color = Vis.λtoRGB(wavelength(first(raygenerator)))
    end
    return Makie.scatter!(ax, xs, ys; color=color, label=label, marker=:cross, markersize=8, kwargs...)
end

Vis.spotdiaggrid(sys::AxisymmetricOpticalSystem{T}, raygenerators::Vector{<:OpticalRayGenerator{T}}; kwargs...) where {T<:Real} = Vis.spotdiaggrid(sys.system, raygenerators; kwargs...)

function Vis.spotdiaggrid(sys::CSGOpticalSystem{T}, raygenerators::Vector{<:OpticalRayGenerator{T}}; colorbynum::Bool=false, kwargs...) where {T<:Real}
    nrows = Int(floor(sqrt(length(raygenerators))))
    ncols = cld(length(raygenerators), nrows)
    cellsize = 2000 ÷ max(nrows, ncols)
    fig = Makie.Figure(size=(ncols * cellsize, nrows * cellsize))
    axes = Makie.Axis[]
    for (i, r) in enumerate(raygenerators)
        row, col = divrem(i - 1, ncols) .+ 1
        ax = Makie.Axis(fig[row, col], xlabel="x (mm)", ylabel="y (mm)", aspect=Makie.DataAspect())
        if colorbynum
            Vis.spotdiag!(ax, sys, r; color=indexedcolor2(i), kwargs...)
        else
            Vis.spotdiag!(ax, sys, r; kwargs...)
        end
        push!(axes, ax)
    end
    # set all plots to the same range
    Makie.linkaxes!(axes...)
    return fig
end

function Vis.surfacesag(object::Union{CSGTree{T},Surface{T}}, resolution::Tuple{Int,Int}, halfsizes::Tuple{T,T}; offset::T=T(10), position::SVector{3,T}=SVector{3,T}(0.0, 0.0, 10.0), direction::SVector{3,T}=SVector{3,T}(0.0, 0.0, -1.0), rotationvec::SVector{3,T}=SVector{3,T}(0.0, 1.0, 0.0)) where {T<:Real}
    direction = normalize(direction)
    if abs(dot(rotationvec, direction)) == one(T)
        rotationvec = SVector{3,T}(1.0, 0.0, 0.0)
    end
    uvec = normalize(cross(normalize(rotationvec), direction))
    vvec = normalize(cross(direction, uvec))
    image = fill(T(NaN), resolution[2], resolution[1]) # image[y, x]
    Threads.@threads for x in 0:(resolution[1]-1)
        for y in 0:(resolution[2]-1)
            u = 2 * (x / (resolution[1] - 1)) - 1.0
            v = 2 * (y / (resolution[2] - 1)) - 1.0
            pos = position + (halfsizes[1] * u * uvec) + (halfsizes[2] * v * vvec)
            r = Ray(pos, direction)
            int = closestintersection(surfaceintersection(object, r), false)
            if int !== nothing
                int = int::Intersection{T,3}
                image[y+1, x+1] = -(α(int) - offset)
            end
        end
    end
    xs = range(-halfsizes[1], halfsizes[1], length=resolution[1])
    ys = range(-halfsizes[2], halfsizes[2], length=resolution[2])
    fig = Makie.Figure(size=(1000, 900))
    ax = Makie.Axis(fig[1, 1], aspect=Makie.DataAspect())
    hm = Makie.heatmap!(ax, xs, ys, permutedims(image); colormap=:jet1, nan_color=:transparent)
    Makie.Colorbar(fig[1, 2], hm, label="sag")
    return fig
end

#############################################################################
## eyebox evaluation

function Vis.eyebox_eval_eye(assembly::LensAssembly{T}, raygen::OpticalRayGenerator{T}, eye_rotation_x::T, eye_rotation_y::T, sample_points_x::Int, sample_points_y::Int; pupil_radius::T=T(2.0), resolution::Int=512, eye_transform::Transform{T}=identitytransform(T)) where {T<:Real}
    out_image = zeros(Float32, resolution, resolution, sample_points_y, sample_points_x)
    xrange = sample_points_x > 1 ? range(-eye_rotation_x, stop=eye_rotation_x, length=sample_points_x) : [0.0]
    yrange = sample_points_y > 1 ? range(-eye_rotation_y, stop=eye_rotation_y, length=sample_points_y) : [0.0]

    for (xi, erx) in enumerate(xrange)
        for (yi, ery) in enumerate(yrange)
            print("Position: ($xi, $yi) / ($sample_points_x, $sample_points_y)\r")
            r = OpticSim.rotmatd(T, ery, erx, 0.0)
            ov = rotation(eye_transform) * SVector{3,T}(0.0, 0.0, 13.0)
            t = r * ov - ov
            sys = OpticSim.Data.ModelEye(assembly, pupil_radius=pupil_radius, detpixels=resolution, transform=eye_transform * Transform(r, t))
            out_image[:, :, yi, xi] = traceMT(sys, raygen, printprog=false)
        end
    end
    # 4D image: shown in 2D with sliders for the other two dimensions
    Vis.showimage(out_image)
    return out_image
end

function Vis.eyebox_eval_planar(assembly::LensAssembly{T}, raygen::OpticalRayGenerator{T}, eyebox::Rectangle{T}, sample_points_x::Int, sample_points_y::Int, vsize::T; pupil_radius::T=T(2.0), resolution::Int=512) where {T<:Real}
    sys = CSGOpticalSystem(assembly, eyebox, 1, 1)
    hits = tracehitsMT(sys, raygen)
    println("Calculating results...")
    det = Rectangle(vsize, vsize, normal(eyebox), centroid(eyebox) - normal(eyebox))
    isempty(hits) && return nothing

    out_image = zeros(Float32, resolution, resolution, sample_points_y, sample_points_x)
    hsu = eyebox.halfsizeu - pupil_radius
    hsv = eyebox.halfsizev - pupil_radius
    xrange = sample_points_x > 1 ? range(-hsu, stop=hsu, length=sample_points_x) : [0.0]
    yrange = sample_points_y > 1 ? range(-hsv, stop=hsv, length=sample_points_y) : [0.0]
    for res in hits
        u, v = uv(res)
        x = u * eyebox.halfsizeu
        y = v * eyebox.halfsizev

        raydirection = normalize((centroid(eyebox) - point(res)) + direction(ray(res)) / dot(-normal(eyebox), direction(ray(res))))
        focal_plane_hit = surfaceintersection(det, Ray(point(res), raydirection))
        if focal_plane_hit isa EmptyInterval
            continue
        end
        px, py = OpticSim.uvtopix(det, uv(OpticSim.halfspaceintersection(focal_plane_hit)), (resolution, resolution))

        for (xi, posx) in enumerate(xrange)
            for (yi, posy) in enumerate(yrange)
                if (x - posx)^2 + (y - posy)^2 < pupil_radius^2
                    out_image[py, px, yi, xi] += OpticSim.sourcepower(ray(res))
                end
            end
        end
    end
    # 4D image: shown in 2D with sliders for the other two dimensions
    Vis.showimage(out_image)
    return out_image
end

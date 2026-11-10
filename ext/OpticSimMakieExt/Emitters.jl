# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# drawing of Emitters (sources and their origin distributions)

const ARROW_LENGTH = 0.5
const ARROW_SIZE = 0.01
const MARKER_SIZE = 1

#-------------------------------------
# draw debug information - local axes and positions
#-------------------------------------
function maybe_draw_debug_info(ax::Makie.AbstractAxis, o::Origins.AbstractOriginDistribution; transform::Geometry.Transform=Transform(), debug::Bool=false, kwargs...)
    debug || return nothing

    dir = forward(transform)
    uv = SVector{3}(right(transform))
    vv = SVector{3}(up(transform))
    pos = origin(transform)
    vsize = visual_size(o)

    # the origin and local axes of the distribution
    Makie.scatter!(ax, [Point3f(pos)], color=:blue, markersize=MARKER_SIZE * vsize)
    arrow_size = ARROW_SIZE * vsize
    drawarrows!(ax, [pos], [dir * ARROW_LENGTH * vsize]; color=:blue, tipsize=arrow_size)
    drawarrows!(ax, [pos], [uv * 0.5 * ARROW_LENGTH * vsize]; color=:red, tipsize=0.5 * arrow_size)
    drawarrows!(ax, [pos], [vv * 0.5 * ARROW_LENGTH * vsize]; color=:green, tipsize=0.5 * arrow_size)

    # all the sample origins
    positions = [Point3f(transform * x) for x in collect(o)]
    Makie.scatter!(ax, positions, color=:green, markersize=MARKER_SIZE * vsize)
    return nothing
end

#-------------------------------------
# point origin
#-------------------------------------
function Vis.draw!(ax::Makie.AbstractAxis, o::Origins.Point; transform::Geometry.Transform=Transform(), kwargs...)
    maybe_draw_debug_info(ax, o; transform=transform, kwargs...)
end

#-------------------------------------
# RectGrid and RectUniform origins
#-------------------------------------
function Vis.draw!(ax::Makie.AbstractAxis, o::Union{Origins.RectGrid,Origins.RectUniform}; transform::Geometry.Transform=Transform(), debug::Bool=false, kwargs...)
    dir = forward(transform)
    uv = SVector{3}(right(transform))
    vv = SVector{3}(up(transform))
    pos = origin(transform)

    plane = OpticSim.Plane(dir, pos)
    rect = OpticSim.Rectangle(plane, o.width / 2, o.height / 2, uv, vv)
    Vis.draw!(ax, rect; kwargs...)

    maybe_draw_debug_info(ax, o; transform=transform, debug=debug)
end

#-------------------------------------
# hexapolar origin
#-------------------------------------
function Vis.draw!(ax::Makie.AbstractAxis, o::Origins.Hexapolar; transform::Geometry.Transform=Transform(), debug::Bool=false, kwargs...)
    dir = forward(transform)
    uv = SVector{3}(right(transform))
    vv = SVector{3}(up(transform))
    pos = origin(transform)

    plane = OpticSim.Plane(dir, pos)
    ellipse = OpticSim.Ellipse(plane, o.halfsizeu, o.halfsizev, uv, vv)
    Vis.draw!(ax, ellipse; kwargs...)

    maybe_draw_debug_info(ax, o; transform=transform, debug=debug)
end

#-------------------------------------
# source
#-------------------------------------
function Vis.draw!(ax::Makie.AbstractAxis, s::Sources.AbstractSource{T}; parent_transform::Geometry.Transform=Transform(), debug::Bool=false, kwargs...) where {T<:Real}
    origins = Sources.origins(s)
    Vis.draw!(ax, origins; transform=parent_transform * Sources.transform(s), debug=debug, kwargs...)

    if debug
        starts = Point3f[]
        dirs = Point3f[]
        for optical_ray in s
            r = parent_transform * OpticSim.ray(optical_ray)
            push!(starts, Point3f(OpticSim.origin(r)))
            push!(dirs, Point3f(OpticSim.direction(r) * OpticSim.power(optical_ray) * ARROW_LENGTH * visual_size(origins)))
        end
        drawarrows!(ax, starts, dirs; color=:yellow, tipsize=ARROW_SIZE * visual_size(origins))
    end
end

#-------------------------------------
# optical rays, drawn as line segments of length proportional to power
#-------------------------------------
function Vis.draw!(ax::Makie.AbstractAxis, rays::AbstractVector{OpticSim.OpticalRay{T,3}};
    debug::Bool=false,  # make sure debug does not end up in kwargs (Makie would error)
    kwargs...
) where {T<:Real}
    points = Vector{Point3f}(undef, 2 * length(rays))
    for (i, optical_ray) in enumerate(rays)
        points[2i-1] = Point3f(OpticSim.origin(optical_ray))
        points[2i] = Point3f(OpticSim.origin(optical_ray) + OpticSim.direction(optical_ray) * OpticSim.power(optical_ray))
    end
    Makie.linesegments!(ax, points; kwargs..., color=:green, linewidth=2)
end

#-------------------------------------
# ray list source: mark the ray origins
#-------------------------------------
function Vis.draw!(ax::Makie.AbstractAxis, s::Sources.RayListSource; parent_transform::Geometry.Transform=Transform(), debug::Bool=false, markersize=8, color=:black, kwargs...)
    points = [Point3f(parent_transform * OpticSim.origin(r)) for r in s.rays]
    Makie.scatter!(ax, points; markersize=markersize, color=color)
end

#-------------------------------------
# composite source
#-------------------------------------
function Vis.draw!(ax::Makie.AbstractAxis, s::Sources.CompositeSource{T}; parent_transform::Geometry.Transform=Transform(), kwargs...) where {T<:Real}
    for source in s.sources
        Vis.draw!(ax, source; parent_transform=parent_transform * Sources.transform(s), kwargs...)
    end
end

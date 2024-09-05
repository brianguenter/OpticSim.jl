# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using ..OpticSim, .Geometry
using LinearAlgebra
using Distributions
using StaticArrays

import Makie

using .Emitters
using .Emitters.Spectrum
using .Emitters.Directions
using .Emitters.Origins
using .Emitters.AngularPower
using .Emitters.Sources

const ARRROW_LENGTH = 0.5
const ARRROW_SIZE = 0.01
const MARKER_SIZE = 1


#-------------------------------------
# draw debug information - local axes and positions
#-------------------------------------
function maybe_draw_debug_info(ax::Makie.AbstractAxis, o::Origins.AbstractOriginDistribution; transform::Geometry.Transform = Transform(), debug::Bool=false, kwargs...) where {T<:Real}

    dir = forward(transform)
    uv = SVector{3}(right(transform))
    vv = SVector{3}(up(transform))
    pos = origin(transform)

    if (debug)
        # draw the origin and normal of the surface
        Makie.scatter!(ax, pos, color=:blue, markersize = MARKER_SIZE * visual_size(o))

        # normal
        arrow_size = ARRROW_SIZE * visual_size(o)
        arrow_start = pos
        arrow_end = dir * ARRROW_LENGTH * visual_size(o) 
        Makie.arrows!(ax, [Makie.Point3f(arrow_start)], [Makie.Point3f(arrow_end)], arrowsize=arrow_size, linewidth=arrow_size * 0.5, linecolor=:blue, arrowcolor=:blue)
        arrow_end = uv * 0.5 * ARRROW_LENGTH * visual_size(o) 
        Makie.arrows!(ax, [Makie.Point3f(arrow_start)], [Makie.Point3f(arrow_end)], arrowsize= 0.5 * arrow_size, linewidth=arrow_size * 0.5, linecolor=:red, arrowcolor=:red)
        arrow_end = vv * 0.5 * ARRROW_LENGTH * visual_size(o) 
        Makie.arrows!(ax, [Makie.Point3f(arrow_start)], [Makie.Point3f(arrow_end)], arrowsize= 0.5 * arrow_size, linewidth=arrow_size * 0.5, linecolor=:green, arrowcolor=:green)

        # draw all the samples origins
        positions = map(x -> transform*x, collect(o))
        positions = collect(Makie.Point3f, positions)
        Makie.scatter!(ax, positions, color=:green, markersize = MARKER_SIZE * visual_size(o))
    end

end


#-------------------------------------
# draw point origin
#-------------------------------------
function OpticSim.Vis.draw!(ax::Makie.AbstractAxis, o::Origins.Point; transform::Geometry.Transform = Transform(), kwargs...) where {T<:Real}
        maybe_draw_debug_info(ax, o; transform=transform, kwargs...)
end

#-------------------------------------
# draw RectGrid and RectUniform origins
#-------------------------------------
function OpticSim.Vis.draw!(ax::Makie.AbstractAxis, o::Union{Origins.RectGrid, Origins.RectUniform}; transform::Geometry.Transform = Transform(), kwargs...) where {T<:Real}
    dir = forward(transform)
    uv = SVector{3}(right(transform))
    vv = SVector{3}(up(transform))
    pos = origin(transform)

    # @info "RECT: transform $(pos)"

    plane = OpticSim.Plane(dir, pos)
    rect = OpticSim.Rectangle(plane, o.width / 2, o.height / 2, uv, vv)
    
    OpticSim.Vis.draw!(ax, rect;  kwargs...)

    maybe_draw_debug_info(ax, o; transform=transform, kwargs...)
end


#-------------------------------------
# draw hexapolar origin
#-------------------------------------
function OpticSim.Vis.draw!(ax::Makie.AbstractAxis, o::Origins.Hexapolar; transform::Geometry.Transform = Transform(), kwargs...) where {T<:Real}
    dir = forward(transform)
    uv = SVector{3}(right(transform))
    vv = SVector{3}(up(transform))
    pos = origin(transform)

    plane = OpticSim.Plane(dir, pos)
    ellipse = OpticSim.Ellipse(plane, o.halfsizeu, o.halfsizev, uv, vv)
    
    OpticSim.Vis.draw!(ax, ellipse;  kwargs...)

    maybe_draw_debug_info(ax, o; transform=transform, kwargs...)
end

#-------------------------------------
# draw source
#-------------------------------------
function OpticSim.Vis.draw!(ax::Makie.AbstractAxis, s::S; parent_transform::Geometry.Transform = Transform(), debug::Bool=false, kwargs...) where {T<:Real,S<:Sources.AbstractSource{T}}
   
    OpticSim.Vis.draw!(ax, Emitters.Sources.origins(s);  transform=parent_transform * Emitters.Sources.transform(s), debug=debug, kwargs...)

    if (debug)
        m = zeros(T, length(s), 7)
        for (index, optical_ray) in enumerate(s)
            ray = OpticSim.ray(optical_ray)
            ray = parent_transform * ray
            m[index, 1:7] = [ray.origin... ray.direction... OpticSim.power(optical_ray)]
        end
        
        m[:, 4:6] .*= m[:, 7] * ARRROW_LENGTH * visual_size(Emitters.Sources.origins(s))  

        color = :yellow
        arrow_size = ARRROW_SIZE * visual_size(Emitters.Sources.origins(s))
        Makie.arrows!(ax, m[:,1], m[:,2], m[:,3], m[:,4], m[:,5], m[:,6]; kwargs...,  arrowcolor=color, linecolor=color, arrowsize=arrow_size, linewidth=arrow_size*0.5)
    end
end

#-------------------------------------
# draw optical rays
#-------------------------------------
function OpticSim.Vis.draw!(ax::Makie.AbstractAxis, rays::AbstractVector{OpticSim.OpticalRay{T, 3}};
    debug::Bool = false,  # make sure debug does not end up in kwargs (Makie would error)
    kwargs...
) where {T<:Real}
    m = zeros(T, length(rays)*2, 3)
    for (index, optical_ray) in enumerate(rays)
        ray = OpticSim.ray(optical_ray)
        m[(index-1)*2+1, 1:3] = [origin(optical_ray)...]
        m[(index-1)*2+2, 1:3] = [(OpticSim.origin(optical_ray) + OpticSim.direction(optical_ray) * 1 * OpticSim.power(optical_ray))... ]
    end
    
    color = :green
    Makie.linesegments!(ax, m[:,1], m[:,2], m[:,3]; kwargs...,  color = color, linewidth = 2, )
end

#-------------------------------------
# draw composite source
#-------------------------------------
function OpticSim.Vis.draw!(ax::Makie.AbstractAxis, s::Sources.CompositeSource{T}; parent_transform::Geometry.Transform = Transform(), kwargs...) where {T<:Real}
    for source in s.sources
        OpticSim.Vis.draw!(ax, source; parent_transform=parent_transform*Emitters.Sources.transform(s), kwargs...)
    end
end

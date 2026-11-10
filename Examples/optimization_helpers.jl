# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# Shared code for the lens optimization problems in optimization_examples.jl: the OptimizationExample type, paraxial
# solves, pupil rays, the merit function and a lens geometry check.

export OptimizationExample, rmsspotradii, surfaceseparations, fieldrays

"""
    OptimizationExample

A lens optimization problem.

- `name`, `source`: what the problem is and where its specification comes from
- `variables`: a description of each free variable, with units
- `x0`: starting values of the free variables
- `bounds`: `(lower, upper)` for each free variable
- `build(x)`: the optical system for variable values `x`
- `merit(x)`: the value to minimize, the mean over the fields of the squared RMS spot radius (mm²)
- `fields`: field angles in degrees, used by [`rmsspotradii`](@ref)

Optimization.jl takes the bounds as separate vectors: `lb = first.(ex.bounds)`, `ub = last.(ex.bounds)`.

The rays for each field and wavelength fill the entrance pupil in a hexapolar grid. The example functions take a
`rings` keyword: `rings` rings give 1 + 3 rings (rings + 1) rays, 127 for the default of 6.
"""
struct OptimizationExample{B,M}
    name::String
    source::String
    variables::Vector{String}
    x0::Vector{Float64}
    bounds::Vector{Tuple{Float64,Float64}}
    build::B
    merit::M
    fields::Vector{Float64}
    rays::Vector{Vector{Tuple{OpticalRay{Float64,3},Float64}}} # per field: (ray, weight)

    function OptimizationExample(name, source, variables, x0, bounds, build::B, merit::M, fields, rays) where {B,M}
        @assert length(variables) == length(x0) == length(bounds) "one name, starting value and bounds per variable"
        for (v, x, (l, u)) in zip(variables, x0, bounds)
            @assert l <= x <= u "$v: starting value $x is outside its bounds ($l, $u)"
        end
        new{B,M}(name, source, variables, x0, bounds, build, merit, fields, rays)
    end
end

function Base.show(io::IO, ex::OptimizationExample)
    println(io, ex.name)
    println(io, "  ", ex.source)
    for (v, x, (l, u)) in zip(ex.variables, ex.x0, ex.bounds)
        println(io, "  ", rpad(v, 44), " start ", rpad(round(x, sigdigits=6), 12), " bounds (", l, ", ", u, ")")
    end
end
#############################################################################
# Paraxial helpers, used to recreate Zemax solves and to find the entrance pupil

# largest curvature magnitude that keeps a surface's radius at least 1% larger than its semi-diameter
maxcurvature(semidiameter) = 1 / (1.01 * semidiameter)

# Paraxial (y, nu = n u) ray trace through surfaces with curvatures `c`, followed by thicknesses `t`, with refractive index
# `n[i]` after surface i and `n0` before the first. Returns the ray at the surface after the last thickness.
function paraxialtrace(y, nu, c, t, n, n0)
    nbefore = n0
    for i in eachindex(c)
        nu = nu - y * c[i] * (n[i] - nbefore)
        y = y + t[i] * nu / n[i]
        nbefore = n[i]
    end
    return y, nu
end

# Curvature of the last surface that gives effective focal length `f` in air, as Zemax's F number solve does.
# `c`, `t`, `n` describe the surfaces before the last one (see paraxialtrace).
function lastcurvature(f, c, t, n, n0=1.0)
    y, nu = paraxialtrace(one(f), zero(f), c, t, n, n0)
    return (nu + 1 / f) / (y * (1 - n[end]))
end

# Distance from the first surface to the entrance pupil, for an aperture stop after the surfaces `c`, `t`, `n`.
function entrancepupildistance(c, t, n, n0=1.0)
    isempty(c) && return 0.0
    a, _ = paraxialtrace(1.0, 0.0, c, t, n, n0) # parallel ray at height 1
    b, _ = paraxialtrace(0.0, n0, c, t, n, n0)  # ray through the first vertex at slope 1
    return b / a                                # the chief ray crosses the axis this far after the first vertex
end

#############################################################################
# Rays and merit function

# Collimated rays filling an entrance pupil of diameter `epd`, `zpupil` after the first lens vertex (which is at z = 0;
# light travels toward -z), at field angle `angle` in degrees. Returns (ray, weight) pairs, `weight` being the
# wavelength's weight.
function pupilrays(angle, epd, zpupil, wavelengths, weights; rings=6)
    d = SVector(0.0, sind(angle), -cosd(angle))
    pupil = collect(Origins.Hexapolar(rings, epd / 2, epd / 2))
    rays = Tuple{OpticalRay{Float64,3},Float64}[]
    for (λ, w) in zip(wavelengths, weights), p in pupil
        start = SVector(p[1], p[2], -zpupil) - 30.0 * d # start 30mm before the pupil, in front of the lens
        push!(rays, (OpticalRay(start, d, 1.0, λ), w))
    end
    return rays
end

# Weighted RMS radius of where `rays` hit the detector, about their weighted centroid, in mm. NaN if none hit.
function rmsspotradius(system, rays, ::Type{T}) where {T}
    points = SVector{3,T}[]
    weights = T[]
    for (r, w) in rays
        rT = OpticalRay(SVector{3,T}(OpticSim.origin(r)), SVector{3,T}(OpticSim.direction(r)), T(OpticSim.power(r)), T(OpticSim.wavelength(r)))
        res = trace(system, rT, test=true)
        if res isa LensTrace
            push!(points, point(res))
            push!(weights, T(w))
        end
    end
    isempty(points) && return T(NaN)
    centroid = sum(weights .* points) / sum(weights)
    return sqrt(sum(w * sum(abs2, p - centroid) for (w, p) in zip(weights, points)) / sum(weights))
end

# mean over the fields of the squared RMS spot radius
#
# NaN if a surface cannot span its clear aperture, which can happen to a radius set by a solve rather than bounded
# directly: such a system cannot be built. Solvers treat NaN as a failed step and backtrack.
function spotmerit(prescription::DataFrame, fieldrays, ::Type{T}) where {T}
    validsurfaces(prescription) || return T(NaN)
    system = AxisymmetricOpticalSystem{T}(prescription, 1, 1)
    return sum(rmsspotradius(system, rays, T)^2 for rays in fieldrays) / length(fieldrays)
end

# true if every lens surface reaches its semi-diameter: 1 - (1 + conic) c² h² > 0, and for spheres |radius| > h
function validsurfaces(p::DataFrame)
    for i in 2:(size(p, 1)-1)
        R, h = p[i, :Radius], p[i, :SemiDiameter]
        isfinite(R) || continue
        k = hasproperty(p, :Conic) && !ismissing(p[i, :Conic]) ? p[i, :Conic] : 0
        1 - (1 + k) * h^2 / R^2 > 0 || return false
    end
    return true
end

"""
    rmsspotradii(example::OptimizationExample, x) -> Vector

RMS spot radius about the centroid for each field of `example` with variables `x`, in μm.
"""
rmsspotradii(ex::OptimizationExample, x) = [1000 * rmsspotradius(ex.build(x), rays, Float64) for rays in ex.rays]

"""
    fieldrays(example::OptimizationExample) -> RayListSource

The rays `example` traces, at its middle wavelength, ready to draw: each field's rays get that field's number as
`sourcenum`, so drawing with `colorbysourcenum = true` colors each field differently. Use a smaller `rings` when
creating the example for a less cluttered drawing:

```julia
ex = Examples.singlet_optimization(rings=3)
Vis.drawtracerays(ex.build(ex.x0); raygenerator=Examples.fieldrays(ex), view=:side, trackallrays=true, test=true, colorbysourcenum=true)
```
"""
function fieldrays(ex::OptimizationExample)
    rays = OpticalRay{Float64,3}[]
    for (field, fieldrays) in enumerate(ex.rays)
        wavelengths = unique(OpticSim.wavelength(r) for (r, _) in fieldrays)
        λ = wavelengths[(length(wavelengths)+1)÷2]
        for (r, _) in fieldrays
            OpticSim.wavelength(r) == λ && push!(rays, OpticalRay(OpticSim.origin(r), OpticSim.direction(r), 1.0, λ; sourcenum=field))
        end
    end
    return Sources.RayListSource(rays)
end


#############################################################################
# Geometry check

# sag of a prescription surface at height h: conic plus even aspheric terms
function surfacesag(row, h)
    R = row.Radius
    k = hasproperty(row, :Conic) && !ismissing(row.Conic) ? row.Conic : 0.0
    z = 0.0
    if isfinite(R)
        c = 1 / R
        t = 1 - (1 + k) * c^2 * h^2
        t < 0 && return NaN # the surface does not exist at this height
        z = c * h^2 / (1 + sqrt(t))
    end
    if row.SurfaceType == "Aspheric" && hasproperty(row, :Parameters) && !ismissing(row.Parameters)
        for (order, a) in row.Parameters
            z += a * h^parse(Int, order)
        end
    end
    return z
end

"""
    surfaceseparations(example::OptimizationExample, x; samples = 200) -> Vector{NamedTuple}

For each pair of neighboring surfaces of `example`'s system with variables `x`: the smallest axial distance between
them over their shared clear aperture, and the height where it occurs. A negative distance means the surfaces cross: a
lens intersects itself (glass thickness) or its neighbor (air gap). NaN means a surface does not exist at some height
(its conic and curvature do not reach the semi-diameter).
"""
function surfaceseparations(ex::OptimizationExample, x; samples=200)
    p = ex.build(Float64.(x)).prescription
    result = NamedTuple[]
    for i in 2:(size(p, 1)-1)
        a, b = p[i, :], p[i+1, :]
        medium = ismissing(a.Material) || a.Material == AGFFileReader.Air ? "air" : "glass"
        semidiameter = min(a.SemiDiameter, b.SemiDiameter)
        heights = range(0.0, semidiameter, length=samples)
        separation = [a.Thickness + surfacesag(b, h) - surfacesag(a, h) for h in heights]
        j = argmin(replace(separation, NaN => -Inf))
        push!(result, (surfaces=(i - 1, i), medium=medium, center=a.Thickness, minimum=separation[j], height=heights[j]))
    end
    return result
end

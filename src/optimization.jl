# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# Lens optimization support: the LensOptimizationProblem type for describing an optimization problem, paraxial solves,
# emitters that fill the entrance pupil, a threaded spot size merit function and a lens geometry check. The example
# problems are in Examples/optimization_examples.jl. The "Optimization functions" section of docs/src/optimization.md
# says what each function is for and when to use it.

export LensOptimizationProblem, optimizelens, rmsspotradii, surfaceseparations, fieldrays, spotmerit_and_gradient

"""
    optimizelens(lens::LensOptimizationProblem; solver = Optim.Fminbox(Optim.LBFGS()), maxiters = 300, gradient = :auto) -> (u, retcode, solution)

Create an Optimization.jl problem for `lens`, with its starting values, bounds and merit function, and solve it. The
solver enforces the bounds. Returns the optimized variables `u`, the solver's return code and its full solution.

`gradient` chooses how gradients are computed:
- `:sequential`: the spot merit and its gradient from compiled sequential ray Jacobians ([`spotmerit_and_gradient`](@ref)).
  Exact, and typically 10 to 20 times faster than `:forwarddiff` for a whole optimization. Needs a `lens` created with a
  `prescription` function.
- `:forwarddiff`: ForwardDiff through `lens.merit`, which works for any merit function.
- `:auto` (the default): `:sequential` if `lens` has a `prescription`, otherwise `:forwarddiff`.

The solver works with each variable scaled to run from 0 to 1 across its bounds. Lens variables have very different
scales: a curvature change of 0.001/mm affects the image about as much as a spacing change of millimeters, so without
scaling the gradient is thousands of times larger for curvatures, and gradient based solvers barely move spacings.

Available when Optimization, OptimizationOptimJL and ForwardDiff are loaded:

```julia
using OpticSim, Optimization, OptimizationOptimJL, ForwardDiff
ex = Examples.singlet_optimization()
result = optimizelens(ex)
```
"""
function optimizelens end

"""
    LensOptimizationProblem

A lens optimization problem.

- `name`, `source`: what the problem is and where its specification comes from
- `variables`: a description of each free variable, with units
- `x0`: starting values of the free variables
- `bounds`: `(lower, upper)` for each free variable
- `build(x)`: the optical system for variable values `x`
- `merit(x)`: the value to minimize, the mean over the fields of the squared RMS spot radius (mm²)
- `fields`: field angles in degrees, used by [`rmsspotradii`](@ref)
- `rays`: for each field, the rays traced at `x0`, as (ray, weight) pairs
- `aimedrays(x)`: the rays traced for variables `x`; the same as `rays` unless the problem aims its rays at the
  entrance pupil of each design
- `prescription(x)`: optional (keyword argument), the lens prescription (DataFrame) for variables `x`, the one `build(x)`
  makes a system of. Give it when `merit(x)` is the spot merit of that prescription traced with `aimedrays(x)`
  (`OpticSim.spotmerit(prescription(x), aimedrays(x), T)`): [`optimizelens`](@ref) then computes the merit and its
  gradient with compiled sequential ray Jacobians ([`spotmerit_and_gradient`](@ref)), many times faster than ForwardDiff
  through the ray trace. The constructor checks that the two merits agree at `x0`.

Optimization.jl takes the bounds as separate vectors: `lb = first.(ex.bounds)`, `ub = last.(ex.bounds)`.

The rays for each field and wavelength come from an emitter ([`OpticSim.pupilsource`](@ref)) that fills the entrance pupil in a
hexapolar grid. The example problems in `Examples` take a `rings` keyword: `rings` rings give 1 + 3 rings (rings + 1) rays, 127 for
6.
"""
struct LensOptimizationProblem{B,M,A,P}
    name::String
    source::String
    variables::Vector{String}
    x0::Vector{Float64}
    bounds::Vector{Tuple{Float64,Float64}}
    build::B
    merit::M
    fields::Vector{Float64}
    rays::Vector{Vector{Tuple{OpticalRay{Float64,3},Float64}}} # per field: (ray, weight)
    aimedrays::A
    prescription::P # nothing, or x -> prescription DataFrame

    function LensOptimizationProblem(name, source, variables, x0, bounds, build::B, merit::M, fields, rays; aimedrays::A=(_ -> rays), prescription::P=nothing) where {B,M,A,P}
        @assert length(variables) == length(x0) == length(bounds) "one name, starting value and bounds per variable"
        for (v, x, (l, u)) in zip(variables, x0, bounds)
            @assert l <= x <= u "$v: starting value $x is outside its bounds ($l, $u)"
        end
        if prescription !== nothing
            m = merit(x0)
            s = first(sequentialspotmerit(prescription(x0), aimedrays(x0); withgradient=false))
            @assert isapprox(m, s, rtol=1e-8) "merit(x0) = $m differs from the sequential spot merit of prescription(x0), $s: with a prescription, merit(x) must be the spot merit of prescription(x) traced with aimedrays(x)"
        end
        new{B,M,A,P}(name, source, variables, x0, bounds, build, merit, fields, rays, aimedrays, prescription)
    end
end

function Base.show(io::IO, ex::LensOptimizationProblem)
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
# Clear apertures that follow the design

"""
    axialrayheights(prescription, height, λ) -> Vector

Real ray trace of a ray parallel to the axis at `height` (mm), wavelength `λ` (μm), through `prescription`, in the
meridional plane. Returns the ray's height at each row from the first lens surface to the last row before the image
(lens surfaces and the stop), NaN from where it misses a surface or is totally internally reflected. Spherical and flat
surfaces only.
"""
function axialrayheights(p::DataFrame, height, λ)
    n = size(p, 1)
    T = float(promote_type(typeof(height), nonmissingtype(eltype(p.Thickness)), nonmissingtype(eltype(p.Radius))))
    nindex(i) = isglass(p, i) ? T(AGFFileReader.index(p[i, :Material], λ)) : one(T)
    y, dz, dy = T(height), one(T), zero(T) # height and direction; the ray starts at the first surface's vertex plane
    nbefore = one(T)
    heights = T[]
    for i in 2:n-1
        row = p[i, :]
        conic = hasproperty(row, :Conic) && !ismissing(row.Conic) ? row.Conic : 0.0
        (row.SurfaceType == "Aspheric" || !iszero(conic)) && throw(ArgumentError("axialrayheights: row $i is not spherical"))
        c = isfinite(row.Radius) ? 1 / row.Radius : zero(T)
        # from the vertex plane to the surface: c (z² + y²) - 2z = 0 along the ray, the root nearest the plane
        b = dz - c * y * dy
        disc = b^2 - c^2 * y^2
        disc < 0 && return [heights; fill(T(NaN), n - 1 - i)] # misses the surface
        t = c * y^2 / (b + sqrt(disc))
        z, y = t * dz, y + t * dy
        push!(heights, y)
        # refract: the unit normal (1 - c z, -c y) points along the ray
        nafter = nindex(i)
        if nafter != nbefore
            nz, ny = 1 - c * z, -c * y
            cosi = dz * nz + dy * ny
            μ = nbefore / nafter
            k = 1 - μ^2 * (1 - cosi^2)
            k < 0 && return [heights; fill(T(NaN), n - 1 - i)] # total internal reflection
            g = sqrt(k) - μ * cosi
            dz, dy = μ * dz + g * nz, μ * dy + g * ny
        end
        nbefore = nafter
        # on to the next vertex plane
        i < n - 1 && (t = (p[i, :Thickness] - z) / dz; y += t * dy)
    end
    return heights
end

"""
    passaxialbeam(prescription, height, wavelengths) -> DataFrame

A copy of `prescription` whose lens surfaces are each at least large enough to pass the axial beam: a semi-diameter at
least the height of its marginal ray, the ray parallel to the axis at `height` (the entrance pupil radius), at each of
`wavelengths`. If the stop is smaller than that beam, the marginal ray is the one through the edge of the stop instead.
The stop keeps its size.

Use it in an optimization's `build(x)` so the clear apertures follow the design, like a lens design program's automatic
semi-diameters: only the stop then limits the axial beam, and the optimizer can't stop the lens down with another
surface. A semi-diameter is never made smaller than in `prescription`, so off-axis vignetting by the rims stays. If the
beam can't pass a surface the semi-diameters become NaN and the system is invalid ([`OpticSim.surfacescross`](@ref)).
"""
function passaxialbeam(p::DataFrame, height, wavelengths)
    stop = findfirst(==("Stop"), p.SurfaceType)
    heights = map(wavelengths) do λ
        h = axialrayheights(p, height, λ)
        if stop !== nothing && !isnan(h[stop-1]) && abs(h[stop-1]) > p[stop, :SemiDiameter]
            # the stop limits the beam: find the entrance height whose ray just reaches the stop's edge
            lo, hi = 0.0, Float64(primalvalue(height))
            for _ in 1:50
                mid = (lo + hi) / 2
                abs(primalvalue(axialrayheights(p, mid, λ)[stop-1])) > p[stop, :SemiDiameter] ? (hi = mid) : (lo = mid)
            end
            h = axialrayheights(p, lo, λ)
        end
        abs.(h)
    end
    q = copy(p)
    S = promote_type(nonmissingtype(eltype(p.SemiDiameter)), eltype(first(heights)))
    sd = S[i == 1 || i == size(p, 1) || i == stop ? p[i, :SemiDiameter] : max(p[i, :SemiDiameter], maximum(h[i-1] for h in heights)) for i in 1:size(p, 1)]
    q[!, :SemiDiameter] = sd
    return q
end

#############################################################################
# Rays and merit function

"""
    pupilsource(angle, epd, zpupil, λ; rings = 6, distance = 30.0) -> Sources.Source

An emitter for a field angle: a collimated beam filling an entrance pupil of diameter `epd`, `zpupil` after the first
lens vertex (which is at z = 0; light travels toward -z), at field angle `angle` in degrees, wavelength `λ` in μm. Its
origins are a hexapolar grid with `rings` rings, `1 + 3 rings (rings + 1)` rays, in a plane parallel to the pupil
`distance` before it along the beam, which must put them in front of the lens.
"""
function pupilsource(angle, epd, zpupil, λ; rings=6, distance=30.0)
    # the number type of the pupil position: with automatic differentiation the pupil can depend on the optimization
    # variables, and the source then generates rays of that type
    T = float(promote_type(typeof(zpupil), typeof(distance)))
    d = SVector{3,T}(0, sind(angle), -cosd(angle))
    # the grid's origins are in the z = 0 plane of the source's frame; translating the frame shifts them to the pupil
    # plane, moved `distance` back along the beam
    Emitters.Sources.Source(T; transform=translation(SVector{3,T}(0, 0, -zpupil) - distance * d),
        origins=Emitters.Origins.Hexapolar(rings, T(epd / 2), T(epd / 2)), directions=Emitters.Directions.Constant(d),
        spectrum=Emitters.Spectrum.DeltaFunction(T(λ)))
end

# The rays of `pupilsource` for each wavelength, as (ray, weight) pairs, `weight` being the wavelength's weight.
function pupilrays(angle, epd, zpupil, wavelengths, weights; rings=6, distance=30.0)
    [(r, w) for (λ, w) in zip(wavelengths, weights) for r in pupilsource(angle, epd, zpupil, λ; rings, distance)]
end

# Weighted RMS radius of where `rays` hit the detector, about their weighted centroid, in mm. NaN if none hit.
#
# Only rays that pass through each lens surface once, in order, count, as in a sequential lens design program. OpticSim
# traces non-sequentially, so a ray can also reach the image by a stray path, e.g. reflecting inside a thick element or
# skimming a lens rim; such light is flare, not part of the image.
rmsspotradius(system, rays, ::Type{T}) where {T} = rmsradius(only(fieldhits(system, [rays], T)))

# weighted RMS radius of (point, weight) pairs about their weighted centroid; NaN if there are none
function rmsradius(hits::AbstractVector{Tuple{SVector{3,T},T}}) where {T}
    isempty(hits) && return T(NaN)
    total = sum(last, hits)
    centroid = sum(w * p for (p, w) in hits) / total
    return sqrt(sum(w * sum(abs2, p - centroid) for (p, w) in hits) / total)
end

"""
    fieldhits(system, fieldrays, T) -> Vector{Vector{Tuple{SVector{3,T},T}}}

Traces each field's rays through `system` and returns, for each field, the image points of the rays that pass through
each lens surface once, in order, with their weights.

The rays are traced in parallel from a work queue (start Julia with several threads, e.g. `julia -t 8`, to use this).
The rays of all fields are split into small chunks; one worker per thread repeatedly takes the next chunk from a shared
atomic counter, so faster cores, or threads whose rays happen to be cheap (vignetted rays stop early), take more chunks
and no thread sits idle at the end. The workers are `:static`, so a worker never moves between threads mid trace, as
OpticSim's interval pools are per thread. Each worker traces through its own copy of the optical system, because tracing
writes to the system's detector image; the copies share the lens geometry. The results are kept per chunk and combined in
chunk order, so the result is the same whichever thread traced which chunk, and for any number of threads.
"""
function fieldhits(system::AxisymmetricOpticalSystem, fieldrays, ::Type{T}; chunksize::Int=16) where {T}
    sequentialhits = lenssurfaces(system.prescription)
    jobs = [(field, r, w) for (field, rays) in enumerate(fieldrays) for (r, w) in rays]
    chunks = collect(Iterators.partition(eachindex(jobs), chunksize))
    results = [Tuple{Int,SVector{3,T},T}[] for _ in chunks]
    next = Threads.Atomic{Int}(1) # the next chunk to trace
    function worker()
        csg = copy(system.system)
        while (k = Threads.atomic_add!(next, 1)) <= length(chunks)
            for j in chunks[k]
                field, r, w = jobs[j]
                rT = OpticalRay(SVector{3,T}(origin(r)), SVector{3,T}(direction(r)), T(power(r)), T(wavelength(r)))
                res = trace(csg, rT, test=true)
                if res isa LensTrace && nhits(res) == sequentialhits
                    push!(results[k], (field, point(res), T(w)))
                end
            end
        end
    end
    nworkers = min(Threads.nthreads(), length(chunks))
    if nworkers <= 1
        worker()
    else
        Threads.@threads :static for _ in 1:nworkers
            worker()
        end
    end
    hits = [Tuple{SVector{3,T},T}[] for _ in fieldrays]
    for result in results, (field, p, w) in result
        push!(hits[field], (p, w))
    end
    return hits
end

# mean over the fields of the squared RMS spot radius
#
# NaN if a surface cannot span its clear aperture, which can happen to a radius set by a solve rather than bounded
# directly: such a system cannot be built. Solvers treat NaN as a failed step and backtrack.
function spotmerit(prescription::DataFrame, fieldrays, ::Type{T}) where {T}
    (validsurfaces(prescription) && !surfacescross(prescription)) || return T(NaN)
    system = AxisymmetricOpticalSystem{T}(prescription, 1, 1)
    return sum(rmsradius(h)^2 for h in fieldhits(system, fieldrays, T)) / length(fieldrays)
end

# number of surfaces a ray refracts at on its sequential path: the rows where the material changes (a stop in air
# isn't one, a stop on a lens surface is)
lenssurfaces(p::DataFrame) = count(i -> !isequal(p[i, :Material], p[i-1, :Material]), 2:(size(p, 1)-1))

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
    rmsspotradii(lens::LensOptimizationProblem, x) -> Vector

RMS spot radius about the centroid for each field of `lens` with variables `x`, in μm.
"""
rmsspotradii(lens::LensOptimizationProblem, x) = [1000 * rmsspotradius(lens.build(x), rays, Float64) for rays in lens.aimedrays(x)]

"""
    fieldrays(lens::LensOptimizationProblem) -> RayListSource

The rays `lens` traces, at its middle wavelength, ready to draw: each field's rays get that field's number as
`sourcenum`, so drawing with `colorbysourcenum = true` colors each field differently. Use a smaller `rings` when
creating the problem for a less cluttered drawing:

```julia
ex = Examples.singlet_optimization(rings=3)
Vis.drawtracerays(ex.build(ex.x0); raygenerator=fieldrays(ex), view=:side, trackallrays=true, test=true, colorbysourcenum=true)
```
"""
function fieldrays(lens::LensOptimizationProblem)
    rays = OpticalRay{Float64,3}[]
    for (field, fieldrays) in enumerate(lens.rays)
        wavelengths = unique(wavelength(r) for (r, _) in fieldrays)
        λ = wavelengths[(length(wavelengths)+1)÷2]
        for (r, _) in fieldrays
            wavelength(r) == λ && push!(rays, OpticalRay(origin(r), direction(r), 1.0, λ; sourcenum=field))
        end
    end
    return Emitters.Sources.RayListSource(rays)
end


#############################################################################
# Geometry check
#
# A design can only be built if its parts are in order along the axis at every height: no lens element passes through its
# own other surface (negative glass thickness) or into the next element, the stop or the image. The check follows how
# AxisymmetricOpticalSystem builds the parts:
#
# - a lens element (a glass row of the prescription, between surfaces i and i + 1) is a solid of radius equal to the larger
#   of its two surfaces' semi-diameters. A spherical surface ends at its own semi-diameter, or at its sphere's equator
#   if the sphere is too small, and a flat ring continues from there to the rim. Conic and aspheric surfaces continue past
#   their own semi-diameter to the rim.
# - the stop is an opaque plane with a hole of its semi-diameter; it extends without limit outside the hole.
# - the image is a plane of its semi-diameter.
#
# Cemented elements share a surface, so their separation there is zero, which is allowed.

# sag of a prescription surface at height h: conic plus even aspheric terms. NaN where the surface doesn't exist.
function prescriptionsag(row, h)
    R = coalesce(row.Radius, Inf)
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

# sag of a surface as OpticSim builds it into a lens element: a spherical surface ends at its own semidiameter, or at its
# sphere's equator if that comes first, and a flat ring at that height continues to the element rim. Conic and aspheric
# surfaces continue to the rim.
function elementsag(row, h)
    R = coalesce(row.Radius, Inf)
    conic = hasproperty(row, :Conic) && !ismissing(row.Conic) ? row.Conic : 0.0
    if isfinite(R) && iszero(conic) && row.SurfaceType != "Aspheric"
        stop = min(coalesce(row.SemiDiameter, Inf), abs(R))
        h >= stop && return stop >= abs(R) ? R : prescriptionsag(row, stop) # the flat ring
    end
    return prescriptionsag(row, h)
end

isglass(p::DataFrame, i) = !(ismissing(p[i, :Material]) || p[i, :Material] == AGFFileReader.Air)

# The solid parts along the axis, in order: lens elements, the stop and the image. Each is (front, back, inner, outer): its
# front and back surface rows, and the heights between which it is present.
function lensparts(p::DataFrame)
    n = size(p, 1)
    parts = Tuple{Int,Int,Float64,Float64}[]
    for i in 2:n
        sd = Float64(primalvalue(coalesce(p[i, :SemiDiameter], Inf)))
        if i == n
            push!(parts, (i, i, 0.0, sd))           # the image
        elseif p[i, :SurfaceType] == "Stop"
            push!(parts, (i, i, sd, Inf))           # opaque outside the hole
        elseif isglass(p, i)
            push!(parts, (i, i + 1, 0.0, max(sd, Float64(primalvalue(p[i+1, :SemiDiameter])))))
        end
    end
    return parts
end

# axial position of surface row j at height h, the first surface's vertex being at 0
function surfacez(p::DataFrame, j, h)
    vertex = sum((p[k, :Thickness] for k in 2:(j-1)); init=0.0)
    return vertex + elementsag(p[j, :], h)
end

# heights to check: evenly spaced, plus every height where a part starts or ends or a flat ring begins, where the closest
# approach often is
function checkheights(p::DataFrame, parts, samples)
    top = maximum(outer for (_, _, _, outer) in parts if isfinite(outer))
    hs = collect(range(0.0, top, length=samples))
    for (_, _, inner, outer) in parts
        isfinite(outer) && push!(hs, outer)
        inner > 0 && push!(hs, inner)
    end
    for i in 2:size(p, 1)
        R = Float64(primalvalue(coalesce(p[i, :Radius], Inf)))
        isfinite(R) && abs(R) < top && push!(hs, abs(R))
    end
    return sort!(unique!(hs))
end

# Calls f(kind, a, b, separation, h) for each height h and each clearance there: the glass thickness of each element
# ("glass", between its surfaces a and b) and the gap from each part to the next part present at that height ("air", from
# surface a to surface b).
function foreachclearance(f, p::DataFrame; samples=32)
    parts = lensparts(p)
    for h in checkheights(p, parts, samples)
        previous = nothing
        for (front, back, inner, outer) in parts
            inner <= h <= outer || continue
            back != front && f("glass", front, back, surfacez(p, back, h) - surfacez(p, front, h), h)
            if previous !== nothing && previous != front # cemented elements share a surface; no gap to check
                f("air", previous, front, surfacez(p, front, h) - surfacez(p, previous, h), h)
            end
            previous = back
        end
    end
end

"""
    surfacescross(prescription; samples = 32) -> Bool

True if a lens element passes through itself (negative glass thickness at some height) or into the next element, the
stop or the image, or if a surface doesn't exist out to its element's rim. Such a system can't be built, so
`spotmerit` returns NaN for it. Follows how OpticSim builds elements: see [`surfaceseparations`](@ref).
"""
function surfacescross(p::DataFrame; samples=32)
    crossed = false
    foreachclearance(p; samples) do _, _, _, separation, _
        (isnan(separation) || separation < 0) && (crossed = true)
    end
    return crossed
end

"""
    surfaceseparations(lens::LensOptimizationProblem, x; samples = 200) -> Vector{NamedTuple}

The clearances of the system of `lens` with variables `x`: the thinnest point of each lens element, and the smallest gap
between each part (lens element or stop) and the part after it, each with the height where it occurs. Negative means
the parts intersect; NaN that a surface doesn't reach its element's rim (its conic and curvature don't allow it).

The parts are what OpticSim builds: an element's radius is the larger semi-diameter of its two surfaces, a spherical
surface too small to reach it ends in a flat ring, and the stop is opaque outside its hole. The gaps are checked at
every height, against whatever part comes next there: outside a small element, a larger element's rim is checked against
the part beyond it.

Each entry has `surfaces` (the two surfaces, numbered from 1 for the first lens surface; the stop has a number too),
`medium` ("glass" or "air"), `center` (their distance on the axis), `minimum` and `height`.
"""
function surfaceseparations(lens::LensOptimizationProblem, x; samples=200)
    p = lens.build(Float64.(x)).prescription
    smallest = Dict{Tuple{String,Int,Int},Tuple{Float64,Float64}}()
    foreachclearance(p; samples) do kind, a, b, separation, h
        key = (kind, a, b)
        s = isnan(separation) ? -Inf : separation
        (!haskey(smallest, key) || s < smallest[key][1]) && (smallest[key] = (s, h))
    end
    result = NamedTuple[]
    for ((kind, a, b), (s, h)) in sort!(collect(smallest), by=e -> (e[1][2], e[1][3], e[1][1] == "air"))
        center = surfacez(p, b, 0.0) - surfacez(p, a, 0.0)
        push!(result, (surfaces=(a - 1, b - 1), medium=kind, center=center, minimum=s == -Inf ? NaN : s, height=h))
    end
    return result
end

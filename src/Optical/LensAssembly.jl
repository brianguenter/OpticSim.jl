# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    LensAssembly{T<:Real}

Structure which contains the elements of the optical system, these can be [`CSGTree`](@ref) or [`Surface`](@ref) objects.

Create one with:
```julia
LensAssembly(elements::Union{Surface{T},CSGTree{T},LensAssembly{T}}...; axis = SVector(0.0, 0.0, 1.0)) where {T<:Real}
```
At least one element is required. Elements which are themselves `LensAssembly`s are unpacked into the new assembly.

The result is a [`TupleLensAssembly`](@ref).
"""
abstract type LensAssembly{T<:Real} end
export LensAssembly

Base.show(io::IO, ::Type{<:LensAssembly}) = print(io, "LensAssembly")
Base.show(io::IO, a::LensAssembly{T}) where {T<:Real} = print(io, "LensAssembly{$T}($(length(elements(a))) elements)")

"""
    TupleLensAssembly{T<:Real,G<:Tuple} <: LensAssembly{T}

The [`LensAssembly`](@ref) created by `LensAssembly(elements...)`.

```julia
struct TupleLensAssembly{T,G<:Tuple} <: LensAssembly{T}
    axis::SVector{3,T}
    groups::G # Tuple{Vector{E1}, Vector{E2}, ...}: one Vector per distinct element type
end
```

Elements are grouped by type, with one `Vector` per distinct element type. Intersection code loops over each `Vector`
knowing the concrete element type, so there is no dynamic dispatch. The size of the type depends only on the number of
distinct element types, not the number of elements, so it stays small (and fast to compile) for assemblies with many
elements.
"""
struct TupleLensAssembly{T<:Real,G<:Tuple} <: LensAssembly{T}
    axis::SVector{3,T}
    groups::G
end
export TupleLensAssembly

# first element is a separate argument so T is always bound (at least one element is required)
function LensAssembly(first::Union{Surface{T},CSGTree{T},LensAssembly{T}}, rest::Vararg{Union{Surface{T},CSGTree{T},LensAssembly{T}}}; axis::SVector{3,T}=SVector{3,T}(0.0, 0.0, 1.0)) where {T<:Real}
    # unpack any nested lens assemblies
    flat = Any[]
    for e in (first, rest...)
        if e isa LensAssembly
            append!(flat, elements(e))
        else
            push!(flat, e)
        end
    end

    # group by type, in order of first appearance
    types = DataType[]
    groups = Vector[]
    for e in flat
        i = findfirst(==(typeof(e)), types)
        if i === nothing
            push!(types, typeof(e))
            push!(groups, typeof(e)[]) # Vector with the concrete element type
            i = length(groups)
        end
        push!(groups[i], e)
    end

    # Tuple(groups) has type Tuple{Vector{E1}, Vector{E2}, ...}, which becomes part of the assembly's type
    return TupleLensAssembly(normalize(axis), Tuple(groups))
end

"""
    elements(ass::LensAssembly)

All elements of `ass`, as a `Tuple`.
"""
elements(ass::TupleLensAssembly) = Tuple(Iterators.flatten(ass.groups))
export elements

# closest intersection of r with a single element, or nothing if it misses
function elementintersection(element, r::AbstractRay{T,N})::Union{Nothing,Intersection{T,N}} where {T<:Real,N}
    intvl = surfaceintersection(element, r)
    if intvl isa EmptyInterval
        return nothing
    end
    return closestintersection(intvl::Union{Interval{T},DisjointUnion{T}})
end

# whichever of two intersections is closer along the ray; the first one wins ties
function closer(a::Union{Nothing,Intersection}, b::Union{Nothing,Intersection})
    if a === nothing
        return b
    elseif b === nothing || α(a) <= α(b)
        return a
    else
        return b
    end
end

function closestintersection(ass::TupleLensAssembly{T}, r::AbstractRay{T,N})::Union{Nothing,Intersection{T,N}} where {T<:Real,N}
    closest::Union{Nothing,Intersection{T,N}} = nothing
    return closestingroups(closest, ass.groups, r)
end

# Visits each group in turn, one statement per group, so each group's loop is compiled for its own concrete element
# type. (Recursing on the tail of the Tuple instead stops being specialized for assemblies with several groups, and then
# allocates a new tail Tuple at every step.)
@generated function closestingroups(closest, groups::G, r::AbstractRay) where {G<:Tuple}
    visits = [:(closest = closestingroup(closest, groups[$i], r)) for i in 1:fieldcount(G)]
    return quote
        $(visits...)
        return closest
    end
end

# The closest of `closest` and the intersections with the elements of one group. The Vector's element type is concrete,
# so surfaceintersection needs no dynamic dispatch.
function closestingroup(closest::Union{Nothing,Intersection{T,N}}, group::AbstractVector, r::AbstractRay{T,N})::Union{Nothing,Intersection{T,N}} where {T<:Real,N}
    for element in group
        # nothing in an element can be hit before the ray enters its bounding box, so once an intersection has been found,
        # elements whose box starts beyond it can be skipped with a cheap box test instead of a full CSG evaluation
        if element isa CSGTree && closest !== nothing && entrydistance(BoundingBox(element), r) > α(closest)
            continue
        end
        closest = closer(closest, elementintersection(element, r))
    end
    return closest
end

function BoundingBox(la::LensAssembly{T}) where {T<:Real}
    es = elements(la)
    bbox = BoundingBox(es[1])
    for b in es[2:end]
        if b isa ParametricSurface{T} || b isa CSGTree{T}
            bbox = union(bbox, b)
        end
    end
    return bbox
end

axis(a::LensAssembly{T}) where {T<:Real} = a.axis
export axis

#############################################################################

"""
    LensTrace{T<:Real,N}

Contains an intersection point and the ray segment leading to it from within an optical trace.
The ray carries the path length, power, wavelength, number of intersections and source number, all of which are accessible directly on this class too.


Has the following accessor methods:
```julia
ray(a::LensTrace{T,N}) -> OpticalRay{T,N}
intersection(a::LensTrace{T,N}) -> Intersection{T,N}
power(a::LensTrace{T,N}) -> T
wavelength(a::LensTrace{T,N}) -> T
pathlength(a::LensTrace{T,N}) -> T
point(a::LensTrace{T,N}) -> SVector{N,T}
uv(a::LensTrace{T,N}) -> SVector{2,T}
sourcenum(a::LensTrace{T,N}) -> Int
nhits(a::LensTrace{T,N}) -> Int
```
"""
struct LensTrace{T<:Real,N}
    ray::OpticalRay{T,N}
    intersection::Intersection{T,N}
end
export LensTrace

ray(a::LensTrace{T,N}) where {T<:Real,N} = a.ray
intersection(a::LensTrace{T,N}) where {T<:Real,N} = a.intersection
power(a::LensTrace{T,N}) where {T<:Real,N} = power(ray(a))
wavelength(a::LensTrace{T,N}) where {T<:Real,N} = wavelength(ray(a))
pathlength(a::LensTrace{T,N}) where {T<:Real,N} = pathlength(ray(a))
point(a::LensTrace{T,N}) where {T<:Real,N} = point(intersection(a))
uv(a::LensTrace{T,N}) where {T<:Real,N} = uv(intersection(a))
sourcenum(a::LensTrace{T,N}) where {T<:Real,N} = sourcenum(ray(a))
nhits(a::LensTrace{T,N}) where {T<:Real,N} = nhits(ray(a))
export intersection

function Base.print(io::IO, a::LensTrace{T,N}) where {T,N}
    println(io, "Ray\n$(ray(a))")
    println(io, "Lensintersection\n$(intersection(a))")
end

struct NoPower end
const nopower = NoPower()

# isnopower(::NoPower) = true
# isnopower(::Any) = false

#############################################################################

"""
    trace(assembly::LensAssembly{T}, r::OpticalRay{T}, temperature::T = 20.0, pressure::T = 1.0; trackrays = nothing, test = false)

Returns the ray as it exits the assembly in the form of a [`LensTrace`](@ref) object if it hits any element in the assembly, otherwise `nothing`.
Recursive rays are offset by a small amount (`RAY_OFFSET`) to prevent it from immediately reintersecting the same lens element.

`trackrays` can be passed an empty vector to accumulate the `LensTrace` objects at each intersection of `ray` with a surface in the assembly.
"""
function trace(assembly::LensAssembly{T}, r::OpticalRay{T,N}, temperature::T=T(AGFFileReader.TEMP_REF), pressure::T=T(AGFFileReader.PRESSURE_REF); trackrays::Union{Nothing,Vector{LensTrace{T,N}}}=nothing, test::Bool=false, recursion::Int=0)::Union{Nothing,NoPower,LensTrace{T,N}} where {T<:Real,N}
    if power(r) < POWER_THRESHOLD || recursion > TRACE_RECURSION_LIMIT
        return nopower
    end

    intsct = closestintersection(assembly, r)
    if intsct === nothing
        return nothing
    else
        surfintsct = point(intsct)
        nml = normal(intsct)
        opticalinterface = interface(intsct)
        λ = wavelength(r)

        # The interface is one of five concrete types. Julia only splits unions of up to four types automatically, so test
        # for the common ones explicitly and each call compiles to a direct call.
        firstray = recursion == 0
        temp = if opticalinterface isa FresnelInterface{T}
            processintersection(opticalinterface, surfintsct, nml, r, temperature, pressure, test, firstray)
        elseif opticalinterface isa ParaxialInterface{T}
            processintersection(opticalinterface, surfintsct, nml, r, temperature, pressure, test, firstray)
        else
            processintersection(opticalinterface, surfintsct, nml, r, temperature, pressure, test, firstray)
        end

        if temp === nothing
            return nopower
        else
            raydirection, raypower, raypathlength = temp
            if trackrays !== nothing
                # we want the power that is hitting the surface, i.e. power on incident ray, but the path length at this intersection
                # i.e. with the path length for this intersection added and power modulated by absorption only
                push!(trackrays, LensTrace(OpticalRay(ray(r), raypower, λ, opl=raypathlength, nhits=nhits(r), sourcenum=sourcenum(r), sourcepower=sourcepower(r)), intsct))
            end
            offsetray = OpticalRay(surfintsct + RAY_OFFSET * raydirection, raydirection, raypower, λ, opl=raypathlength, nhits=nhits(r) + 1, sourcenum=sourcenum(r), sourcepower=sourcepower(r))
            res = trace(assembly, offsetray, temperature, pressure, trackrays=trackrays, test=test, recursion=recursion + 1)
            if res === nothing
                return LensTrace(OpticalRay(surfintsct, raydirection, raypower, λ, opl=raypathlength, nhits=nhits(r), sourcenum=sourcenum(r), sourcepower=sourcepower(r)), intsct)
            else
                return res
            end
        end
    end
end

# """approximate focal length of simple double convex lens. For plano convex set one of the radii to Inf"""
focallength(index::Float64, rf::Float64, rb::Float64) = 1 / ((index - 1) * (1 / rf - 1 / rb))

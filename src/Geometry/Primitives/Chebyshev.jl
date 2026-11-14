# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
Module to enclose [Chebyshev polynomial](https://en.wikipedia.org/wiki/Chebyshev_polynomials) specific functionality.
"""
module Chebyshev
using OpticSim: NaNsafeacos

"""
    T(n::Int, q::R, fast::Bool = true) -> R

Evaluate Chebyshev polynomial of the first kind ``T_n(q)``.

`fast` will use trigonometric definition, rather than the recursive definition which is much faster but slightly less precise.
"""
@inline function T(n::Int, q::R, fast::Bool=true)::R where {R<:Real}
    @assert n >= 0
    if n === 0
        return one(R)
    elseif n === 1
        return q
    elseif fast && n > 3
        if abs(q) < one(R)
            return cos(n * NaNsafeacos(q))
        elseif q >= one(R)
            return cosh(n * acosh(q))
        else
            return (-1)^n * cosh(n * acosh(-q))
        end
    else
        return 2q * T(n - 1, q, fast) - T(n - 2, q, fast)
    end
end

"""
    U(n::Int, q::R, fast::Bool = true) -> R

Evaluate Chebyshev polynomial of the second kind ``U_n(q)``.

`fast` will use trigonometric definition, rather than the recursive definition which is much faster but slightly less precise.
"""
@inline function U(n::Int, q::R, fast::Bool=true)::R where {R<:Real}
    @assert n >= 0
    if n === 0
        return one(R)
    elseif n === 1
        return 2q
    elseif abs(q) < one(R) && fast && q > 3
        # much faster but not stable at |q| = 1
        aq = NaNsafeacos(q)
        return sin((n + 1) * aq) / sin(aq)
    else
        return 2q * U(n - 1, q, fast) - U(n - 2, q, fast)
    end
end

"""
    dTdq(n::Int, q::R, fast::Bool = true) -> R

Evaluate derivative of Chebyshev polynomial of the first kind ``\\frac{dT_n}{dq}(q)``.

`fast` will use trigonometric definition, rather than the recursive definition which is much faster but slightly less precise.
"""
@inline function dTdq(n::Int, q::R, fast::Bool=true)::R where {R<:Real}
    @assert n >= 0
    if n === 0
        return zero(R)
    elseif n === 1
        return one(R)
    elseif fast && n > 4
        if abs(q) == one(R)
            return q^(n + 1) * n^2
        elseif abs(q) < one(R)
            return n * sin(n * acos(q)) / sqrt(1 - q^2)
        elseif q > one(R)
            return n * sinh(n * acosh(q)) / sqrt(q^2 - 1)
        else
            return -n * (-1)^n * sinh(n * acosh(-q)) / sqrt(q^2 - 1)
        end
    else
        return n * U(n - 1, q, fast)
    end
end

# dUdq(n::Int, q::R)::R where {R<:Real} = ((n + 1) * T(n + 1, q) - q * U(n, q)) / (q^2 - one(T))
# d2Tq2(n::Int, q::R)::R where {R<:Real} = n * ((n + 1) * T(n, q) - U(n, q)) / (q^2 - 1)

end # module Chebyshev

"""
    ChebyshevSurface{T,N,P,IJ} <: ParametricSurface{T,N}

Rectangular surface incorporating Chebyshev polynomials as well as radius and conic terms.
`T` is the datatype, `N` is the dimensionality, `P` is the number of Chebyshev terms and `IJ` is the tuple of their `(i, j)` indices.

The surface is centered at the origin and treated as being the cap of an infinite rectangular prism, thus creating a true half-space.
**Note that the surface is vertically offset so that the center (i.e., `(u,v) == (0,0)`) lies at 0 on the z-axis.**

```julia
ChebyshevSurface(halfsizeu, halfsizev, chebycoeff; radius = Inf, conic = 0)
```

`chebycoeff` is a vector containing tuples of the form `(i, j, v)` where `v` is the value of the coefficient ``c_{ij}``.

The sag is defined by the equation

```math
z(u,v) = \\frac{c(u^2 + v^2)^2}{1 + \\sqrt{1 - (1+k)c^2(u^2 + v^2)}} + \\sum_{i}^{P}\\sum_{j}^{Q}c_{ij}T_i(u)T_j(v)
```

where ``c = \\frac{1}{\\texttt{radius}}``, ``k = \\texttt{conic}`` and ``T_n`` is the nᵗʰ Chebyshev polynomial of the first kind.
"""
struct ChebyshevSurface{T,N,P,IJ} <: ParametricSurface{T,N}
    halfsizeu::T
    halfsizev::T
    curvature::T
    conic::T
    boundingprism::BoundingBox{T}
    chebycoeff::SVector{P,Tuple{Int,Int,T}}
    offset::T

    function ChebyshevSurface(halfsizeu::T, halfsizev::T, chebycoeff::Union{Nothing,Vector{Tuple{Int,Int,T}}}; radius::T=typemax(T), conic::T=zero(T)) where {T<:Real}
        @assert !isnan(halfsizeu) && !isnan(halfsizev) && !isnan(radius) && !isnan(conic)
        @assert halfsizeu > zero(T) && halfsizev > zero(T)
        @assert one(T) - (1 / radius)^2 * (conic + one(T)) * (halfsizeu^2 + halfsizev^2) > 0 "Invalid surface (conic/radius combination: $radius, $conic)"
        offset = zero(T)
        if chebycoeff === nothing
            P = 0
        else
            for (i, j, c) in chebycoeff
                @assert i >= 0 && j >= 0 "i and j must be non-negative"
                if i % 2 == 0 && j % 2 == 0
                    offset += c * (-1)^(i ÷ 2) * (-1)^(j ÷ 2)
                end
            end
            chebycoeff = filter(k -> abs(k[3]) > zero(T), chebycoeff)
            P = length(chebycoeff)
        end
        bounding_prism = BoundingBox(-halfsizeu, halfsizeu, -halfsizev, halfsizev, typemin(T), typemax(T))
        coeffs = SVector{P,Tuple{Int,Int,T}}(P === 0 ? [] : chebycoeff)
        IJ = Tuple((i, j) for (i, j, _) in coeffs) #(i, j) indices of the terms, used to generate the derivative code of `sag`
        return new{T,3,P,IJ}(halfsizeu, halfsizev, 1 / radius, conic, bounding_prism, coeffs, offset)
    end
end
export ChebyshevSurface

uvrange(::Type{<:ChebyshevSurface{T,N,P}}) where {T<:Real,N,P} = ((-one(T), one(T)), (-one(T), one(T)))

boundingobj(z::ChebyshevSurface{T}) where {T<:Real} = z.boundingprism
halfsizeu(z::ChebyshevSurface{T}) where {T<:Real} = z.halfsizeu
halfsizev(z::ChebyshevSurface{T}) where {T<:Real} = z.halfsizev

#=
The shape of a ChebyshevSurface is defined once, by `sag` below, in Cartesian coordinates. `point`, `partials`, `normal`
and `surfacecurvature` are derived from it (see SagCurvature.jl). `IJ` is the tuple of (i, j) indices of the Chebyshev
terms. It is a type parameter so that derivative code can be generated once per term structure. Coefficients, in the
order `sag` reads them: [curvature, conic, 1/halfsizeu, 1/halfsizev, offset, c₁, …, c_P]
=#
function sag(::Type{<:ChebyshevSurface{T,N,P,IJ}}, x, y, c) where {T,N,P,IJ}
    r2 = x^2 + y^2
    curvature, conic = c[1], c[2]
    h = curvature * r2 / (1 + sqrt(1 - (1 + conic) * curvature^2 * r2)) - c[5]
    u = x * c[3]
    v = y * c[4]
    for k in 1:P
        i, j = IJ[k]
        h += c[5+k] * chebyshevT(i, u) * chebyshevT(j, v)
    end
    return h
end

#Chebyshev polynomial of the first kind by its three term recurrence; branches only on the degree n, never on the value of q
function chebyshevT(n::Int, q)
    n == 0 && return one(q)
    tₙ₋₁, tₙ = one(q), q
    for _ in 2:n
        tₙ₋₁, tₙ = tₙ, 2q * tₙ - tₙ₋₁
    end
    return tₙ
end

sagcoefficients(s::ChebyshevSurface{T,N,P}) where {T,N,P} = SVector{5 + P,T}(s.curvature, s.conic, inv(s.halfsizeu), inv(s.halfsizev), s.offset, map(t -> t[3], s.chebycoeff)...)
sagdefined(s::ChebyshevSurface, r2) = conicdefined(s.curvature, s.conic, r2)

function uv(s::ChebyshevSurface{T,3,P}, p::SVector{3,T}) where {T<:Real,P}
    return SVector{2,T}(p[1] / s.halfsizeu, p[2] / s.halfsizev)
end

function onsurface(surf::ChebyshevSurface{T,3,P}, p::SVector{3,T}) where {T<:Real,P}
    u, v = uv(surf, p)
    if abs(u) > one(T) || abs(v) > one(T)
        return false
    else
        surfpoint = point(surf, u, v)
        return samepoint(p[3], surfpoint[3])
    end
end

function inside(surf::ChebyshevSurface{T,3,P}, p::SVector{3,T}) where {T<:Real,P}
    u, v = uv(surf, p)
    if abs(u) > one(T) || abs(v) > one(T)
        return false
    else
        surfpoint = point(surf, u, v)
        return p[3] < surfpoint[3]
    end
end

#########################################################################################################

# Assumes the ray has been transformed into the canonical coordinate frame which has the vertical axis passing through (0,0,0) and aligned with the z axis.
function surfaceintersection(surf::AcceleratedParametricSurface{T,3,<:ChebyshevSurface{T,3,P}}, r::AbstractRay{T,3}) where {T<:Real,P}
    bboxint = surfaceintersection(surf.surface.boundingprism, r)
    if isemptyinterval(bboxint)
        return EmptyInterval(T)
    else
        if doesintersect(surf.triangles_bbox, r) || inside(surf.triangles_bbox, origin(r))
            surfint = triangulatedintersection(surf, r)
            if !(surfint isa EmptyInterval{T})
                return intervalintersection(bboxint, surfint)
            end
        end
        # hasn't hit the surface
        if lower(bboxint) isa RayOrigin{T} && upper(bboxint) isa Infinity{T}
            if inside(surf.surface, origin(r))
                return Interval(RayOrigin(T), Infinity(T))
            else
                return EmptyInterval(T)
            end
            # otherwise check that the intersection is underneath the surface
        else
            p = point(closestintersection(bboxint, false))
            ρ, ϕ = uv(surf, p)
            surfpoint = point(surf.surface, ρ, ϕ)
            if p[3] < surfpoint[3]
                return bboxint # TODO!! UV (and interface) issues?
            else
                return EmptyInterval(T)
            end
        end
    end
end

function BoundingBox(surf::ChebyshevSurface{T,3,P}) where {T<:Real,P}
    xmin = -surf.halfsizeu
    xmax = surf.halfsizeu
    ymin = -surf.halfsizev
    ymax = surf.halfsizev
    # polynomials range between -1 and 1 so we have to sum the absolute value of every coefficient to get the theoretical max
    zmax = P > 0 ? sum(abs(c) for (_, _, c) in surf.chebycoeff) : zero(T)
    zmin = -zmax
    q = one(T) - (one(T) + surf.conic) * surf.curvature^2 * (surf.halfsizeu^2 + surf.halfsizev^2)
    if q < zero(T)
        throw(ErrorException("The surface is invalid, no bounding box can be constructed"))
    end
    hmax = surf.curvature * (surf.halfsizeu^2 + surf.halfsizev^2) / (one(T) + sqrt(q))
    if hmax > zero(T)
        zmax += hmax
    else
        zmin += hmax
    end
    return BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)
end
nsagcoefficients(::Type{<:ChebyshevSurface{T,N,P}}) where {T,N,P} = 5 + P

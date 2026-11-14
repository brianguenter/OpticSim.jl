# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
Module to enclose [Zernike polynomial](https://en.wikipedia.org/wiki/Zernike_polynomials) specific functionality.
"""
module Zernike

"""
    OSAtoNM(j::Int) -> Tuple{Int, Int}

Convert OSA zernike index `j` to `(N,M)` form according to formula `J = N * (N + 2) + M`.
"""
function OSAtoNM(j::Int)::Tuple{Int,Int}
    n = Int(ceil((-3 + sqrt(9 + 8j)) / 2))
    m = 2j - n * (n + 2)
    return (Int(n), Int(m))
end

"""
    NolltoNM(j::Int) -> Tuple{Int, Int}

Convert Noll zernike index `j` to `(N,M)` form.
"""
function NolltoNM(j::Int)
    n = Int(ceil((-3 + sqrt(1 + 8j)) / 2))
    jr = j - Int(n * (n + 1) / 2)
    if mod(n, 4) ∈ (0, 1)
        m1 = jr
        m2 = -(jr - 1)
        if iseven(n - m1)
            m = m1
        else
            m = m2
        end
    else # mod(n,4) ∈ (2,3)
        m1 = jr - 1
        m2 = -(jr)
        if iseven(n - m1)
            m = m1
        else
            m = m2
        end
    end
    return (Int(n), Int(m))
end

"""
    normalisation(::Type{T}, N::Int, M::Int) -> T

Normalisation coefficient for Zernike polynomial term ``Z_{n}^{m}``.
"""
@inline normalisation(::Type{T}, N::Int, M::Int) where {T<:Real} = T(sqrt((2 * (N + 1)) / (1 + (M == 0 ? 1 : 0))))

"""
    ζ(N::Int, M::Int, ρ::T, ϕ::T) -> Tuple{T,T}

Evaluate Zernike polynomial term ``Z_{n}^{m}(\\rho, \\phi)``.
"""
@inline function ζ(N::Int, M::Int, ρ::T, ϕ::T)::T where {T<:Real}
    aM = abs(M)
    if M < 0
        return normalisation(T, N, M) * R(N, aM, ρ) * sin(aM * ϕ)
    else
        return normalisation(T, N, M) * R(N, aM, ρ) * cos(aM * ϕ)
    end
end

"""
    δζ(N::Int, M::Int, ρ::T, ϕ::T) -> Tuple{T,T}

Evaluate partial derivatives of Zernike polynomial term ``Z_{n}^{m}(\\rho, \\phi)``.
"""
@inline function δζ(N::Int, M::Int, ρ::T, ϕ::T)::Tuple{T,T} where {T<:Real}
    n = normalisation(T, N, M)
    aM = abs(M)
    RNM = R(N, aM, ρ)
    ρ2 = ρ^2
    δRNM = ((2 * N * aM * (ρ2 - 1) + (N - aM) * (aM + N * (2 * ρ2 - 1))) * RNM - (N + aM) * (N - aM) * R(N - 2, aM, ρ)) / (2 * N * ρ * (ρ2 - 1))
    Mϕ = aM * ϕ
    if M < 0
        δρ = n * δRNM * sin(Mϕ)
        δϕ = n * RNM * aM * cos(Mϕ)
    else
        δρ = n * δRNM * cos(Mϕ)
        δϕ = n * RNM * aM * -sin(Mϕ)
    end
    return δρ, δϕ
end

"""
    R(N::Int, M::Int, ρ::T) -> T

Evaluate radial polynomial ``R_{n}^{m}(\\rho)``.
"""
@inline function R(N::Int, M::Int, ρ::T)::T where {T<:Real}
    if (N - M) % 2 === 1
        return zero(T)
    end
    total = zero(T)
    @simd for k in 0:((N-M)÷2)
        total += ((-1)^k * factorial(N - k)) / (factorial(k) * factorial((N + M) ÷ 2 - k) * factorial((N - M) ÷ 2 - k)) * ρ^(N - 2 * k)
    end
    return total
end

end # module Zernike

#########################################################################################################

"""
Either `ZernikeIndexingOSA` or `ZernikeIndexingNoll`, see [Zernike polynomials wikipedia entry](https://en.wikipedia.org/wiki/Zernike_polynomials) for details.
"""
@enum ZernikeIndexType ZernikeIndexingOSA ZernikeIndexingNoll
export ZernikeIndexType, ZernikeIndexingOSA, ZernikeIndexingNoll

"""
    ZernikeSurface{T,N,P,Q,M,Z} <: ParametricSurface{T,N}

Surface incorporating the Zernike polynomials - radius, conic and aspherics are defined relative to absolute semi-diameter, Zernike terms are normalized according to the `normradius` parameter.
`T` is the datatype, `N` is the dimensionality, `P` is the number of Zernike terms, `Q` is the number of aspheric terms, `M` is the Aspheric Type and `Z` is the tuple of `(n, m)` indices of the Zernike terms.

The surface is centered at the origin and treated as being the cap of an infinite cylinder, thus creating a true half-space.
Outside of `0 <= ρ <= 1` the height of the surface is not necessarily well defined, so NaN may be returned.

For convenience the input `zcoeff` can be indexed using either OSA or Noll convention, indicated using the `indexing` argument as either `ZernikeIndexingOSA` or `ZernikeIndexingNoll`.

```julia
ZernikeSurface(semidiameter, radius = Inf, conic = 0, zcoeff = nothing, aspherics = nothing, normradius = semidiameter, indexing = ZernikeIndexingOSA)
```

`zcoeff` and `aspherics` should be vectors containing tuples of the form `(i, v)` where `i` is either the index of the Zernike term
for the corresponding `indexing`, or the polynomial power of the aspheric term (may be even or odd) and 
`v` is the corresponding coefficient ``A_i`` or ``\\alpha_i`` respectively..  `M` will be determined from the terms entered to optimize
the evaluation of the aspheric polynomial.

The sag is defined by the equation

```math
z(r,\\phi) = \\frac{cr^2}{1 + \\sqrt{1 - (1+k)c^2r^2}} + \\sum_{i}^{Q}\\alpha_ir^{2i} + \\sum_{i}^PA_iZ_i(\\rho, \\phi)
```

where ``\\rho = \\frac{r}{\\texttt{normradius}}``, ``c = \\frac{1}{\\texttt{radius}}``, ``k = \\texttt{conic}`` and ``Z_n`` is the nᵗʰ Zernike polynomial.
"""
struct ZernikeSurface{T,N,P,Q,M,Z} <: ParametricSurface{T,N}
    asp::AsphericSurface{T,N,Q,M}
    coeffs::SVector{P,Tuple{Int,Int,T}}
    boundingcylinder::Cylinder{T,N}

    function ZernikeSurface(semidiameter::T; radius::T=typemax(T), conic::T=zero(T), zcoeff::Union{Nothing,Vector{Tuple{Int,T}}}=nothing, aspherics::Union{Nothing,Vector{Tuple{Int,T}}}=nothing, normradius::T=semidiameter, indexing::ZernikeIndexType=ZernikeIndexingOSA) where {T<:Real}
        asp = AsphericSurface(semidiameter; radius, conic, aspherics, normradius)
        Q = length(asp.aspherics) #this is not the same as the aspherics variable passed to the function!

        zcs = []
        if zcoeff !== nothing
            for (i, k) in zcoeff
                if abs(k) > zero(T)
                    if indexing === ZernikeIndexingOSA
                        R, S = Zernike.OSAtoNM(i)
                    else
                        R, S = Zernike.NolltoNM(i)
                    end
                    push!(zcs, (R, S, k))
                end
            end
        end
        P = length(zcs)
        M = asphericType(asp)
        Z = Tuple((R, S) for (R, S, _) in zcs) #(n, m) indices of the terms, used to generate the derivative code of `sag`
        new{T,3,P,Q,M,Z}(asp::AsphericSurface{T,3,Q,M}, SVector{P,Tuple{Int,Int,T}}(zcs), Cylinder(semidiameter, interface=opaqueinterface(T))) # TODO!! incorrect interface on cylinder
    end

end
export ZernikeSurface

uvrange(::Type{<:ZernikeSurface{T,N,P,Q,M}}) where {T<:Real,N,P,Q,M} = ((zero(T), one(T)), (-T(π), T(π))) # ρ and ϕ

semidiameter(z::ZernikeSurface{T}) where {T<:Real} = z.asp.semidiameter
halfsizeu(z::ZernikeSurface{T}) where {T<:Real} = semidiameter(z)
halfsizev(z::ZernikeSurface{T}) where {T<:Real} = semidiameter(z)

boundingobj(z::ZernikeSurface{T}) where {T<:Real} = z.boundingcylinder

#=
The shape of a ZernikeSurface is defined once, by `sag` below, in Cartesian coordinates. `point`, `partials`, `normal`
and `surfacecurvature` are derived from it (see SagCurvature.jl). `Z` is the tuple of (n, m) indices of the Zernike terms.
It is a type parameter so that derivative code can be generated once per term structure.
=#
function sag(::Type{<:ZernikeSurface{T,N,P,Q,M,Z}}, x, y, c) where {T,N,P,Q,M,Z}
    h = sag(AsphericSurface{T,N,Q,M}, x, y, c) #reads c[1:2+Q]
    invnormradius = c[3+Q]
    ux = x * invnormradius
    uy = y * invnormradius
    u2 = ux^2 + uy^2
    for i in 1:P
        n, m = Z[i]
        h += c[3+Q+i] * zernikecartesian(n, m, ux, uy, u2)
    end
    return h
end

"""
    zernikecartesian(n, m, ux, uy, u2)

Normalized Zernike term ``Z_n^m`` at normalized Cartesian coordinates `(ux, uy)`, `u2 = ux² + uy²`. Written as a polynomial
in `ux`, `uy` (using ``u^{|m|}\\cos(|m|\\phi) = \\Re((u_x + i u_y)^{|m|})``) so it is smooth at the origin and has no branches
on the coordinate values.
"""
function zernikecartesian(n::Int, m::Int, ux, uy, u2)
    am = abs(m)
    isodd(n - am) && return zero(u2)
    re, im = one(ux), zero(ux) #(ux + i uy)^|m|
    for _ in 1:am
        re, im = re * ux - im * uy, re * uy + im * ux
    end
    angular = m < 0 ? im : re
    radial = zero(u2)
    for k in 0:((n-am)÷2)
        radial += ((-1)^k * factorial(n - k)) / (factorial(k) * factorial((n + am) ÷ 2 - k) * factorial((n - am) ÷ 2 - k)) * u2^((n - am - 2k) ÷ 2)
    end
    return sqrt((2 * (n + 1)) / (1 + (m == 0 ? 1 : 0))) * radial * angular
end

sagcoefficients(z::ZernikeSurface{T,N,P,Q}) where {T,N,P,Q} = SVector{3 + Q + P,T}(sagcoefficients(z.asp)..., inv(z.asp.normradius), map(t -> t[3], z.coeffs)...)
sagsmooth(::Type{<:ZernikeSurface{T,N,P,Q,M}}) where {T,N,P,Q,M} = sagsmooth(AsphericSurface{T,N,Q,M})
sagdefined(z::ZernikeSurface, r2) = sagdefined(z.asp, r2)

function uv(z::ZernikeSurface{T,3,P,Q,M}, p::SVector{3,T}) where {T<:Real,P,Q,M}
    # avoid divide by zero for ForwardDiff
    ϕ = NaNsafeatan(p[2], p[1])
    if p[1] == zero(T) && p[2] == zero(T)
        ρ = zero(T)
    else
        ρ = sqrt(p[1]^2 + p[2]^2) / semidiameter(z)
    end
    return SVector{2,T}(ρ, ϕ)
end

function onsurface(surf::ZernikeSurface{T,3,P,Q,M}, p::SVector{3,T}) where {T<:Real,P,Q,M}
    ρ, ϕ = uv(surf, p)
    if ρ > one(T)
        return false
    else
        surfpoint = point(surf, ρ, ϕ)
        return samepoint(p[3], surfpoint[3])
    end
end

function inside(surf::ZernikeSurface{T,3,P,Q,M}, p::SVector{3,T}) where {T<:Real,P,Q,M}
    ρ, ϕ = uv(surf, p)
    if ρ > one(T)
        return false
    else
        surfpoint = point(surf, ρ, ϕ)
        return p[3] < surfpoint[3]
    end
end

#########################################################################################################

# Assumes the ray has been transformed into the canonical coordinate frame which has the vertical axis passing through (0,0,0) and aligned with the z axis.
# `<:` matters: ZernikeSurface{T,3,P,Q,M} leaves the type parameter Z free, and without `<:` this method would never match a concrete surface
function surfaceintersection(surf::AcceleratedParametricSurface{T,3,<:ZernikeSurface{T,3,P,Q,M}}, r::AbstractRay{T,3}) where {T<:Real,P,Q,M}
    cylint = surfaceintersection(surf.surface.boundingcylinder, r)
    if isemptyinterval(cylint)
        return EmptyInterval(T)
    else
        if doesintersect(surf.triangles_bbox, r) || inside(surf.triangles_bbox, origin(r))
            surfint = triangulatedintersection(surf, r)
            if !(surfint isa EmptyInterval{T})
                return intervalintersection(cylint, surfint)
            end
        end
        # hasn't hit the surface
        if lower(cylint) isa RayOrigin{T} && upper(cylint) isa Infinity{T}
            if inside(surf.surface, origin(r))
                return Interval(RayOrigin(T), Infinity(T))
            else
                return EmptyInterval(T)
            end
            # otherwise check that the intersection is underneath the surface
        else
            p = point(closestintersection(cylint, false))
            ρ, ϕ = uv(surf, p)
            surfpoint = point(surf.surface, ρ, ϕ)
            if p[3] < surfpoint[3]
                return cylint # TODO!! UV (and interface) issues?
            else
                return EmptyInterval(T)
            end
        end
    end
end

function AcceleratedParametricSurface(surf::T, numsamples::Int=17; interface::NullOrFresnel{S}=NullInterface(S)) where {S<:Real,N,T<:ZernikeSurface{S,N}}
    # Zernike uses ρ, ϕ uv space so need to modify extension of triangulation
    a = AcceleratedParametricSurface(surf, triangulate(surf, numsamples, true, false, true, false), interface=interface)
    emptytrianglepool!(S)
    return a
end

function BoundingBox(surf::ZernikeSurface{T,3,P,Q,M}) where {T<:Real,P,Q,M}
    bb = BoundingBox(surf.asp)
    # zernike terms have condition than |Zᵢ| <= 1
    # so this gives us a (loose) bounding box
    ak = P > 0 ? sum(abs.(Zernike.normalisation(T, n, m) * k for (n, m, k) in surf.coeffs)) : zero(T)

    return BoundingBox(bb.xmin, bb.xmax, bb.ymin, bb.ymax, bb.zmin - ak, bb.zmax + ak) #could just return bb, but this way is safer
end
nsagcoefficients(::Type{<:ZernikeSurface{T,N,P,Q}}) where {T,N,P,Q} = 3 + Q + P

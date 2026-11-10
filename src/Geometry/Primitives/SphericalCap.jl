# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    SphericalCap{T} <: ParametricSurface{T}

Spherical cap surface, creates a half-space which is essentially the subtraction of a sphere from an infinite plane.
Only the spherical cap itself is visualized, not the plane.
The positive normal side is outside the surface.

**Can be used as a detector in [`AbstractOpticalSystem`](@ref)s.**

```julia
SphericalCap(radius::T, ϕmax::T, [surfacenormal::SVector{3,T}, centrepoint::SVector{3,T}]; interface::NullOrFresnel{T} = nullinterface(T))
```

The minimal case returns a spherical cap centered at the origin with `surfacenormal = [0, 0, 1]`.
"""
struct SphericalCap{T} <: ParametricSurface{T,3}
    radius::T
    ϕmax::T
    zmax::T
    centrenormal::SVector{3,T}
    centrepoint::SVector{3,T}
    interface::NullOrFresnel{T}
    uvec::SVector{3,T}
    vvec::SVector{3,T}

    function SphericalCap(radius::T, ϕmax::T; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
        @assert !isnan(radius)
        @assert radius > zero(T) && zero(T) < ϕmax < T(π)
        zmax = radius * (one(T) + sin(-(π / 2 - ϕmax)))
        new{T}(radius, ϕmax, zmax, SVector{3,T}(0.0, 0.0, 1.0), SVector{3,T}(0.0, 0.0, 0.0), interface, SVector{3,T}(1.0, 0.0, 0.0), SVector{3,T}(0.0, 1.0, 0.0))
    end

    function SphericalCap(radius::T, ϕmax::T, centrenormal::SVector{3,T}, centrepoint::SVector{3,T}; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
        @assert radius > zero(T) && zero(T) < ϕmax < T(π)
        n̂ = normalize(centrenormal)
        rotationvec = SVector{3,T}(0.0, 1.0, 0.0)
        if abs(dot(rotationvec, n̂)) == one(T)
            rotationvec = SVector{3,T}(1.0, 0.0, 0.0)
        end
        uvec = normalize(cross(normalize(rotationvec), n̂))
        vvec = normalize(cross(n̂, uvec))
        zmax = radius * (one(T) + sin(-(π / 2 - ϕmax)))
        new{T}(radius, ϕmax, zmax, n̂, centrepoint, interface, uvec, vvec)
    end
end
export SphericalCap

Base.show(io::IO, a::SphericalCap{T}) where {T<:Real} = print(io, "SphericalCap{$T}($(a.radius), $(a.centrepoint), $(a.centrenormal), $(a.ϕmax), $(interface(a)))")

interface(a::SphericalCap{T}) where {T<:Real} = a.interface
radius(a::SphericalCap{T}) where {T<:Real} = a.radius
centroid(r::SphericalCap{T}) where {T<:Real} = r.centrepoint

uvrange(::Type{SphericalCap{T}}) where {T<:Real} = ((-T(π), T(π)), (zero(T), one(T))) # θ and ρ

function normal(r::SphericalCap{T}, θ::T, ρ::T) where {T<:Real}
    ϕ = π / 2 - ρ * r.ϕmax
    return -normalize(cos(-ϕ) * (cos(-θ) * r.uvec + sin(-θ) * r.vvec) + sin(-ϕ) * r.centrenormal)
end

function point(r::SphericalCap{T}, θ::T, ρ::T) where {T<:Real}
    ϕ = π / 2 - ρ * r.ϕmax
    return centroid(r) + radius(r) * (cos(-ϕ) * (cos(-θ) * r.uvec + sin(-θ) * r.vvec) + sin(-ϕ) * r.centrenormal + r.centrenormal)
end

uv(r::SphericalCap{T}, x::T, y::T, z::T) where {T<:Real} = uv(r, SVector{3,T}(x, y, z))
function uv(r::SphericalCap{T}, p::SVector{3,T}) where {T<:Real}
    prel = p - centroid(r) - radius(r) * r.centrenormal
    v = dot(prel, r.vvec) / radius(r)
    u = dot(prel, r.uvec) / radius(r)
    θ = -NaNsafeatan(v, u)
    z = clamp(dot(prel, -r.centrenormal) / radius(r), -one(T), one(T))
    ϕ = NaNsafeasin(z)
    ρ = (π / 2 - ϕ) / r.ϕmax
    return θ, ρ
end

# The solid is the half space below the plane through the rim of the cap, minus the ball the cap is part of. The cap is
# the part of the ball's surface that lies in that half space.
function inside(r::SphericalCap{T}, p::SVector{3,T}) where {T<:Real}
    inball = norm(p - centroid(r) - radius(r) * r.centrenormal) < radius(r)
    return !inball && dot(r.centrenormal, p - centroid(r)) < r.zmax
end

function onsurface(r::SphericalCap{T}, p::SVector{3,T}) where {T<:Real}
    prel = (p - centroid(r) - radius(r) * r.centrenormal) / radius(r)
    v = dot(prel, r.vvec)
    u = dot(prel, r.uvec)
    if (u == zero(T) && v == zero(T)) || sqrt(u^2 + v^2) < sin(r.ϕmax)
        return norm(prel) == one(T) && dot(r.centrenormal, p - centroid(r)) < r.zmax
    else
        false
    end
end

function uvtopix(::SphericalCap{T}, uv::SVector{2,T}, imsize::Tuple{Int,Int}) where {T<:Real}
    θ, ρ = uv
    h, w = imsize
    u = (cos(θ) * ρ + one(T)) / 2
    v = (sin(θ) * ρ + one(T)) / 2
    pixu = Int(floor((w - 1) * u)) + 1
    pixv = Int(floor((h - 1) * v)) + 1
    return pixu, pixv
end

function surfaceintersection(sph::SphericalCap{T}, r::AbstractRay{T,3}) where {T<:Real}
    rad = radius(sph)
    orel = origin(r) - centroid(sph) - rad * sph.centrenormal

    ox, oy, oz = orel
    dx, dy, dz = direction(r)

    a = dx^2 + dy^2 + dz^2
    b = 2 * (ox * dx + oy * dy + oz * dz)
    c = (ox^2 + oy^2 + oz^2) - rad^2

    # the solid is the half space below the rim plane minus the ball. Crossings of the rim plane are not part of the cap,
    # so they get a null interface: they bound the solid for CSG but rays don't interact with them.
    halfspace = surfaceintersection(Plane(sph.centrenormal, centroid(sph) + sph.zmax * sph.centrenormal), r)
    halfspace isa EmptyInterval{T} && return halfspace

    temp = quadraticroots(a, b, c)
    if temp === nothing || isapprox(temp[1], temp[2], rtol = 1e-12, atol = 2 * eps(T))
        return halfspace # the ray misses the ball or just touches it
    end
    t1, t2 = minmax(temp...)
    t2 <= zero(T) && return halfspace # the ball is behind the ray

    # the ray is inside the ball between t1 and t2. The normals point into the ball, out of the solid.
    function ballintersection(t)
        pt = point(r, t)
        θ, ρ = uv(sph, pt)
        return Intersection(t, pt, normal(sph, θ, ρ), θ, ρ, interface(sph))
    end
    exitball = positivehalfspace(ballintersection(t2))
    outsideball = t1 > zero(T) ? DisjointUnion(rayorigininterval(ballintersection(t1)), exitball) : exitball
    return intervalintersection(halfspace, outsideball)
end

BoundingBox(a::SphericalCap{T}) where {T<:Real} = BoundingBox(Plane(a.centrenormal, a.centrepoint + a.radius * a.centrenormal))

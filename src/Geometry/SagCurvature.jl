# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

#########################################################################################################
# Sag surfaces: z = sag(x, y) in the surface's local frame.
#
# Each sag surface type `S` defines its shape once, as a branch free function of Cartesian x, y:
#
#     sag(::Type{S}, x, y, c)         height of the surface at (x, y)
#     sagcoefficients(s::S)           SVector of the numbers `c` that `sag` reads (curvature, conic,
#                                     polynomial coefficients, ...)
#
# `sag` must accept any number type: Float64, ForwardDiff.Dual (optimization) and FastDifferentiation.Node
# (symbolic differentiation). It must not branch on the values of its arguments; domain checks (for
# example the conic square root becoming imaginary) belong in the callers. The structure of the sum
# (number of terms, term indices) must be determined by the type `S` alone, so the derivative code below
# is generated once per surface type and order, with the coefficients as runtime inputs.
#
# `point`, `partials`, `normal` and `surfacecurvature` of every sag surface are built from `sag`, so the
# formula for each surface exists in exactly one place.
#########################################################################################################

"""
    sag(::Type{S}, x, y, c)

Height of a sag surface of type `S` at local coordinates `(x, y)`, with coefficients `c` as returned by
[`sagcoefficients`](@ref). Branch free and generic in the number type.
"""
function sag end

"""
    sagcoefficients(surface)

`SVector` of the coefficients that [`sag`](@ref) reads for `surface`.
"""
function sagcoefficients end

"""
    sagsmooth(::Type{S}) -> Bool

`true` (the default) if the `sag` expression of `S` can be differentiated at `x = y = 0`. Odd aspheric
terms contain `sqrt(x² + y²)`, whose symbolic derivative is 0/0 on the axis even when the surface is
smooth there. For those types the derivatives on the axis are taken a vanishing distance out along the
requested direction, which gives the one-sided radial derivative.
"""
sagsmooth(::Type) = true

#distance used to step off the axis for surfaces that are not smooth there. Squared it is still a normal Float64.
const ONAXIS_OFFSET = 1e-150

#Expression that evaluates `outputs(h, hx, hy, hxx, hxy, hyy)` of sag(S, x, y, c), as the body of a generated function with arguments x, y, c.
function _sagderivative_body(S::Type, K::Int, order::Int)
    vars = FastDifferentiation.make_variables(:sag_xy, 2)
    cvars = FastDifferentiation.make_variables(:sag_c, K)
    h = FastDifferentiation.Node(sag(S, vars[1], vars[2], cvars))
    hx = FastDifferentiation.derivative(h, vars[1])
    hy = FastDifferentiation.derivative(h, vars[2])
    outputs = if order == 1
        SVector(hx, hy)
    else
        SVector(h, hx, hy, FastDifferentiation.derivative(hx, vars[1]), FastDifferentiation.derivative(hx, vars[2]), FastDifferentiation.derivative(hy, vars[2]))
    end
    return _inline_derivative_code(outputs, vars, cvars, :(SVector(x, y)), :c)
end

#Code that evaluates the FastDifferentiation expressions `outputs` of the variables `vars` and `cvars`, given expressions for the values of `vars` and `cvars`; used as the body of a generated function.
function _inline_derivative_code(outputs, vars, cvars, varsvalue, cvarsvalue)
    lambda = FastDifferentiation.make_Expr(outputs, vars, cvars)
    #make_Expr returns `(input_variables1, input_variables2) -> begin <input size check>; @inbounds <body> end`. Generated functions can not contain closures, so bind the input names and splice in the body.
    #The input sizes are fixed by the method signature, so the size check is dropped: it loops over inputs of different types, which is not type stable and allocates.
    argnames = lambda.args[1].args
    body = filter(lambda.args[2].args) do arg
        arg isa LineNumberNode || (arg isa Expr && arg.head === :macrocall && arg.args[1] === Symbol("@inbounds"))
    end
    @assert any(arg -> arg isa Expr, body) "unexpected code layout from FastDifferentiation.make_Expr"
    return quote
        $(argnames[1]) = $varsvalue
        $(argnames[2]) = $cvarsvalue
        $(body...)
    end
end

"""
    saggradient(::Type{S}, x, y, c) -> SVector(∂h/∂x, ∂h/∂y)

Gradient of the sag of surface type `S`. Compiled once per surface type from [`sag`](@ref).
"""
@generated function saggradient(::Type{S}, x, y, c::SVector{K}) where {S,K}
    return _sagderivative_body(S, K, 1)
end

"""
    saghessian(::Type{S}, x, y, c) -> SVector(h, hx, hy, hxx, hxy, hyy)

Height, gradient and second derivatives of the sag of surface type `S`. Compiled once per surface type
from [`sag`](@ref).
"""
@generated function saghessian(::Type{S}, x, y, c::SVector{K}) where {S,K}
    return _sagderivative_body(S, K, 2)
end

"""
    sagkernel(surface, x, y) -> (K, c)

The kernel type `K` and coefficients `c` for evaluating the sag of `surface` near local coordinates `(x, y)`:
`sag(K, x, y, c)`. For most surfaces this is `(typeof(surface), sagcoefficients(surface))`. Surfaces whose sag is piecewise,
like a [`GridSagSurface`](@ref), return the kernel and coefficients of the piece that contains `(x, y)`. `K` must depend only
on the type of `surface`, never on `(x, y)` or on field values, so that callers stay type stable.
"""
@inline sagkernel(surf, x, y) = (typeof(surf), sagcoefficients(surf))

"""
    sagdefined(surface, r2) -> Bool

`false` where the sag of `surface` is undefined at squared radius `r2` (for example where the conic square root becomes
imaginary). `point`, `partials` and `normal` return NaN there.
"""
function sagdefined end

#coordinates at which to evaluate derivatives for local point (x, y); steps off the axis along (cosϕ, sinϕ) when the sag is not differentiable there
@inline function _derivative_point(::Type{K}, x::T, y::T, cosϕ, sinϕ) where {K,T}
    if !sagsmooth(K) && iszero(x) && iszero(y)
        return T(ONAXIS_OFFSET) * cosϕ, T(ONAXIS_OFFSET) * sinϕ
    end
    return x, y
end

"""
    sagpartials_polar(surface, x, y, cosϕ, sinϕ, scale) -> (∂P/∂ρ, ∂P/∂ϕ)

Partial derivatives of P(ρ, ϕ) = (x, y, sag(x, y)) with x = scale·ρ·cosϕ and y = scale·ρ·sinϕ, for sag
surfaces that are parameterized by normalized polar coordinates.
"""
@inline function sagpartials_polar(surf, x::T, y::T, cosϕ::T, sinϕ::T, scale::T) where {T}
    K, c = sagkernel(surf, x, y)
    dx, dy = _derivative_point(K, x, y, cosϕ, sinϕ)
    hx, hy = saggradient(K, dx, dy, c)
    return SVector{3,T}(scale * cosϕ, scale * sinϕ, scale * (cosϕ * hx + sinϕ * hy)), SVector{3,T}(-y, x, -y * hx + x * hy)
end

"""
    sagnormal(surface, x, y) -> SVector{3}

Unit normal (−∂h/∂x, −∂h/∂y, 1)/‖·‖ of a sag surface at local coordinates `(x, y)`. On the axis of a surface
that is not differentiable there it returns (0, 0, 1).
"""
@inline function sagnormal(surf, x::T, y::T) where {T}
    K, c = sagkernel(surf, x, y)
    if !sagsmooth(K) && iszero(x) && iszero(y)
        return SVector{3,T}(0, 0, 1)
    end
    hx, hy = saggradient(K, x, y, c)
    return normalize(SVector{3,T}(-hx, -hy, one(T)))
end

"""
    curvature_from_sag(derivs) -> (normal, d1, d2, C)

Normal, orthonormal tangent basis and curvature matrix (see [`surfacecurvature`](@ref)) of the surface
(x, y, h(x, y)) from `derivs = (h, hx, hy, hxx, hxy, hyy)`. The normal is (−hx, −hy, 1)/‖·‖, pointing
toward +z, and `d1` is along ∂/∂x.
"""
function curvature_from_sag(derivs::SVector{6,T}) where {T}
    _, hx, hy, hxx, hxy, hyy = derivs
    Su = SVector{3,T}(1, 0, hx)
    Sv = SVector{3,T}(0, 1, hy)
    D = sqrt(1 + hx^2 + hy^2)
    n = SVector{3,T}(-hx / D, -hy / D, 1 / D)
    d1 = normalize(Su)
    d2 = cross(n, d1)
    # second fundamental form in the (x, y) parameterization: n⋅S_xx = hxx/D etc.
    II = SMatrix{2,2,T,4}(hxx, hxy, hxy, hyy) / D
    # coordinates of d1, d2 in the (Su, Sv) basis: solve [Su Sv] g = d using the first fundamental form
    I1 = SMatrix{2,2,T,4}(dot(Su, Su), dot(Sv, Su), dot(Su, Sv), dot(Sv, Sv))
    G = inv(I1) * SMatrix{2,2,T,4}(dot(Su, d1), dot(Sv, d1), dot(Su, d2), dot(Sv, d2))
    C = transpose(G) * II * G
    C = (C + transpose(C)) / 2
    return (normal=n, d1=d1, d2=d2, C=SMatrix{2,2,T,4}(C))
end

"""
    sagsurfacecurvature(surface, x, y) -> (normal, d1, d2, C)

[`surfacecurvature`](@ref) of a sag surface at local coordinates `(x, y)`.
"""
function sagsurfacecurvature(surf, x::T, y::T) where {T}
    K, c = sagkernel(surf, x, y)
    dx, dy = _derivative_point(K, x, y, one(T), zero(T))
    return curvature_from_sag(saghessian(K, dx, dy, c))
end

#########################################################################################################
# Sag surfaces parameterized by normalized polar coordinates (ρ, ϕ): the point is at radius r = ρ ⋅ semidiameter.

const PolarSagSurface{T} = Union{AsphericSurface{T,3},ZernikeSurface{T,3},QTypeSurface{T,3},GridSagSurface{T,3,<:ZernikeSurface}}

function point(z::PolarSagSurface{T}, ρ::T, ϕ::T)::SVector{3,T} where {T<:Real}
    r = ρ * semidiameter(z)
    if !sagdefined(z, r^2)
        return SVector{3,T}(NaN, NaN, NaN)
    end
    sinϕ, cosϕ = sincos(ϕ)
    x = r * cosϕ
    y = r * sinϕ
    K, c = sagkernel(z, x, y)
    return SVector{3,T}(x, y, sag(K, x, y, c))
end

function partials(z::PolarSagSurface{T}, ρ::T, ϕ::T)::Tuple{SVector{3,T},SVector{3,T}} where {T<:Real}
    rad = semidiameter(z)
    r = ρ * rad
    if !sagdefined(z, r^2)
        return SVector{3,T}(NaN, NaN, NaN), SVector{3,T}(NaN, NaN, NaN)
    end
    sinϕ, cosϕ = sincos(ϕ)
    return sagpartials_polar(z, r * cosϕ, r * sinϕ, cosϕ, sinϕ, rad)
end

function normal(z::PolarSagSurface{T}, ρ::T, ϕ::T)::SVector{3,T} where {T<:Real}
    r = ρ * semidiameter(z)
    if !sagdefined(z, r^2)
        return SVector{3,T}(NaN, NaN, NaN)
    end
    sinϕ, cosϕ = sincos(ϕ)
    return sagnormal(z, r * cosϕ, r * sinϕ)
end

#########################################################################################################
# Sag surfaces parameterized by normalized rectangular coordinates (u, v): the point is at x = u ⋅ halfsizeu, y = v ⋅ halfsizev.

const RectangularSagSurface{T} = Union{ChebyshevSurface{T,3},GridSagSurface{T,3,<:ChebyshevSurface}}

function point(s::RectangularSagSurface{T}, u::T, v::T)::SVector{3,T} where {T<:Real}
    x = u * halfsizeu(s)
    y = v * halfsizev(s)
    if !sagdefined(s, x^2 + y^2)
        return SVector{3,T}(NaN, NaN, NaN)
    end
    K, c = sagkernel(s, x, y)
    return SVector{3,T}(x, y, sag(K, x, y, c))
end

function partials(s::RectangularSagSurface{T}, u::T, v::T)::Tuple{SVector{3,T},SVector{3,T}} where {T<:Real}
    hu = halfsizeu(s)
    hv = halfsizev(s)
    x = u * hu
    y = v * hv
    if !sagdefined(s, x^2 + y^2)
        return SVector{3,T}(NaN, NaN, NaN), SVector{3,T}(NaN, NaN, NaN)
    end
    K, c = sagkernel(s, x, y)
    hx, hy = saggradient(K, x, y, c)
    return SVector{3,T}(hu, 0, hu * hx), SVector{3,T}(0, hv, hv * hy)
end

function normal(s::RectangularSagSurface{T}, u::T, v::T)::SVector{3,T} where {T<:Real}
    x = u * halfsizeu(s)
    y = v * halfsizev(s)
    if !sagdefined(s, x^2 + y^2)
        return SVector{3,T}(NaN, NaN, NaN)
    end
    return sagnormal(s, x, y)
end

#########################################################################################################

surfacecurvature(s::Union{PolarSagSurface{T},RectangularSagSurface{T}}, p::SVector{3,T}) where {T<:Real} = sagsurfacecurvature(s, p[1], p[2])

#########################################################################################################
# Implicit surfaces: f(p) = 0 on the surface, with ∇f pointing the same way as the surface's `normal`.
#
#     implicit(::Type{S}, p, c)       value of f at the 3D local point p, branch free and generic in the number type
#     implicitcoefficients(s::S)      SVector of the numbers `c` that `implicit` reads
#
# Along a unit tangent t the normal curvature is n⋅∂²S/∂s² = -tᵀ(∇∇f)t/|∇f|, which follows from differentiating f(S(s)) = 0
# twice. As for sag surfaces the derivative code is generated once per surface type from `implicit`.

"""
    implicit(::Type{S}, p, c)

Implicit function of surface type `S`: zero on the surface, gradient pointing along the surface's [`normal`](@ref).
"""
function implicit end

"""
    implicitcoefficients(surface)

`SVector` of the coefficients that [`implicit`](@ref) reads for `surface`.
"""
function implicitcoefficients end

"""
    implicitderivatives(::Type{S}, p, c) -> SVector(fx, fy, fz, fxx, fxy, fxz, fyy, fyz, fzz)

Gradient and second derivatives of the implicit function of surface type `S`. Compiled once per surface type from
[`implicit`](@ref).
"""
@generated function implicitderivatives(::Type{S}, p::SVector{3}, c::SVector{K}) where {S,K}
    vars = FastDifferentiation.make_variables(:implicit_p, 3)
    cvars = FastDifferentiation.make_variables(:implicit_c, K)
    f = FastDifferentiation.Node(implicit(S, vars, cvars))
    g = [FastDifferentiation.derivative(f, v) for v in vars]
    d(i, j) = FastDifferentiation.derivative(g[i], vars[j])
    outputs = SVector(g[1], g[2], g[3], d(1, 1), d(1, 2), d(1, 3), d(2, 2), d(2, 3), d(3, 3))
    return _inline_derivative_code(outputs, vars, cvars, :p, :c)
end

"""
    curvature_from_implicit(derivs) -> (normal, d1, d2, C)

Normal, orthonormal tangent basis and curvature matrix (see [`surfacecurvature`](@ref)) from the gradient and second
derivatives of an implicit function, as returned by [`implicitderivatives`](@ref).
"""
function curvature_from_implicit(derivs::SVector{9,T}) where {T}
    g = SVector{3,T}(derivs[1], derivs[2], derivs[3])
    H = SMatrix{3,3,T,9}(derivs[4], derivs[5], derivs[6], derivs[5], derivs[7], derivs[8], derivs[6], derivs[8], derivs[9])
    gnorm = norm(g)
    n = g / gnorm
    # any unit tangent: project the coordinate axis least aligned with n
    a = abs(n[1]) < T(0.9) ? SVector{3,T}(1, 0, 0) : SVector{3,T}(0, 1, 0)
    d1 = normalize(a - dot(a, n) * n)
    d2 = cross(n, d1)
    C = -SMatrix{2,2,T,4}(dot(d1, H * d1), dot(d2, H * d1), dot(d1, H * d2), dot(d2, H * d2)) / gnorm
    return (normal=n, d1=d1, d2=d2, C=C)
end

implicitsurfacecurvature(s::S, p::SVector{3,T}) where {S,T} = curvature_from_implicit(implicitderivatives(S, p, implicitcoefficients(s)))

#sphere centered at the origin, outward normal
implicit(::Type{<:Sphere}, p, c) = p[1]^2 + p[2]^2 + p[3]^2 - c[1]^2
implicitcoefficients(s::Sphere{T}) where {T} = SVector{1,T}(s.radius)

#cylinder along the z axis, outward normal
implicit(::Type{<:Cylinder}, p, c) = p[1]^2 + p[2]^2 - c[1]^2
implicitcoefficients(s::Cylinder{T}) where {T} = SVector{1,T}(s.radius)

#plane n⋅p = d
implicit(::Type{<:Plane}, p, c) = c[1] * p[1] + c[2] * p[2] + c[3] * p[3] - c[4]
implicitcoefficients(s::Plane{T}) where {T} = SVector{4,T}(s.normal..., dot(s.normal, s.pointonplane))

#spherical cap: sphere of the cap's radius centered at centrepoint + radius⋅centrenormal; the cap's normal points toward the center
implicit(::Type{<:SphericalCap}, p, c) = c[1]^2 - ((p[1] - c[2])^2 + (p[2] - c[3])^2 + (p[3] - c[4])^2)
implicitcoefficients(s::SphericalCap{T}) where {T} = SVector{4,T}(s.radius, (s.centrepoint + s.radius * s.centrenormal)...)

surfacecurvature(s::Union{Sphere{T,3},Cylinder{T,3},Plane{T,3},SphericalCap{T}}, p::SVector{3,T}) where {T<:Real} = implicitsurfacecurvature(s, p)
surfacecurvature(s::Union{PlanarShape{T},StopSurface{T}}, p::SVector{3,T}) where {T<:Real} = surfacecurvature(s.plane, p)
surfacecurvature(s::AcceleratedParametricSurface{T}, p::SVector{3,T}) where {T<:Real} = surfacecurvature(s.surface, p)

#########################################################################################################
# Curvature at a ray hit in world space

"""
    leaves(tree::CSGTree) -> Vector{LeafNode}

All leaf nodes of an evaluated CSG tree. Their `id`s are the [`surfaceid`](@ref)s of the intersections they produce.
"""
leaves(a::LeafNode) = LeafNode[a]
leaves(a::Union{UnionNode,IntersectionNode}) = vcat(leaves(a.leftchild), leaves(a.rightchild))
leaves(a::ComplementNode) = leaves(a.child)
export leaves

"""
    surfacecurvature(leaf::LeafNode, p::SVector{3}) -> (normal, d1, d2, C)

[`surfacecurvature`](@ref) of the leaf's surface at the world space point `p`, with `normal`, `d1` and `d2` in world
space. `C` does not change under the leaf's rigid transform. The normal is the surface's own outward normal, not the
possibly reversed normal of an intersection inside a complemented CSG node; `C` is measured along the returned normal.
"""
function surfacecurvature(a::LeafNode{T}, p::SVector{3,T}) where {T<:Real}
    loc = surfacecurvature(a.geometry, a.invtransform * p)
    m = Geometry.matrix(a.transform)
    R = SMatrix{3,3,T,9}(m[1, 1], m[2, 1], m[3, 1], m[1, 2], m[2, 2], m[3, 2], m[1, 3], m[2, 3], m[3, 3])
    return (normal=R * loc.normal, d1=R * loc.d1, d2=R * loc.d2, C=loc.C)
end

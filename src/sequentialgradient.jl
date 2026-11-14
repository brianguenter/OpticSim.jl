# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# Fast gradients of the spot merit function for sequential lenses, with FastDifferentiation.
#
# A ray that passes every surface of a sequential lens in order goes through the same chain of operations (intersect a
# surface, refract, intersect the next ...) whatever the ray and whatever the lens parameters. So the map
#     (lens parameters θ, ray start o and d, the ray's intersection parameters tⱼ) -> (x, y) on the image
# is one fixed expression per lens layout (the number and kind of its surfaces). It is built symbolically once with
# FastDifferentiation, together with its Jacobian, and compiled. A ray's image point and derivatives then take two steps:
#   1. trace the ray numerically, with a Newton iteration for each intersection, and record each tⱼ;
#   2. evaluate the compiled function with θ, o, d and the tⱼ.
#
# The tⱼ are inputs, but they depend on θ: an intersection moves when its surface does. Inside the compiled function
# each tⱼ goes through one Newton step, t̂ = t - F(t)/F'(t), where F(t) = 0 is the intersection condition. At the root
# t̂ = t and ∂t̂/∂t = F F''/F'² = 0, so the derivative of t̂ with respect to θ is -F_θ/F_t, the implicit function theorem:
# the Jacobian is exact although the tⱼ are passed in as numbers. The comment at the Newton step in
# buildsequentialjacobian has the derivation; the Optimization page of the documentation explains it at more length.
#
# The gradient with respect to an optimization problem's variables x then follows by the chain rule: the merit's gradient
# with respect to θ (and to each ray's start) from the compiled Jacobians, then one forward mode pass through the map
# x -> (θ, ray starts), which traces no rays. That last step needs ForwardDiff; see `spotmerit_and_gradient`.

#############################################################################
# Lens layout and parameters

# Kind of a surface: :curved (conic with even aspheric terms), :plane, :stop (an aperture only) or :image, and the
# number of its even aspheric terms (orders 4, 6, ...). A lens's layout is the list of these; it fixes the compiled
# function.
struct SequentialSurface
    kind::Symbol
    naspherics::Int
end

"""
    SequentialLens{T}

A sequential lens read from a prescription by [`OpticSim.sequentiallens`](@ref): its layout, its parameters `θ`, the
clear semi-diameter of each surface and the refractive index of each medium at each wavelength.

`θ` holds, surface by surface in order: for a curved surface c, k, a₄, a₆, ..., z; for a plane, a stop or the image, z.
z is the vertex position (light travels toward -z) and a curved surface is z = vertex - s(r) with
s = c r²/(1 + √(1 - (1 + k) c² r²)) + a₄ r⁴ + a₆ r⁶ + ..., as in OpticSim.
"""
struct SequentialLens{T}
    layout::Vector{SequentialSurface}
    θ::Vector{T}
    semidiameters::Vector{Float64}
    indices::Dict{Float64,Vector{Float64}} # wavelength => index of the object space, then after each refracting surface
end

# refractive index of a material at λ (μm) relative to air, as OpticSim uses it
sequentialindex(m, λ) = (m === missing || AGFFileReader.isair(m)) ? 1.0 : Float64(AGFFileReader.index(m, λ))

# the even aspheric coefficients of orders 4, 6, ... of a prescription row (OpticSim's unnormalized "Parameters")
function sequentialaspherics(row, ::Type{T}) where {T}
    hasproperty(row, :Parameters) || return T[]
    params = row.Parameters
    (params === missing || isempty(params)) && return T[]
    orders = Dict(parse(Int, k) => v for (k, v) in params)
    all(j -> iseven(j) && j >= 4, keys(orders)) || throw(ArgumentError("the sequential gradient supports even aspheric terms of order 4 and up, got orders $(sort(collect(keys(orders))))"))
    return T[get(orders, j, zero(T)) for j in 4:2:maximum(keys(orders))]
end

"""
    sequentiallens(prescription, wavelengths) -> SequentialLens

The sequential lens described by a prescription (a DataFrame as used by [`AxisymmetricOpticalSystem`](@ref)), with the
refractive indices at `wavelengths` (μm). Supports Standard and Aspheric surfaces with even aspheric terms, stops and a
flat image. Works with any number type in the prescription, so it can be differentiated with ForwardDiff.
"""
function sequentiallens(p::DataFrame, wavelengths)
    T = promote_type(Float64, nonmissingtype(eltype(p.Thickness)), nonmissingtype(eltype(p.Radius)))
    thickness = [ismissing(t) || isinf(primalvalue(t)) ? zero(T) : T(t) for t in p.Thickness]
    vertices = -cumsum(thickness)
    layout = SequentialSurface[]
    θ = T[]
    semidiameters = Float64[]
    materials = Any[p[1, :Material]]
    for i in 2:size(p, 1)
        row = p[i, :]
        z = vertices[i-1]
        if row.SurfaceType == "Image"
            isinf(primalvalue(row.Radius)) || throw(ArgumentError("the sequential gradient supports a flat image only"))
            push!(layout, SequentialSurface(:image, 0))
            push!(θ, z)
            push!(semidiameters, Inf)
            break
        end
        if isequal(p[i, :Material], p[i-1, :Material])
            row.SurfaceType == "Stop" || continue # a dummy surface: nothing changes there
            push!(layout, SequentialSurface(:stop, 0))
            push!(θ, z)
            push!(semidiameters, Float64(primalvalue(row.SemiDiameter)))
            continue
        end
        row.SurfaceType in ("Standard", "Aspheric", "Stop") || throw(ArgumentError("the sequential gradient doesn't support $(row.SurfaceType) surfaces"))
        R = T(row.Radius)
        c = isinf(primalvalue(R)) ? zero(T) : 1 / R
        k = hasproperty(row, :Conic) && !ismissing(row.Conic) ? T(row.Conic) : zero(T)
        a = sequentialaspherics(row, T)
        if iszero(primalvalue(c)) && isempty(a)
            push!(layout, SequentialSurface(:plane, 0))
            push!(θ, z)
        else
            push!(layout, SequentialSurface(:curved, length(a)))
            append!(θ, [c, k])
            append!(θ, a)
            push!(θ, z)
        end
        push!(semidiameters, Float64(primalvalue(row.SemiDiameter)))
        push!(materials, p[i, :Material])
    end
    indices = Dict(Float64(λ) => [sequentialindex(m, λ) for m in materials] for λ in wavelengths)
    return SequentialLens{T}(layout, θ, semidiameters, indices)
end

#############################################################################
# Sag, shared by the numeric trace and the compiled function

# sag s and ds/d(r²) at r² of a curved surface with curvature c, conic k and aspheric coefficients a (orders 4, 6, ...)
@inline function sequentialsag(r2, c, k, a)
    root = sqrt(1 - (1 + k) * c^2 * r2)
    s = c * r2 / (1 + root)
    ds = c / (2root)
    p = r2
    for (j, aj) in enumerate(a) # a[j] multiplies r^(2j + 2) = (r²)^(j + 1)
        s += aj * p * r2
        ds += (j + 1) * aj * p
        p *= r2
    end
    return s, ds
end

#############################################################################
# Numeric trace

"""
    sequentialtrace(lens, θ, o, d, λ) -> ((x, y), ts) or nothing

Traces the ray from `o` along unit `d` at wavelength `λ` through `lens` with parameters `θ`: (x, y) is where it meets the
image and `ts` its intersection parameter at each curved surface. `nothing` if it misses a surface or falls outside a
clear semi-diameter, or on total internal reflection: the ray doesn't take the sequential path to the image.
"""
function sequentialtrace(lens::SequentialLens, θ::AbstractVector{Float64}, o::SVector{3,Float64}, d::SVector{3,Float64}, λ::Float64)
    n = lens.indices[λ]
    ts = Float64[]
    j = 1 # next parameter in θ
    m = 1 # medium the ray is in
    for (i, s) in enumerate(lens.layout)
        if s.kind === :curved
            c, k = θ[j], θ[j+1]
            a = view(θ, (j+2):(j+1+s.naspherics))
            zv = θ[j+2+s.naspherics]
            j += 3 + s.naspherics
            t = (zv - o[3]) / d[3] # start at the vertex plane
            converged = false
            for _ in 1:30
                P = o + t * d
                r2 = P[1]^2 + P[2]^2
                1 - (1 + k) * c^2 * r2 < 0 && return nothing
                sg, ds = sequentialsag(r2, c, k, a)
                Δ = (P[3] - zv + sg) / (d[3] + 2ds * (P[1] * d[1] + P[2] * d[2]))
                t -= Δ
                if abs(Δ) < 1e-13 * (1 + abs(t))
                    converged = true
                    break
                end
            end
            converged || return nothing
            P = o + t * d
            r2 = P[1]^2 + P[2]^2
            (r2 > lens.semidiameters[i]^2 || 1 - (1 + k) * c^2 * r2 < 0) && return nothing
            _, ds = sequentialsag(r2, c, k, a)
            N = normalize(SVector(2ds * P[1], 2ds * P[2], 1.0))
            push!(ts, t)
        else
            zv = θ[j]
            j += 1
            P = o + ((zv - o[3]) / d[3]) * d
            s.kind === :image && return SVector(P[1], P[2]), ts
            P[1]^2 + P[2]^2 > lens.semidiameters[i]^2 && return nothing
            if s.kind === :stop
                o = P
                continue
            end
            N = SVector(0.0, 0.0, 1.0)
        end
        # Snell's law; the normal N points toward +z, against the ray
        μ = n[m] / n[m+1]
        m += 1
        cs = -dot(N, d)
        q = 1 - μ^2 * (1 - cs^2)
        q < 0 && return nothing
        d = normalize(μ * d + (μ * cs - sqrt(q)) * N)
        o = P
    end
    throw(ArgumentError("the lens has no image surface"))
end

#############################################################################
# Compiled Jacobian

const SEQUENTIAL_JACOBIANS = Dict{Vector{SequentialSurface},Any}()
const SEQUENTIAL_JACOBIANS_LOCK = ReentrantLock()

"""
    sequentialjacobian(layout) -> (evaluate!, nθ)

The compiled function of a lens layout: `evaluate!(out, θ, n, o, d, ts)` writes [x, y, vec(J)...] to `out`, where (x, y) is the
image point of the ray from `o` along `d` whose intersection parameters are `ts`, and J = ∂(x, y)/∂[θ; o; d] is a
2 × (nθ + 6) matrix. `n` are the media's refractive indices. Built with FastDifferentiation the first time a layout is
used, then cached.
"""
function sequentialjacobian(layout::Vector{SequentialSurface})
    lock(SEQUENTIAL_JACOBIANS_LOCK) do
        get!(SEQUENTIAL_JACOBIANS, layout) do
            buildsequentialjacobian(layout)
        end
    end
end

function buildsequentialjacobian(layout::Vector{SequentialSurface})
    FD = FastDifferentiation
    nθ = sum(s.kind === :curved ? 3 + s.naspherics : 1 for s in layout)
    θ = FD.make_variables(:θ, nθ)
    n = FD.make_variables(:n, count(s -> s.kind in (:curved, :plane), layout) + 1)
    o0 = FD.make_variables(:o, 3)
    d0 = FD.make_variables(:d, 3)
    tvars = FD.make_variables(:t, max(1, count(s -> s.kind === :curved, layout)))
    o = collect(o0)
    d = collect(d0)
    j = 1
    m = 1
    ci = 1
    image = nothing
    for s in layout
        if s.kind === :curved
            c, k = θ[j], θ[j+1]
            a = θ[(j+2):(j+1+s.naspherics)]
            zv = θ[j+2+s.naspherics]
            j += 3 + s.naspherics
            t = tvars[ci]
            ci += 1
            # One Newton step from the recorded root makes the derivatives with respect to θ exact.
            #
            # The ray meets the surface at P = o + t d where t solves F(t, θ) = P_z - zv + s(r²) = 0. No formula gives t;
            # the numeric trace finds it by iteration and passes it in. Used directly, t would be a constant to
            # FastDifferentiation, which would then miss dt/dθ: the intersection sliding along the ray as the surface
            # moves. The implicit function theorem gives that derivative from F alone: differentiating F(t(θ), θ) = 0,
            #     F_t dt/dθ + F_θ = 0,   so   dt/dθ = -F_θ / F_t.
            # The Newton step t̂ = t - F/F_t produces it. At the root (F = 0):
            #   - its value is t̂ = t, so nothing changes numerically;
            #   - ∂t̂/∂t = 1 - (F_t² - F F_tt)/F_t² = F F_tt/F_t² = 0, so it doesn't matter that t is a constant input;
            #   - ∂t̂/∂θ = -F_θ/F_t + F ∂F_t/∂θ / F_t² = -F_θ/F_t, the implicit function theorem's dt/dθ.
            # Everything after uses t̂, so the chain rule carries dt/dθ to the image point, through the later surfaces
            # too. This holds to first order, which a gradient needs; second derivatives would not be exact. It relies on
            # t being a converged root: sequentialtrace iterates until the Newton update is below 1e-13 (1 + |t|).
            P = o .+ t .* d
            sg, ds = sequentialsag(P[1]^2 + P[2]^2, c, k, a)
            t̂ = t - (P[3] - zv + sg) / (d[3] + 2ds * (P[1] * d[1] + P[2] * d[2]))
            P = o .+ t̂ .* d
            _, ds = sequentialsag(P[1]^2 + P[2]^2, c, k, a)
            N = [2ds * P[1], 2ds * P[2], 1]
            N = N ./ sqrt(N[1]^2 + N[2]^2 + N[3]^2)
        else
            zv = θ[j]
            j += 1
            P = o .+ ((zv - o[3]) / d[3]) .* d
            if s.kind === :image
                image = [P[1], P[2]]
                break
            end
            if s.kind === :stop
                o = P
                continue
            end
            N = [0, 0, 1]
        end
        μ = n[m] / n[m+1]
        m += 1
        cs = -(N[1] * d[1] + N[2] * d[2] + N[3] * d[3])
        d = μ .* d .+ (μ * cs - sqrt(1 - μ^2 * (1 - cs^2))) .* N
        d = d ./ sqrt(d[1]^2 + d[2]^2 + d[3]^2)
        o = P
    end
    image = FD.Node.(image)
    J = FD.jacobian(image, [θ; o0; d0])
    f! = FD.make_function([image; vec(J)], θ, n, o0, d0, tvars; in_place=true)
    return (evaluate! = f!, nθ = nθ)
end

#############################################################################
# Spot merit with its gradient with respect to θ and the ray starts

# The spot merit (mean over the fields of the squared RMS spot radius about the centroid) of `fieldrays` traced through the
# sequential lens of `prescription`, and, if `withgradient`, its gradient with respect to the lens parameters θ and to
# each ray's start [o; d] (a 6 × number of rays matrix, zero for rays that don't reach the image). NaN if the lens can't
# be built or a field's rays all miss the image.
function sequentialspotmerit(prescription::DataFrame, fieldrays; withgradient::Bool=true)
    nothinggradient = (nothing, nothing)
    (validsurfaces(prescription) && !surfacescross(prescription)) || return NaN, nothinggradient...
    wavelengths = unique(Float64(primalvalue(wavelength(r))) for rays in fieldrays for (r, _) in rays)
    lens = sequentiallens(prescription, wavelengths)
    θ = Float64.(primalvalue.(lens.θ))
    jacobian = withgradient ? sequentialjacobian(lens.layout) : nothing
    nθ = withgradient ? jacobian.nθ : 0
    jobs = [(f, r, Float64(primalvalue(w))) for (f, rays) in enumerate(fieldrays) for (r, w) in rays]
    points = Vector{SVector{2,Float64}}(undef, length(jobs))
    jacobians = Matrix{Float64}(undef, withgradient ? 2 * (nθ + 6) : 0, length(jobs))
    reached = falses(length(jobs))
    chunks = collect(Iterators.partition(eachindex(jobs), 64))
    Threads.@threads for chunk in chunks
        out = withgradient ? zeros(2 + 2 * (nθ + 6)) : Float64[]
        for i in chunk
            _, r, _ = jobs[i]
            o = SVector{3,Float64}(primalvalue.(origin(r)))
            d = SVector{3,Float64}(primalvalue.(direction(r)))
            λ = Float64(primalvalue(wavelength(r)))
            hit = sequentialtrace(lens, θ, o, d, λ)
            hit === nothing && continue
            point, ts = hit
            reached[i] = true
            if withgradient
                jacobian.evaluate!(out, θ, lens.indices[λ], o, d, isempty(ts) ? [0.0] : ts)
                points[i] = SVector(out[1], out[2])
                jacobians[:, i] .= view(out, 3:length(out))
            else
                points[i] = point
            end
        end
    end
    nfields = length(fieldrays)
    merit = 0.0
    gθ = zeros(nθ)
    gstart = zeros(6, withgradient ? length(jobs) : 0)
    for f in 1:nfields
        rays = [i for i in eachindex(jobs) if reached[i] && jobs[i][1] == f]
        isempty(rays) && return NaN, nothinggradient...
        W = sum(jobs[i][3] for i in rays)
        centroid = sum(jobs[i][3] * points[i] for i in rays) / W
        merit += sum(jobs[i][3] * sum(abs2, points[i] - centroid) for i in rays) / (W * nfields)
        withgradient || continue
        # d/dv Σ w |p - c|²/W = 2/W Σ w (p - c)ᵀ ∂p/∂v; the centroid's own change drops out, as Σ w (p - c) = 0
        for i in rays
            g = (2 * jobs[i][3] / (W * nfields)) .* (points[i] - centroid)
            J = reshape(view(jacobians, :, i), 2, nθ + 6)
            gθ .+= g[1] .* view(J, 1, 1:nθ) .+ g[2] .* view(J, 2, 1:nθ)
            gstart[:, i] .= g[1] .* view(J, 1, (nθ+1):(nθ+6)) .+ g[2] .* view(J, 2, (nθ+1):(nθ+6))
        end
    end
    return merit, (withgradient ? (gθ, gstart) : nothinggradient)...
end

"""
    spotmerit_and_gradient(lens::LensOptimizationProblem, x; withgradient = true) -> (merit, gradient)

The spot merit of `lens` at `x` (the mean over the fields of the squared RMS spot radius about the centroid, as
`OpticSim.spotmerit`) and its gradient with respect to the variables, computed with compiled sequential ray
Jacobians (FastDifferentiation) instead of ForwardDiff through the ray trace: exact, and many times faster. `lens` must
have been created with a `prescription` function. With `withgradient = false` only the merit is computed and the gradient
is `nothing`. Needs ForwardDiff to be loaded, for the derivatives of the prescription with respect to `x`.
"""
function spotmerit_and_gradient end

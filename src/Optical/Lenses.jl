# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

opticinterface(::Type{S}, insidematerial, outsidematerial, reflectance=zero(S), interfacemode=ReflectOrTransmit) where {S<:Real} = FresnelInterface{S}(insidematerial, outsidematerial, reflectance=reflectance, transmission=one(S) - reflectance, interfacemode=interfacemode)

"""
    SphericalLens(insidematerial, frontvertex, frontradius, backradius, thickness, semidiameter;  lastmaterial = AGFFileReader.Air, nextmaterial = AGFFileReader.Air, frontsurfacereflectance = 0.0, backsurfacereflectance = 0.0, frontdecenter = (0, 0), backdecenter = (0, 0), frontsemidiameter = semidiameter, backsemidiameter = semidiameter, interfacemode = ReflectOrTransmit)

Constructs a simple cylindrical lens with spherical front and back surfaces. The side walls of the lens are absorbing.

The lens has radius `semidiameter`. A surface with a smaller `frontsemidiameter` or `backsemidiameter` ends there, and a flat,
absorbing ring runs from its edge to the lens rim, as on a real lens whose surfaces have different clear apertures. The
same ring is added where a sphere is too small to reach the rim.
"""
function SphericalLens(insidematerial::T, frontvertex::S, frontradius::S, backradius::S, thickness::S, semidiameter::S; lastmaterial::Q=AGFFileReader.Air, nextmaterial::R=AGFFileReader.Air, frontsurfacereflectance::S=zero(S), backsurfacereflectance::S=zero(S), frontdecenter::Tuple{S,S}=(zero(S), zero(S)), backdecenter::Tuple{S,S}=(zero(S), zero(S)), frontsemidiameter::S=semidiameter, backsemidiameter::S=semidiameter, interfacemode=ReflectOrTransmit) where {R<:AGFFileReader.AbstractGlass,Q<:AGFFileReader.AbstractGlass,T<:AGFFileReader.AbstractGlass,S<:Real}
    @assert !isnan(frontradius)
    @assert !isnan(backradius)

    return ConicLens(insidematerial, frontvertex, frontradius, zero(S), backradius, zero(S), thickness, semidiameter, lastmaterial=lastmaterial, nextmaterial=nextmaterial, frontsurfacereflectance=frontsurfacereflectance, backsurfacereflectance=backsurfacereflectance, frontdecenter=frontdecenter, backdecenter=backdecenter, frontsemidiameter=frontsemidiameter, backsemidiameter=backsemidiameter, interfacemode=interfacemode)
end
export SphericalLens

"""
    ConicLens(insidematerial, frontvertex, frontradius, frontconic, backradius, backconic, thickness, semidiameter;  lastmaterial = AGFFileReader.Air, nextmaterial = AGFFileReader.Air, frontsurfacereflectance = 0.0, backsurfacereflectance = 0.0, frontdecenter = (0, 0), backdecenter = (0, 0), frontsemidiameter = semidiameter, backsemidiameter = semidiameter, interfacemode = ReflectOrTransmit)

Constructs a simple cylindrical lens with front and back surfaces with a radius and conic term. The side walls of the lens are absorbing.

`frontsemidiameter` and `backsemidiameter` end a spherical surface early, as for [`SphericalLens`](@ref); conic surfaces extend to
the lens rim.
"""
function ConicLens(insidematerial::T, frontvertex::S, frontradius::S, frontconic::S, backradius::S, backconic::S, thickness::S, semidiameter::S; lastmaterial::Q=AGFFileReader.Air, nextmaterial::R=AGFFileReader.Air, frontsurfacereflectance::S=zero(S), backsurfacereflectance::S=zero(S), nsamples::Int=17, frontdecenter::Tuple{S,S}=(zero(S), zero(S)), backdecenter::Tuple{S,S}=(zero(S), zero(S)), frontsemidiameter::S=semidiameter, backsemidiameter::S=semidiameter, interfacemode=ReflectOrTransmit) where {R<:AGFFileReader.AbstractGlass,Q<:AGFFileReader.AbstractGlass,T<:AGFFileReader.AbstractGlass,S<:Real}
    @assert !isnan(frontradius)
    @assert !isnan(frontconic)

    return AsphericLens(insidematerial, frontvertex, frontradius, frontconic, nothing, backradius, backconic, nothing, thickness, semidiameter, lastmaterial=lastmaterial, nextmaterial=nextmaterial, frontsurfacereflectance=frontsurfacereflectance, backsurfacereflectance=backsurfacereflectance, nsamples=nsamples, frontdecenter=frontdecenter, backdecenter=backdecenter, frontsemidiameter=frontsemidiameter, backsemidiameter=backsemidiameter, interfacemode=interfacemode)
end
export ConicLens

# Magnitude of the sag of a conic surface with vertex radius `radius` at height `r` from the axis. Throws an
# ArgumentError if the surface doesn't extend out to `r`.
function conicsag(radius::S, conic::S, r::S) where {S<:Real}
    c = 1 / abs(radius)
    q = 1 - (1 + conic) * c^2 * r^2
    q < zero(S) && throw(ArgumentError("a conic surface with radius $radius and conic $conic does not extend out to the semidiameter $r"))
    return c * r^2 / (1 + sqrt(q))
end

# Magnitude of the sag of a sphere of radius `radius` where it ends: at height `r`, or at its equator if `r` is beyond it.
spheresag(radius::S, r::S) where {S<:Real} = r >= abs(radius) ? abs(radius) : conicsag(radius, zero(S), r)

"""
    AsphericLens(insidematerial, frontvertex, frontradius, frontconic, frontaspherics, backradius, backconic, backaspherics, thickness, semidiameter;  lastmaterial = AGFFileReader.Air, nextmaterial = AGFFileReader.Air, frontsurfacereflectance = 0.0, backsurfacereflectance = 0.0, frontdecenter = (0, 0), backdecenter = (0, 0), frontsemidiameter = semidiameter, backsemidiameter = semidiameter, interfacemode = ReflectOrTransmit)

Constructs a simple cylindrical lens with front and back surfaces with a radius, conic and aspheric terms. The side walls of the lens are absorbing.

`frontsemidiameter` and `backsemidiameter` end a spherical surface early, as for [`SphericalLens`](@ref); conic and aspheric
surfaces extend to the lens rim.
"""
function AsphericLens(insidematerial::T, frontvertex::S, frontradius::S, frontconic::S, frontaspherics::Union{Nothing,Vector{Pair{Int,S}}}, backradius::S, backconic::S, backaspherics::Union{Nothing,Vector{Pair{Int,S}}}, thickness::S, semidiameter::S; lastmaterial::Q=AGFFileReader.Air, nextmaterial::R=AGFFileReader.Air, frontsurfacereflectance::S=zero(S), backsurfacereflectance::S=zero(S), nsamples::Int=17, frontdecenter::Tuple{S,S}=(zero(S), zero(S)), backdecenter::Tuple{S,S}=(zero(S), zero(S)), frontsemidiameter::S=semidiameter, backsemidiameter::S=semidiameter, interfacemode=ReflectOrTransmit) where {R<:AGFFileReader.AbstractGlass,Q<:AGFFileReader.AbstractGlass,T<:AGFFileReader.AbstractGlass,S<:Real}
    @assert semidiameter > zero(S)
    @assert !isnan(frontradius)

    frontdecenter_l = abs(frontdecenter[1]) > eps(S) || abs(frontdecenter[2]) > eps(S) ? sqrt(frontdecenter[1]^2 + frontdecenter[2]^2) : zero(S)
    backdecenter_l = abs(backdecenter[1]) > eps(S) || abs(backdecenter[2]) > eps(S) ? sqrt(backdecenter[1]^2 + backdecenter[2]^2) : zero(S)

    if isinf(frontradius) && (frontaspherics === nothing)
        #planar
        lens_front = leaf(Plane(SVector{3,S}(0, 0, 1), SVector{3,S}(0, 0, frontvertex), vishalfsizeu=semidiameter, vishalfsizev=semidiameter, interface=interface = opticinterface(S, insidematerial, lastmaterial, frontsurfacereflectance, interfacemode)))
    elseif frontconic == zero(S) && (frontaspherics === nothing)
        #spherical
        if frontradius < zero(S)
            ϕmax = NaNsafeasin(min(abs((min(semidiameter, abs(frontradius)) + frontdecenter_l) / frontradius), one(S))) + S(π / 50)
            lens_front = leaf(SphericalCap(abs(frontradius), ϕmax, SVector{3,S}(0, 0, 1), SVector{3,S}(frontdecenter[1], frontdecenter[2], frontvertex), interface=opticinterface(S, insidematerial, lastmaterial, frontsurfacereflectance, interfacemode)))
        else
            offset = translation(S, frontdecenter[1], frontdecenter[2], frontvertex - frontradius)
            surf = Sphere(abs(frontradius), interface=opticinterface(S, insidematerial, lastmaterial, frontsurfacereflectance, interfacemode))
            lens_front = leaf(surf, offset)
        end
        if frontsemidiameter < semidiameter || abs(frontradius) - frontdecenter_l <= semidiameter
            # The surface ends before the lens rim: at its own semidiameter, smaller than the lens's, or where the sphere
            # can't go further (its radius is smaller than the lens's semidiameter). A flat ring at that height fills the
            # gap to the rim, as on a real lens. It lies outside the optical aperture; like the barrel it absorbs, as a real
            # lens has it blackened or covered by its mount.
            plane = leaf(Plane(SVector{3,S}(0, 0, 1), SVector{3,S}(0, 0, frontvertex - sign(frontradius) * spheresag(frontradius, frontsemidiameter)), vishalfsizeu=semidiameter, vishalfsizev=semidiameter, interface=FresnelInterface{S}(insidematerial, lastmaterial, reflectance=zero(S), transmission=zero(S))))
            if frontradius > zero(S)
                lens_front = lens_front ∪ plane
            else
                lens_front = lens_front ∩ plane
            end
        end
    else
        #conic or aspheric
        if frontaspherics !== nothing
            frontaspherics = [(i, -k) for (i, k) in frontaspherics]
        end
        surf = AcceleratedParametricSurface(AsphericSurface(semidiameter + frontdecenter_l + S(0.01), radius=-frontradius, conic=frontconic, aspherics=frontaspherics), nsamples, interface=opticinterface(S, insidematerial, lastmaterial, frontsurfacereflectance, interfacemode))
        lens_front = leaf(surf, translation(S, frontdecenter[1], frontdecenter[2], frontvertex))
    end
    if isinf(backradius) && (backaspherics === nothing)
        #planar
        lens_rear = leaf(Plane(SVector{3,S}(0, 0, -1), SVector{3,S}(0, 0, frontvertex - thickness), vishalfsizeu=semidiameter, vishalfsizev=semidiameter, interface=interface = opticinterface(S, insidematerial, nextmaterial, backsurfacereflectance, interfacemode)))
    elseif backconic == zero(S) && (backaspherics === nothing)
        #spherical
        if backradius < zero(S)
            offset = translation(S, backdecenter[1], backdecenter[2], (frontvertex - thickness) - backradius)
            surf = Sphere(abs(backradius), interface=opticinterface(S, insidematerial, nextmaterial, backsurfacereflectance, interfacemode))
            lens_rear = leaf(surf, offset)
            if backradius == semidiameter
                # in this case we need to clip the sphere
                plane = leaf(Plane(SVector{3,S}(0, 0, 1), SVector{3,S}(0, 0, frontvertex - thickness + backradius)))
                lens_rear = lens_rear ∩ plane
            end
        else
            ϕmax = NaNsafeasin(min(abs((min(semidiameter, abs(backradius)) + backdecenter_l) / backradius), one(S))) + S(π / 50)
            lens_rear = leaf(SphericalCap(abs(backradius), ϕmax, SVector{3,S}(0, 0, -1), SVector{3,S}(backdecenter[1], backdecenter[2], frontvertex - thickness), interface=opticinterface(S, insidematerial, nextmaterial, backsurfacereflectance, interfacemode)))
        end
        if backsemidiameter < semidiameter || abs(backradius) - backdecenter_l <= semidiameter
            # a flat ring fills the gap between the end of the surface and the lens rim; it absorbs, like the front one
            plane = leaf(Plane(SVector{3,S}(0, 0, -1), SVector{3,S}(0, 0, frontvertex - thickness - sign(backradius) * spheresag(backradius, backsemidiameter)), vishalfsizeu=semidiameter, vishalfsizev=semidiameter, interface=FresnelInterface{S}(insidematerial, nextmaterial, reflectance=zero(S), transmission=zero(S))))
            if backradius < zero(S)
                lens_rear = lens_rear ∪ plane
            else
                lens_rear = lens_rear ∩ plane
            end
        end
    else
        #conic or aspheric
        if backaspherics !== nothing
            backaspherics = Tuple{Int,S}.(backaspherics)
        end
        surf = AcceleratedParametricSurface(AsphericSurface(semidiameter + backdecenter_l + S(0.01), radius=backradius, conic=backconic, aspherics=backaspherics), nsamples, interface=opticinterface(S, insidematerial, nextmaterial, backsurfacereflectance, interfacemode))
        lens_rear = leaf(surf, Transform(zero(S), S(π), zero(S), backdecenter[1], backdecenter[2], frontvertex - thickness))
    end
    # how far a concave surface rises above its vertex at the lens edge, so the barrel is long enough to enclose it. A
    # sphere rises only to where it ends: its own semidiameter, or its equator, where a flat ring takes over.
    concavesag(radius, conic, aspherics, ownsemidiameter) = conic == zero(S) && aspherics === nothing ? spheresag(radius, ownsemidiameter) : conicsag(radius, conic, semidiameter)
    extra_front = frontradius >= zero(S) || isinf(frontradius) ? zero(S) : concavesag(frontradius, frontconic, frontaspherics, frontsemidiameter)
    extra_back = backradius >= zero(S) || isinf(backradius) ? zero(S) : concavesag(backradius, backconic, backaspherics, backsemidiameter)
    barrel_center = ((frontvertex + extra_front) + (frontvertex - thickness - extra_back)) / 2
    lens_barrel = leaf(Cylinder(semidiameter, (thickness + extra_back + extra_front) * 2, interface=FresnelInterface{S}(insidematerial, AGFFileReader.Air, reflectance=zero(S), transmission=zero(S))), translation(S, zero(S), zero(S), barrel_center))
    lens_csg = lens_front ∩ lens_rear ∩ lens_barrel
    return lens_csg
end
export AsphericLens

"""
    FresnelLens(insidematerial, frontvertex, radius, thickness, semidiameter, groovedepth; conic = 0.0, aspherics = nothing, outsidematerial = AGFFileReader.Air, reverse = false)

Create a Fresnel lens as a CSG generator, can be concave or convex. Groove positions are found iteratively based on `groovedepth`. For negative radii the vertex on the central surface is at `frontvertex`, so the total thickness of the lens is `thickness` + `groovedepth`.
The side wall absorbs. By default the lens is flipped so the flat back faces +z; `reverse = true` keeps the grooved side facing +z.
**Aspherics currently not supported**.
"""
function FresnelLens(insidematerial::G, frontvertex::T, radius::T, thickness::T, semidiameter::T, groovedepth::T; conic::T=0.0, aspherics::Union{Nothing,Vector{Pair{Int,T}}}=nothing, outsidematerial::H=AGFFileReader.Air, reverse::Bool=false) where {T<:Real,G<:AGFFileReader.AbstractGlass,H<:AGFFileReader.AbstractGlass}
    @assert abs(radius) > semidiameter
    # TODO - CSG gets more complicated as the surface can be both convex and concave and we need to find the groove positions iteratively
    aspherics === nothing || throw(ArgumentError("FresnelLens: aspherics are not supported"))
    interface = FresnelInterface{T}(insidematerial, outsidematerial)

    # make the fresnel
    if conic == zero(T)
        if radius == zero(T) || isinf(radius)
            throw(ArgumentError("FresnelLens: radius must be finite and nonzero, got $radius"))
        end
        # spherical lens which makes things much easier, we can just calculate groove positions directly
        n = 1
        if radius > zero(T)
            sphere = Sphere(radius, interface=interface)
            fresnel = leaf(sphere, translation(T, zero(T), zero(T), frontvertex - radius))
        else
            ϕmax = NaNsafeasin(abs(semidiameter / radius)) + T(π / 50)
            sphere = SphericalCap(abs(radius), ϕmax, interface=interface)
            fresnel = leaf(sphere, translation(T, zero(T), zero(T), frontvertex))
        end
        while true
            offset = n * groovedepth
            cylrad = sqrt(offset * (2 * abs(radius) - offset))
            if cylrad >= semidiameter
                break
            end
            if radius > zero(T)
                cylinder = leaf(Cylinder(cylrad, groovedepth * 1.5, interface=interface), translation(T, zero(T), zero(T), frontvertex - groovedepth / 2))
                newsurf = leaf(sphere, translation(T, zero(T), zero(T), frontvertex - radius + offset)) - cylinder
                fresnel = newsurf ∪ fresnel
            else
                cylinder = leaf(Cylinder(cylrad, groovedepth * 1.5, interface=interface), translation(T, zero(T), zero(T), frontvertex + groovedepth / 2))
                newsurf = leaf(sphere, translation(T, zero(T), zero(T), frontvertex - offset)) ∪ cylinder
                fresnel = newsurf ∩ fresnel
            end
            n += 1
        end
    else # conic
        if radius == zero(T) || isinf(radius)
            throw(ArgumentError("FresnelLens: radius must be finite and nonzero, got $radius"))
        end
        n = 1
        surface = AcceleratedParametricSurface(ZernikeSurface(semidiameter + T(0.01), radius=-radius, conic=conic), interface=interface)
        fresnel = leaf(surface, translation(T, zero(T), zero(T), frontvertex))
        while true
            offset = n * groovedepth
            if radius < zero(T)
                offset *= -1
            end
            cylrad = sqrt(2 * radius * offset - (conic + one(T)) * offset^2)
            if cylrad >= semidiameter
                break
            end
            if radius > zero(T)
                cylinder = leaf(Cylinder(cylrad, groovedepth * 1.5, interface=interface), translation(T, zero(T), zero(T), frontvertex - groovedepth / 2))
                newsurf = leaf(surface, translation(T, zero(T), zero(T), frontvertex + offset)) - cylinder
                fresnel = newsurf ∪ fresnel
            else
                cylinder = leaf(Cylinder(cylrad, groovedepth * 1.5, interface=interface), translation(T, zero(T), zero(T), frontvertex + groovedepth / 2))
                newsurf = leaf(surface, translation(T, zero(T), zero(T), frontvertex + offset)) ∪ cylinder
                fresnel = newsurf ∩ fresnel
            end
            n += 1
        end
    end

    # the side wall absorbs, like the other lens types
    barrelinterface = FresnelInterface{T}(insidematerial, outsidematerial, reflectance=zero(T), transmission=zero(T))
    outer_barrel = leaf(Cylinder(semidiameter, thickness * 2, interface=barrelinterface), translation(T, zero(T), zero(T), frontvertex - thickness / 2))
    fresnel = outer_barrel ∩ fresnel
    backplane = leaf(Plane(zero(T), zero(T), -one(T), zero(T), zero(T), frontvertex - thickness, vishalfsizeu=semidiameter, vishalfsizev=semidiameter, interface=interface))
    fresnel = backplane ∩ fresnel
    if !reverse
        fresnel = transform(fresnel, rotationd(T, zero(T), T(180.0), zero(T)))
    end
    return fresnel
end
export FresnelLens

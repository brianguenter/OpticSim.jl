# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

export cooketriplet, doubleconvexlensonly


"""
    hemisphere()

Create a geometric hemisphere
"""
function hemisphere()::CSGTree
    sph = Sphere(10.0)
    pln = Plane(0.0, 0.0, -1.0, 0.0, 0.0, 0.0)
    # CSV operations create a csggenerator which instantiates the csg tree after applying a rigid body transformation.
    # This allows you to make as many instances of the object as you want with different transformations. We just want
    # the CSGTree object rather than a generator.
    return (sph ∩ pln)()
end

"""
    opticalhemisphere()

Create an optical hemisphere that has optical material properties so it will reflect and refract light. In the previous
example the hemisphere object had optical properties of Air, which is the default optical interface, so it won't refract
or reflect light.
"""
function opticalhemisphere()::CSGOpticalSystem
    sph = Sphere(10.0, interface=FresnelInterface{Float64}(Examples_N_BK7, Air))
    pln = Plane(0.0, 0.0, -1.0, 0.0, 0.0, 0.0, interface=FresnelInterface{Float64}(Examples_N_BK7, Air))
    assy = LensAssembly((sph ∩ pln)())
    return CSGOpticalSystem(assy, Rectangle(1.0, 1.0, SVector{3,Float64}(0.0, 0.0, 1.0), SVector{3,Float64}(0.0, 0.0, -11.0), interface=opaqueinterface()))
end

function cooketriplet(::Type{T}=Float64, detpix::Int=1000) where {T<:Real}
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType=["Object", "Standard", "Standard", "Standard", "Stop", "Standard", "Standard", "Image"],
            Radius=[Inf, 26.777, 66.604, -35.571, 35.571, 35.571, -26.777, Inf],
            # OptimizeRadius = [false, true, true, true, true, true, true, false],
            Thickness=[Inf, 4.0, 2.0, 4.0, 2.0, 4.0, 44.748, missing],
            # OptimizeThickness = [false, true, true, true, true, true, true, false],
            Material=[Air, Examples_N_SK16, Air, Examples_N_SF2, Air, Examples_N_SK16, Air, missing],
            SemiDiameter=[Inf, 8.580, 7.513, 7.054, 6.033, 7.003, 7.506, 15.0]
        ),
        detpix,
        detpix
    )
end

function cooketripletfirstelement(::Type{T}=Float64) where {T<:Real}
    AxisymmetricOpticalSystem(
        DataFrame(
            SurfaceType=["Object", "Standard", "Standard", "Image"],
            Radius=[Inf, -35.571, 35.571, Inf],
            Thickness=[Inf, 4.0, 44.748, missing],
            Material=[Air, Examples_N_SK16, Air, missing],
            SemiDiameter=[Inf, 7.054, 6.033, 15.0]
        )
    )
end

function convexplano(::Type{T}=Float64) where {T<:Real}
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType=["Object", "Standard", "Standard", "Image"],
            Radius=[Inf, 60.0, Inf, Inf],
            Thickness=[Inf, 10.0, 57.8, missing],
            Material=[Air, Examples_N_BK7, Air, missing],
            SemiDiameter=[Inf, 9.0, 9.0, 15.0]
        )
    )
end



function doubleconvex(frontradius::T, rearradius::T) where {T<:Real}
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType=["Object", "Standard", "Standard", "Image"],
            Radius=[convert(T, Inf64), frontradius, rearradius, convert(T, Inf64)],
            # OptimizeRadius = [false, true, true, false],
            Thickness=[convert(T, Inf64), convert(T, 10.0), convert(T, 57.8), missing],
            # OptimizeThickness = [false, false, false, false],

            Material=[Air, Examples_N_BK7, Air, missing],
            SemiDiameter=[convert(T, Inf64), convert(T, 9.0), convert(T, 9.0), convert(T, 15.0)]
        )
    )
end

function doubleconvexconic(::Type{T}=Float64) where {T<:Real}
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType=["Object", "Standard", "Standard", "Image"],
            Radius=[Inf64, 60, -60, Inf64],
            # OptimizeRadius = [false, true, true, false],
            Thickness=[Inf64, 10.0, 57.8, missing],
            # OptimizeThickness = [false, false, false, false],
            Conic=[missing, 0.01, 0.01, missing],
            # OptimizeConic = [false, true, true, false],
            Material=[Air, Examples_N_BK7, Air, missing],
            SemiDiameter=[Inf64, 9.0, 9.0, 15.0]
        )
    )
end

function doubleconvexlensonly(frontradius::T, rearradius::T) where {T<:Real}
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType=["Object", "Standard", "Standard", "Image"],
            Radius=[convert(T, Inf64), frontradius, rearradius, convert(T, Inf64)],
            # OptimizeRadius = [false, true, true, false],
            Thickness=[convert(T, Inf64), convert(T, 10.0), convert(T, 57.8), missing],
            # OptimizeThickness = [false, false, false, false],
            Material=[Air, Examples_N_BK7, Air, missing],
            SemiDiameter=[convert(T, Inf64), convert(T, 9.0), convert(T, 9.0), convert(T, 15.0)]
        )
    )
end

function doubleconvex(
    ::Type{T}=Float64;
    temperature::Unitful.Temperature=AGFFileReader.TEMP_REF_UNITFUL,
    pressure::T=convert(T, PRESSURE_REF)
) where {T<:Real}
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType=["Object", "Standard", "Standard", "Image"],
            Radius=[Inf64, 60, -60, Inf64],
            # OptimizeRadius = [false, true, true, false],
            Thickness=[Inf64, 10.0, 57.8, missing],
            # OptimizeThickness = [false, true, true, false],
            Material=[Air, Examples_N_BK7, Air, missing],
            SemiDiameter=[Inf64, 9.0, 9.0, 15.0]
        );
        temperature,
        pressure
    )
end

function doubleconcave(::Type{T}=Float64) where {T<:Real}
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType=["Object", "Standard", "Standard", "Image"],
            Radius=[Inf64, -41.0, 41.0, Inf64],
            Thickness=[Inf64, 10.0, 57.8, missing],
            Material=[Air, Examples_N_BK7, Air, missing],
            SemiDiameter=[Inf64, 9.0, 9.0, 15.0]
        )
    )
end

function planoconcaverefl(::Type{T}=Float64) where {T<:Real}
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType=["Object", "Standard", "Standard", "Image"],
            Radius=[Inf64, Inf64, -41.0, Inf64],
            Thickness=[Inf64, 10.0, -57.8, missing],
            Material=[Air, Examples_N_BK7, Air, missing],
            SemiDiameter=[Inf64, 9.0, 9.0, 25.0],
            Reflectance=[missing, missing, 1.0, missing]
        )
    )
end

function concaveplano(::Type{T}=Float64) where {T<:Real}
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType=["Object", "Standard", "Standard", "Image"],
            Radius=[Inf64, -41.0, Inf64, Inf64],
            Thickness=[Inf64, 10.0, 57.8, missing],
            Material=[Air, Examples_N_BK7, Air, missing],
            SemiDiameter=[Inf64, 9.0, 9.0, 15.0]
        )
    )
end

function planoplano(::Type{T}=Float64) where {T<:Real}
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType=["Object", "Standard", "Standard", "Image"],
            Radius=[Inf64, Inf64, Inf64, Inf64],
            Thickness=[Inf64, 10.0, 57.8, missing],
            Material=[Air, Examples_N_BK7, Air, missing],
            SemiDiameter=[Inf64, 9.0, 9.0, 15.0]
        )
    )
end

"""This example no longer works correctly. The visualization code needs to be updated to support RayListSource"""
function prism_refraction()
    # build the triangular prism
    int = FresnelInterface{Float64}(Examples_N_SF14, AGFFileReader.Air)
    s = 2.0
    prism = (
        Plane(
            SVector(0.0, -1.0, 0.0),
            SVector(0.0, -s, 0.0),
            interface=int,
            vishalfsizeu=2 * s,
            vishalfsizev=2 * s
        ) ∩
        Plane(
            SVector(0.0, sind(30), cosd(30)),
            SVector(0.0, s * sind(30), s * cosd(30)),
            interface=int,
            vishalfsizeu=2 * s,
            vishalfsizev=2 * s
        ) ∩
        Plane(
            SVector(0.0, sind(30), -cosd(30)),
            SVector(0.0, s * sind(30), -s * cosd(30)),
            interface=int,
            vishalfsizeu=2 * s,
            vishalfsizev=2 * s
        )
    )
    sys = CSGOpticalSystem(LensAssembly(prism()), Rectangle(15.0, 15.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -20.0), interface=opaqueinterface()))
    # create some 'white' light
    rays = Vector{OpticalRay{Float64,3}}(undef, 0)
    for i in 0:7
        λ = ((i / 7) * 200 + 450) / 1000
        r = OpticalRay(SVector(0.0, -3.0, 10.0), SVector(0.0, 0.5, -1.0), 1.0, λ)
        push!(rays, r)
    end
    raygen = Emitters.Sources.RayListSource(rays)
end

"""FresnelLens constructor doesn't work. Don't use"""
function fresnel(convex=true; kwargs...)
    lens = FresnelLens(Examples_N_BK7, 0.0, convex ? 15.0 : -15.0, 1.0, 8.0, 0.8, conic=0.1)
    sys = CSGOpticalSystem(LensAssembly(lens()), Rectangle(15.0, 15.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -25.0), interface=opaqueinterface()))
end

function all_examples()
    examples = [hemisphere(),
        opticalhemisphere(),
        cooketriplet(),
        cooketripletfirstelement(),
        convexplano(),
        doubleconvex(15.0, -15.0),
        doubleconvexconic(),
        doubleconvexlensonly(15.0, -15.0),
        doubleconcave(),
        planoconcaverefl(),
        concaveplano(),
        planoplano(),
        prism_refraction()
    ]
end
export all_examples
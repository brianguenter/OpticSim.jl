# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

export cooketriplet, doublegauss, biotar, doubleconvexlensonly


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

"""
    doublegauss(T = Float64, detpix = 1000)

A six element double Gauss lens, f/1.7 with a 100.2mm focal length and a 19° field of view: two cemented doublets with
the stop between them, and a singlet on each side. From US patent 2,784,643 (Lautenbacher & Brendel, Agfa, 1957), as
tabulated in COMSOL's double Gauss lens model. The glasses are OHARA S-LAM3 and S-BAH11 and Schott N-SF5, modeled
from their catalog nd and Vd.
"""
function doublegauss(::Type{T}=Float64, detpix::Int=1000) where {T<:Real}
    S_LAM3 = AGFFileReader.modelglass(1.71700, 47.92, 0.0)
    S_BAH11 = AGFFileReader.modelglass(1.66672, 48.32, 0.0)
    N_SF5 = AGFFileReader.modelglass(1.67270, 32.25, 0.0)
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType=["Object", "Standard", "Standard", "Standard", "Standard", "Standard", "Stop",
                "Standard", "Standard", "Standard", "Standard", "Standard", "Image"],
            Radius=[Inf, 75.050, 270.700, 39.270, Inf, 25.650, Inf, -31.870, Inf, -43.510, 221.140, -88.790, Inf],
            Thickness=[Inf, 9.000, 0.100, 16.510, 2.000, 10.990, 13.000, 7.030, 8.980, 0.100, 7.980, 61.418, missing],
            Material=[Air, S_LAM3, Air, S_BAH11, N_SF5, Air, Air, N_SF5, S_LAM3, Air, S_BAH11, Air, missing],
            SemiDiameter=[Inf, 33.0, 33.0, 27.5, 24.5, 19.5, 18.6, 18.5, 21.0, 21.0, 23.0, 23.0, 42.5]
        ),
        detpix,
        detpix
    )
end

# The f/1.4 Biotar, Example 2 of Willy Merté's patent DE 485,798 (US 1,786,916), Carl Zeiss Jena, 1927, normalized to a
# 100mm focal length. r4 reads "-115.0" in the scanned patent; every other radius of Example 2 is Example 1's scaled by
# 0.967-0.970, and Example 1's r4 is -1118.78, so r4 is taken as -1082. This gives the patent's 100mm focal length (99.9),
# while -115 gives 98.95 and makes the second element's edge thickness negative.
const BIOTAR_RADII = [83.6, 321.0, 44.8, -1082.0, 28.3, -38.5, 50.5, -53.2, 106.0, -120.0]
# thickness after each surface
const BIOTAR_THICKNESSES = [10.75, 1.65, 15.55, 5.05, 18.9, 5.05, 21.22, 0.97, 13.9]
# glasses after surfaces 1 to 9 (nothing after surface 2, 5 and 8: air), modelled from the patent's nd and ν
const BIOTAR_GLASSES = (
    AGFFileReader.modelglass(1.64238, 48.0, 0.0), # I, V and VI
    AGFFileReader.modelglass(1.62306, 56.9, 0.0), # II
    AGFFileReader.modelglass(1.57566, 41.2, 0.0), # III
    AGFFileReader.modelglass(1.67270, 32.2, 0.0), # IV
)
# The patent gives no clear apertures or stop position. The stop is the aperture stop: 22.2mm in radius, which passes
# an f/1.43 axial beam, essentially the patent's f/1.4. The deep inner surfaces 5 and 6 meet at a height of about 22.8mm,
# which limits them to 22.5mm, and a 22.2mm stop fits between them only 10.75 to 11.86mm after surface 5; it is placed
# at 11.3mm. The other clear apertures pass the stop's beam. At the edge of an f/1.4 beam the second element (the front
# doublet's crown) comes to a knife edge, so the patent lens loses its outermost axial rays there.
const BIOTAR_SEMIDIAMETERS = [36.0, 34.5, 33.0, 32.0, 22.5, 22.5, 28.0, 28.0, 28.0, 27.0]
const BIOTAR_STOP_SEMIDIAMETER = 22.2
const BIOTAR_STOP_FRACTION = 11.3 / 18.9 # where the stop is in the air space after surface 5

# `stop` is the distance from surface 5 to the stop, within the air space after surface 5
function biotarprescription(radii::AbstractVector{T}, thicknesses::AbstractVector{T}, image::T; stop=BIOTAR_STOP_FRACTION * thicknesses[5]) where {T}
    gI, gII, gIII, gIV = BIOTAR_GLASSES
    d = thicknesses
    DataFrame(
        SurfaceType=["Object", "Standard", "Standard", "Standard", "Standard", "Standard", "Stop",
            "Standard", "Standard", "Standard", "Standard", "Standard", "Image"],
        Radius=[T(Inf); radii[1:5]; T(Inf); radii[6:10]; T(Inf)],
        Thickness=[T(Inf), d[1], d[2], d[3], d[4], stop, d[5] - stop, d[6], d[7], d[8], d[9], image, missing],
        Material=[Air, gI, Air, gII, gIII, Air, Air, gIV, gI, Air, gI, Air, missing],
        SemiDiameter=[T(Inf); T.(BIOTAR_SEMIDIAMETERS[1:5]); T(BIOTAR_STOP_SEMIDIAMETER); T.(BIOTAR_SEMIDIAMETERS[6:10]); T(33.0)],
    )
end

"""
    biotar(T = Float64, detpix = 1000; image = 65.54)

Zeiss's f/1.4 Biotar, the first double Gauss lens: Example 2 of Willy Merté's 1927 patent DE 485,798 (US 1,786,916),
normalized to a 100mm focal length. Six elements in four groups: a singlet, a cemented doublet, the stop, a cemented
doublet and a singlet. `image` is the distance from the last surface to the image; the default is the best focus for
the 5 fields and 3 wavelengths of [`biotar_optimization`](@ref).

The prescription is read from a scan of the patent, with one radius corrected (see `BIOTAR_RADII`); the clear
apertures and the exact stop position are not in the patent and are chosen here.
"""
biotar(::Type{T}=Float64, detpix::Int=1000; image=65.54) where {T<:Real} =
    AxisymmetricOpticalSystem{T}(biotarprescription(T.(BIOTAR_RADII), T.(BIOTAR_THICKNESSES), T(image)), detpix, detpix)

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
    # draw the result; requires a Makie backend, e.g. `using GLMakie`
    Vis.drawtracerays(sys, raygenerator=raygen, test=true, trackallrays=true)
end

"""Fresnel lens with a conic profile"""
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
        prism_refraction(),
        fresnel(),
        fresnel(false)
    ]
end
export all_examples
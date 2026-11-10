# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# Lens optimization problems, recreated from published Zemax examples where one with a complete specification exists.
#
# Each problem is an OptimizationExample: the free variables with starting values and bounds, a function that builds
# the optical system from the variables, and a merit function to minimize. The merit function accepts any number type,
# so ForwardDiff can differentiate it. The problems do not depend on any optimization package; see the Optimization
# page of the documentation for solving them with Optimization.jl.
# The OptimizationExample type and the shared helpers are in optimization_helpers.jl. Runnable examples that solve
# these problems and print the results are in run_optimization_examples.jl.

export singlet_optimization, aspheric_singlet_optimization, cooke40_optimization, cooke_triplet_optimization

#############################################################################
# Glass

# Schott F2, used by the Zemax Cooke triplet sample. It is not one of AGFFileReader's built-in example glasses, so it is
# defined here: F2's Sellmeier coefficients, nd and Vd (Schott catalog), with N-SF2's remaining data (thermal
# coefficients, transmission), which only matter away from 20°C or for absorption.
const Examples_F2 = AGFFileReader.Glass(2, 1.34533359, 0.00997743871, 0.209073176, 0.0470450767, 0.937357162, 111.886764, 0.0, 0.0, NaN, NaN, 0.365, 2.5, 3.1e-6, 1.75e-8, 6.62e-11, 7.51e-7, 8.99e-10, 0.277, 20.0, 0.0081, 1.0, 1.4, 6.68, 1.0, 1, 1.0, [(0.365, 0.007, 25.0), (0.37, 0.06, 25.0), (0.38, 0.4, 25.0), (0.39, 0.68, 25.0), (0.4, 0.83, 25.0), (0.405, 0.865, 25.0), (0.42, 0.926, 25.0), (0.436, 0.949, 25.0), (0.46, 0.961, 25.0), (0.5, 0.975, 25.0), (0.546, 0.986, 25.0), (0.58, 0.987, 25.0), (0.62, 0.984, 25.0), (0.66, 0.984, 25.0), (0.7, 0.987, 25.0), (1.06, 0.997, 25.0), (1.53, 0.984, 25.0), (1.97, 0.93, 25.0), (2.325, 0.76, 25.0), (2.5, 0.67, 25.0)], 1.62004, 1.2, 0.0, 0, 36.37, 0, 2.718, 0.0)

#############################################################################
# Singlet: Zemax knowledge base, "How To Design a Singlet Lens" (Dan Hill, 2005)

const SINGLET_SOURCE = "Zemax knowledge base, \"How To Design a Singlet Lens\" (Dan Hill, 2005): F/4 N-BK7 singlet, focal length 100mm, entrance pupil 25mm, fields 0, 3.5 and 5 degrees, 632.8nm, center thickness 2 to 12mm, minimize RMS spot radius about the centroid"

# Singlet prescription. `front` = (surface type, conic, aspheric coefficients) of the front surface.
function singletprescription(c1::T, t1::T, t2::T; conic=zero(T), aspherics=nothing) where {T}
    n = T(AGFFileReader.index(AGFFileReader.Examples_N_BK7, 0.6328))
    c2 = lastcurvature(T(100), [c1], [t1], [n]) # F/4 solve: keeps the focal length at 100mm
    aspheric = aspherics !== nothing
    DataFrame(
        SurfaceType=["Object", aspheric ? "Aspheric" : "Standard", "Standard", "Image"],
        Radius=[T(Inf), 1 / c1, 1 / c2, T(Inf)],
        Conic=[missing, conic, zero(T), missing],
        Parameters=[missing, aspheric ? ["4" => aspherics[1], "6" => aspherics[2]] : Pair{String,T}[], Pair{String,T}[], missing],
        Thickness=[T(Inf), t1, t2, missing],
        Material=[AGFFileReader.Air, AGFFileReader.Examples_N_BK7, AGFFileReader.Air, missing],
        SemiDiameter=[T(Inf), T(SINGLET_SEMIDIAMETER), T(SINGLET_SEMIDIAMETER), T(15.0)],
    )
end

const SINGLET_FIELDS = [0.0, 3.5, 5.0]
const SINGLET_SEMIDIAMETER = 14.0
singletrays(rings) = [pupilrays(θ, 25.0, 0.0, [0.6328], [1.0]; rings) for θ in SINGLET_FIELDS] # stop at the front surface

"""
    singlet_optimization(; rings = 6) -> OptimizationExample

Zemax's singlet tutorial: an F/4 N-BK7 singlet with 100mm focal length, entrance pupil 25mm, fields 0, 3.5 and 5
degrees, at 632.8nm. The front curvature, center thickness and distance to the image are free; the back radius is set
by an F number solve that keeps the focal length at 100mm. The tutorial's 2 to 12mm center thickness limit is a bound.

The tutorial starts with a flat front surface; here the front curvature starts at 0.001/mm (radius 1000mm) and is kept
positive, since a curvature of exactly zero is an infinite radius. The tutorial reports final RMS spot radii of
50.1, 46.4 and 82.0μm (relative to the chief ray). It also requires an edge thickness above 2mm, which bounds alone
cannot express; check it on the result.
"""
function singlet_optimization(; rings=6)
    rays = singletrays(rings)
    presc(x) = singletprescription(x[1], x[2], x[3])
    build(x::AbstractVector{T}) where {T} = AxisymmetricOpticalSystem{T}(presc(x), 1, 1)
    merit(x::AbstractVector{T}) where {T} = spotmerit(presc(x), rays, T)
    OptimizationExample(
        "Singlet lens (Zemax tutorial)", SINGLET_SOURCE,
        ["front curvature (1/mm)", "center thickness (mm)", "back surface to image (mm)"],
        [0.001, 4.0, 100.0], [(0.001, 0.05), (2.0, 12.0), (50.0, 150.0)],
        build, merit, SINGLET_FIELDS, rays)
end

"""
    aspheric_singlet_optimization(; rings = 6) -> OptimizationExample

The singlet of [`singlet_optimization`](@ref) with an even aspheric front surface: its conic and 4th and 6th order
terms are free as well. The aspheric variables are the sag each term adds at the lens edge (semi-diameter 14mm), in mm,
which keeps them on the same scale as the other variables; they are converted to OpticSim's (unnormalized) coefficients.

This problem is derived from Zemax's singlet tutorial, not taken from a published Zemax example: no asphere example
with a complete specification was found. Compare its result with the spherical singlet's.
"""
function aspheric_singlet_optimization(; rings=6)
    rays = singletrays(rings)
    presc(x) = singletprescription(x[1], x[2], x[3]; conic=x[4], aspherics=(x[5] / SINGLET_SEMIDIAMETER^4, x[6] / SINGLET_SEMIDIAMETER^6))
    build(x::AbstractVector{T}) where {T} = AxisymmetricOpticalSystem{T}(presc(x), 1, 1)
    merit(x::AbstractVector{T}) where {T} = spotmerit(presc(x), rays, T)
    OptimizationExample(
        "Aspheric singlet lens (derived from the Zemax tutorial)", SINGLET_SOURCE * "; front surface made an even asphere",
        ["front curvature (1/mm)", "center thickness (mm)", "back surface to image (mm)", "front conic",
            "front 4th order term (sag at edge, mm)", "front 6th order term (sag at edge, mm)"],
        [0.001, 4.0, 100.0, 0.0, 0.0, 0.0], [(0.001, 0.05), (2.0, 12.0), (50.0, 150.0), (-5.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)],
        build, merit, SINGLET_FIELDS, rays)
end

#############################################################################
# Cooke triplet: Zemax sample file "Cooke 40 degree field.zmx"

const COOKE40_SOURCE = "Zemax sample \"Cooke 40 degree field.zmx\" (\"A simple Cooke triplet\"): entrance pupil 10mm, fields 0, 14 and 20 degrees, wavelengths 0.48, 0.55 and 0.65μm weighted 0.8, 1 and 0.8, SK16/F2/SK16, stop on the back of the middle lens, last radius set by an F/5 solve (focal length 50mm)"

const COOKE40_FIELDS = [0.0, 14.0, 20.0]
const COOKE40_WAVELENGTHS = [0.48, 0.55, 0.65]
const COOKE40_WEIGHTS = [0.8, 1.0, 0.8]
const COOKE40_GLASS = [AGFFileReader.Examples_N_SK16, Examples_F2, AGFFileReader.Examples_N_SK16]
const COOKE40_THICKNESS = [3.25895583, 0.999974567, 2.95207564] # glass center thicknesses, fixed

# x = curvatures of surfaces 1 to 5, the two air gaps and the distance from the last surface to the image
function cooke40prescription(x::AbstractVector{T}) where {T}
    c1, c2, c3, c4, c5, gap1, gap2, back = x
    nsk, nf = (T(AGFFileReader.index(g, 0.55)) for g in COOKE40_GLASS[1:2])
    tsk1, tf, tsk2 = T.(COOKE40_THICKNESS)
    c6 = lastcurvature(T(50), [c1, c2, c3, c4, c5], [tsk1, gap1, tf, gap2, tsk2], [nsk, one(T), nf, one(T), nsk])
    DataFrame(
        SurfaceType=["Object", "Standard", "Standard", "Standard", "Stop", "Standard", "Standard", "Image"],
        Radius=[T(Inf), 1 / c1, 1 / c2, 1 / c3, 1 / c4, 1 / c5, 1 / c6, T(Inf)],
        Thickness=[T(Inf), tsk1, gap1, tf, gap2, tsk2, back, missing],
        Material=[AGFFileReader.Air, COOKE40_GLASS[1], AGFFileReader.Air, COOKE40_GLASS[2], AGFFileReader.Air, COOKE40_GLASS[3], AGFFileReader.Air, missing],
        SemiDiameter=[T(Inf), T(9.5), T(9.5), T(5.0), T(5.0), T(7.5), T(7.5), T(20.0)],
    )
end

"""
    cooke40_optimization(; rings = 6) -> OptimizationExample

The Cooke triplet of Zemax's sample file "Cooke 40 degree field.zmx": entrance pupil 10mm, fields 0, 14 and 20
degrees, wavelengths 0.48, 0.55 and 0.65μm weighted 0.8, 1 and 0.8, glasses SK16, F2 and SK16 (SK16 is modelled with
the equivalent N-SK16). The stop is on the back of the middle lens.

Free: the curvatures of the first five surfaces, both air gaps and the distance to the image. Glass thicknesses are
fixed. Each curvature is bounded so its radius stays at least 1% larger than the surface's semi-diameter. The last
curvature is set by the sample's F/5 solve, which keeps the focal length at 50mm; a box bound cannot hold a solved
value, so the merit function returns NaN if the solve gives a radius smaller than the last surface's semi-diameter.

The starting values are the sample's own, already optimized, design, so the optimizer should find little to improve:
this example checks that OpticSim and Zemax agree on what is optimal.
"""
function cooke40_optimization(; rings=6)
    nsk, nf = (AGFFileReader.index(g, 0.55) for g in COOKE40_GLASS[1:2])
    x0 = [4.54264791e-2, -2.2948389e-3, -4.50181210e-2, 4.92806889e-2, 1.25496334e-2, 6.0075511, 4.75040893, 42.2077801]
    # entrance pupil: the image of the stop (after surface 4) through surfaces 1 to 3, at the starting design
    zpupil = entrancepupildistance(x0[1:3], [COOKE40_THICKNESS[1], x0[6], COOKE40_THICKNESS[2]], [nsk, 1.0, nf])
    rays = [pupilrays(θ, 10.0, zpupil, COOKE40_WAVELENGTHS, COOKE40_WEIGHTS; rings) for θ in COOKE40_FIELDS]
    cmax = maxcurvature.([9.5, 9.5, 5.0, 5.0, 7.5]) # semi-diameters of surfaces 1 to 5
    presc(x) = cooke40prescription(x)
    build(x::AbstractVector{T}) where {T} = AxisymmetricOpticalSystem{T}(presc(x), 1, 1)
    merit(x::AbstractVector{T}) where {T} = spotmerit(presc(x), rays, T)
    OptimizationExample(
        "Cooke triplet, 40 degree field (Zemax sample)", COOKE40_SOURCE,
        vcat(["curvature surface $i (1/mm)" for i in 1:5], ["air gap lens 1 to 2 (mm)", "air gap lens 2 to 3 (mm)", "last surface to image (mm)"]),
        x0, [[(-c, c) for c in cmax]; [(1.0, 15.0), (1.0, 15.0), (30.0, 60.0)]],
        build, merit, COOKE40_FIELDS, rays)
end

#############################################################################
# Cooke triplet: OpticSim's cooketriplet()

"""
    cooke_triplet_optimization(; rings = 6) -> OptimizationExample

OpticSim's [`cooketriplet`](@ref) as an optimization problem: the five lens radii (the stop surface excluded, as
curvatures), the two air gaps between the lenses and the distance to the image are free. Radii may change by 25% and stay
at least 1% larger than their surface's semi-diameter, gaps
range from 0.5 to 5mm and the image distance from 35 to 55mm. Fields 0 and 5 degrees at 0.55μm, entrance pupil 12mm.

This is not a published Zemax problem: `cooketriplet` comes without a field or aperture specification, so these are
chosen here.
"""
function cooke_triplet_optimization(; rings=6)
    radii = [26.777, 66.604, -35.571, 35.571, -26.777]
    x0 = [1 ./ radii; 2.0; 2.0; 44.748]
    bounds = [[extrema((1 / (0.75r), 1 / (1.25r))) for r in radii]; [(0.5, 5.0), (0.5, 5.0), (35.0, 55.0)]]
    # also keep each radius at least 1% larger than its semi-diameter (the 25% limits already do for these radii)
    cmax = maxcurvature.([8.580, 7.513, 7.054, 7.003, 7.506])
    bounds[1:5] = [(max(l, -c), min(u, c)) for ((l, u), c) in zip(bounds[1:5], cmax)]
    function prescription(x::AbstractVector{T}) where {T}
        c1, c2, c3, c4, c5, gap1, gap2, back = x
        DataFrame(
            SurfaceType=["Object", "Standard", "Standard", "Standard", "Stop", "Standard", "Standard", "Image"],
            Radius=[T(Inf), 1 / c1, 1 / c2, 1 / c3, T(35.571), 1 / c4, 1 / c5, T(Inf)],
            Thickness=[T(Inf), T(4.0), gap1, T(4.0), gap2, T(4.0), back, missing],
            Material=[Air, Examples_N_SK16, Air, Examples_N_SF2, Air, Examples_N_SK16, Air, missing],
            SemiDiameter=[T(Inf), T(8.580), T(7.513), T(7.054), T(6.033), T(7.003), T(7.506), T(15.0)],
        )
    end
    nsk, nsf = AGFFileReader.index(Examples_N_SK16, 0.55), AGFFileReader.index(Examples_N_SF2, 0.55)
    zpupil = entrancepupildistance(x0[1:3], [4.0, 2.0, 4.0], [nsk, 1.0, nsf]) # stop at surface 4
    rays = [pupilrays(θ, 12.0, zpupil, [0.55], [1.0]; rings) for θ in (0.0, 5.0)]
    presc(x) = prescription(x)
    build(x::AbstractVector{T}) where {T} = AxisymmetricOpticalSystem{T}(presc(x), 1, 1)
    merit(x::AbstractVector{T}) where {T} = spotmerit(presc(x), rays, T)
    OptimizationExample(
        "Cooke triplet (OpticSim cooketriplet)", "OpticSim's cooketriplet(); field, aperture and bounds chosen for this example",
        vcat(["curvature surface $i (1/mm)" for i in (1, 2, 3, 5, 6)], ["air gap lens 1 to 2 (mm)", "air gap lens 2 to 3 (mm)", "last surface to image (mm)"]),
        x0, bounds, build, merit, [0.0, 5.0], rays)
end

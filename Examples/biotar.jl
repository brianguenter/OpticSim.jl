# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# Zeiss's f/1.4 Biotar (Merté, 1927): the lens data and its specification, in one place. Used by biotar() below and by
# the optimization problems in optimization_examples.jl.

export biotar

#############################################################################
# The lens, from the patent

# Example 2 of Willy Merté's patent DE 485,798 (US 1,786,916), Carl Zeiss Jena, 1927, normalized to a 100mm focal
# length. r4 reads "-115.0" in the scanned patent; every other radius of Example 2 is Example 1's scaled by 0.967-0.970,
# and Example 1's r4 is -1118.78, so r4 is taken as -1082. This gives the patent's 100mm focal length (99.9), while -115
# gives 98.95 and makes the second element's edge thickness negative.
const BIOTAR_RADII = [83.6, 321.0, 44.8, -1082.0, 28.3, -38.5, 50.5, -53.2, 106.0, -120.0]
const BIOTAR_CURVATURES = 1 ./ BIOTAR_RADII
# thickness after each of surfaces 1 to 9
const BIOTAR_THICKNESSES = [10.75, 1.65, 15.55, 5.05, 18.9, 5.05, 21.22, 0.97, 13.9]
# glasses modelled from the patent's nd and ν
const BIOTAR_GLASS_I = AGFFileReader.modelglass(1.64238, 48.0, 0.0)   # elements I, V and VI
const BIOTAR_GLASS_II = AGFFileReader.modelglass(1.62306, 56.9, 0.0)  # element II
const BIOTAR_GLASS_III = AGFFileReader.modelglass(1.57566, 41.2, 0.0) # element III
const BIOTAR_GLASS_IV = AGFFileReader.modelglass(1.67270, 32.2, 0.0)  # element IV
# the material after each of surfaces 1 to 9: a singlet, a cemented doublet, the air space with the stop, a cemented
# doublet and a singlet
const BIOTAR_MATERIALS = [BIOTAR_GLASS_I, Air, BIOTAR_GLASS_II, BIOTAR_GLASS_III, Air, BIOTAR_GLASS_IV, BIOTAR_GLASS_I, Air, BIOTAR_GLASS_I]
const BIOTAR_GLASS_THICKNESSES = Tuple(findall(!=(Air), BIOTAR_MATERIALS)) # the BIOTAR_THICKNESSES that are glass

const BIOTAR_FOCAL_LENGTH_NOMINAL = 100.0 # the patent's
const BIOTAR_FNUMBER = 1.4
const BIOTAR_ENTRANCE_PUPIL_DIAMETER = BIOTAR_FOCAL_LENGTH_NOMINAL / BIOTAR_FNUMBER

#############################################################################
# Chosen here: the patent gives no clear apertures, stop position or image distance

# The stop is the aperture stop: 22.2mm in radius, which passes an f/1.43 axial beam, essentially the patent's f/1.4. The
# deep inner surfaces 5 and 6 meet at a height of about 22.8mm, which limits them to 22.5mm, and a 22.2mm stop fits
# between them only 10.75 to 11.86mm after surface 5; it is placed at 11.3mm. The other clear apertures pass the stop's
# beam. At the edge of an f/1.4 beam the second element (the front doublet's crown) comes to a knife edge, so the patent
# lens loses its outermost axial rays there.
const BIOTAR_SEMIDIAMETERS = [36.0, 34.5, 33.0, 32.0, 22.5, 22.5, 28.0, 28.0, 28.0, 27.0]
const BIOTAR_STOP_SEMIDIAMETER = 22.2
const BIOTAR_STOP_DISTANCE = 11.3 # from surface 5 to the stop
const BIOTAR_IMAGE_SEMIDIAMETER = 33.0
# the last surface to the image: the patent lens refocused for the fields and wavelengths below
const BIOTAR_IMAGE_DISTANCE = 65.54

#############################################################################
# The specification the lens is evaluated and optimized for

const BIOTAR_SOURCE = "Willy Merté, patent DE 485,798 / US 1,786,916 (Carl Zeiss Jena, 1927), Example 2, f/1.4, normalized to a 100mm focal length; 5 fields over a 17.3° half field (silent-era 35mm cine full aperture), F, d and C lines"
const BIOTAR_WAVELENGTHS = [0.4861, 0.5876, 0.6563] # F, d and C lines
const BIOTAR_REFERENCE_WAVELENGTH = BIOTAR_WAVELENGTHS[2] # d line, for paraxial solves
# relative field heights 0, 0.3, 0.5, 0.7 and 1 of a 17.3° half field: the first f/1.4 Biotar was a 5cm lens for 35mm cine
# film. Its silent-era full aperture (before optical soundtracks) was 0.980 x 0.735in, 24.89 x 18.67mm, 4:3, with a 31.1mm
# diagonal, and atan(15.56/50) = 17.3°. (The patent itself gives no field.)
const BIOTAR_FIELDS = [atand(k * hypot(24.89, 18.67) / BIOTAR_FOCAL_LENGTH_NOMINAL) for k in (0.0, 0.3, 0.5, 0.7, 1.0)]

#############################################################################
# Functions

# refractive index after each of surfaces 1 to 9 at λ
biotarindices(λ) = [m == Air ? 1.0 : AGFFileReader.index(m, λ) for m in BIOTAR_MATERIALS]

# focal length of the patent's lens at the reference wavelength: 99.91mm, the patent's 100 within the precision of its
# data
function biotarfocallength()
    n = biotarindices(BIOTAR_REFERENCE_WAVELENGTH)
    _, nu = paraxialtrace(1.0, 0.0, BIOTAR_CURVATURES, [BIOTAR_THICKNESSES; 0.0], [n; 1.0], 1.0)
    return -1 / nu
end

"""
    biotarprescription(radii, thicknesses, image; stop = BIOTAR_STOP_DISTANCE) -> DataFrame

The Biotar's prescription with the given `radii` (10), `thicknesses` (after surfaces 1 to 9), distance from the last
surface to the image `image`, and distance `stop` from surface 5 to the stop, within the air space after surface 5.

Each surface's clear aperture is the patent lens's, enlarged where needed to pass the whole f/1.4 axial beam (see
[`OpticSim.passaxialbeam`](@ref)), so only the stop limits that beam whatever the radii and spacings. For the patent's
lens nothing is enlarged.
"""
function biotarprescription(radii::AbstractVector{T}, thicknesses::AbstractVector{T}, image::T; stop=BIOTAR_STOP_DISTANCE) where {T}
    d = thicknesses
    m = BIOTAR_MATERIALS
    sd = BIOTAR_SEMIDIAMETERS
    p = DataFrame(
        SurfaceType=["Object", "Standard", "Standard", "Standard", "Standard", "Standard", "Stop",
            "Standard", "Standard", "Standard", "Standard", "Standard", "Image"],
        Radius=[T(Inf); radii[1:5]; T(Inf); radii[6:10]; T(Inf)],
        Thickness=[T(Inf); d[1:4]; stop; d[5] - stop; d[6:9]; image; missing],
        Material=[Air; m[1:4]; m[5]; m[5]; m[6:9]; Air; missing],
        SemiDiameter=[T(Inf); T.(sd[1:5]); T(BIOTAR_STOP_SEMIDIAMETER); T.(sd[6:10]); T(BIOTAR_IMAGE_SEMIDIAMETER)],
    )
    return passaxialbeam(p, BIOTAR_ENTRANCE_PUPIL_DIAMETER / 2, BIOTAR_WAVELENGTHS)
end

# the patent's lens, with the distance from the last surface to the image `image`
biotarpatentprescription(image::T) where {T} = biotarprescription(T.(BIOTAR_RADII), T.(BIOTAR_THICKNESSES), image)

"""
    biotar(T = Float64, detpix = 1000; image = BIOTAR_IMAGE_DISTANCE)

Zeiss's f/1.4 Biotar, the first double Gauss lens: Example 2 of Willy Merté's 1927 patent DE 485,798 (US 1,786,916),
normalized to a 100mm focal length. Six elements in four groups: a singlet, a cemented doublet, the stop, a cemented
doublet and a singlet. `image` is the distance from the last surface to the image; the default, 65.54mm, is the best
focus for the 5 fields and 3 wavelengths of [`biotar_optimization`](@ref).

The prescription is read from a scan of the patent, with one radius corrected (see `BIOTAR_RADII`); the clear
apertures and the exact stop position are not in the patent and are chosen here. All the lens data is in
`Examples/biotar.jl`.
"""
biotar(::Type{T}=Float64, detpix::Int=1000; image=BIOTAR_IMAGE_DISTANCE) where {T<:Real} =
    AxisymmetricOpticalSystem{T}(biotarpatentprescription(T(image)), detpix, detpix)

# Rays filling the f/1.4 entrance pupil for each field, at the F, d and C lines, from one emitter per field and
# wavelength. The rays are aimed: the entrance pupil, the image of the stop through surfaces 1 to 5, is computed
# paraxially for the given design (curvatures `c` of surfaces 1 to 5, thicknesses `d` after surfaces 1 to 4, and the
# distance `stop` from surface 5 to the stop), so the beams follow the stop as it and the front groups change.
function biotarrays(rings; c=BIOTAR_CURVATURES[1:5], d=BIOTAR_THICKNESSES[1:4], stop=BIOTAR_STOP_DISTANCE)
    n = biotarindices(BIOTAR_REFERENCE_WAVELENGTH)
    zpupil = entrancepupildistance(c, [d; stop], n[1:5])
    # the pupil is inside the lens: start the rays 20mm in front of the first surface
    [pupilrays(θ, BIOTAR_ENTRANCE_PUPIL_DIAMETER, zpupil, BIOTAR_WAVELENGTHS, ones(3); rings, distance=zpupil + 20.0) for θ in BIOTAR_FIELDS]
end

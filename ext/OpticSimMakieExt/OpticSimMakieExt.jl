# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# Makie implementation of the OpticSim.Vis drawing functions. Loaded automatically when a Makie backend
# (GLMakie, CairoMakie, WGLMakie) is loaded alongside OpticSim.
module OpticSimMakieExt

using OpticSim
using OpticSim: Vis, euclideancontrolpoints, vertex, makiemesh, detector, centroid, lower, upper, intervals, α, point, Spline
using OpticSim.Geometry
using OpticSim.Emitters
using OpticSim.Emitters: Sources, Origins, Directions, Spectrum, visual_size
using StaticArrays
using LinearAlgebra
import Makie
using Makie: RGB, RGBA, RGBf, Point2f, Point3f
import Makie.GeometryBasics
import Makie.FileIO

include("Visualization.jl")
include("Emitters.jl")

# writes the image files for the outpath keyword of trace and traceMT
function OpticSim.savedetectorimage(path::AbstractString, image::AbstractMatrix{<:Real})
    brightest = maximum(image)
    scaled = brightest > 0 ? image ./ brightest : image
    FileIO.save(path, Makie.Colors.Gray.(scaled))
end

end # module OpticSimMakieExt

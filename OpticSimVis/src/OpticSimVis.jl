# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module OpticSimVis

using OpticSimCore
using OpticSimCore.Geometry
using OpticSimCore.Repeat
using OpticSimCore.Emitters
using OpticSimCore.Emitters: Spectrum, Directions, Origins, AngularPower, Sources
using OpticSimCore: euclideancontrolpoints, evalcsg, vertex, makiemesh, detector, centroid, lower, upper, intervals, α

using Unitful
using ImageView
using Images
using ColorTypes
using ColorSchemes
using Distributions
using StaticArrays
using LinearAlgebra
import Makie
import GeometryBasics
import Plots
import Luxor
using FileIO

include("Visualization.jl")
include("Emitters.jl")
include("VisRepeatingStructures.jl")

end # module OpticSimVis

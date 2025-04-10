# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""Contains example usage of the features in the OpticSim.jl package."""
module Examples
using ..OpticSim
using Vis
using ..OpticSim.Geometry
using ..OpticSim.Emitters

using AGFFileReader
using StaticArrays
using DataFrames: DataFrame
using Images
using Unitful
using Plots
using LinearAlgebra
import Luxor

include("docs_examples.jl")
include("other_examples.jl")
include("eyemodels.jl")

end #module Examples
export Examples

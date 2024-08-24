# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module OpticSim

using Reexport
@reexport using OpticSimCore
import OpticSimCore as Core
import OpticSimExamples as Examples
import OpticSimExamples.Data
import OpticSimVis as Vis
import OpticSimOptimization as Optimization
import OpticSimCloud as Cloud

export Examples
export Vis
export Optimization
export Cloud

# define the NotebooksUtils module
include("NotebooksUtils/NotebooksUtils.jl")

end # module

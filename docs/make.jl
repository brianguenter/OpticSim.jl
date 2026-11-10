# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using Documenter
using OpticSim
using AGFFileReader
import CairoMakie # loads the OpticSim Vis extension; renders static images without a GPU

# return figures from Vis functions so Documenter embeds them in the page
OpticSim.Vis.set_current_mode(:docs)
CairoMakie.activate!(type="png", px_per_unit=1)

makedocs(
    sitename="OpticSim.jl",
    format=Documenter.HTML(
        # prettyurls = get(ENV, "CI", nothing) == "true",
        assets=[asset("assets/logo.svg", class=:ico, islocal=true)],
        size_threshold=300 * 2^10,
        example_size_threshold=0, # always write figures to separate files rather than inlining them in the html
    ),
    modules=[OpticSim],
    pages=[
        "Home" => "index.md",
        "Geometry" => [
            "Primitives" => "primitives.md",
            "Transforms and Vectors" => "transforms_and_vectors.md",
            "CSG" => "csg.md",
        ],
        "Optical" => [
            "Systems" => "systems.md",
            "Emitters" => "emitters.md",
            "Interfaces" => "interfaces.md",
            "Lenses" => "lenses.md"
        ],
        "Visualization" => "vis.md",
        "Optimization" => "optimization.md",
        "Reference" => "ref.md"
    ],
    expandfirst=["systems.md"],
    checkdocs=:public
)

deploydocs(
    repo="github.com/brianguenter/OpticSim.jl.git",
    devbranch="main",
    push_preview=true,
)

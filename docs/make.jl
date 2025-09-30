# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using Documenter
using OpticSim
using AGFFileReader


makedocs(
    sitename="OpticSim.jl",
    format=Documenter.HTML(
        # prettyurls = get(ENV, "CI", nothing) == "true",
        assets=[asset("assets/logo.svg", class=:ico, islocal=true)],
        size_threshold=300 * 2^10
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

<p align="center">
  <a href="https://microsoft.github.io/OpticSim.jl/dev/">
    <img src=docs/src/assets/logo.svg height=128px style="text-align:center">
  </a>
</p>

# OpticSim.jl

<table>
<thead>
  <tr>
    <th>Documentation</th>
    <th>Build Status</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>
      <a href="https://brianguenter.github.io/OpticSim.jl/stable/">
        <img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="docs stable">
      </a>
      <a href="https://brianguenter.github.io/OpticSim.jl/dev/">
        <img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="docs dev">
      </a>
    </td>
    <td>
      <a href="https://github.com/brianguenter/OpticSim.jl/actions/workflows/CI.yml">
        <img src="https://github.com/brianguenter/OpticSim.jl/workflows/CI/badge.svg" alt="CI action">
      </a>
      <a href="https://codecov.io/gh/brianguenterOpticSim.jl">
        <img src="https://codecov.io/gh/brianguenter/OpticSim.jl/branch/main/graph/badge.svg?token=9QxvIHt5F5" alt="codecov">
      </a>
    </td>
  </tr>
</tbody>
</table>

OpticSim.jl is a [Julia](https://julialang.org/) package for geometric optics. One of the design goals of OpticSim is to make it easy to create optical systems procedurally. Unlike Zemax, Code V, or other interactive optical design systems OpticSim.jl has limited support for interactivity.

A large variety of surface types are supported, and these can be composed into complex 3D objects through the use of constructive solid geometry (CSG). A substantial catalog of optical materials is provided through the GlassCat submodule.

# Installation

Before you can use the software you will need to download glass files. See the documentation for detailed information about how to do this.

*Warning*: During installation OpticSim automatically downloads glass catalogs from a variety of public sources. The Schott website keeps moving their catalog on their website so our software can't find it to download. This caused all our examples to fail because they use Schott glasses. We have replaced the glasses in the examples with hard coded glass files so all examples now work (if they don't work for you please file an issue). 

**If you want to use Schott glasses you will have to manually download and install the Schott catalog.**

## Contributing

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

This project welcomes contributions and suggestions.


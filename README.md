
<div align="right">
  <details>
    <summary >🌐 Language</summary>
    <div>
      <div align="center">
        <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=en">English</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=zh-CN">简体中文</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=zh-TW">繁體中文</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=ja">日本語</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=ko">한국어</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=hi">हिन्दी</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=th">ไทย</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=fr">Français</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=de">Deutsch</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=es">Español</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=it">Italiano</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=ru">Русский</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=pt">Português</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=nl">Nederlands</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=pl">Polski</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=ar">العربية</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=fa">فارسی</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=tr">Türkçe</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=vi">Tiếng Việt</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=id">Bahasa Indonesia</a>
        | <a href="https://openaitx.github.io/view.html?user=brianguenter&project=OpticSim.jl&lang=as">অসমীয়া</
      </div>
    </div>
  </details>
</div>

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

A large variety of surface types are supported, and these can be composed into complex 3D objects through the use of constructive solid geometry (CSG). A substantial catalog of optical materials is provided through the AGFFileReader package.

### Package status
This package is currently undergoing a significant rewrite. The latest versions of the package do not have full functionality yet. The core ray tracing works (in the package `OpticSim`) but the glass catalog download, visualization, and repeating structures code has been moved into separate packages: `AGFFileReader`,`OpticSimVisualization`,`OpticSimRepeatingStructures`. The last two packages are not yet full functional.
## Contributing

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

This project welcomes contributions and suggestions.


[![License: Apache 2.0](https://img.shields.io/badge/License-Apache-blue.svg)](https://github.com/COMODO-research/Comodo.jl/blob/main/LICENSE)
[![example workflow](https://github.com/COMODO-research/Comodo.jl/actions/workflows/test.yml/badge.svg)](https://github.com/COMODO-research/Comodo.jl/blob/main/.github/workflows/test.yml) [![Doc](https://img.shields.io/badge/docs-dev-blue.svg)](https://comodo-research.github.io/Comodo.jl/dev/)  [![codecov](https://codecov.io/gh/COMODO-research/Comodo.jl/graph/badge.svg?token=2ZOXAXXX1I)](https://codecov.io/gh/COMODO-research/Comodo.jl)  [![Join the chat at https://gitter.im/Comodo.jl](https://badges.gitter.im/Comodo.jl.svg)](https://app.gitter.im/#/room/#comodo:gitter.im?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
![](assets/img/COMODO.png)  

# About Comodo
Comodo is a [Julia](https://julialang.org/) package for **computational (bio)mechanics and computational design**, and offers functionality for geometry processing, meshing, finite element analysis, automated design, topology optimisation, and image-based modelling. 

Loosely Comodo could stand for **Com**putational **Mo**delling for **D**esign **O**ptimization. A more philosophical angle would be to say that **DO** is like *-do* in the Japanese art *Judo* (ju=柔=gentle, do=道=way), so in this sense Comodo stands for *"the way of computational modelling"*. Comodo is perhaps best defined by its scope. Comodo aims to be a "one-stop-shop" for researchers in computational (bio)mechanics and computational design. It will feature tools for geometry processing, meshing, automated design / topology optimization, finite element analysis, as well as (e.g. medical) image processing and segmentation. 

Comodo.jl started out as a modern re-implementation in Julia of the MATLAB toolbox [GIBBON](https://github.com/gibbonCode/GIBBON). However, rather than literally porting each functional unit, it instead aims to follow a similar philosophy and cover similar but more advanced core functionaly.

# Installation
```julia
julia> ]
(@v1.xx) pkg> add https://github.com/COMODO-research/Comodo.jl
```

or 

```julia
julia> using Pkg
julia> Pkg.add(url = "https://github.com/COMODO-research/Comodo.jl")
```

# Getting started
To get started install the package, study the documentation, and test some of the demos provided in the [`examples`](https://github.com/COMODO-research/Comodo.jl/tree/main/examples) folder. 

<img src="https://github.com/COMODO-research/Comodo_data_docs/blob/main/img_anim/comodo_snippets.gif" alt="Comodo snippets" width="50%"/>


# Documentation 
[Functional Documentation](https://comodo-research.github.io/Comodo.jl/dev/)

# Combining with finite element analysis
For finite element analysis users are encouraged to combine the Comodo capabilities with the open source C++ solver [FEBio](https://febio.org/), e.g. based on the Julia wrapper [FEBio.jl](https://github.com/febiosoftware/FEBio.jl). In addition, users may want to explore the Julia packages [Gridap.jl](https://github.com/gridap/Gridap.jl) and [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl).

<img src="https://github.com/COMODO-research/Comodo_data_docs/blob/main/img_anim/febio_example_01.gif" alt="febio functionality" width="50%"/>

# Testing 
You can test Comodo by running
```julia
pkg> test Comodo
```
The source for the tests is [`runtests.jl`, found in the `test` folder](https://github.com/COMODO-research/Comodo.jl/blob/main/test/runtests.jl)

# Roadmap
New functionality to add:
- [ ] Levelset methods
- [ ] Lattice structure creation: 
	- Element based conversions 
	- Triply periodic minimal and spinodoid surfaces
	- Boundary conforming lattices
- [ ] Surface stitching method

# How to contribute? 
Your help would be greatly appreciated! If you can contribute please do so by posting a pull-request. I am very much open to fully acknowledging your contributions e.g. by listing you as a contributor properly whereever possible, by welcoming you on board as a long term contributor, or by inviting you to be a co-author on publications featuring Comodo functionality. 

# License 
Comodo.jl is released open source under the [Apache 2.0 license](https://github.com/COMODO-research/Comodo.jl/blob/main/LICENSE).

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache-blue.svg)](https://github.com/COMODO-research/Comodo.jl/blob/main/LICENSE)
[![example workflow](https://github.com/COMODO-research/Comodo.jl/actions/workflows/test.yml/badge.svg)](https://github.com/COMODO-research/Comodo.jl/blob/main/.github/workflows/test.yml) 
[![cov](https://github.com/COMODO-research/Comodo.jl/badges/coverage.svg)](https://github.com/COMODO-research/Comodo.jl/actions)
[![Join the chat at https://gitter.im/Comodo.jl](https://badges.gitter.im/Comodo.jl.svg)](https://app.gitter.im/#/room/#comodo:gitter.im?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
![](assets/img/COMODO.png) 

[![Mastodon](https://img.shields.io/badge/-MASTODON-%232B90D9?style=for-the-badge&logo=mastodon&logoColor=white)](https://fosstodon.org/@kevinmoerman) 




# About Comodo
Comodo is a [Julia](https://julialang.org/) package **computational (bio)mechanics and computational design**, and offers functionality for geometry processing, meshing, finite element analysis, automated design, topology optimisation, and image-based modelling. 

Loosely Comodo could stand for **Com**putational **Mo**delling for **D**esign **O**ptimization. A more philosophical angle would be to say that **DO** is like *-do* in the Japanese art *Judo* (ju=柔=gentle, do=道=way), so in this sense Comodo stands for *"the way of computational modelling"*. Comodo is perhaps best defined by its scope. Comodo aims to be a "one-stop-shop" for researchers in computational (bio)mechanics and computational design. It will feature tools for geometry processing, meshing, automated design / topology optimization, finite element analysis, as well as (e.g. medical) image processing and segmentation. 

Comodo.jl started out as a modern re-implementation in Julia of the MATLAB toolbox [GIBBON](https://github.com/gibbonCode/GIBBON). However, rather than literally porting each functional unit, it instead aims to follow a similar philosophy and cover similar but more advanced core functionaly.

## Installation
```julia
pkg> add https://github.com/COMODO-research/Comodo.jl
```
# Getting started
The `examples` folder contains examples on the use of Comodo's functionality. 

# Documentation 
Under construction

# Combinding with finite element analysis
For finite element analysis users are encouraged to combine the Comodo capabilities with the open source C++ solver [FEBio](https://febio.org/), e.g. based on the Julia wrapper [FEBio.jl](https://github.com/febiosoftware/FEBio.jl). In addition, users may want to explore the Julia packages [Gridap.jl](https://github.com/gridap/Gridap.jl) and [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl).

# Testing 
Under construction.

# Roadmap
A detailed roadmap is under construction. 

# How to contribute? 
Your help would be greatly appreciated! If you can contribute please do so by posting a pull-request. I am very much open to fully acknowledging your contributions e.g. by listing you as a contributor properly whereever possible, by welcoming you on board as a long term contributor, or by inviting you to be a co-author on publications featuring Comodo functionality. 

# License 
Comodo.jl is released open source under the [Apache 2.0 license](https://github.com/COMODO-research/Comodo.jl/blob/main/LICENSE).

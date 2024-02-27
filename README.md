# COMODO.jl

![](assets/img/COMODO.png)

Loosely Comodo could stand for **Com**putational **Mo**delling for **D**esign **O**ptimization. But really this Julia package aims to be a "one-stop-shop" for researchers in computational (bio)mechanics and computational design. It will feature tools for geometry processing, meshing, automated design / topology optimization, finite element analysis, as well as (e.g. medical) image processing and segmentation. 
For finite element analysis this project currently features the use of the open source C++ solver FEBio. In the future the use of the Julia packages Gridap.jl and Ferrite.jl for finite element analysis will be explored. 

Comodo.jl started out as a modern re-implementation in Julia of the MATLAB toolbox [GIBBON (The Geometry and Image-Based Bioengineering add-ON)](https://github.com/gibbonCode/GIBBON). However, rather than literally porting each functional unit, it instead aims to follow a similar philosophy and cover similar but more advanced core functionaly.

# Getting started
The `docs` folder contains examples on the use of Comodo's functionality. 

# License 

## Installation
```julia
pkg> add https://github.com/COMODO-research/Comodo.jl
```

## Getting started
The `docs` folder contains usage examples. See for instance: [`demo_febio_0001_cube_uniaxial_hyperelastic.jl`](https://github.com/COMODO-research/Comodo.jl/blob/main/docs/demo_febio_0001_cube_uniaxial_hyperelastic.jl). 

![](assets/img/febio_example_01.gif) 

## Documentation 
Under construction

## Testing 
Under construction

## Roadmap
A detailed roadmap is under construction. 

## How to contribute? 
Your help would be greatly appreciated! If you can contribute please do so by posting a pull-request. 

## License 
Comodo.jl is released open source under the [Apache 2.0 license](https://github.com/COMODO-research/Comodo.jl/blob/main/LICENSE).
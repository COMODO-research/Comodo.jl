# GIBBON.jl

![](assets/img/gibbonLogo.png)

This Julia package aims to be a "one-stop-shop" for researchers in computational (bio)mechanics and computational design. It will feature tools for geometry processing, meshing, automated design / topology optimization, finite element analysis, as well as image processing and segmentation. In terms of finite element analysis linking with the open source solver FEBio is one of the targets. In addition, an interface for Gridap.jl will be added. 

In some sense GIBBON.jl is a modern re-implementation in Julia of the MATLAB toolbox [GIBBON (The Geometry and Image-Based Bioengineering add-ON)](https://github.com/gibbonCode/GIBBON). However, rather than literally porting each functional unit, it instead aims to follow a similar philosophy and cover the same core functionaly.

# Getting started
The `docs` folder contains examples on the use of GIBBON's functions. 

# License 

## Installation
```julia
pkg> add https://github.com/gibbonCode/Gibbon.jl
```

## Getting started
The `docs` folder contains usage examples. See for instance: `demo_febio_input_file_01.jl`. This example should produce:

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
GIBBON.jl is released open source under the [Apache 2.0 license](https://github.com/gibbonCode/Gibbon.jl/blob/master/LICENSE).
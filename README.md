# simple-rasterizer
Software (CPU) 3D Rasterizer capable of drawing 3D objects with lighting

<img src="assets/rasterizer.gif" width="500">

## Features

* Simple C code
* Load OBJ/MTL files
* Perspective-correct texture support
* Point/Ambient/Directional lighting
* Materials with specular and diffuse properties
* Camera controls
* Camera and backface culling
* z-Buffering 
* GUI with controls for lighting, movement etc

## Current limitations

* Diffuse colors on MTL files are used for as ambient colors as well
* No multi-threading/SIMD
* Support only a subset of the OBJ specification

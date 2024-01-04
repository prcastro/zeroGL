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

<!-- ## TODO
* Create the concept of a canvas (pixels, depthbuffer, width, height)
* Make game->framebuffer a canvas
* Make all drawing functions accept a canvas instead of a framebuffer
* Remove WIDTH and HEIGHT from simplerenderer.h. Put it on main.c and use it only when instatiating the canvas.
* Look how olive.c uses a canvas as a texture too
-->
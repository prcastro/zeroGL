# simple-rasterizer
Software (CPU) 3D Rasterizer capable of drawing 3D objects with lighting

<img src="assets/rasterizer.gif" width="500">

## Getting started

To run the demo, download the latest .zip file from the releases page, uncompress it and run

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

## Building the project

### Windows

#### build.bat
Install Visual Studio Community 2022 (the C/C++ compilers may suffice) and then double-click the `build.bat` file. Alternatively, from the command-line run:

```console
> build.bat
```

#### Using cl.exe
Install Visual Studio Community 2022 (the C/C++ compilers may suffice), then run:

```console
> "c:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat" x86_amd64

> cl.exe /Zi /EHsc /O2 /nologo /Fotarget\obj\ /Fetarget\main.exe main.c external\upng\upng.c /I external\sdl2\include /link /LIBPATH:external\sdl2\lib\x64 SDL2.lib SDL2main.lib
```

Change the `vcvarsall.bat` path according to your VS installation.

#### VS Code
If you're using Visual Studio Code, then you can install Visual Studio Community 2022 (the C/C++ compilers may suffice) and also install the C/C++ extensions from Microsoft.

Open `main.c` on the editor and then click on `Run C/C++ File` on the top right corner. Remember to open VS Code from within a developer terminal (or run `vcvarsall.bat` on the terminal before opening VS Code from within it).

#### Zig

You can also bypass the all Visual Studio bullcrap and use zig instead. Just [install zig](https://ziglang.org/learn/getting-started/#installing-zig) (it's basically downloading a binary and adding it to your path) and then run:

```console
> zig build-exe main.c external\upng\upng.c -O ReleaseFast --library c -I.\external\sdl2\include -Lexternal\sdl2\lib\x64 -lSDL2 -lSDL2main -femit-bin=target\main.exe
```

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

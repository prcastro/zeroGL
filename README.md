# zeroGL
A zero-dependencies, single-header 3D graphics library

<img src="./docs/img/rasterizer.gif" width="500">

## Features

* Single-header library written in C with no external dependencies
* Programmable shaders
* Perspective-correct texture support
* Point/Ambient/Directional lighting
* Materials with specular and diffuse properties
* Camera and backface culling
* z-Buffering
* Load OBJ/MTL files (in separate library)

## Getting started

To run the demo, download the latest .zip file from the releases page, uncompress it and run

## Building the project

### All platforms using Zig

The project comes with a distribution for sdl2 under `external/` . On other systems, just install SDL v2.28.2 (using `apt` or `brew`) and find the compiled libraries using `sdl2-config --libs` . To compile, just [install zig](https://ziglang.org/learn/getting-started/#installing-zig) (it's basically downloading a binary and adding it to your path) and then run:

**Windows**
```console
> zig build-exe main.c -O ReleaseFast --library c -I.\external\sdl2\include -Lexternal\sdl2\lib\x64 -lSDL2 -lSDL2main -femit-bin=target\main.exe
```

**linux, macOS** (replace `/usr/local/lib` with the path you found with `sdl2-config --libs`):
```console
$ zig build-exe main.c -O ReleaseFast --library c -I./external/sdl2/include -L/usr/local/lib -lSDL2main -lSDL2 -femit-bin=target/main
```

### Windows

#### CMake
Install Visual Studio Community 2022 (the C/C++ compilers may suffice) and CMake. Then, from the command-line run:
```console
> cmake -B build/ -G "Visual Studio 17 2022"

> cmake --build build/ --config Release
```

The binary will be written to target/Release/main.exe

#### build.bat
Install Visual Studio Community 2022 (the C/C++ compilers may suffice) and then double-click the `build.bat` file. Alternatively, from the command-line run:

```console
> build.bat
```

The binary will be written to target/main.exe


#### Using cl.exe
Install Visual Studio Community 2022 (the C/C++ compilers may suffice), then run:

```console
> "c:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat" x86_amd64

> cl.exe /Zi /EHsc /O2 /nologo /Fotarget\obj\ /Fetarget\main.exe main.c /I external\sdl2\include /link /LIBPATH:external\sdl2\lib\x64 SDL2.lib SDL2main.lib
```

Change the `vcvarsall.bat` path according to your VS installation. The binary will be written to target/Release/main.exe.

#### VS Code
If you're using Visual Studio Code, then you can install Visual Studio Community 2022 (the C/C++ compilers may suffice) and also install the C/C++ extensions from Microsoft.

Open `main.c` on the editor and then click on `Run C/C++ File` on the top right corner. Remember to open VS Code from within a developer terminal (or run `vcvarsall.bat` on the terminal before opening VS Code from within it).

### Linux

On Linux sytems where SDL2 is installed, the exemple main program can also be built by running (again, replace `/usr/local/lib` with the path you found with `sdl2-config --libs`):

```console
gcc main.c -I./external/sdl2/include -L/usr/local/lib -lSDL2main -lSDL2 -lm -o target/main
```

## Current limitations

* No alpha compositing
* No multi-threading/SIMD
* Support only a subset of the OBJ specification

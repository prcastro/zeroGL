@echo off

rem Change this to your Visual Studio 2022 path.
set "VS_PATH=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build"
call "%VS_PATH%\vcvarsall.bat" x86_amd64
cl.exe /Zi /EHsc /O2 /nologo /Fobuild\obj\ /Febuild\main.exe main.c /I external\sdl3\include /link /LIBPATH:external\sdl3\build\RelWithDebInfo SDL3.lib
copy external\sdl3\build\RelWithDebInfo\SDL3.dll build\SDL3.dll

@echo off

rem Change this to your Visual Studio 2022 path.
set "VS_PATH=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build"
call "%VS_PATH%\vcvarsall.bat" x86_amd64
cl.exe /Zi /EHsc /O2 /nologo /Fotarget\obj\ /Fetarget\main.exe main.c /I external\sdl3\include /link /LIBPATH:external\sdl3\lib\x64 SDL3.lib

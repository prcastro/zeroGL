cmake_minimum_required(VERSION 3.10)
project(main C)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build/obj)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build/obj)

add_executable(main main.c)
if (MSVC)
    target_compile_options(main PRIVATE /Zi /EHsc /nologo)
else()
    target_compile_options(main PRIVATE -g -Wall -Wextra)
endif()

# Link Nuklear
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/external/nuklear)

# Link SDL3 (optionally use vendored version)
option(USE_VENDORED_SDL3 "Use vendored SDL3" ON)

if(USE_VENDORED_SDL3)
  add_subdirectory(${CMAKE_SOURCE_DIR}/external/sdl3 EXCLUDE_FROM_ALL)
else()
  find_package(SDL3 REQUIRED CONFIG REQUIRED COMPONENTS SDL3-shared)
endif()

target_link_libraries(main PRIVATE SDL3::SDL3)

# Link STB
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/external/stb)

# Link my libraries in root folder
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR})

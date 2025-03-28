cmake_minimum_required(VERSION 3.10)
project(main C)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/target)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/target/obj)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/target/obj)

add_executable(main main.c)
target_compile_options(main PRIVATE /Zi /EHsc /nologo)

# Link Nuklear
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/external/nuklear)

# Link SDL3
add_subdirectory(${CMAKE_SOURCE_DIR}/external/sdl3 EXCLUDE_FROM_ALL)
target_link_libraries(main PRIVATE SDL3::SDL3)

# Link STB
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/external/stb)

# Link my libraries in root folder
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR})

cmake_minimum_required(VERSION 3.10)
project(main C)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/target)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/target/obj)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/target/obj)

add_executable(main main.c)
target_compile_options(main PRIVATE /Zi /EHsc /nologo)

# Link Nuklear
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/external/nuklear)

# Link SDL2
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/external/sdl2/include)
target_link_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/external/sdl2/lib/x64)
set(SDL2_DLL_PATH "${CMAKE_SOURCE_DIR}/external/sdl2/lib/x64/SDL2.dll")
target_link_libraries(main PRIVATE SDL2 SDL2main)

add_custom_command(TARGET main POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${SDL2_DLL_PATH}"          # Source DLL
    "$<TARGET_FILE_DIR:main>/"  # Destination: same directory as executable
)

# Link UPNG
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/external/stb)

# Link my libraries in root folder
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR})

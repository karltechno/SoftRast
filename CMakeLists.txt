cmake_minimum_required(VERSION 3.2)
project(SoftRast)

add_subdirectory("kt")

set(SOFT_RAST_FILES
    "SoftRast/Platform/Window_Win32.cpp"
    "SoftRast/Platform/Window_Win32.h"
    "SoftRast/Main.cpp"
    "SoftRast/Rasterizer.cpp"
    "SoftRast/Rasterizer.h"
)

include_directories("External")

add_executable(SoftRast ${SOFT_RAST_FILES})
target_link_libraries(SoftRast kt)
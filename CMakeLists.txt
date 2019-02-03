cmake_minimum_required(VERSION 3.2)
project(SoftRast)

if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

add_subdirectory("kt")
add_subdirectory("SoftRast")
add_subdirectory("Viewer")
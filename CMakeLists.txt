cmake_minimum_required(VERSION 3.2)
project(SoftRast)

if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS /Zi)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /DEBUG")
endif()

add_subdirectory("kt")
add_subdirectory("SoftRast")
add_subdirectory("Viewer")
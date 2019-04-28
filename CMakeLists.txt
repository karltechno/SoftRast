cmake_minimum_required(VERSION 3.2)
project(SoftRast)

if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS /Zi)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /DEBUG")
endif()

if((CMAKE_SYSTEM_NAME MATCHES "Windows") AND (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    set(CLANG_ON_WINDOWS TRUE)
else()
    set(CLANG_ON_WINDOWS FALSE)
endif()


if(MSVC OR CLANG_ON_WINDOWS)
    add_compile_options(/arch:AVX2 /fp:fast /Oi)
elseif((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    add_compile_options(-mavx2 -ffast-math)
endif()


add_subdirectory("kt")
add_subdirectory("SoftRast")
add_subdirectory("Viewer")
add_subdirectory("microprofile")
cmake_minimum_required(VERSION 3.2)
project(SoftRast)

add_subdirectory("kt")

if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

set(SOFT_RAST_FILES
    "SoftRast/Platform/Window_Win32.cpp"
    "SoftRast/Platform/Window_Win32.h"
    "SoftRast/Platform/Input_Win32.cpp"
    "SoftRast/Platform/Input_Win32.h"
    "SoftRast/Main.cpp"
    "SoftRast/Rasterizer.cpp"
    "SoftRast/Rasterizer.h"
    "SoftRast/Obj.h"
    "SoftRast/Obj.cpp"
    "SoftRast/Camera.h"
    "SoftRast/Camera.cpp"
    "SoftRast/Input.h"
    "SoftRast/Input.cpp"
	"SoftRast/SoftRastTypes.h"
)

include_directories("External")

add_executable(SoftRast ${SOFT_RAST_FILES})
target_link_libraries(SoftRast kt)
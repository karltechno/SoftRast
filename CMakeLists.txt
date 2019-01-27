cmake_minimum_required(VERSION 3.2)
project(SoftRast)

add_subdirectory("kt")

if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS /Zi)
    SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} /DEBUG")
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
	"SoftRast/Texture.h"
	"SoftRast/Texture.cpp"
	"SoftRast/stb_image.h"
    "SoftRast/stb_image.cpp"
    "SoftRast/Renderer.h"
    "SoftRast/Renderer.cpp"
)

file(COPY "SoftRast/Models/" DESTINATION "${CMAKE_BINARY_DIR}/Models")

add_executable(SoftRast ${SOFT_RAST_FILES})
target_link_libraries(SoftRast kt)
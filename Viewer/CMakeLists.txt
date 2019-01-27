cmake_minimum_required(VERSION 3.2)

set(VIEWER_FILES
    "Platform/Window_Win32.cpp"
    "Platform/Window_Win32.h"
    "Platform/Input_Win32.cpp"
    "Platform/Input_Win32.h"
    "Main.cpp"
    "Obj.h"
    "Obj.cpp"
    "Camera.h"
    "Camera.cpp"
    "Input.h"
    "Input.cpp"
)

file(COPY "Models/" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/Models")

add_executable(Viewer ${VIEWER_FILES})
target_link_libraries(Viewer kt)
target_link_libraries(Viewer SoftRast)
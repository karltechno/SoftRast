set(MICROPROFILE_SOURCES
    "microprofile.h"
    "microprofile.cpp"
    "microprofile_html.h"
    "microprofile.config.h"
)

add_library(microprofile ${MICROPROFILE_SOURCES})
target_include_directories(microprofile PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(microprofile PUBLIC MICROPROFILE_USE_CONFIG MICROPROFILE_GPU_TIMERS=0)

if(MSVC)
target_link_libraries(microprofile ws2_32)
endif()
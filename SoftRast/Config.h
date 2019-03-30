#pragma once
#include <stdint.h>

#define SR_USE_REVERSE_Z (1)

#define SR_DEBUG_SINGLE_THREADED (1)

namespace sr
{


namespace Config
{
// Todo: not hardcoded.
constexpr uint32_t c_screenWidth = 1280;
constexpr uint32_t c_screenHeight = 720;

constexpr uint32_t c_binHeightLog2 = 7;
constexpr uint32_t c_binWidthLog2 = 7;
constexpr uint32_t c_binHeight = 1 << c_binHeightLog2;
constexpr uint32_t c_binWidth = 1 << c_binWidthLog2;

constexpr int32_t c_subPixelBits = 8;
constexpr int32_t c_subPixelStep = 1 << c_subPixelBits;
constexpr int32_t c_subPixelMask = c_subPixelStep - 1;

constexpr uint32_t c_maxVaryings = 16;

constexpr uint32_t c_simdWidth = 8;

constexpr uint32_t c_maxTexDimLog2 = 14; // 16k

#if SR_USE_REVERSE_Z
constexpr float c_depthMin = 1.0f;
constexpr float c_depthMax = 0.0f;
#else
constexpr float c_depthMin = 0.0f;
constexpr float c_depthMax = 1.0f;
#endif

}

}
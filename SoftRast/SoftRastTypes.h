#pragma once
#include <stdint.h>

namespace sr
{

enum class WindingOrder : uint32_t
{
	CCW,
	CW,
	Default = CCW
};


enum class IndexType
{
	u16,
	u32
};

// Todo: move into constants file
namespace Config
{
constexpr uint32_t c_binHeightLog2 = 7;
constexpr uint32_t c_binWidthLog2 = 7;
constexpr uint32_t c_binHeight = 1 << c_binHeightLog2;
constexpr uint32_t c_binWidth = 1 << c_binWidthLog2;

constexpr int32_t c_subPixelBits = 8;
constexpr int32_t c_subPixelStep = 1 << c_subPixelBits;
constexpr int32_t c_subPixelMask = c_subPixelStep - 1;

constexpr uint32_t c_maxVaryings = 16;

constexpr uint32_t c_simdWidth = 8;

// Todo: use these values, also switch to reverse-z
constexpr float c_depthMin = 0.0f;
constexpr float c_depthMax = 1.0f;

}

}
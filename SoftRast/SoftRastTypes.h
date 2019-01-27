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

namespace Config
{
constexpr uint32_t c_tileHeightLog2 = 6;
constexpr uint32_t c_tileWidthLog2 = 6;
constexpr uint32_t c_tileHeight = 1 << c_tileHeightLog2;
constexpr uint32_t c_tileWidth = 1 << c_tileHeightLog2;

constexpr uint32_t c_maxVaryings = 16;

}

}
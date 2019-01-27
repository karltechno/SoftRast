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

constexpr uint32_t c_maxVaryings = 16;


}
#pragma once
#include <stdint.h>
#include <stdio.h>

#include <kt/Array.h>
#include <kt/Serialization.h>

#include "SIMDUtil.h"
#include "SoftRastTypes.h"
#include "Config.h"


namespace sr
{

namespace Tex
{

void CalcMipDims2D(uint32_t _x, uint32_t _y, uint32_t _level, uint32_t o_dims[2]);

struct TextureData
{
	KT_NO_COPY(TextureData);

	TextureData() = default;
	TextureData(TextureData&&) = default;
	TextureData& operator=(TextureData&&) = default;
		
	void CreateFromFile(char const* _file);
	void Clear();

	kt::Array<uint8_t> m_texels;
	uint32_t m_mipOffsets[Config::c_maxTexDimLog2];

	uint32_t m_widthLog2 = 0;
	uint32_t m_heightLog2 = 0;
	uint32_t m_numMips = 0;

	uint32_t m_bytesPerPixel = 0;
};

void SampleWrap(TextureData const& _tex, __m256 _u, __m256 _v, __m256 dudx, __m256 dudy, __m256 dvdx, __m256 dvdy, float o_colour[4 * 8]);

}
}

namespace kt
{

template <>
void Serialize(ISerializer* _s, sr::Tex::TextureData& _tex);

}
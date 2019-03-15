#pragma once
#include <stdint.h>
#include <stdio.h>

#include <kt/Array.h>
#include <kt/Serialization.h>

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

	uint32_t m_width = 0;
	uint32_t m_height = 0;
	uint32_t m_numMips = 0;

	uint32_t m_bitsPerPixel = 0;
	
	uint32_t m_rowStride = 0;
};

void SampleClamp_Slow(TextureData const& _tex, float const _u, float const _v, float o_colour[4]);

void SampleWrap_Slow(TextureData const& _tex, float const _u, float const _v, float o_colour[4]);


}
}

namespace kt
{

template <>
void Serialize(ISerializer* _s, sr::Tex::TextureData& _tex);

}
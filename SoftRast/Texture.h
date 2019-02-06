#pragma once
#include <stdint.h>
#include <stdio.h>


namespace sr
{

namespace Tex
{

struct TextureData
{
	void CreateFromFile(char const* _file);
	void Clear();

	uint8_t* m_data = nullptr;
	uint32_t m_width = 0;
	uint32_t m_height = 0;
	uint32_t m_bitsPerPixel = 0;
	
	uint32_t m_rowStride = 0;
};

void SampleClamp_Slow(TextureData const& _tex, float const _u, float const _v, float o_colour[4]);

void SampleWrap_Slow(TextureData const& _tex, float const _u, float const _v, float o_colour[4]);


}
}
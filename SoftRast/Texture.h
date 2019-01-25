#pragma once
#include <stdint.h>
#include <stdio.h>


namespace sr
{

namespace tex
{

struct TextureData
{
	void CreateFromFile(char const* _file);
	void Clear();

	uint8_t* m_data = nullptr;
	uint32_t m_width = 0;
	uint32_t m_height = 0;
	uint32_t m_bitsPerPixel = 0;
};

}
}
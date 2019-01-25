#include "Texture.h"
#include "kt/Memory.h"
#include "stb_image.h"
#include "kt/Logging.h"

namespace sr
{

namespace Tex
{


void TextureData::CreateFromFile(char const* _file)
{
	Clear();
	static uint32_t const req_comp = 4;
	int x, y, comp;
	uint8_t* mem = stbi_load(_file, &x, &y, &comp, req_comp);
	if (!mem)
	{
		KT_LOG_ERROR("Failed to load texture: %s", _file);
	}

	m_width = x;
	m_height = y;
	m_data = mem;
	m_bitsPerPixel = req_comp * 8;
}

void TextureData::Clear()
{
	kt::Free(m_data);
}

}
}
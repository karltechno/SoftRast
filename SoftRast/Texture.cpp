#include <string.h>

#include <kt/Memory.h>
#include <kt/Logging.h>

#include "Texture.h"
#include "stb_image.h"
#include "stb_image_resize.h"

namespace kt
{

template <>
void Serialize(ISerializer* _s, sr::Tex::TextureData& _tex)
{
	Serialize(_s, _tex.m_texels);
	Serialize(_s, _tex.m_widthLog2);
	Serialize(_s, _tex.m_heightLog2);
	Serialize(_s, _tex.m_bytesPerPixel);
	Serialize(_s, _tex.m_mipOffsets);
	Serialize(_s, _tex.m_numMips);
}

}

namespace sr
{

namespace Tex
{


void TextureData::CreateFromFile(char const* _file)
{
	Clear();
	static uint32_t const req_comp = 4;
	int x, y, comp;
	uint8_t* srcImageData = stbi_load(_file, &x, &y, &comp, req_comp);
	if (!srcImageData)
	{
		KT_LOG_ERROR("Failed to load texture: %s", _file);
	}

	KT_SCOPE_EXIT(stbi_image_free(srcImageData));

	KT_ASSERT(kt::IsPow2(x) && kt::IsPow2(y));

	m_widthLog2 = kt::FloorLog2(uint32_t(x));
	m_heightLog2 = kt::FloorLog2(uint32_t(y));

	KT_ASSERT(m_heightLog2 < Config::c_maxTexDimLog2);
	KT_ASSERT(m_widthLog2 < Config::c_maxTexDimLog2);
	
	m_bytesPerPixel = req_comp;

	// Calculate mips.
	
	uint32_t const fullMipChainLen = kt::FloorLog2(kt::Max(uint32_t(x), uint32_t(y))) + 1; // +1 for base tex.
	m_numMips = fullMipChainLen;
	
	struct MipInfo
	{
		uint32_t m_offs;
		uint32_t m_dims[2];
	};

	uint32_t curMipDataOffset = x * y * req_comp;

	MipInfo* mipInfos = (MipInfo*)KT_ALLOCA(sizeof(MipInfo) * (fullMipChainLen - 1));
	m_mipOffsets[0] = 0;

	for (uint32_t mipIdx = 0; mipIdx < fullMipChainLen - 1; ++mipIdx)
	{
		CalcMipDims2D(uint32_t(x), uint32_t(y), mipIdx + 1, mipInfos[mipIdx].m_dims);
		mipInfos[mipIdx].m_offs = curMipDataOffset;
		m_mipOffsets[mipIdx + 1] = curMipDataOffset;
		curMipDataOffset += mipInfos[mipIdx].m_dims[0] * mipInfos[mipIdx].m_dims[1] * req_comp;
	}

	m_texels.Resize(curMipDataOffset);
	uint8_t* texWritePointer = m_texels.Data();

	memcpy(texWritePointer, srcImageData, req_comp * x * y);

	for (uint32_t mipIdx = 0; mipIdx < fullMipChainLen - 1; ++mipIdx)
	{
		MipInfo const& mipInfo = mipInfos[mipIdx];
		stbir_resize_uint8(srcImageData, x, y, 0, texWritePointer + mipInfo.m_offs, mipInfo.m_dims[0], mipInfo.m_dims[1], 0, req_comp);
	}
}

void TextureData::Clear()
{
	m_texels.ClearAndFree();
}

void CalcMipDims2D(uint32_t _x, uint32_t _y, uint32_t _level, uint32_t o_dims[2])
{
	o_dims[0] = kt::Max<uint32_t>(1u, _x >> _level);
	o_dims[1] = kt::Max<uint32_t>(1u, _y >> _level);
}

void SampleClamp_Slow(TextureData const& _tex, uint32_t const _mipIdx, float const _u, float const _v, float o_colour[4])
{
	uint32_t const mipClamped = kt::Min(_mipIdx, _tex.m_numMips);

	uint32_t const width	= 1u << (kt::Min(_tex.m_widthLog2, mipClamped) - mipClamped);
	uint32_t const height	= 1u << (kt::Min(_tex.m_heightLog2, mipClamped) - mipClamped);

	uint32_t const pitch = width * _tex.m_bytesPerPixel;

	uint32_t const clampU = uint32_t(kt::Clamp<int32_t>(int32_t(_u * width), 0, width - 1));
	uint32_t const clampV = uint32_t(kt::Clamp<int32_t>(int32_t(_v * height), 0, height - 1));

	uint32_t const offs = clampV * pitch + clampU * _tex.m_bytesPerPixel;
	uint8_t const* pix = _tex.m_texels.Data() + (_tex.m_mipOffsets[mipClamped] + offs);
	static const float recip255 = 1.0f / 255.0f;
	o_colour[0] = pix[0] * recip255;
	o_colour[1] = pix[1] * recip255;
	o_colour[2] = pix[2] * recip255;
	o_colour[3] = pix[3] * recip255;
}

void SampleWrap_Slow(TextureData const& _tex, uint32_t const _mipIdx, float const _u, float const _v, float o_colour[4])
{
	uint32_t const mipClamped = kt::Min(_mipIdx, _tex.m_numMips);

	uint32_t const width = 1u << (_tex.m_widthLog2 - mipClamped);
	uint32_t const height = 1u << (_tex.m_heightLog2 - mipClamped);

	uint32_t const pitch = width * _tex.m_bytesPerPixel;

	float const uSign = _u < 0.0f ? -1.0f : 1.0f;
	float const vSign = _v < 0.0f ? -1.0f : 1.0f;

	float const absU = uSign * _u;
	float const absV = vSign * _v;

	float const fracU = absU - int32_t(absU);
	float const fracV = absV - int32_t(absV);

	float const uWrap = uSign < 0.0f ? (1.0f - fracU) : fracU;
	float const vWrap = vSign < 0.0f ? (1.0f - fracV) : fracV;

	uint32_t const clampU = uint32_t(kt::Clamp<int32_t>(int32_t(uWrap * width), 0, width - 1));
	uint32_t const clampV = uint32_t(kt::Clamp<int32_t>(int32_t(vWrap * height), 0, height - 1));

	uint32_t const offs = clampV * pitch + clampU * _tex.m_bytesPerPixel;
	uint8_t const* pix = _tex.m_texels.Data() + (_tex.m_mipOffsets[mipClamped] + offs);
	static const float recip255 = 1.0f / 255.0f;
	o_colour[0] = pix[0] * recip255;
	o_colour[1] = pix[1] * recip255;
	o_colour[2] = pix[2] * recip255;
	o_colour[3] = pix[3] * recip255;
}

}
}
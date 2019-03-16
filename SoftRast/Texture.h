#pragma once
#include <stdint.h>
#include <stdio.h>

#include <kt/Array.h>
#include <kt/Serialization.h>

#include "SIMDUtil.h"
#include "SoftRastTypes.h"

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

void SampleClamp_Slow(TextureData const& _tex, uint32_t const _mipIdx, float const _u, float const _v, float o_colour[4]);

void SampleWrap_Slow(TextureData const& _tex, uint32_t const _mipIdx, float const _u, float const _v, float o_colour[4]);

inline void SampleWrap(TextureData const& _tex, __m256 _u, __m256 _v,  float o_colour[4 * 8])
{
	uint32_t const mipClamped = 0;

	__m256i const mipFloor = _mm256_min_epi32(_mm256_set1_epi32(_tex.m_numMips - 1), _mm256_set1_epi32(mipClamped));

	__m256i const one = _mm256_set1_epi32(1);

	__m256i const widthLog2 = _mm256_set1_epi32(_tex.m_widthLog2);
	__m256i const heightLog2 = _mm256_set1_epi32(_tex.m_heightLog2);

	__m256i const width = _mm256_sllv_epi32(one, _mm256_sub_epi32(widthLog2, _mm256_min_epi32(widthLog2, mipFloor)));
	__m256i const height = _mm256_sllv_epi32(one, _mm256_sub_epi32(heightLog2, _mm256_min_epi32(heightLog2, mipFloor)));

	__m256i const pitch = _mm256_mullo_epi32(width, _mm256_set1_epi32(_tex.m_bytesPerPixel)); // todo: always four?

	__m256 const signBit = SR_AVX_LOAD_CONST_FLOAT(simdutil::c_avxSignBit);

	__m256 const uSign = _mm256_and_ps(signBit, _u);
	__m256 const vSign = _mm256_and_ps(signBit, _v);

	__m256 const absU = _mm256_xor_ps(uSign, _u);
	__m256 const absV = _mm256_xor_ps(vSign, _v);

	// Todo: naive casting faster than floor (roundps) ?
	__m256 const fracU = _mm256_sub_ps(absU, _mm256_floor_ps(absU));
	__m256 const fracV = _mm256_sub_ps(absV, _mm256_floor_ps(absV));

	__m256 const uWrap = _mm256_blendv_ps(fracU, _mm256_sub_ps(_mm256_set1_ps(1.0f), fracU), uSign);
	__m256 const vWrap = _mm256_blendv_ps(fracV, _mm256_sub_ps(_mm256_set1_ps(1.0f), fracV), vSign);

	__m256 const widthF = _mm256_cvtepi32_ps(width);
	__m256i const widthMinusOne = _mm256_sub_epi32(width, one);

	__m256 const heightF = _mm256_cvtepi32_ps(height);
	__m256i const heightMinusOne = _mm256_sub_epi32(height, one);

	__m256i const clampU = _mm256_min_epi32(widthMinusOne, _mm256_max_epi32(
		_mm256_setzero_si256(), _mm256_cvtps_epi32(_mm256_mul_ps(widthF, uWrap))));

	__m256i const clampV = _mm256_min_epi32(heightMinusOne, _mm256_max_epi32(
		_mm256_setzero_si256(), _mm256_cvtps_epi32(_mm256_mul_ps(heightF, vWrap))));

	__m256i const offs = _mm256_add_epi32(_mm256_mullo_epi32(pitch, clampV), _mm256_mullo_epi32(clampU, _mm256_set1_epi32(_tex.m_bytesPerPixel)));

	KT_ALIGNAS(32) uint32_t offsArr[8];
	KT_ALIGNAS(32) uint32_t mips[8];

	_mm256_store_si256((__m256i*)offsArr, offs);
	_mm256_store_si256((__m256i*)mips, mipFloor);
	static const float recip255 = 1.0f / 255.0f;
	
	for (uint32_t i = 0; i < 8; ++i)
	{
		uint8_t const* pix = _tex.m_texels.Data() + (offsArr[i] + _tex.m_mipOffsets[mips[i]]); // todo: mip offset
		o_colour[0 + i * 4] = pix[0] * recip255;
		o_colour[1 + i * 4] = pix[1] * recip255;
		o_colour[2 + i * 4] = pix[2] * recip255;
		o_colour[3 + i * 4] = pix[3] * recip255;
	}
}

}
}

namespace kt
{

template <>
void Serialize(ISerializer* _s, sr::Tex::TextureData& _tex);

}
#include <string.h>
#include <intrin.h>

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

static uint32_t MortonEncode(uint32_t _x, uint32_t _y)
{
	constexpr uint32_t pdep_x_mask = 0x55555555; // 0b010101 ...
	constexpr uint32_t pdep_y_mask = 0xAAAAAAAA; // 0b101010 ...
	return _pdep_u32(_x, pdep_x_mask) | _pdep_u32(_y, pdep_y_mask);
}

constexpr uint32_t c_texTileSizeLog2 = 5;
constexpr uint32_t c_texTileSize = 1 << c_texTileSizeLog2;
constexpr uint32_t c_texTileMask = c_texTileSize - 1;

static void TileTexture(uint8_t const* _src, uint8_t* _dest, uint32_t const dimX_noPad, uint32_t const dimY_noPad)
{
	uint32_t const mipTileWidth = uint32_t(kt::AlignUp(dimX_noPad, c_texTileSize)) >> c_texTileSizeLog2;

	for (uint32_t yy = 0; yy < dimY_noPad; ++yy)
	{
		for (uint32_t xx = 0; xx < dimX_noPad; ++xx)
		{
			uint32_t const linearOffs = yy * dimX_noPad + xx;

			uint32_t const tileX = xx >> c_texTileSizeLog2;
			uint32_t const tileY = yy >> c_texTileSizeLog2;

			uint32_t const inTileAddressX = xx & c_texTileMask;
			uint32_t const inTileAddressY = yy & c_texTileMask;

			uint32_t const morton = MortonEncode(inTileAddressX, inTileAddressY);

			uint32_t const tiledOffs = (tileY * mipTileWidth + tileX) * (c_texTileSize * c_texTileSize) + morton;
			memcpy(_dest + tiledOffs * 4, _src + linearOffs * 4, 4);
		}
	}
}

void TextureData::CreateFromFile(char const* _file)
{
	Clear();
	static uint32_t const req_comp = 4;
	int x, y, comp;
	uint8_t* srcImageData = stbi_load(_file, &x, &y, &comp, req_comp);
	if (!srcImageData)
	{
		KT_LOG_ERROR("Failed to load texture: %s", _file);
		return;
	}

	KT_SCOPE_EXIT(stbi_image_free(srcImageData));

	KT_ASSERT(kt::IsPow2(x) && kt::IsPow2(y));

	m_widthLog2 = kt::FloorLog2(uint32_t(x));
	m_heightLog2 = kt::FloorLog2(uint32_t(y));

	KT_ASSERT(m_heightLog2 < Config::c_maxTexDimLog2);
	KT_ASSERT(m_widthLog2 < Config::c_maxTexDimLog2);
	KT_ASSERT((x % c_texTileSize) == 0);
	KT_ASSERT((y % c_texTileSize) == 0); // TODO: Pad if this isn't true (but will be with 2^x x>=5)
	
	m_bytesPerPixel = req_comp;

	// Calculate mips.
	
	uint32_t const fullMipChainLen = kt::FloorLog2(kt::Max(uint32_t(x), uint32_t(y))) + 1; // +1 for base tex.
	m_numMips = fullMipChainLen;
	
	struct MipInfo
	{
		uint32_t m_offs;
		uint32_t m_dims[2];
	};

	uint32_t curMipDataOffset = 0;

	MipInfo* mipInfos = (MipInfo*)KT_ALLOCA(sizeof(MipInfo) * fullMipChainLen);

	for (uint32_t mipIdx = 0; mipIdx < fullMipChainLen; ++mipIdx)
	{
		CalcMipDims2D(uint32_t(x), uint32_t(y), mipIdx, mipInfos[mipIdx].m_dims);
		mipInfos[mipIdx].m_offs = curMipDataOffset;
		m_mipOffsets[mipIdx] = curMipDataOffset;

		// Align the offset to account for tiling
		uint32_t const mipDimX_tilePad = uint32_t(kt::AlignUp(mipInfos[mipIdx].m_dims[0], c_texTileSize));
		uint32_t const mipDimY_tilePad = uint32_t(kt::AlignUp(mipInfos[mipIdx].m_dims[1], c_texTileSize));

		curMipDataOffset += mipDimX_tilePad * mipDimY_tilePad * req_comp;
	}


	m_texels.Resize(curMipDataOffset);
	uint8_t* texWritePointer = m_texels.Data();

	// tile mip 0
	TileTexture(srcImageData, texWritePointer, mipInfos[0].m_dims[0], mipInfos[0].m_dims[1]);

	uint32_t const largestMipSize = x * y * 4;
	uint8_t* tempResizeBuff = (uint8_t*)kt::Malloc(largestMipSize);
	KT_SCOPE_EXIT(kt::Free(tempResizeBuff));

	for (uint32_t mipIdx = 1; mipIdx < fullMipChainLen; ++mipIdx)
	{
		MipInfo const& mipInfo = mipInfos[mipIdx];
		uint8_t* mipPtr = texWritePointer + mipInfo.m_offs;

		uint32_t const mipDimX = mipInfo.m_dims[0];
		uint32_t const mipDimY = mipInfo.m_dims[1];

		stbir_resize_uint8(srcImageData, x, y, 0, tempResizeBuff, mipDimX, mipDimY, 0, req_comp);
		TileTexture(tempResizeBuff, mipPtr, mipDimX, mipDimY);
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

__m256i CalcMipLevels(TextureData const& _tex, __m256 _dudx, __m256 _dudy, __m256 _dvdx, __m256 _dvdy)
{
	__m256 const height = _mm256_set1_ps(float(1u << _tex.m_heightLog2));
	__m256 const width = _mm256_set1_ps(float(1u << _tex.m_widthLog2));
	
	__m256 const dudx_tex = _mm256_mul_ps(_dudx, width);
	__m256 const dudy_tex = _mm256_mul_ps(_dudy, height);

	__m256 const dvdx_tex = _mm256_mul_ps(_dvdx, width);
	__m256 const dvdy_tex = _mm256_mul_ps(_dvdy, height);

	// inner product
	__m256 const du_dot2 = _mm256_fmadd_ps(dudx_tex, dudx_tex, _mm256_mul_ps(dudy_tex, dudy_tex));
	__m256 const dv_dot2 = _mm256_fmadd_ps(dvdx_tex, dvdx_tex, _mm256_mul_ps(dvdy_tex, dvdy_tex));

	// Todo: with proper log2 we can use identity log2(x^(1/2)) == 0.5 * log2(x) and remove sqrt.
	__m256 const maxCoord = _mm256_sqrt_ps(_mm256_max_ps(du_dot2, dv_dot2));

	return _mm256_min_epi32(_mm256_set1_epi32(_tex.m_numMips - 1), _mm256_max_epi32(_mm256_setzero_si256(), simdutil::ExtractExponent(maxCoord)));
}


static __m256i BoundCoordsWrap(__m256i _coord, __m256i _bound)
{
	// Assuming width and height are powers of two.
	__m256i const one = _mm256_set1_epi32(1);
	return _mm256_and_si256(_coord, _mm256_sub_epi32(_bound, one));
}

static void GatherQuads
(
	TextureData const& _tex, 
	uint32_t _mips[8], 
	__m256i _mipWidth, 
	__m256i _x0, 
	__m256i _y0, 
	__m256i _x1,
	__m256i _y1,
	float o_x0y0[8 * 4],
	float o_x1y0[8 * 4],
	float o_x0y1[8 * 4],
	float o_x1y1[8 * 4]
)
{
	__m256i const tileX0 = _mm256_srli_epi32(_x0, c_texTileSizeLog2);
	__m256i const tileY0 = _mm256_srli_epi32(_y0, c_texTileSizeLog2);
	__m256i const tileX1 = _mm256_srli_epi32(_x1, c_texTileSizeLog2);
	__m256i const tileY1 = _mm256_srli_epi32(_y1, c_texTileSizeLog2);

	__m256i const c_texTileMaskAvx = _mm256_set1_epi32(c_texTileMask);

	__m256i const inTileAddressX0 = _mm256_and_si256(_x0, c_texTileMaskAvx);
	__m256i const inTileAddressY0 = _mm256_and_si256(_y0, c_texTileMaskAvx);
	__m256i const inTileAddressX1 = _mm256_and_si256(_x1, c_texTileMaskAvx);
	__m256i const inTileAddressY1 = _mm256_and_si256(_y1, c_texTileMaskAvx);

	// TODO: Broken for non pow2 textures (we assert not supporting those though!)
	__m256i const texTileSize = _mm256_set1_epi32(c_texTileSize);
	__m256i const mipTileWidth = _mm256_srli_epi32(_mm256_max_epi32(_mipWidth, texTileSize), c_texTileSizeLog2);

	KT_ALIGNAS(32) uint32_t offs_x0y0[8];
	KT_ALIGNAS(32) uint32_t offs_x1y0[8];
	KT_ALIGNAS(32) uint32_t offs_x0y1[8];
	KT_ALIGNAS(32) uint32_t offs_x1y1[8];

	_mm256_store_si256((__m256i*)offs_x0y0, _mm256_mullo_epi32(_mm256_set1_epi32(c_texTileSize * c_texTileSize), _mm256_add_epi32(_mm256_mullo_epi32(tileY0, mipTileWidth), tileX0)));
	_mm256_store_si256((__m256i*)offs_x1y0, _mm256_mullo_epi32(_mm256_set1_epi32(c_texTileSize * c_texTileSize), _mm256_add_epi32(_mm256_mullo_epi32(tileY0, mipTileWidth), tileX1)));
	_mm256_store_si256((__m256i*)offs_x0y1, _mm256_mullo_epi32(_mm256_set1_epi32(c_texTileSize * c_texTileSize), _mm256_add_epi32(_mm256_mullo_epi32(tileY1, mipTileWidth), tileX0)));
	_mm256_store_si256((__m256i*)offs_x1y1, _mm256_mullo_epi32(_mm256_set1_epi32(c_texTileSize * c_texTileSize), _mm256_add_epi32(_mm256_mullo_epi32(tileY1, mipTileWidth), tileX1)));

	KT_ALIGNAS(32) uint32_t linear_x0[8];
	KT_ALIGNAS(32) uint32_t linear_y0[8];
	KT_ALIGNAS(32) uint32_t linear_x1[8];
	KT_ALIGNAS(32) uint32_t linear_y1[8];
	_mm256_store_si256((__m256i*)linear_x0, inTileAddressX0);
	_mm256_store_si256((__m256i*)linear_y0, inTileAddressY0);
	_mm256_store_si256((__m256i*)linear_x1, inTileAddressX1);
	_mm256_store_si256((__m256i*)linear_y1, inTileAddressY1);

	for (uint32_t i = 0; i < 8; ++i)
	{
		uint8_t const* mipPtr = _tex.m_texels.Data() + _tex.m_mipOffsets[_mips[i]];

		{
			uint8_t const* pix_x0y0 = mipPtr + 4 * (offs_x0y0[i] + MortonEncode(linear_x0[i], linear_y0[i]));
			__m128i const x0y0 = _mm_set1_epi32(*(uint32_t*)pix_x0y0);
			_mm_storeu_ps(o_x0y0 + i * 4, _mm_mul_ps(_mm_set1_ps(1.0f / 255.0f), _mm_cvtepi32_ps(_mm_cvtepu8_epi32(x0y0))));
		}

		{
			uint8_t const* pix_x1y0 = mipPtr + 4 * (offs_x1y0[i] + MortonEncode(linear_x1[i], linear_y0[i]));
			__m128i const x1y0 = _mm_set1_epi32(*(uint32_t*)pix_x1y0);
			_mm_storeu_ps(o_x1y0 + i * 4, _mm_mul_ps(_mm_set1_ps(1.0f / 255.0f), _mm_cvtepi32_ps(_mm_cvtepu8_epi32(x1y0))));
		}

		{
			uint8_t const* pix_x1y1 = mipPtr + 4 * (offs_x1y1[i] + MortonEncode(linear_x1[i], linear_y1[i]));
			__m128i const x1y1 = _mm_set1_epi32(*(uint32_t*)pix_x1y1);
			_mm_storeu_ps(o_x1y1 + i * 4, _mm_mul_ps(_mm_set1_ps(1.0f / 255.0f), _mm_cvtepi32_ps(_mm_cvtepu8_epi32(x1y1))));
		}

		{
			uint8_t const* pix_x0y1 = mipPtr + 4 * (offs_x0y1[i] + MortonEncode(linear_x0[i], linear_y1[i]));
			__m128i const x0y1 = _mm_set1_epi32(*(uint32_t*)pix_x0y1);
			_mm_storeu_ps(o_x0y1 + i * 4, _mm_mul_ps(_mm_set1_ps(1.0f / 255.0f), _mm_cvtepi32_ps(_mm_cvtepu8_epi32(x0y1))));
		}
	}
}

void SampleWrap
(
	TextureData const& _tex,
	__m256 _u,
	__m256 _v,
	__m256 _dudx,
	__m256 _dudy,
	__m256 _dvdx,
	__m256 _dvdy,
	float o_colour[4 * 8]
)
{
	__m256i const mipFloor = CalcMipLevels(_tex, _dudx, _dudy, _dvdx, _dvdy);

	__m256i const one = _mm256_set1_epi32(1);

	__m256i const widthLog2 = _mm256_set1_epi32(_tex.m_widthLog2);
	__m256i const heightLog2 = _mm256_set1_epi32(_tex.m_heightLog2);

	__m256i const width = _mm256_sllv_epi32(one, _mm256_sub_epi32(widthLog2, _mm256_min_epi32(widthLog2, mipFloor)));
	__m256i const height = _mm256_sllv_epi32(one, _mm256_sub_epi32(heightLog2, _mm256_min_epi32(heightLog2, mipFloor)));

	__m256 const signBit = SR_AVX_LOAD_CONST_FLOAT(simdutil::c_avxSignBit);

	__m256 const uSign = _mm256_and_ps(signBit, _u);
	__m256 const vSign = _mm256_and_ps(signBit, _v);

	__m256 const absU = _mm256_xor_ps(uSign, _u);
	__m256 const absV = _mm256_xor_ps(vSign, _v);

	__m256 const fracU = _mm256_sub_ps(absU, _mm256_floor_ps(absU));
	__m256 const fracV = _mm256_sub_ps(absV, _mm256_floor_ps(absV));

	__m256 const fracU_wrap = _mm256_blendv_ps(fracU, _mm256_sub_ps(_mm256_set1_ps(1.0f), fracU), uSign);
	__m256 const fracV_wrap = _mm256_blendv_ps(fracV, _mm256_sub_ps(_mm256_set1_ps(1.0f), fracV), vSign);

	__m256 const widthF = _mm256_cvtepi32_ps(width);
	__m256 const heightF = _mm256_cvtepi32_ps(height);

	__m256 const u_texSpace = _mm256_mul_ps(widthF, fracU_wrap);
	__m256 const v_texSpace = _mm256_mul_ps(heightF, fracV_wrap);

	__m256 const u_texSpace_floor = _mm256_floor_ps(u_texSpace);
	__m256 const v_texSpace_floor = _mm256_floor_ps(v_texSpace);

	__m256 const u_interp = _mm256_sub_ps(u_texSpace, u_texSpace_floor);
	__m256 const v_interp = _mm256_sub_ps(v_texSpace, v_texSpace_floor);

	__m256i const x0 = BoundCoordsWrap(_mm256_cvtps_epi32(u_texSpace_floor), width);
	__m256i const y0 = BoundCoordsWrap(_mm256_cvtps_epi32(v_texSpace_floor), height);

	__m256i const x1 = BoundCoordsWrap(_mm256_add_epi32(x0, one), width);
	__m256i const y1 = BoundCoordsWrap(_mm256_add_epi32(y0, one), height);

	KT_ALIGNAS(32) uint32_t mips[8];
	_mm256_store_si256((__m256i*)mips, mipFloor);

	KT_ALIGNAS(32) float x0y0_gather[8 * 4];
	KT_ALIGNAS(32) float x0y1_gather[8 * 4];
	KT_ALIGNAS(32) float x1y0_gather[8 * 4];
	KT_ALIGNAS(32) float x1y1_gather[8 * 4];

	GatherQuads(_tex, mips, width, x0, y0, x1, y1, x0y0_gather, x1y0_gather, x0y1_gather, x1y1_gather);

	static const __m256i permuteMaskBase = _mm256_setr_epi32(0, 0, 0, 0, 1, 1, 1, 1);

	__m256i permuteMask = permuteMaskBase;
	__m256i const two = _mm256_set1_epi32(2);

	for (uint32_t i = 0; i < 4; ++i)
	{
		__m256 const interpV_perm = _mm256_permutevar8x32_ps(v_interp, permuteMask);
		__m256 const interpU_perm = _mm256_permutevar8x32_ps(u_interp, permuteMask);

		permuteMask = _mm256_add_epi32(two, permuteMask);

		__m256 const x0y0 = _mm256_load_ps(x0y0_gather + i * 8);
		__m256 const x0y1 = _mm256_load_ps(x0y1_gather + i * 8);

		__m256 const left = _mm256_fmadd_ps(interpV_perm, x0y1, _mm256_fnmadd_ps(interpV_perm, x0y0, x0y0));

		__m256 const x1y0 = _mm256_load_ps(x1y0_gather + i * 8);
		__m256 const x1y1 = _mm256_load_ps(x1y1_gather + i * 8);

		__m256 const right = _mm256_fmadd_ps(interpV_perm, x1y1, _mm256_fnmadd_ps(interpV_perm, x1y0, x1y0));

		_mm256_store_ps(o_colour + i * 8, _mm256_fmadd_ps(interpU_perm, right, _mm256_fnmadd_ps(interpU_perm, left, left)));
	}
}

}
}
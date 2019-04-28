#pragma once
#include <kt/kt.h>
#include "SIMDUtil.h"
#include "Obj.h"

#include <immintrin.h>
#include <stdint.h>

namespace sr
{
namespace shader
{

uint32_t constexpr c_objVertexStride = sizeof(sr::Obj::Vertex) / sizeof(float) + 4;

struct Derivatives
{
	__m256 dudx;
	__m256 dudy;
	__m256 dvdx;
	__m256 dvdy;
};

struct OBJVaryings
{
	__m256 pos_x;
	__m256 pos_y;
	__m256 pos_z;
	__m256 norm_x;
	__m256 norm_y;
	__m256 norm_z;
	__m256 u;
	__m256 v;
};

KT_FORCEINLINE Derivatives UnpackDerivatives(float const* _varyings, uint32_t _stride = c_objVertexStride)
{
	Derivatives ret;
	// load as rows and transpose sub 4x4 matricies.
	// [dudx0 dudy0 dvdx0 dvdy0 | dudx4 dudy4 dvdx4 dvdy4] 
	// [dudx1 dudy1 dvdx1 dvdy1 | dudx5 dudy5 dvdx5 dvdy5] 
	// [dudx2 dudy2 dvdx2 dvdy2 | dudx6 dudy6 dvdx6 dvdy6] 
	// [dudx3 dudy3 dvdx3 dvdy3 | dudx7 dudy7 dvdx7 dvdy7] 
	ret.dudx = _mm256_loadu2_m128(_varyings + _stride * 4, _varyings + _stride * 0);
	ret.dudy = _mm256_loadu2_m128(_varyings + _stride * 5, _varyings + _stride * 1);
	ret.dvdx = _mm256_loadu2_m128(_varyings + _stride * 6, _varyings + _stride * 2);
	ret.dvdy = _mm256_loadu2_m128(_varyings + _stride * 7, _varyings + _stride * 3);

	sr::simdutil::Transpose4x4SubMatricies(ret.dudx, ret.dudy, ret.dvdx, ret.dvdy);
	return ret;
}

KT_FORCEINLINE OBJVaryings UnpackOBJVaryings(float const* _varyings)
{
	OBJVaryings ret;
	// offset for derivatives
	float const* vertexVaryingBegin = _varyings + 4;
	ret.pos_x = _mm256_loadu_ps(vertexVaryingBegin);
	ret.pos_y = _mm256_loadu_ps(vertexVaryingBegin + c_objVertexStride);
	ret.pos_z = _mm256_loadu_ps(vertexVaryingBegin + c_objVertexStride * 2);
	ret.norm_x = _mm256_loadu_ps(vertexVaryingBegin + c_objVertexStride * 3);
	ret.norm_y = _mm256_loadu_ps(vertexVaryingBegin + c_objVertexStride * 4);
	ret.norm_z = _mm256_loadu_ps(vertexVaryingBegin + c_objVertexStride * 5);
	ret.u = _mm256_loadu_ps(vertexVaryingBegin + c_objVertexStride * 6);
	ret.v = _mm256_loadu_ps(vertexVaryingBegin + c_objVertexStride * 7);

	sr::simdutil::Transpose8x8(ret.pos_x, ret.pos_y, ret.pos_z, ret.norm_x, ret.norm_y, ret.norm_z, ret.u, ret.v);
	return ret;
}

KT_FORCEINLINE void UnlitDiffuseShader(void const* _uniforms, float const* _varyings, uint32_t o_texels[8], uint32_t _execMask)
{
	sr::Tex::TextureData* tex = (sr::Tex::TextureData*)_uniforms;

	if (!tex || tex->m_texels.Size() == 0)
	{
		memset(o_texels, 0xFFFFFFFF, 8 * sizeof(uint32_t));
		return;
	}

	OBJVaryings const objVaryings = UnpackOBJVaryings(_varyings);
	Derivatives const derivs = UnpackDerivatives(_varyings);

	__m256 r;
	__m256 g;
	__m256 b;
	__m256 a;

	sr::Tex::SampleWrap(*tex, objVaryings.u, objVaryings.v, derivs.dudx, derivs.dudy, derivs.dvdx, derivs.dvdy, r, g, b, a, _execMask);
	sr::simdutil::RGBA32SoA_To_RGBA8AoS(r, g, b, a, o_texels);
}

KT_FORCEINLINE void VisualizeNormalsShader(void const* _uniforms, float const* _varyings, uint32_t o_texels[8], uint32_t _execMask)
{
	OBJVaryings const objVaryings = UnpackOBJVaryings(_varyings);

	__m256 const mulAndAdd = _mm256_set1_ps(0.5f);

	__m256 const r = _mm256_fmadd_ps(objVaryings.norm_x, mulAndAdd, mulAndAdd);
	__m256 const g = _mm256_fmadd_ps(objVaryings.norm_y, mulAndAdd, mulAndAdd);
	__m256 const b = _mm256_fmadd_ps(objVaryings.norm_z, mulAndAdd, mulAndAdd);
	__m256 const a = _mm256_set1_ps(1.0f);
	sr::simdutil::RGBA32SoA_To_RGBA8AoS(r, g, b, a, o_texels);
}

}

}
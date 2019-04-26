#pragma once
#include <kt/kt.h>
#include <immintrin.h>

#define SR_AVX_CONST1_INT(_name, _x)	KT_ALIGNAS(32) static int32_t const _name[8] = {_x, _x, _x, _x, _x, _x, _x, _x};
#define SR_AVX_CONST1_UINT(_name, _x)	KT_ALIGNAS(32) static uint32_t const _name[8] = {_x, _x, _x, _x, _x, _x, _x, _x};
#define SR_AVX_CONST1_FLOAT(_name, _x)	KT_ALIGNAS(32) static float const _name[8] = {_x, _x, _x, _x, _x, _x, _x, _x};


#define SR_AVX_LOAD_CONST_FLOAT(_const) (_mm256_load_ps((float*)_const))
#define SR_AVX_LOAD_CONST_INT(_const) (_mm256_load_si256((__m256i*)_const))


namespace sr
{

namespace simdutil
{

SR_AVX_CONST1_UINT(c_avxSignBit, 0x80000000);
SR_AVX_CONST1_UINT(c_avxSignMask, 0x7FFFFFFF);
SR_AVX_CONST1_UINT(c_avxExponentMask, 0x000000FF);

KT_FORCEINLINE __m256i ExtractExponent(__m256 _v)
{
	__m256i asuint = _mm256_castps_si256(_v);
	asuint = _mm256_srai_epi32(asuint, 23);
	asuint = _mm256_and_si256(asuint, SR_AVX_LOAD_CONST_INT(c_avxExponentMask));
	return _mm256_sub_epi32(asuint, _mm256_set1_epi32(127));
}

KT_FORCEINLINE void Transpose8x8(__m256& _r0, __m256& _r1, __m256& _r2, __m256& _r3, __m256& _r4, __m256& _r5, __m256& _r6, __m256& _r7)
{
	// Reference:
	// pg 423+ https://www.intel.com/content/dam/www/public/us/en/documents/manuals/64-ia-32-architectures-optimization-manual.pdf
	// https://stackoverflow.com/questions/25622745/transpose-an-8x8-float-using-avx-avx2
	__m256 t0 = _mm256_unpacklo_ps(_r0, _r1);
	__m256 t1 = _mm256_unpackhi_ps(_r0, _r1);
	__m256 t2 = _mm256_unpacklo_ps(_r2, _r3);
	__m256 t3 = _mm256_unpackhi_ps(_r2, _r3);
	__m256 t4 = _mm256_unpacklo_ps(_r4, _r5);
	__m256 t5 = _mm256_unpackhi_ps(_r4, _r5);
	__m256 t6 = _mm256_unpacklo_ps(_r6, _r7);
	__m256 t7 = _mm256_unpackhi_ps(_r6, _r7);

	_r0 = _mm256_shuffle_ps(t0, t2, _MM_SHUFFLE(1, 0, 1, 0));
	_r1 = _mm256_shuffle_ps(t0, t2, _MM_SHUFFLE(3, 2, 3, 2));
	_r2 = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(1, 0, 1, 0));
	_r3 = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(3, 2, 3, 2));
	_r4 = _mm256_shuffle_ps(t4, t6, _MM_SHUFFLE(1, 0, 1, 0));
	_r5 = _mm256_shuffle_ps(t4, t6, _MM_SHUFFLE(3, 2, 3, 2));
	_r6 = _mm256_shuffle_ps(t5, t7, _MM_SHUFFLE(1, 0, 1, 0));
	_r7 = _mm256_shuffle_ps(t5, t7, _MM_SHUFFLE(3, 2, 3, 2));

	t0 = _mm256_permute2f128_ps(_r0, _r4, 0x20);
	t1 = _mm256_permute2f128_ps(_r1, _r5, 0x20);
	t2 = _mm256_permute2f128_ps(_r2, _r6, 0x20);
	t3 = _mm256_permute2f128_ps(_r3, _r7, 0x20);
	t4 = _mm256_permute2f128_ps(_r0, _r4, 0x31);
	t5 = _mm256_permute2f128_ps(_r1, _r5, 0x31);
	t6 = _mm256_permute2f128_ps(_r2, _r6, 0x31);
	t7 = _mm256_permute2f128_ps(_r3, _r7, 0x31);

	_r0 = t0;
	_r1 = t1;
	_r2 = t2;
	_r3 = t3;
	_r4 = t4;
	_r5 = t5;
	_r6 = t6;
	_r7 = t7;
}

KT_FORCEINLINE void Transpose4x4SubMatricies(__m256& _r0, __m256& _r1, __m256& _r2, __m256& _r3)
{
	__m256 const x0x1y0y1 = _mm256_unpacklo_ps(_r0, _r1);
	__m256 const z0z1w0w1 = _mm256_unpackhi_ps(_r0, _r1);
	__m256 const x2x3y2y3 = _mm256_unpacklo_ps(_r2, _r3);
	__m256 const z2z3w2w3 = _mm256_unpackhi_ps(_r2, _r3);

	_r0 = _mm256_shuffle_ps(x0x1y0y1, x2x3y2y3, _MM_SHUFFLE(1, 0, 1, 0));
	_r1 = _mm256_shuffle_ps(x0x1y0y1, x2x3y2y3, _MM_SHUFFLE(3, 2, 3, 2));
	_r2 = _mm256_shuffle_ps(z0z1w0w1, z2z3w2w3, _MM_SHUFFLE(1, 0, 1, 0));
	_r3 = _mm256_shuffle_ps(z0z1w0w1, z2z3w2w3, _MM_SHUFFLE(3, 2, 3, 2));
}

KT_FORCEINLINE void RGBA32SoA_To_RGBA8AoS(__m256 _r, __m256 _g, __m256 _b, __m256 _a, uint32_t* o_texels)
{
	__m256 const mul = _mm256_set1_ps(255.0f);
	__m256 const add = _mm256_set1_ps(0.5f);

	__m256i const ri32 = _mm256_cvtps_epi32(_mm256_fmadd_ps(_r, mul, add));
	__m256i const gi32 = _mm256_cvtps_epi32(_mm256_fmadd_ps(_g, mul, add));
	__m256i const bi32 = _mm256_cvtps_epi32(_mm256_fmadd_ps(_b, mul, add));
	__m256i const ai32 = _mm256_cvtps_epi32(_mm256_fmadd_ps(_a, mul, add));
	
	__m128i const r16 = _mm_packus_epi32(_mm256_castsi256_si128(ri32), _mm256_extractf128_si256(ri32, 1));
	__m128i const g16 = _mm_packus_epi32(_mm256_castsi256_si128(gi32), _mm256_extractf128_si256(gi32, 1));
	__m128i const b16 = _mm_packus_epi32(_mm256_castsi256_si128(bi32), _mm256_extractf128_si256(bi32, 1));
	__m128i const a16 = _mm_packus_epi32(_mm256_castsi256_si128(ai32), _mm256_extractf128_si256(ai32, 1));
	

	// rrrrrrrrbbbbbbbb
	__m128i const r8b8 = _mm_packus_epi16(r16, b16);
	// ggggggggaaaaaaaa
	__m128i const g8a8 = _mm_packus_epi16(g16, a16);

	// rgrgrgrgrgrgrgrg
	__m128i const rg = _mm_unpacklo_epi8(r8b8, g8a8);
	// babababababababa
	__m128i const ba = _mm_unpackhi_epi8(r8b8, g8a8);

	__m128i const rgba_lo = _mm_unpacklo_epi16(rg, ba);
	__m128i const rgba_hi = _mm_unpackhi_epi16(rg, ba);

	__m128i* texelStoreLo = (__m128i*)o_texels;
	__m128i* texelStoreHi = (__m128i*)(&o_texels[4]);

	_mm_store_si128(texelStoreLo, rgba_lo);
	_mm_store_si128(texelStoreHi, rgba_hi);
}

}

}
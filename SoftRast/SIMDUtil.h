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


}

}
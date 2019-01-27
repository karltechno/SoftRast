#pragma once
#include "SoftRastTypes.h"
#include <immintrin.h>

namespace sr
{

namespace Renderer
{

using PixelShaderFn = void(*)(void const* _uniforms, __m256 const _varyings[sr::c_maxVaryings], float o_colour[4 * 8], __m256 const& _execMask);

struct GenericDrawBuffer
{
	void const* m_ptr;
	uint32_t m_num;
	uint32_t m_stride;
};

struct DrawCall
{
	PixelShaderFn m_pixelShader = nullptr;
	void const* m_pixelUniforms = nullptr;

	GenericDrawBuffer m_indexBuffer;
	GenericDrawBuffer m_positionBuffer;
	GenericDrawBuffer m_attributeBuffer;
};

}
}
#include "Rasterizer.h"

#include <kt/Vec4.h>
#include <kt/Vec3.h>
#include <kt/Mat4.h>
#include <kt/Memory.h>
#include <kt/Logging.h>
#include <float.h>

namespace sr
{

namespace Raster
{


void FillScreenTest(FrameBuffer& _buffer, uint8_t const color[3])
{
	uint8_t* p = _buffer.ptr;
	uint8_t* pEnd = _buffer.ptr + _buffer.width * _buffer.height * 4;

	while (p != pEnd)
	{
		*p++ = color[0];
		*p++ = color[1];
		*p++ = color[2];
		*p++ = 0xFF;
	}
}

void ClearDepthBufferTest(DepthBuffer& _buffer)
{
	float* p = _buffer.ptr;
	float* pEnd = _buffer.ptr + _buffer.width * _buffer.height;

	while (p != pEnd)
	{
		*p++ = 1.0f;
	}
}

constexpr int32_t c_subPixelBits = 4;
constexpr int32_t c_subPixelStep = 1 << c_subPixelBits;
constexpr int32_t c_subPixelMask = c_subPixelStep - 1;

constexpr int32_t c_blockSize = 8;

struct EdgeConstants
{
	void Set(int32_t const _v0[2], int32_t const _v1[2], int32_t _xmin, int32_t _ymin)
	{
#if 1
		c = (-_v0[0] * (_v1[1] - _v0[1]) + _v0[1] * (_v1[0] - _v0[0]));
#else
		int64_t const c_lhs = -_v0[0] * (_v1[1] - _v0[1]);
 		int64_t const c_rhs = _v0[1] * (_v1[0] - _v0[0]);
		c = int32_t((c_lhs + c_rhs) >> (c_subPixelBits + 1));
#endif
		dy_rasterCoord = (_v1[1] - _v0[1]);
		dx_rasterCoord = (_v0[0] - _v1[0]);

		// Left/horizontal fill rule
		if (dy_rasterCoord < 0 || dx_rasterCoord == 0 && dy_rasterCoord > 0)
		{
			c += 1;
		}

		cur_base = c + dy_rasterCoord * (_xmin << c_subPixelBits) + dx_rasterCoord * (_ymin << c_subPixelBits);
		eval = cur_base;

		dy_subpix = dy_rasterCoord << c_subPixelBits;
		dx_subpix = dx_rasterCoord << c_subPixelBits;
	}

	void IncX()
	{
		eval += dy_subpix;
	}

	void NextLine()
	{
		cur_base += dx_subpix;
		eval = cur_base;
	}

	int32_t c;

	int32_t dy_rasterCoord;
	int32_t dx_rasterCoord;

	int32_t dy_subpix;
	int32_t dx_subpix;

	int32_t cur_base;

	int32_t eval;
};

struct RasterEdgeConstants
{
	void Setup(int32_t const _v0[2], int32_t const _v1[2], int32_t _xminFp, int32_t _yminFp)
	{
#if 1
		int32_t c = -_v0[0] * (_v1[1] - _v0[1]) + _v0[1] * (_v1[0] - _v0[0]);
#else
		int64_t const c_lhs = -_v0[0] * (_v1[1] - _v0[1]);
		int64_t const c_rhs = _v0[1] * (_v1[0] - _v0[0]);
		c = int32_t((c_lhs + c_rhs) >> (c_subPixelBits + 1));
#endif
		dy = (_v1[1] - _v0[1]);
		dx = (_v0[0] - _v1[0]);

		// Left/horizontal fill rule
		if (dy < 0 || dx == 0 && dy > 0)
		{
			c += 1;
		}

		dy <<= c_subPixelBits;
		dx <<= c_subPixelBits;

		cBboxOrigin = c + dx * _yminFp + dy * _xminFp;
		dx_nextBlockRow = dy - dx * c_blockSize;
	}

	// Constant term at top left bbox vertex.
	int32_t cBboxOrigin;

	int32_t dx; 
	int32_t dy;

	// Skip to next row from block (-1*dx*y)
	int32_t dx_nextBlockRow;

	__m256i dx_step8;
	__m256i dy_step8;
};

struct BlockConstantsSIMD_8x8
{
	void Setup(EdgeConstants const& _e0, EdgeConstants const& _e1, EdgeConstants const& _e2)
	{
		e0_dx = _mm256_set1_epi32(_e0.dx_subpix);
		e1_dx = _mm256_set1_epi32(_e1.dx_subpix);
		e2_dx = _mm256_set1_epi32(_e2.dx_subpix);

		e0_dy = _e0.dy_subpix;
		e1_dy = _e1.dy_subpix;
		e2_dy = _e2.dy_subpix;
	}

	void InitRows(int32_t _e0eval, int32_t _e1eval, int32_t _e2eval, __m256i& o_e0, __m256i& o_e1, __m256i& o_e2) const
	{
		__m256i const ramp = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
		o_e0 = _mm256_add_epi32(_mm256_set1_epi32(_e0eval), _mm256_mullo_epi32(ramp, _mm256_set1_epi32(e0_dy)));
		o_e1 = _mm256_add_epi32(_mm256_set1_epi32(_e1eval), _mm256_mullo_epi32(ramp, _mm256_set1_epi32(e1_dy)));
		o_e2 = _mm256_add_epi32(_mm256_set1_epi32(_e2eval), _mm256_mullo_epi32(ramp, _mm256_set1_epi32(e2_dy)));
	}

	int32_t e0_dy, e1_dy, e2_dy;

	__m256i e0_dx;
	__m256i e1_dx;
	__m256i e2_dx;
};

static void ShadeWholeBlock
(
	RasterContext const& _ctx,
	PipelineTri const& _tri, 
	EdgeConstants const& _e0, 
	EdgeConstants const& _e1, 
	EdgeConstants const& _e2,
	int32_t _xmin,
	int32_t _ymin
)
{
	int32_t const x0_block = _xmin << c_subPixelBits;
	int32_t const y0_block = _ymin << c_subPixelBits;

	int32_t const x1_block = (_xmin + c_blockSize) << c_subPixelBits;
	int32_t const y1_block = (_ymin + c_blockSize) << c_subPixelBits;

	for (int32_t y_fp = y0_block; y_fp < y1_block; y_fp += c_subPixelStep)
	{
		for (int32_t x_fp = x0_block; x_fp < x1_block; x_fp += c_subPixelStep)
		{
			int32_t const x = x_fp >> c_subPixelBits;
			int32_t const y = y_fp >> c_subPixelBits;

			uint8_t* pixel = _ctx.frameBuffer->ptr + y * _ctx.vpWidth * 4 + x * 4; // assumes 32bit framebuffer

			int32_t const e0Eval = _e0.c + _e0.dy_rasterCoord * x_fp + _e0.dx_rasterCoord * y_fp;
			int32_t const e1Eval = _e1.c + _e1.dy_rasterCoord * x_fp + _e1.dx_rasterCoord * y_fp;
			int32_t const e2Eval = _e2.c + _e2.dy_rasterCoord * x_fp + _e2.dx_rasterCoord * y_fp;

			KT_ASSERT((e0Eval | e1Eval | e2Eval) > 0);

			float const baryV0 = e1Eval * _tri.halfRecipTriArea_fp;
			float const baryV1 = e2Eval * _tri.halfRecipTriArea_fp;
			float const baryV2 = 1.0f - baryV0 - baryV1;

			float const recipZ = baryV0 * _tri.verts[0].zOverW + baryV1 * _tri.verts[1].zOverW + baryV2 * _tri.verts[2].zOverW;

			float const e0_persp = e0Eval * _tri.verts[2].recipW;
			float const e1_persp = e1Eval * _tri.verts[0].recipW;
			float const e2_persp = e2Eval * _tri.verts[1].recipW;

			float const recip_persp_bary = 1.0f / (e0_persp + e1_persp + e2_persp);

			float const v0_persp = e1_persp * recip_persp_bary;
			float const v1_persp = e2_persp * recip_persp_bary;
			float const v2_persp = 1.0f - v1_persp - v0_persp;

			float* depthPtr = _ctx.depthBuffer->At(x, y);

			if (*depthPtr > recipZ && recipZ > 0.0f)
			{
				*depthPtr = recipZ;
				uint8_t thecol = (uint8_t)(255 * (1.0f - fmodf(kt::Min(v2_persp, kt::Min(v0_persp, v1_persp)), 0.15f)));

				float const* vary0 = _tri.verts[0].varyings;
				float const* vary1 = _tri.verts[1].varyings;
				float const* vary2 = _tri.verts[2].varyings;

#if 1
				pixel[0] = uint8_t((vary0[0] * v0_persp + vary1[0] * v1_persp + vary2[0] * v2_persp) * 255.0f);
				pixel[1] = uint8_t((vary0[1] * v0_persp + vary1[1] * v1_persp + vary2[1] * v2_persp) * 255.0f);
				pixel[2] = uint8_t((vary0[2] * v0_persp + vary1[2] * v1_persp + vary2[2] * v2_persp) * 255.0f);
				pixel[3] = 0xFF;
#else
				pixel[0] = 255;
				pixel[1] = 0;
				pixel[2] = 0;
				pixel[3] = 0xFF;
#endif
			}

		}
	}
}

static void RasterizePartialBlockSIMD_8x8
(
	RasterContext const& _ctx,
	PipelineTri const& _tri,
	BlockConstantsSIMD_8x8 const& _wholeBlock,
	int32_t _e0evalBlockCorner,
	int32_t _e1evalBlockCorner,
	int32_t _e2evalBlockCorner,
	int32_t _xminRaster,
	int32_t _yminRaster
)
{
	__m256i e0row, e1row, e2row;

	_wholeBlock.InitRows(_e0evalBlockCorner, _e1evalBlockCorner, _e2evalBlockCorner, e0row, e1row, e2row);

	// Should store these before in simd state?
	__m256 const recipTriArea = _mm256_broadcast_ss(&_tri.halfRecipTriArea_fp);

	__m256 const zOverW_v0 = _mm256_broadcast_ss(&_tri.verts[0].zOverW);
	__m256 const zOverW_v1 = _mm256_broadcast_ss(&_tri.verts[1].zOverW);
	__m256 const zOverW_v2 = _mm256_broadcast_ss(&_tri.verts[2].zOverW);

	__m256 const recipW_v0 = _mm256_broadcast_ss(&_tri.verts[0].recipW);
	__m256 const recipW_v1 = _mm256_broadcast_ss(&_tri.verts[1].recipW);
	__m256 const recipW_v2 = _mm256_broadcast_ss(&_tri.verts[2].recipW);

	__m256 const one256 = _mm256_set1_ps(1.0f);

	__m256i const minusOne256 = _mm256_set1_epi32(-1);

	for (int32_t y = _yminRaster; y < (_yminRaster + c_blockSize); ++y)
	{
		uint8_t* pixelBegin = _ctx.frameBuffer->ptr + y * _ctx.vpWidth * 4 + _xminRaster * 4; // assumes 32bit framebuffer

		__m256i edgeMask = _mm256_or_si256(_mm256_or_si256(e0row, e1row), e2row);
		edgeMask = _mm256_srai_epi32(edgeMask, 31); // or sign bits, arith shift right so sign bit = ~0 and no sign bit = 0
		edgeMask = _mm256_xor_si256(edgeMask, _mm256_cmpeq_epi32(edgeMask, edgeMask)); // flip bits so negative (~0) is 0.

		if (_mm256_movemask_ps(_mm256_castsi256_ps(edgeMask)))
		{
			__m256 const e0_ps = _mm256_cvtepi32_ps(e0row);
			__m256 const e1_ps = _mm256_cvtepi32_ps(e1row);
			__m256 const e2_ps = _mm256_cvtepi32_ps(e2row);

			__m256 const baryv0 = _mm256_mul_ps(recipTriArea, e1_ps);
			__m256 const baryv1 = _mm256_mul_ps(recipTriArea, e2_ps);
			__m256 const baryv2 = _mm256_sub_ps(_mm256_sub_ps(one256, baryv0), baryv1);

			__m256 const recipZTemp1 = _mm256_mul_ps(baryv2, zOverW_v2);
			__m256 const recipZTemp2 = _mm256_fmadd_ps(baryv1, zOverW_v1, recipZTemp1);
			__m256 const recipZ = _mm256_fmadd_ps(baryv0, zOverW_v0, recipZTemp2);

			__m256 const e0_persp = _mm256_mul_ps(e0_ps, recipW_v2);
			__m256 const e1_persp = _mm256_mul_ps(e1_ps, recipW_v0);
			__m256 const e2_persp = _mm256_mul_ps(e2_ps, recipW_v1);

			__m256 const recip_persp_bary = _mm256_div_ps(one256, _mm256_add_ps(_mm256_add_ps(e0_persp, e1_persp), e2_persp));

			__m256 const v0_persp = _mm256_mul_ps(e1_persp, recip_persp_bary);
			__m256 const v1_persp = _mm256_mul_ps(e2_persp, recip_persp_bary);
			__m256 const v2_persp = _mm256_sub_ps(one256, _mm256_add_ps(v1_persp, v0_persp));


			float* depthPtr = _ctx.depthBuffer->At(_xminRaster, y);
			__m256 const depthGather = _mm256_loadu_ps(depthPtr);

			__m256 depthCmpMask = _mm256_and_ps(_mm256_cmp_ps(recipZ, _mm256_setzero_ps(), _CMP_GT_OQ), _mm256_cmp_ps(depthGather, recipZ, _CMP_GT_OQ));
			depthCmpMask = _mm256_and_ps(_mm256_castsi256_ps(edgeMask), depthCmpMask);

			if (int32_t maskReg = _mm256_movemask_ps(depthCmpMask))
			{
				// blend new and old depth
				__m256 const newDepth = _mm256_blendv_ps(depthGather, recipZ, depthCmpMask);
				_mm256_storeu_ps(depthPtr, newDepth);

				__m256 const r0 = _mm256_broadcast_ss(&_tri.verts[0].varyings[0]);
				__m256 const r1 = _mm256_broadcast_ss(&_tri.verts[1].varyings[0]);
				__m256 const r2 = _mm256_broadcast_ss(&_tri.verts[2].varyings[0]);

				__m256 const g0 = _mm256_broadcast_ss(&_tri.verts[0].varyings[1]);
				__m256 const g1 = _mm256_broadcast_ss(&_tri.verts[1].varyings[1]);
				__m256 const g2 = _mm256_broadcast_ss(&_tri.verts[2].varyings[1]);

				__m256 const b0 = _mm256_broadcast_ss(&_tri.verts[0].varyings[2]);
				__m256 const b1 = _mm256_broadcast_ss(&_tri.verts[1].varyings[2]);
				__m256 const b2 = _mm256_broadcast_ss(&_tri.verts[2].varyings[2]);

				__m256 const c_255 = _mm256_set1_ps(255.0f);

				__m256 const red = _mm256_mul_ps(c_255, _mm256_fmadd_ps(v2_persp, r2, (_mm256_fmadd_ps(v1_persp, r1, _mm256_mul_ps(v0_persp, r0)))));
				__m256 const green = _mm256_mul_ps(c_255, _mm256_fmadd_ps(v2_persp, g2, (_mm256_fmadd_ps(v1_persp, g1, _mm256_mul_ps(v0_persp, g0)))));
				__m256 const blue = _mm256_mul_ps(c_255, _mm256_fmadd_ps(v2_persp, b2, (_mm256_fmadd_ps(v1_persp, b1, _mm256_mul_ps(v0_persp, b0)))));
				
				float red_store[8];
				float green_store[8];
				float blue_store[8];
				_mm256_store_ps(red_store, red);
				_mm256_store_ps(green_store, green);
				_mm256_store_ps(blue_store, blue);

				for (uint32_t pixIdx = 0; pixIdx < c_blockSize; ++pixIdx)
				{
					if (maskReg & (1 << pixIdx))
					{
						pixelBegin[pixIdx * 4 + 0] = uint8_t(red_store[pixIdx]);
						pixelBegin[pixIdx * 4 + 1] = uint8_t(green_store[pixIdx]);
						pixelBegin[pixIdx * 4 + 2] = uint8_t(blue_store[pixIdx]);
						pixelBegin[pixIdx * 4 + 3] = 0xFF;
					}
				}
			}
		}

		e0row = _mm256_add_epi32(e0row, _wholeBlock.e0_dx);
		e1row = _mm256_add_epi32(e1row, _wholeBlock.e1_dx);
		e2row = _mm256_add_epi32(e2row, _wholeBlock.e2_dx);
	}
}

static void RasterizeWholeBlockSIMD_8x8
(
	RasterContext const& _ctx,
	PipelineTri const& _tri,
	BlockConstantsSIMD_8x8 const& _wholeBlock,
	int32_t _e0evalBlockCorner,
	int32_t _e1evalBlockCorner,
	int32_t _e2evalBlockCorner,
	int32_t _xminRaster,
	int32_t _yminRaster
)
{
	__m256i e0row, e1row, e2row;

	_wholeBlock.InitRows(_e0evalBlockCorner, _e1evalBlockCorner, _e2evalBlockCorner, e0row, e1row, e2row);

	// Should store these before in simd state?
	__m256 const recipTriArea = _mm256_broadcast_ss(&_tri.halfRecipTriArea_fp);

	__m256 const zOverW_v0 = _mm256_broadcast_ss(&_tri.verts[0].zOverW);
	__m256 const zOverW_v1 = _mm256_broadcast_ss(&_tri.verts[1].zOverW);
	__m256 const zOverW_v2 = _mm256_broadcast_ss(&_tri.verts[2].zOverW);

	__m256 const recipW_v0 = _mm256_broadcast_ss(&_tri.verts[0].recipW);
	__m256 const recipW_v1 = _mm256_broadcast_ss(&_tri.verts[1].recipW);
	__m256 const recipW_v2 = _mm256_broadcast_ss(&_tri.verts[2].recipW);

	__m256 const one256 = _mm256_set1_ps(1.0f);

	for (int32_t y = _yminRaster; y < (_yminRaster + c_blockSize); ++y)
	{
		uint8_t* pixelBegin = _ctx.frameBuffer->ptr + y * _ctx.vpWidth * 4 + _xminRaster * 4; // assumes 32bit framebuffer

		__m256 const e0_ps = _mm256_cvtepi32_ps(e0row);
		__m256 const e1_ps = _mm256_cvtepi32_ps(e1row);
		__m256 const e2_ps = _mm256_cvtepi32_ps(e2row);

		__m256 const baryv0 = _mm256_mul_ps(recipTriArea, e1_ps);
		__m256 const baryv1 = _mm256_mul_ps(recipTriArea, e2_ps);
		__m256 const baryv2 = _mm256_sub_ps(_mm256_sub_ps(one256, baryv0), baryv1);

		__m256 const recipZTemp1 = _mm256_mul_ps(baryv2, zOverW_v2);
		__m256 const recipZTemp2 = _mm256_fmadd_ps(baryv1, zOverW_v1, recipZTemp1);
		__m256 const recipZ = _mm256_fmadd_ps(baryv0, zOverW_v0, recipZTemp2);

		__m256 const e0_persp = _mm256_mul_ps(e0_ps, recipW_v2);
		__m256 const e1_persp = _mm256_mul_ps(e1_ps, recipW_v0);
		__m256 const e2_persp = _mm256_mul_ps(e2_ps, recipW_v1);

		__m256 const recip_persp_bary = _mm256_div_ps(one256, _mm256_add_ps(_mm256_add_ps(e0_persp, e1_persp), e2_persp));

		__m256 const v0_persp = _mm256_mul_ps(e1_persp, recip_persp_bary);
		__m256 const v1_persp = _mm256_mul_ps(e2_persp, recip_persp_bary);
		__m256 const v2_persp = _mm256_sub_ps(one256, _mm256_add_ps(v1_persp, v0_persp));


		float* depthPtr = _ctx.depthBuffer->At(_xminRaster, y);
		__m256 const depthGather = _mm256_loadu_ps(depthPtr);

		__m256 const depthCmpMask = _mm256_and_ps(_mm256_cmp_ps(recipZ, _mm256_setzero_ps(), _CMP_GT_OQ), _mm256_cmp_ps(depthGather, recipZ, _CMP_GT_OQ));

		if (int32_t maskReg = _mm256_movemask_ps(depthCmpMask))
		{
			// blend new and old depth
			__m256 const newDepth = _mm256_blendv_ps(recipZ, depthGather, depthCmpMask);
			_mm256_storeu_ps(depthPtr, newDepth);
				
			__m256 const r0 = _mm256_broadcast_ss(&_tri.verts[0].varyings[0]);
			__m256 const r1 = _mm256_broadcast_ss(&_tri.verts[1].varyings[0]);
			__m256 const r2 = _mm256_broadcast_ss(&_tri.verts[2].varyings[0]);

			__m256 const g0 = _mm256_broadcast_ss(&_tri.verts[0].varyings[1]);
			__m256 const g1 = _mm256_broadcast_ss(&_tri.verts[1].varyings[1]);
			__m256 const g2 = _mm256_broadcast_ss(&_tri.verts[2].varyings[1]);

			__m256 const b0 = _mm256_broadcast_ss(&_tri.verts[0].varyings[2]);
			__m256 const b1 = _mm256_broadcast_ss(&_tri.verts[1].varyings[2]);
			__m256 const b2 = _mm256_broadcast_ss(&_tri.verts[2].varyings[2]);

			__m256 const c_255 = _mm256_set1_ps(255.0f);

			__m256 const red	= _mm256_mul_ps(c_255, _mm256_fmadd_ps(v2_persp, r2, (_mm256_fmadd_ps(v1_persp, r1, _mm256_mul_ps(v0_persp, r0)))));
			__m256 const green	= _mm256_mul_ps(c_255, _mm256_fmadd_ps(v2_persp, g2, (_mm256_fmadd_ps(v1_persp, g1, _mm256_mul_ps(v0_persp, g0)))));
			__m256 const blue	= _mm256_mul_ps(c_255, _mm256_fmadd_ps(v2_persp, b2, (_mm256_fmadd_ps(v1_persp, b1, _mm256_mul_ps(v0_persp, b0)))));

			float red_store[8];
			float green_store[8];
			float blue_store[8];
			_mm256_store_ps(red_store, red);
			_mm256_store_ps(green_store, green);
			_mm256_store_ps(blue_store, blue);

			for (uint32_t pixIdx = 0; pixIdx < c_blockSize; ++pixIdx)
			{
				if (maskReg & (1 << pixIdx))
				{
					pixelBegin[pixIdx * 4 + 0] = uint8_t(red_store[pixIdx]);
					pixelBegin[pixIdx * 4 + 1] = uint8_t(green_store[pixIdx]);
					pixelBegin[pixIdx * 4 + 2] = uint8_t(blue_store[pixIdx]);
					pixelBegin[pixIdx * 4 + 3] = 0;
				}
			}
		}

		e0row = _mm256_add_epi32(e0row, _wholeBlock.e0_dx);
		e1row = _mm256_add_epi32(e1row, _wholeBlock.e1_dx);
		e2row = _mm256_add_epi32(e2row, _wholeBlock.e2_dx);
	}
}

static void RasterizeWholeBlockSIMD_8x8_WithShader
(
	RasterContext const& _ctx,
	PipelineTri const& _tri,
	BlockConstantsSIMD_8x8 const& _wholeBlock,
	int32_t _e0evalBlockCorner,
	int32_t _e1evalBlockCorner,
	int32_t _e2evalBlockCorner,
	int32_t _xminRaster,
	int32_t _yminRaster,
	Renderer::DrawCall const& _call
)
{
	__m256i e0row, e1row, e2row;

	_wholeBlock.InitRows(_e0evalBlockCorner, _e1evalBlockCorner, _e2evalBlockCorner, e0row, e1row, e2row);

	// Should store these before in simd state?
	__m256 const recipTriArea = _mm256_broadcast_ss(&_tri.halfRecipTriArea_fp);

	__m256 const zOverW_v0 = _mm256_broadcast_ss(&_tri.verts[0].zOverW);
	__m256 const zOverW_v1 = _mm256_broadcast_ss(&_tri.verts[1].zOverW);
	__m256 const zOverW_v2 = _mm256_broadcast_ss(&_tri.verts[2].zOverW);

	__m256 const recipW_v0 = _mm256_broadcast_ss(&_tri.verts[0].recipW);
	__m256 const recipW_v1 = _mm256_broadcast_ss(&_tri.verts[1].recipW);
	__m256 const recipW_v2 = _mm256_broadcast_ss(&_tri.verts[2].recipW);

	__m256 const one256 = _mm256_set1_ps(1.0f);


	for (int32_t y = _yminRaster; y < (_yminRaster + c_blockSize); ++y)
	{
		uint8_t* pixelBegin = _ctx.frameBuffer->ptr + y * _ctx.vpWidth * 4 + _xminRaster * 4; // assumes 32bit framebuffer

		__m256 const e0_ps = _mm256_cvtepi32_ps(e0row);
		__m256 const e1_ps = _mm256_cvtepi32_ps(e1row);
		__m256 const e2_ps = _mm256_cvtepi32_ps(e2row);

		__m256 const baryv0 = _mm256_mul_ps(recipTriArea, e1_ps);
		__m256 const baryv1 = _mm256_mul_ps(recipTriArea, e2_ps);
		__m256 const baryv2 = _mm256_sub_ps(_mm256_sub_ps(one256, baryv0), baryv1);

		__m256 const recipZTemp1 = _mm256_mul_ps(baryv2, zOverW_v2);
		__m256 const recipZTemp2 = _mm256_fmadd_ps(baryv1, zOverW_v1, recipZTemp1);
		__m256 const recipZ = _mm256_fmadd_ps(baryv0, zOverW_v0, recipZTemp2);

		__m256 const e0_persp = _mm256_mul_ps(e0_ps, recipW_v2);
		__m256 const e1_persp = _mm256_mul_ps(e1_ps, recipW_v0);
		__m256 const e2_persp = _mm256_mul_ps(e2_ps, recipW_v1);

		__m256 const recip_persp_bary = _mm256_div_ps(one256, _mm256_add_ps(_mm256_add_ps(e0_persp, e1_persp), e2_persp));

		__m256 const v0_persp = _mm256_mul_ps(e1_persp, recip_persp_bary);
		__m256 const v1_persp = _mm256_mul_ps(e2_persp, recip_persp_bary);
		__m256 const v2_persp = _mm256_sub_ps(one256, _mm256_add_ps(v1_persp, v0_persp));


		float* depthPtr = _ctx.depthBuffer->At(_xminRaster, y);
		__m256 const depthGather = _mm256_loadu_ps(depthPtr);

		__m256 const depthCmpMask = _mm256_and_ps(_mm256_cmp_ps(recipZ, _mm256_setzero_ps(), _CMP_GT_OQ), _mm256_cmp_ps(depthGather, recipZ, _CMP_GT_OQ));

		if (int32_t maskReg = _mm256_movemask_ps(depthCmpMask))
		{
			__m256 rowVaryings[c_maxVaryings];

			for (uint32_t i = 0; i < c_maxVaryings; ++i)
			{
				__m256 const v0_vary = _mm256_broadcast_ss(&_tri.verts[0].varyings[i]);
				__m256 const v1_vary = _mm256_broadcast_ss(&_tri.verts[1].varyings[i]);
				__m256 const v2_vary = _mm256_broadcast_ss(&_tri.verts[2].varyings[i]);
				rowVaryings[i] = _mm256_fmadd_ps(v0_persp, v0_vary, _mm256_fmadd_ps(v1_persp, v1_vary, _mm256_mul_ps(v2_persp, v2_vary)));
			}
			
			KT_ALIGNAS(32) float colourRGBA[4 * 8];

			_call.m_pixelShader(_call.m_pixelUniforms, rowVaryings, colourRGBA, depthCmpMask);

			for (uint32_t pixIdx = 0; pixIdx < c_blockSize; ++pixIdx)
			{
				if (maskReg & (1 << pixIdx))
				{
					pixelBegin[pixIdx * 4 + 0] = uint8_t(kt::Min(1.0f, colourRGBA[pixIdx * 4]) * 255.0f);
					pixelBegin[pixIdx * 4 + 1] = uint8_t(kt::Min(1.0f, colourRGBA[pixIdx * 4 + 1]) * 255.0f);
					pixelBegin[pixIdx * 4 + 2] = uint8_t(kt::Min(1.0f, colourRGBA[pixIdx * 4 + 2]) * 255.0f);
					pixelBegin[pixIdx * 4 + 3] = 0;
				}
			}
		}

		e0row = _mm256_add_epi32(e0row, _wholeBlock.e0_dx);
		e1row = _mm256_add_epi32(e1row, _wholeBlock.e1_dx);
		e2row = _mm256_add_epi32(e2row, _wholeBlock.e2_dx);
	}
}



static void ShadePartialBlock
(
	RasterContext const& _ctx,
	PipelineTri const& _tri,
	EdgeConstants const& _e0,
	EdgeConstants const& _e1,
	EdgeConstants const& _e2,
	int32_t _xmin,
	int32_t _ymin
)
{
	int32_t const x0_block = _xmin << c_subPixelBits;
	int32_t const y0_block = _ymin << c_subPixelBits;

	int32_t const x1_block = (_xmin + c_blockSize) << c_subPixelBits;
	int32_t const y1_block = (_ymin + c_blockSize) << c_subPixelBits;

	int32_t e0Base = _e0.c + _e0.dy_rasterCoord * x0_block + _e0.dx_rasterCoord * y0_block;
	int32_t e1Base = _e1.c + _e1.dy_rasterCoord * x0_block + _e1.dx_rasterCoord * y0_block;
	int32_t e2Base = _e2.c + _e2.dy_rasterCoord * x0_block + _e2.dx_rasterCoord * y0_block;

	int32_t e0Eval = e0Base;
	int32_t e1Eval = e1Base;
	int32_t e2Eval = e2Base;


	for (int32_t y_fp = y0_block; y_fp < y1_block; y_fp += c_subPixelStep)
	{
		for (int32_t x_fp = x0_block; x_fp < x1_block; x_fp += c_subPixelStep)
		{
			int32_t const x = x_fp >> c_subPixelBits;
			int32_t const y = y_fp >> c_subPixelBits;

			bool inTri = true;

			inTri &= (e0Eval >= 0);
			inTri &= (e1Eval >= 0);
			inTri &= (e2Eval >= 0);

			uint8_t* pixel = _ctx.frameBuffer->ptr + y * _ctx.vpWidth * 4 + x * 4;

			if (inTri)
			{
				float const baryV0 = e1Eval * _tri.halfRecipTriArea_fp;
				float const baryV1 = e2Eval * _tri.halfRecipTriArea_fp;
				float const baryV2 = 1.0f - baryV0 - baryV1;

				float const recipZ = baryV0 * _tri.verts[0].zOverW + baryV1 * _tri.verts[1].zOverW + baryV2 * _tri.verts[2].zOverW;

				float const e0_persp = e0Eval * _tri.verts[2].recipW;
				float const e1_persp = e1Eval * _tri.verts[0].recipW;
				float const e2_persp = e2Eval * _tri.verts[1].recipW;

				float const recip_persp_bary = 1.0f / (e0_persp + e1_persp + e2_persp);

				float const v0_persp = e1_persp * recip_persp_bary;
				float const v1_persp = e2_persp * recip_persp_bary;
				float const v2_persp = 1.0f - v1_persp - v0_persp;

				float* depthPtr = _ctx.depthBuffer->At(x, y);

				if (*depthPtr > recipZ && recipZ > 0.0f)
				{
					*depthPtr = recipZ;

					float const* vary0 = _tri.verts[0].varyings;
					float const* vary1 = _tri.verts[1].varyings;
					float const* vary2 = _tri.verts[2].varyings;

					pixel[0] = uint8_t((vary0[0] * v0_persp + vary1[0] * v1_persp + vary2[0] * v2_persp) * 255.0f);
					pixel[1] = uint8_t((vary0[1] * v0_persp + vary1[1] * v1_persp + vary2[1] * v2_persp) * 255.0f);
					pixel[2] = uint8_t((vary0[2] * v0_persp + vary1[2] * v1_persp + vary2[2] * v2_persp) * 255.0f);
					pixel[3] = 0xFF;
				}
			}

			e0Eval += _e0.dy_subpix;
			e1Eval += _e1.dy_subpix;
			e2Eval += _e2.dy_subpix;
		}
		e0Base += _e0.dx_subpix;
		e1Base += _e1.dx_subpix;
		e2Base += _e2.dx_subpix;
		
		e0Eval = e0Base;
		e1Eval = e1Base;
		e2Eval = e2Base;
	}
}

static void RasterTransformedTri_WithBlocks(RasterContext const& _ctx, PipelineTri const& _tri)
{
	int32_t const xmin = _tri.xmin & ~(c_blockSize - 1);
	int32_t const ymin = _tri.ymin & ~(c_blockSize - 1);
	int32_t const xmax = _tri.xmax;
	int32_t const ymax = _tri.ymax;

	EdgeConstants e0_edge_eq, e1_edge_eq, e2_edge_eq;
	e0_edge_eq.Set(_tri.v0_fp, _tri.v1_fp, xmin, ymin);
	e1_edge_eq.Set(_tri.v1_fp, _tri.v2_fp, xmin, ymin);
	e2_edge_eq.Set(_tri.v2_fp, _tri.v0_fp, xmin, ymin);

	BlockConstantsSIMD_8x8 wholeBlock;
	wholeBlock.Setup(e0_edge_eq, e1_edge_eq, e2_edge_eq);

	for (int32_t y = ymin; y < ymax; y += c_blockSize)
	{
		for (int32_t x = xmin; x < xmax; x += c_blockSize)
		{
			// Block corners
			int32_t const x0_block = x << c_subPixelBits;
			int32_t const y0_block = y << c_subPixelBits;

			int32_t const x1_block = (x + c_blockSize - 1) << c_subPixelBits;
			int32_t const y1_block = (y + c_blockSize - 1) << c_subPixelBits;

			// Test whole block.
			int32_t const e0_x0y0 = e0_edge_eq.c + e0_edge_eq.dy_rasterCoord * x0_block + e0_edge_eq.dx_rasterCoord * y0_block;
			int32_t const e0_x0y1 = e0_edge_eq.c + e0_edge_eq.dy_rasterCoord * x0_block + e0_edge_eq.dx_rasterCoord * y1_block;
			int32_t const e0_x1y0 = e0_edge_eq.c + e0_edge_eq.dy_rasterCoord * x1_block + e0_edge_eq.dx_rasterCoord * y0_block;
			int32_t const e0_x1y1 = e0_edge_eq.c + e0_edge_eq.dy_rasterCoord * x1_block + e0_edge_eq.dx_rasterCoord * y1_block;

			int32_t const e0_blockeval = (e0_x0y0 | e0_x0y1 | e0_x1y0 | e0_x1y1);

			int32_t const e1_x0y0 = e1_edge_eq.c + e1_edge_eq.dy_rasterCoord * x0_block + e1_edge_eq.dx_rasterCoord * y0_block;
			int32_t const e1_x0y1 = e1_edge_eq.c + e1_edge_eq.dy_rasterCoord * x0_block + e1_edge_eq.dx_rasterCoord * y1_block;
			int32_t const e1_x1y0 = e1_edge_eq.c + e1_edge_eq.dy_rasterCoord * x1_block + e1_edge_eq.dx_rasterCoord * y0_block;
			int32_t const e1_x1y1 = e1_edge_eq.c + e1_edge_eq.dy_rasterCoord * x1_block + e1_edge_eq.dx_rasterCoord * y1_block;

			int32_t const e1_blockeval = (e1_x0y0 | e1_x0y1 | e1_x1y0 | e1_x1y1);

			int32_t const e2_x0y0 = e2_edge_eq.c + e2_edge_eq.dy_rasterCoord * x0_block + e2_edge_eq.dx_rasterCoord * y0_block;
			int32_t const e2_x0y1 = e2_edge_eq.c + e2_edge_eq.dy_rasterCoord * x0_block + e2_edge_eq.dx_rasterCoord * y1_block;
			int32_t const e2_x1y0 = e2_edge_eq.c + e2_edge_eq.dy_rasterCoord * x1_block + e2_edge_eq.dx_rasterCoord * y0_block;
			int32_t const e2_x1y1 = e2_edge_eq.c + e2_edge_eq.dy_rasterCoord * x1_block + e2_edge_eq.dx_rasterCoord * y1_block;

			int32_t const e2_blockeval = (e2_x0y0 | e2_x0y1 | e2_x1y0 | e2_x1y1);

			static bool s_doSimd = true;

			if ((e0_blockeval | e1_blockeval | e2_blockeval) > 0)
			{
				if (!s_doSimd)
				{
					ShadeWholeBlock(_ctx, _tri, e0_edge_eq, e1_edge_eq, e2_edge_eq, x, y); // todo sort out where we have fixed point coords
				}
				else
				{
					RasterizeWholeBlockSIMD_8x8(_ctx, _tri, wholeBlock, e0_x0y0, e1_x0y0, e2_x0y0, x, y);
				}
			}
			else /*if(e0_blockeval > 0 || e1_blockeval > 0 || e2_blockeval > 0)*/
			{
				if (!s_doSimd)
				{
					ShadePartialBlock(_ctx, _tri, e0_edge_eq, e1_edge_eq, e2_edge_eq, x, y);
				}
				else
				{
					RasterizePartialBlockSIMD_8x8(_ctx, _tri, wholeBlock, e0_x0y0, e1_x0y0, e2_x0y0, x, y);
				}
			}
		}
	}
}

static void RasterTransformedTri_WithBlocksAndShader(RasterContext const& _ctx, PipelineTri const& _tri, Renderer::DrawCall const& _call)
{
	int32_t const xmin = _tri.xmin & ~(c_blockSize - 1);
	int32_t const ymin = _tri.ymin & ~(c_blockSize - 1);
	int32_t const xmax = _tri.xmax;
	int32_t const ymax = _tri.ymax;

	EdgeConstants e0_edge_eq, e1_edge_eq, e2_edge_eq;
	e0_edge_eq.Set(_tri.v0_fp, _tri.v1_fp, xmin, ymin);
	e1_edge_eq.Set(_tri.v1_fp, _tri.v2_fp, xmin, ymin);
	e2_edge_eq.Set(_tri.v2_fp, _tri.v0_fp, xmin, ymin);

	BlockConstantsSIMD_8x8 wholeBlock;
	wholeBlock.Setup(e0_edge_eq, e1_edge_eq, e2_edge_eq);

	for (int32_t y = ymin; y < ymax; y += c_blockSize)
	{
		for (int32_t x = xmin; x < xmax; x += c_blockSize)
		{
			// Block corners
			int32_t const x0_block = x << c_subPixelBits;
			int32_t const y0_block = y << c_subPixelBits;

			int32_t const x1_block = (x + c_blockSize - 1) << c_subPixelBits;
			int32_t const y1_block = (y + c_blockSize - 1) << c_subPixelBits;

			// Test whole block.
			int32_t const e0_x0y0 = e0_edge_eq.c + e0_edge_eq.dy_rasterCoord * x0_block + e0_edge_eq.dx_rasterCoord * y0_block;
			int32_t const e0_x0y1 = e0_edge_eq.c + e0_edge_eq.dy_rasterCoord * x0_block + e0_edge_eq.dx_rasterCoord * y1_block;
			int32_t const e0_x1y0 = e0_edge_eq.c + e0_edge_eq.dy_rasterCoord * x1_block + e0_edge_eq.dx_rasterCoord * y0_block;
			int32_t const e0_x1y1 = e0_edge_eq.c + e0_edge_eq.dy_rasterCoord * x1_block + e0_edge_eq.dx_rasterCoord * y1_block;

			int32_t const e0_blockeval = (e0_x0y0 | e0_x0y1 | e0_x1y0 | e0_x1y1);

			int32_t const e1_x0y0 = e1_edge_eq.c + e1_edge_eq.dy_rasterCoord * x0_block + e1_edge_eq.dx_rasterCoord * y0_block;
			int32_t const e1_x0y1 = e1_edge_eq.c + e1_edge_eq.dy_rasterCoord * x0_block + e1_edge_eq.dx_rasterCoord * y1_block;
			int32_t const e1_x1y0 = e1_edge_eq.c + e1_edge_eq.dy_rasterCoord * x1_block + e1_edge_eq.dx_rasterCoord * y0_block;
			int32_t const e1_x1y1 = e1_edge_eq.c + e1_edge_eq.dy_rasterCoord * x1_block + e1_edge_eq.dx_rasterCoord * y1_block;

			int32_t const e1_blockeval = (e1_x0y0 | e1_x0y1 | e1_x1y0 | e1_x1y1);

			int32_t const e2_x0y0 = e2_edge_eq.c + e2_edge_eq.dy_rasterCoord * x0_block + e2_edge_eq.dx_rasterCoord * y0_block;
			int32_t const e2_x0y1 = e2_edge_eq.c + e2_edge_eq.dy_rasterCoord * x0_block + e2_edge_eq.dx_rasterCoord * y1_block;
			int32_t const e2_x1y0 = e2_edge_eq.c + e2_edge_eq.dy_rasterCoord * x1_block + e2_edge_eq.dx_rasterCoord * y0_block;
			int32_t const e2_x1y1 = e2_edge_eq.c + e2_edge_eq.dy_rasterCoord * x1_block + e2_edge_eq.dx_rasterCoord * y1_block;

			int32_t const e2_blockeval = (e2_x0y0 | e2_x0y1 | e2_x1y0 | e2_x1y1);

			static bool s_doSimd = true;

			if ((e0_blockeval | e1_blockeval | e2_blockeval) > 0)
			{
				RasterizeWholeBlockSIMD_8x8_WithShader(_ctx, _tri, wholeBlock, e0_x0y0, e1_x0y0, e2_x0y0, x, y, _call);
			}
			else /*if(e0_blockeval > 0 || e1_blockeval > 0 || e2_blockeval > 0)*/
			{
				RasterizePartialBlockSIMD_8x8(_ctx, _tri, wholeBlock, e0_x0y0, e1_x0y0, e2_x0y0, x, y);
			}
		}
	}
}


static void RasterTransformedTri(FrameBuffer& _buffer, DepthBuffer& _depthBuffer, PipelineTri const& _tri)
{
	EdgeConstants e0_edge_eq, e1_edge_eq, e2_edge_eq;
	e0_edge_eq.Set(_tri.v0_fp, _tri.v1_fp, _tri.xmin, _tri.ymin);
	e1_edge_eq.Set(_tri.v1_fp, _tri.v2_fp, _tri.xmin, _tri.ymin);
	e2_edge_eq.Set(_tri.v2_fp, _tri.v0_fp, _tri.xmin, _tri.ymin);


	for (int32_t y = _tri.ymin; y < _tri.ymax; ++y)
	{
		for (int32_t x = _tri.xmin; x < _tri.xmax; ++x)
		{
			bool inTri = true;
#if 1
			int32_t const e0Eval = e0_edge_eq.eval;
			int32_t const e1Eval = e1_edge_eq.eval;
			int32_t const e2Eval = e2_edge_eq.eval;

			inTri &= (e0Eval >= 0);
			inTri &= (e1Eval >= 0);
			inTri &= (e2Eval >= 0);

#else
			float const eval0 = EvalEdge(v0raster, v1raster, x, y);
			float const eval1 = EvalEdge(v1raster, v2raster, x, y);
			float const eval2 = EvalEdge(v2raster, v0raster, x, y);
			inTri = eval0 >= 0.0f && eval1 >= 0.0f && eval2 >= 0.0f;
#endif
			uint8_t* pixel = _buffer.ptr + y * _buffer.width * 4 + x * 4;

			if (inTri)
			{
				float const baryV0 = e1Eval * _tri.halfRecipTriArea_fp;
				float const baryV1 = e2Eval * _tri.halfRecipTriArea_fp;
				float const baryV2 = 1.0f - baryV0 - baryV1;

				float const recipZ = baryV0 * _tri.verts[0].zOverW + baryV1 * _tri.verts[1].zOverW + baryV2 * _tri.verts[2].zOverW;

				float const e0_persp = e0Eval * _tri.verts[2].recipW;
				float const e1_persp = e1Eval * _tri.verts[0].recipW;
				float const e2_persp = e2Eval * _tri.verts[1].recipW;

				float const recip_persp_bary = 1.0f / (e0_persp + e1_persp + e2_persp);

				float const v0_persp = e1_persp * recip_persp_bary;
				float const v1_persp = e2_persp * recip_persp_bary;
				float const v2_persp = 1.0f - v1_persp - v0_persp;

				float* depthPtr = _depthBuffer.At(x, y);

				if (*depthPtr > recipZ && recipZ > 0.0f)
				{
					*depthPtr = recipZ;
					uint8_t thecol = (uint8_t)(255 * (1.0f - fmodf(kt::Min(v2_persp, kt::Min(v0_persp, v1_persp)), 0.15f)));

					float const* vary0 = _tri.verts[0].varyings;
					float const* vary1 = _tri.verts[1].varyings;
					float const* vary2 = _tri.verts[2].varyings;

#if 0
					pixel[0] = (uint8_t)(baryV0 * 255);
					pixel[1] = (uint8_t)(baryV1 * 255);
					pixel[2] = (uint8_t)(baryV2 * 255);
					pixel[3] = 0xFF;
#elif 1
					//pixel[0] = (uint8_t)(v0_persp * 255);
					//pixel[1] = (uint8_t)(v1_persp * 255);
					//pixel[2] = (uint8_t)(v2_persp * 255);
					pixel[0] = uint8_t((vary0[0] * v0_persp + vary1[0] * v1_persp + vary2[0] * v2_persp) * 255.0f);
					pixel[1] = uint8_t((vary0[1] * v0_persp + vary1[1] * v1_persp + vary2[1] * v2_persp) * 255.0f);
					pixel[2] = uint8_t((vary0[2] * v0_persp + vary1[2] * v1_persp + vary2[2] * v2_persp) * 255.0f);
					pixel[3] = 0xFF;
#else
					pixel[0] = thecol;
					pixel[1] = thecol;
					pixel[2] = thecol;
					pixel[3] = 0xFF;
#endif


				}
			}

			e0_edge_eq.IncX();
			e1_edge_eq.IncX();
			e2_edge_eq.IncX();
		}
		e0_edge_eq.NextLine();
		e1_edge_eq.NextLine();
		e2_edge_eq.NextLine();
	}
}

enum VertexClipCode
{
	X_Neg = 0x1,
	X_Pos = 0x2,

	Y_Neg = 0x4,
	Y_Pos = 0x8,

	Z_Near = 0x10,
	Z_Far = 0x20
};

static uint8_t ComputeClipMask(kt::Vec4 const& _v)
{
	uint8_t mask = 0;

	if (_v.x + _v.w < 0.0f) mask |= VertexClipCode::X_Neg;
	if (_v.x - _v.w > 0.0f) mask |= VertexClipCode::X_Pos;
	if (_v.y + _v.w < 0.0f) mask |= VertexClipCode::Y_Neg;
	if (_v.y - _v.w > 0.0f) mask |= VertexClipCode::Y_Pos;
	if (_v.z < 0.0f)		mask |= VertexClipCode::Z_Near;
	if (_v.z - _v.w > 0.0f) mask |= VertexClipCode::Z_Far;

	return mask;
}

// Each frustrum plane can turn a point into an edge (1 vert -> 2 verts). Therefore each clip plane can add 1 vertex. 9 total (including 3 from initial tri)
constexpr uint32_t CLIP_BUFFER_SIZE = 3 + 6;

constexpr uint32_t CLIP_OUT_TRI_BUFF_SIZE = (CLIP_BUFFER_SIZE - 2);

static bool SetupTri(PipelineVert const& _v0, PipelineVert const& _v1, PipelineVert const& _v2, PipelineTri& o_tri, uint32_t const _screenWidth, uint32_t const _screenHeight)
{
	kt::Vec2 const halfScreenCoords(_screenWidth*0.5f, _screenHeight*0.5f);

	kt::Vec4 v0_clip = _v0.transformedPos;
	kt::Vec4 v1_clip = _v1.transformedPos;
	kt::Vec4 v2_clip = _v2.transformedPos;

	// w = recip(w) (1/z)
	v0_clip.w = 1.0f / v0_clip.w;
	v1_clip.w = 1.0f / v1_clip.w;
	v2_clip.w = 1.0f / v2_clip.w;

	v0_clip.z *= v0_clip.w;
	v1_clip.z *= v1_clip.w;
	v2_clip.z *= v2_clip.w;

	kt::Vec2 const v0raster = kt::Vec2(v0_clip.w * v0_clip.x * halfScreenCoords.x + halfScreenCoords.x, v0_clip.w * v0_clip.y * -halfScreenCoords.y + halfScreenCoords.y);
	kt::Vec2 const v1raster = kt::Vec2(v1_clip.w * v1_clip.x * halfScreenCoords.x + halfScreenCoords.x, v1_clip.w * v1_clip.y * -halfScreenCoords.y + halfScreenCoords.y);
	kt::Vec2 const v2raster = kt::Vec2(v2_clip.w * v2_clip.x * halfScreenCoords.x + halfScreenCoords.x, v2_clip.w * v2_clip.y * -halfScreenCoords.y + halfScreenCoords.y);

	int32_t const v0_fp[2] = { int32_t(v0raster.x * c_subPixelStep + 0.5f), int32_t(v0raster.y * c_subPixelStep + 0.5f) };
	int32_t const v1_fp[2] = { int32_t(v1raster.x * c_subPixelStep + 0.5f), int32_t(v1raster.y * c_subPixelStep + 0.5f) };
	int32_t const v2_fp[2] = { int32_t(v2raster.x * c_subPixelStep + 0.5f), int32_t(v2raster.y * c_subPixelStep + 0.5f) };

	int32_t const e0_fp[2] = { v1_fp[0] - v0_fp[0], v1_fp[1] - v0_fp[1] };
	int32_t const e1_fp[2] = { v2_fp[0] - v1_fp[0], v2_fp[1] - v1_fp[1] };
	int32_t const e2_fp[2] = { v0_fp[0] - v2_fp[0], v0_fp[1] - v2_fp[1] };

	int32_t xmin, ymin;
	int32_t xmax, ymax;

	xmin = kt::Max(0, ((kt::Min(kt::Min(v0_fp[0], v1_fp[0]), v2_fp[0]) + c_subPixelMask) >> c_subPixelBits));
	ymin = kt::Max(0, ((kt::Min(kt::Min(v0_fp[1], v1_fp[1]), v2_fp[1]) + c_subPixelMask) >> c_subPixelBits));

	xmax = kt::Min((int32_t)_screenWidth, ((kt::Max(kt::Max(v0_fp[0], v1_fp[0]), v2_fp[0]) + c_subPixelMask) >> c_subPixelBits));
	ymax = kt::Min((int32_t)_screenHeight, ((kt::Max(kt::Max(v0_fp[1], v1_fp[1]), v2_fp[1]) + c_subPixelMask) >> c_subPixelBits));

	int64_t triArea2_fp = (v2_fp[0] - v0_fp[0]) * (v1_fp[1] - v0_fp[1]) - (v2_fp[1] - v0_fp[1]) * (v1_fp[0] - v0_fp[0]);

	if (triArea2_fp <= 0)
	{
		return false;
	}


	o_tri.xmin = xmin;
	o_tri.xmax = xmax;
	o_tri.ymin = ymin;
	o_tri.ymax = ymax;

	o_tri.halfRecipTriArea_fp = 1.0f / triArea2_fp;

	o_tri.v0_fp[0] = v0_fp[0];
	o_tri.v0_fp[1] = v0_fp[1];

	o_tri.v1_fp[0] = v1_fp[0];
	o_tri.v1_fp[1] = v1_fp[1];

	o_tri.v2_fp[0] = v2_fp[0];
	o_tri.v2_fp[1] = v2_fp[1];

	o_tri.verts[0].zOverW = v0_clip.z;
	o_tri.verts[1].zOverW = v1_clip.z;
	o_tri.verts[2].zOverW = v2_clip.z;

	o_tri.verts[0].recipW = v0_clip.w;
	o_tri.verts[1].recipW = v1_clip.w;
	o_tri.verts[2].recipW = v2_clip.w;

	memcpy(o_tri.verts[0].varyings, _v0.varyings, sizeof(_v0.varyings));
	memcpy(o_tri.verts[1].varyings, _v1.varyings, sizeof(_v1.varyings));
	memcpy(o_tri.verts[2].varyings, _v2.varyings, sizeof(_v2.varyings));
	return true;
}

void ClipAndSetup(PipelineVert clipSpaceVerts[3], PipelineTri o_tris[CLIP_OUT_TRI_BUFF_SIZE], uint32_t& o_numOutTris)
{
	uint32_t maskOr = 0;

	uint8_t clipv0 = ComputeClipMask(clipSpaceVerts[0].transformedPos);
	uint8_t clipv1 = ComputeClipMask(clipSpaceVerts[1].transformedPos);
	uint8_t clipv2 = ComputeClipMask(clipSpaceVerts[2].transformedPos);

	maskOr = clipv0 | clipv1 | clipv2;

	if (maskOr == 0)
	{
		if (SetupTri(clipSpaceVerts[0], clipSpaceVerts[1], clipSpaceVerts[2], o_tris[0], 1280, 720)) // hey dont hardcode that!
		{
			o_numOutTris = 1;
		}
		else
		{
			o_numOutTris = 0;
		}
		return;
	}

	if (clipv0 & clipv1 & clipv2)
	{
		// If clip AND mask has any bits set, all verts are the wrong side of a clip plane, so the whole triangle can be culled.
		o_numOutTris = 0;
		return;
	}

	static kt::Vec4 const c_clipPlanes[6] =
	{
		{ 1.0f, 0.0f, 0.0f, 1.0f }, // X_Neg
		{ -1.0f, 0.0f, 0.0f, 1.0f  }, // X_Pos

		{ 0.0f, 1.0f, 0.0f, 1.0f }, // Y_Neg
		{ 0.0f, -1.0f, 0.0f, 1.0f  }, // Y_Pos

		{ 0.0f, 0.0f, 1.0f, 1.0f }, // Z_Near
		{ 0.0f, 0.0f, -1.0f, 1.0f }, // Z_Far
	};

	// Todo: don't copy all this stuff around.

	auto copyVert = [](PipelineVert& out, PipelineVert const& _in) 
	{
		memcpy(&out, &_in, sizeof(PipelineVert));
	};

	PipelineVert input[CLIP_BUFFER_SIZE];
	copyVert(input[0], clipSpaceVerts[0]);
	copyVert(input[1], clipSpaceVerts[1]);
	copyVert(input[2], clipSpaceVerts[2]);

	PipelineVert output[CLIP_BUFFER_SIZE];

	uint32_t numInputVerts = 3;
	uint32_t numOutputVerts = 0;

	do
	{
		uint32_t clipIdx = kt::Cnttz(maskOr);
		maskOr ^= (1 << clipIdx);

		kt::Vec4 const& clipPlane = c_clipPlanes[clipIdx];

		PipelineVert const* v0 = &input[numInputVerts - 1];
		float v0_dot = kt::Dot(clipPlane, v0->transformedPos);

		for (uint32_t nextClipVertIdx = 0; nextClipVertIdx < numInputVerts; ++nextClipVertIdx)
		{
			PipelineVert const*const v1 = &input[nextClipVertIdx];
			float const v1_dot = kt::Dot(clipPlane, v1->transformedPos);
			bool const v0_inside = v0_dot >= 0.0f;
			bool const v1_inside = v1_dot >= 0.0f;

			if (v0_inside)
			{
				KT_ASSERT(numOutputVerts < CLIP_BUFFER_SIZE);
				copyVert(output[numOutputVerts++], *v0);
			}

			if (v0_inside ^ v1_inside)
			{
				if (v1_inside)
				{
					KT_ASSERT(numOutputVerts < CLIP_BUFFER_SIZE);
					float const tInterp = v1_dot / (v1_dot - v0_dot);
					PipelineVert* out = &output[numOutputVerts++];
					out->transformedPos = kt::Lerp(v1->transformedPos, v0->transformedPos, tInterp);

					for (uint32_t varI = 0; varI < c_maxVaryings; ++varI)
					{
						out->varyings[varI] = kt::Lerp(v1->varyings[varI], v0->varyings[varI], tInterp);
					}
				}
				else
				{
					KT_ASSERT(numOutputVerts < CLIP_BUFFER_SIZE);
					float const tInterp = v0_dot / (v0_dot - v1_dot);
					PipelineVert* out = &output[numOutputVerts++];
					out->transformedPos = kt::Lerp(v0->transformedPos, v1->transformedPos, tInterp);
					for (uint32_t varI = 0; varI < c_maxVaryings; ++varI)
					{
						out->varyings[varI] = kt::Lerp(v0->varyings[varI], v1->varyings[varI], tInterp);
					}
				}
			}

			v0_dot = v1_dot;
			v0 = v1;
		}

		if (!numOutputVerts)
		{
			o_numOutTris = 0;
			return;
		}

		memcpy(input, output, sizeof(PipelineVert) * numOutputVerts); // this is slow - should have pointer indirection into buffers so no copy is needed.
		numInputVerts = numOutputVerts;
		numOutputVerts = 0;
	} while (maskOr);

	// Fan triangulation
	o_numOutTris = 0;
	for (uint32_t i = 2; i < numInputVerts; ++i)
	{
		KT_ASSERT(o_numOutTris <= CLIP_OUT_TRI_BUFF_SIZE);
		PipelineTri& tri = o_tris[o_numOutTris];
		if (SetupTri(input[0], input[i - 1], input[i], o_tris[o_numOutTris], 1280, 720)) // ahhh hardcoded
		{
			++o_numOutTris;
		}
	}
}


void RasterClippedTri(FrameBuffer& _buffer, DepthBuffer& _depthBuffer, kt::Vec4 const& _v0, kt::Vec4 const& _v1, kt::Vec4 const& _v2)
{
	kt::Vec2 const halfScreenCoords(_buffer.width*0.5f, _buffer.height*0.5f);

	kt::Vec4 v0_clip = _v0;
	kt::Vec4 v1_clip = _v1;
	kt::Vec4 v2_clip = _v2;

	// w = recip(w) (1/z)
	v0_clip.w = 1.0f / v0_clip.w;
	v1_clip.w = 1.0f / v1_clip.w;
	v2_clip.w = 1.0f / v2_clip.w;

	v0_clip.z *= v0_clip.w;
	v1_clip.z *= v1_clip.w;
	v2_clip.z *= v2_clip.w;

	kt::Vec2 const v0raster = kt::Vec2(v0_clip.w * v0_clip.x * halfScreenCoords.x + halfScreenCoords.x, v0_clip.w * v0_clip.y * -halfScreenCoords.y + halfScreenCoords.y);
	kt::Vec2 const v1raster = kt::Vec2(v1_clip.w * v1_clip.x * halfScreenCoords.x + halfScreenCoords.x, v1_clip.w * v1_clip.y * -halfScreenCoords.y + halfScreenCoords.y);
	kt::Vec2 const v2raster = kt::Vec2(v2_clip.w * v2_clip.x * halfScreenCoords.x + halfScreenCoords.x, v2_clip.w * v2_clip.y * -halfScreenCoords.y + halfScreenCoords.y);

	int32_t const v0_fp[2] = { int32_t(v0raster.x * c_subPixelStep + 0.5f), int32_t(v0raster.y * c_subPixelStep + 0.5f) };
	int32_t const v1_fp[2] = { int32_t(v1raster.x * c_subPixelStep + 0.5f), int32_t(v1raster.y * c_subPixelStep + 0.5f) };
	int32_t const v2_fp[2] = { int32_t(v2raster.x * c_subPixelStep + 0.5f), int32_t(v2raster.y * c_subPixelStep + 0.5f) };

	int32_t const e0_fp[2] = { v1_fp[0] - v0_fp[0], v1_fp[1] - v0_fp[1] };
	int32_t const e1_fp[2] = { v2_fp[0] - v1_fp[0], v2_fp[1] - v1_fp[1] };
	int32_t const e2_fp[2] = { v0_fp[0] - v2_fp[0], v0_fp[1] - v2_fp[1] };

	int32_t xmin, ymin;
	int32_t xmax, ymax;

	xmin = kt::Max(0, ((kt::Min(kt::Min(v0_fp[0], v1_fp[0]), v2_fp[0]) + c_subPixelMask) >> c_subPixelBits));
	ymin = kt::Max(0, ((kt::Min(kt::Min(v0_fp[1], v1_fp[1]), v2_fp[1]) + c_subPixelMask) >> c_subPixelBits));

	xmax = kt::Min((int32_t)_buffer.width, ((kt::Max(kt::Max(v0_fp[0], v1_fp[0]), v2_fp[0]) + c_subPixelMask) >> c_subPixelBits));
	ymax = kt::Min((int32_t)_buffer.height, ((kt::Max(kt::Max(v0_fp[1], v1_fp[1]), v2_fp[1]) + c_subPixelMask) >> c_subPixelBits));

	int64_t triArea2_fp = (v2_fp[0] - v0_fp[0]) * (v1_fp[1] - v0_fp[1]) - (v2_fp[1] - v0_fp[1]) * (v1_fp[0] - v0_fp[0]);

	if (triArea2_fp <= 0.0f)
	{
		return;
	}

	PipelineTri tri;
	tri.xmin = xmin;
	tri.xmax = xmax;
	tri.ymin = ymin;
	tri.ymax = ymax;

	tri.halfRecipTriArea_fp = 1.0f / triArea2_fp;

	tri.v0_fp[0] = v0_fp[0];
	tri.v0_fp[1] = v0_fp[1];

	tri.v1_fp[0] = v1_fp[0];
	tri.v1_fp[1] = v1_fp[1];

	tri.v2_fp[0] = v2_fp[0];
	tri.v2_fp[1] = v2_fp[1];

	tri.verts[0].zOverW = v0_clip.z;
	tri.verts[1].zOverW = v1_clip.z;
	tri.verts[2].zOverW = v2_clip.z;

	tri.verts[0].recipW = v0_clip.w;
	tri.verts[1].recipW = v1_clip.w;
	tri.verts[2].recipW = v2_clip.w;

	RasterTransformedTri(_buffer, _depthBuffer, tri);
}

void SetupAndRasterTriTest(FrameBuffer& _buffer, DepthBuffer& _depthBuffer, kt::Mat4 const& _mtx, kt::Vec3 const& _v0, kt::Vec3 const& _v1, kt::Vec3 const& _v2)
{
	RasterContext ctx;
	ctx.depthBuffer = &_depthBuffer;
	ctx.frameBuffer = &_buffer;
	ctx.vpHeight = _buffer.height;
	ctx.vpWidth = _buffer.width;

	PipelineTri outTris[CLIP_OUT_TRI_BUFF_SIZE];
	uint32_t numClipTris;
	
	PipelineVert verts[3];
	verts[0].transformedPos = _mtx * kt::Vec4(_v0, 1.0f);
	verts[1].transformedPos = _mtx * kt::Vec4(_v1, 1.0f);
	verts[2].transformedPos = _mtx * kt::Vec4(_v2, 1.0f);

	verts[0].varyings[0] = 1.0f;
	verts[0].varyings[1] = 0.0f;
	verts[0].varyings[2] = 0.0f;

	verts[1].varyings[0] = 0.0f;
	verts[1].varyings[1] = 1.0f;
	verts[1].varyings[2] = 0.0f;

	verts[2].varyings[0] = 0.0f;
	verts[2].varyings[1] = 0.0f;
	verts[2].varyings[2] = 1.0f;

	ClipAndSetup(verts, outTris, numClipTris);

	for (uint32_t i = 0; i < numClipTris; ++i)
	{
		//RasterTransformedTri(*ctx.frameBuffer, *ctx.depthBuffer, outTris[i]);
		RasterTransformedTri_WithBlocks(ctx, outTris[i]);
	}

#if 0
	// w = recip(w) (1/z)
	v0_clip.w = 1.0f / v0_clip.w;
	v1_clip.w = 1.0f / v1_clip.w;
	v2_clip.w = 1.0f / v2_clip.w;

	v0_clip.z *= v0_clip.w;
	v1_clip.z *= v1_clip.w;
	v2_clip.z *= v2_clip.w;

	kt::Vec2 const v0raster = kt::Vec2(v0_clip.w * v0_clip.x * halfScreenCoords.x + halfScreenCoords.x, v0_clip.w * v0_clip.y * -halfScreenCoords.y + halfScreenCoords.y);
	kt::Vec2 const v1raster = kt::Vec2(v1_clip.w * v1_clip.x * halfScreenCoords.x + halfScreenCoords.x, v1_clip.w * v1_clip.y * -halfScreenCoords.y + halfScreenCoords.y);
	kt::Vec2 const v2raster = kt::Vec2(v2_clip.w * v2_clip.x * halfScreenCoords.x + halfScreenCoords.x, v2_clip.w * v2_clip.y * -halfScreenCoords.y + halfScreenCoords.y);

	int32_t const v0_fp[2] = { int32_t(v0raster.x * c_subPixelStep + 0.5f), int32_t(v0raster.y * c_subPixelStep + 0.5f) };
	int32_t const v1_fp[2] = { int32_t(v1raster.x * c_subPixelStep + 0.5f), int32_t(v1raster.y * c_subPixelStep + 0.5f) };
	int32_t const v2_fp[2] = { int32_t(v2raster.x * c_subPixelStep + 0.5f), int32_t(v2raster.y * c_subPixelStep + 0.5f) };

	int32_t const e0_fp[2] = { v1_fp[0] - v0_fp[0], v1_fp[1] - v0_fp[1] };
	int32_t const e1_fp[2] = { v2_fp[0] - v1_fp[0], v2_fp[1] - v1_fp[1] };
	int32_t const e2_fp[2] = { v0_fp[0] - v2_fp[0], v0_fp[1] - v2_fp[1] };

	int32_t xmin, ymin;
	int32_t xmax, ymax;

	xmin = kt::Max(0, ((kt::Min(kt::Min(v0_fp[0], v1_fp[0]), v2_fp[0]) + c_subPixelMask) >> c_subPixelBits));
	ymin = kt::Max(0, ((kt::Min(kt::Min(v0_fp[1], v1_fp[1]), v2_fp[1]) + c_subPixelMask) >> c_subPixelBits));

	xmax = kt::Min((int32_t)_buffer.width, ((kt::Max(kt::Max(v0_fp[0], v1_fp[0]), v2_fp[0]) + c_subPixelMask) >> c_subPixelBits));
	ymax = kt::Min((int32_t)_buffer.height, ((kt::Max(kt::Max(v0_fp[1], v1_fp[1]), v2_fp[1]) + c_subPixelMask) >> c_subPixelBits));

	int64_t const triArea2_fp = (v2_fp[0] - v0_fp[0]) * (v1_fp[1] - v0_fp[1]) - (v2_fp[1] - v0_fp[1]) * (v1_fp[0] - v0_fp[0]);

	if (triArea2_fp <= 0.0f)
	{
		return;
	}

	Transformed_Tri tri;
	tri.xmin = xmin;
	tri.xmax = xmax;
	tri.ymin = ymin;
	tri.ymax = ymax;

	tri.halfRecipTriArea_fp = 1.0f / triArea2_fp;

	tri.v0_fp[0] = v0_fp[0];
	tri.v0_fp[1] = v0_fp[1];

	tri.v1_fp[0] = v1_fp[0];
	tri.v1_fp[1] = v1_fp[1];

	tri.v2_fp[0] = v2_fp[0];
	tri.v2_fp[1] = v2_fp[1];

	tri.verts[0].zOverW = v0_clip.z;
	tri.verts[1].zOverW = v1_clip.z;
	tri.verts[2].zOverW = v2_clip.z;

	tri.verts[0].recipW = v0_clip.w;
	tri.verts[1].recipW = v1_clip.w;
	tri.verts[2].recipW = v2_clip.w;

	RasterTransformedTri(_buffer, _depthBuffer, tri);
#endif
}


void DrawSerial_Test(FrameBuffer& _buffer, DepthBuffer& _depthBuffer, kt::Mat4 const& _mtx, Renderer::DrawCall const& _call)
{
	RasterContext ctx;
	ctx.depthBuffer = &_depthBuffer;
	ctx.frameBuffer = &_buffer;
	ctx.vpHeight = _buffer.height;
	ctx.vpWidth = _buffer.width;

	KT_ASSERT(_call.m_indexBuffer.m_num % 3 == 0);
	for (uint32_t i = 0; i < _call.m_indexBuffer.m_num; i += 3)
	{
		PipelineTri outTris[CLIP_OUT_TRI_BUFF_SIZE];
		PipelineVert verts[3];
		uint32_t numClipTris;

		uint32_t idx0, idx1, idx2;

		if (_call.m_indexBuffer.m_stride == 4)
		{
			idx0 = ((uint32_t*)_call.m_indexBuffer.m_ptr)[i];
			idx1 = ((uint32_t*)_call.m_indexBuffer.m_ptr)[i + 1];
			idx2 = ((uint32_t*)_call.m_indexBuffer.m_ptr)[i + 2];
		}
		else if (_call.m_indexBuffer.m_stride == 2)
		{
			idx0 = (uint32_t)((uint16_t*)_call.m_indexBuffer.m_ptr)[i];
			idx1 = (uint32_t)((uint16_t*)_call.m_indexBuffer.m_ptr)[i + 1];
			idx2 = (uint32_t)((uint16_t*)_call.m_indexBuffer.m_ptr)[i + 2];
		}

		KT_ASSERT(idx0 < _call.m_positionBuffer.m_num);
		KT_ASSERT(idx1 < _call.m_positionBuffer.m_num);
		KT_ASSERT(idx2 < _call.m_positionBuffer.m_num);

		uint8_t const* pos0 = (uint8_t*)_call.m_positionBuffer.m_ptr + _call.m_positionBuffer.m_stride * idx0;
		uint8_t const* pos1 = (uint8_t*)_call.m_positionBuffer.m_ptr + _call.m_positionBuffer.m_stride * idx1;
		uint8_t const* pos2 = (uint8_t*)_call.m_positionBuffer.m_ptr + _call.m_positionBuffer.m_stride * idx2;

		uint8_t const* attr0 = (uint8_t*)_call.m_attributeBuffer.m_ptr + _call.m_attributeBuffer.m_stride * idx0;
		uint8_t const* attr1 = (uint8_t*)_call.m_attributeBuffer.m_ptr + _call.m_attributeBuffer.m_stride * idx1;
		uint8_t const* attr2 = (uint8_t*)_call.m_attributeBuffer.m_ptr + _call.m_attributeBuffer.m_stride * idx2;


		verts[0].transformedPos = _mtx * kt::Vec4(*(kt::Vec3*)pos0, 1.0f);
		verts[1].transformedPos = _mtx * kt::Vec4(*(kt::Vec3*)pos1, 1.0f);
		verts[2].transformedPos = _mtx * kt::Vec4(*(kt::Vec3*)pos2, 1.0f);

		memcpy(verts[0].varyings, attr0, _call.m_attributeBuffer.m_stride);
		memcpy(verts[1].varyings, attr1, _call.m_attributeBuffer.m_stride);
		memcpy(verts[2].varyings, attr2, _call.m_attributeBuffer.m_stride);

		ClipAndSetup(verts, outTris, numClipTris);

		for (uint32_t i = 0; i < numClipTris; ++i)
		{
			//RasterTransformedTri(*ctx.frameBuffer, *ctx.depthBuffer, outTris[i]);
			RasterTransformedTri_WithBlocksAndShader(ctx, outTris[i], _call);
		}
	}
}


void DepthBuffer::Init(kt::IAllocator* _allocator, uint32_t _width, uint32_t _height)
{
	if (allocator && ptr)
	{
		allocator->Free(ptr);
	}

	ptr = (float*)_allocator->Alloc(sizeof(float) * _width * _height, 16);
	allocator = _allocator;
	width = _width;
	height = _height;
}

DepthBuffer::~DepthBuffer()
{
	if (allocator && ptr)
	{
		allocator->Free(ptr);
	}
}


}

}
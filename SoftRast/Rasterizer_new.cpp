#include "Rasterizer_new.h"
#include "Binning.h"
#include "Renderer.h"
#include "kt/Sort.h"

namespace sr
{

struct FragmentStream
{

};

struct FragmentBuffer
{
	enum Flags : uint16_t
	{
		None,
		Full8x8Block 
	};

	struct Frag
	{
		uint16_t x : 7;
		uint16_t y : 7;
		uint16_t flags : 2;

		union 
		{
			struct  
			{
				uint16_t chunkIdx;
				uint16_t triIdx;
			};
			uint32_t packedChunkTriIdx;
		};

	};

	static_assert((1 << 7) <= Config::c_binHeight, "Cant fit height in 7 bits");
	static_assert((1 << 7) <= Config::c_binWidth, "Cant fit width in 7 bits");

	// Todo: chunk the buffer
	static uint32_t constexpr c_maxFragments = 2048 * 16;

	KT_ALIGNAS(32) Frag m_fragments[c_maxFragments];
	uint32_t m_numFragments = 0;

	float m_interpolants[Config::c_maxVaryings * c_maxFragments];

	void Reset()
	{
		m_numFragments = 0;
	}
};

struct BlockInterpolants8x8
{
	__m256 attrib_eval[Config::c_maxVaryings];
	__m256 attrib_dy[Config::c_maxVaryings];

	__m256 w_eval;
	__m256 dwdy;

	__m256 zw_eval;
	__m256 dzwdy;
};

struct EdgeEquations8x8
{
	__m256i eval[3];

	__m256i dx[3];
};


void GatherPartialBlockState
(
	BlockInterpolants8x8& o_block,
	EdgeEquations8x8& o_edges,
	BinChunk const& _chunk,
	uint32_t _xIdxTileRelative,
	uint32_t _yIdxTileRelative,
	uint32_t _triIdx
)
{
	__m256 const ramp = _mm256_setr_ps(0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f);
	__m256i const rampi = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);

	__m256i const x_offs = _mm256_set1_epi32(_xIdxTileRelative);
	__m256i const y_offs = _mm256_set1_epi32(_yIdxTileRelative);

	__m256 const fx_offs = _mm256_cvtepi32_ps(x_offs);
	__m256 const fy_offs = _mm256_cvtepi32_ps(y_offs);

	o_block.w_eval = _mm256_broadcast_ss(&_chunk.m_recipW[_triIdx].c0);
	__m256 const dwdx = _mm256_broadcast_ss(&_chunk.m_recipW[_triIdx].dx);
	o_block.dwdy = _mm256_broadcast_ss(&_chunk.m_recipW[_triIdx].dy);

	o_block.zw_eval = _mm256_broadcast_ss(&_chunk.m_zOverW[_triIdx].c0);
	__m256 const dzwdx = _mm256_broadcast_ss(&_chunk.m_zOverW[_triIdx].dx);
	o_block.dzwdy = _mm256_broadcast_ss(&_chunk.m_zOverW[_triIdx].dy);

	// Offset to this 8x8 block
	o_block.w_eval = _mm256_fmadd_ps(dwdx, fx_offs, _mm256_fmadd_ps(o_block.dwdy, fy_offs, o_block.w_eval));
	o_block.zw_eval = _mm256_fmadd_ps(dzwdx, fx_offs, _mm256_fmadd_ps(o_block.dzwdy, fy_offs, o_block.zw_eval));

	o_block.w_eval = _mm256_fmadd_ps(dwdx, ramp, o_block.w_eval);
	o_block.zw_eval = _mm256_fmadd_ps(dzwdx, ramp, o_block.zw_eval);

	//for (uint32_t i = 0; i < Config::c_maxVaryings; ++i)
	//{
	//	o_block.attrib_eval[i]	=	_mm256_broadcast_ss(&_chunk.m_attribsC[_triIdx * _chunk.m_attribsPerTri + i]);
	//	__m256 const attribdx	=	_mm256_broadcast_ss(&_chunk.m_attribsDx[_triIdx * _chunk.m_attribsPerTri + i]);
	//	o_block.attrib_dy[i]	=	_mm256_broadcast_ss(&_chunk.m_attribsDy[_triIdx * _chunk.m_attribsPerTri + i]);
	//	o_block.attrib_eval[i]	=	_mm256_fmadd_ps(attribdx, fx_offs, _mm256_fmadd_ps(o_block.attrib_dy[i], fy_offs, o_block.attrib_eval[i]));
	//	o_block.attrib_eval[i]	=	_mm256_fmadd_ps(attribdx, ramp, o_block.attrib_eval[i]);
	//}

	for (uint32_t i = 0; i < Config::c_maxVaryings; ++i)
	{
		o_block.attrib_eval[i] = _mm256_broadcast_ss(&_chunk.m_attribsC[_triIdx * _chunk.m_attribsPerTri + i]);
		__m256 const attribdx = _mm256_broadcast_ss(&_chunk.m_attribsDx[_triIdx * _chunk.m_attribsPerTri + i]);
		o_block.attrib_dy[i] = _mm256_broadcast_ss(&_chunk.m_attribsDy[_triIdx * _chunk.m_attribsPerTri + i]);
		o_block.attrib_eval[i] = _mm256_fmadd_ps(attribdx, fx_offs, _mm256_fmadd_ps(o_block.attrib_dy[i], fy_offs, o_block.attrib_eval[i]));
		o_block.attrib_eval[i] = _mm256_fmadd_ps(attribdx, ramp, o_block.attrib_eval[i]);
	}


	BinChunk::EdgeEq const& chunkEq = _chunk.m_edgeEq[_triIdx];

	for (uint32_t i = 0; i < 3; ++i)
	{
		int32_t const initial = chunkEq.c[i] + chunkEq.dx[i] * _yIdxTileRelative + chunkEq.dy[i] * _xIdxTileRelative;
		__m256i initial8 = _mm256_set1_epi32(initial);
		o_edges.eval[i] = _mm256_add_epi32(initial8, _mm256_mullo_epi32(rampi, _mm256_set1_epi32(chunkEq.dy[i])));
		o_edges.dx[i] = _mm256_set1_epi32(chunkEq.dx[i]);
	}
}

#define TEST_IT 1

static void ShadePartialBlockAVX_8x8
(
	DrawCall const& _call,
	ColourTile* _colour,
	DepthTile* _depth,
	int32_t _xTileRelative,
	int32_t _yTileRelative,
	BlockInterpolants8x8& _block,
	EdgeEquations8x8& _edges
)
{
	__m256i e0 = _edges.eval[0];
	__m256i e1 = _edges.eval[1];
	__m256i e2 = _edges.eval[2];

	__m256 zOverW = _block.zw_eval;
	__m256 recipW = _block.w_eval;

	for (uint32_t i = 0; i < 8; ++i)
	{
		__m256i edgeMask = _mm256_or_si256(_mm256_or_si256(e0, e1), e2);
		edgeMask = _mm256_srai_epi32(edgeMask, 31); // or sign bits, arith shift right so sign bit = ~0 and no sign bit = 0
		edgeMask = _mm256_xor_si256(edgeMask, _mm256_cmpeq_epi32(edgeMask, edgeMask)); // flip bits so negative (~0) is 0.
		
#if TEST_IT
		if (_mm256_movemask_ps(_mm256_castsi256_ps(edgeMask)))
#endif
		{
			float* depthPtr = _depth->m_depth + (_xTileRelative + (_yTileRelative + i) * Config::c_binWidth);
			__m256 const depthGather = _mm256_loadu_ps(depthPtr);

			__m256 const w = _mm256_div_ps(_mm256_set1_ps(1.0f), recipW);

			__m256 depthCmpMask = _mm256_and_ps(_mm256_cmp_ps(zOverW, _mm256_setzero_ps(), _CMP_GT_OQ), _mm256_cmp_ps(depthGather, zOverW, _CMP_GT_OQ));
			depthCmpMask = _mm256_and_ps(_mm256_castsi256_ps(edgeMask), depthCmpMask);
			
#if TEST_IT
			if (int32_t maskReg = _mm256_movemask_ps(depthCmpMask))
#endif
			{
				__m256 const newDepth = _mm256_blendv_ps(depthGather, zOverW, depthCmpMask);
				_mm256_storeu_ps(depthPtr, newDepth);

				KT_ALIGNAS(32) float colourRGBA[4 * 8];

				__m256 interpolants[Config::c_maxVaryings];
				for (uint32_t i = 0; i < Config::c_maxVaryings; ++i)
				{
					interpolants[i] = _mm256_mul_ps(_block.attrib_eval[i], w);
				}

				// BROKEN
				//_call.m_pixelShader(_call.m_pixelUniforms, interpolants, colourRGBA, depthCmpMask);

				uint8_t* pixelBegin = &_colour->m_colour[(4 * _xTileRelative) + (_yTileRelative + i) * Config::c_binWidth * 4];

				for (uint32_t pixIdx = 0; pixIdx < 8; ++pixIdx)
				{
#if TEST_IT
					if (maskReg & (1 << pixIdx))
#endif
					{
#if 1
						pixelBegin[pixIdx * 4 + 0] = uint8_t(kt::Min(1.0f, colourRGBA[pixIdx * 4]) * 255.0f);
						pixelBegin[pixIdx * 4 + 1] = uint8_t(kt::Min(1.0f, colourRGBA[pixIdx * 4 + 1]) * 255.0f);
						pixelBegin[pixIdx * 4 + 2] = uint8_t(kt::Min(1.0f, colourRGBA[pixIdx * 4 + 2]) * 255.0f);
						pixelBegin[pixIdx * 4 + 3] = 255;
#else
						pixelBegin[pixIdx * 4 + 0] = 255;
						pixelBegin[pixIdx * 4 + 1] = 255;
						pixelBegin[pixIdx * 4 + 2] = 255;
						pixelBegin[pixIdx * 4 + 3] = 255;
#endif
					}
				}
			}
		}

		e0 = _mm256_add_epi32(e0, _edges.dx[0]);
		e1 = _mm256_add_epi32(e1, _edges.dx[1]);
		e2 = _mm256_add_epi32(e2, _edges.dx[2]);

		zOverW = _mm256_add_ps(zOverW, _block.dzwdy);
		recipW = _mm256_add_ps(recipW, _block.dwdy);

		// TOdo: slow
		for (uint32_t k = 0; k < Config::c_maxVaryings; ++k)
		{
			_block.attrib_eval[k] = _mm256_add_ps(_block.attrib_eval[k], _block.attrib_dy[k]);
		}
	}
}

// Todo: naming is really confusing as edge eq dx*y and dy*x are used but the opposite for planes. Should maybe clarify naming.
struct EdgeEquations8x8_New
{
	__m256i tileTopLeftEdge[3];
	__m256i dx[3];
	__m256i dy[3];
};

struct ZOverW8x8
{
	__m256 tileTopLeft;
	__m256 dx;
	__m256 dy;
};

static uint64_t ComputeBlockMask8x8
(
	EdgeEquations8x8_New const& _edges,
	ZOverW8x8 const& _zOverW,
	DrawCall const& _call, 
	DepthTile* _depth, 
	int32_t _xTileRelative, 
	int32_t _yTileRelative
)
{
	uint64_t mask8x8 = 0;

	__m256i edges[3];

	__m256i const xTileSimd = _mm256_set1_epi32(_xTileRelative);
	__m256i const yTileSimd = _mm256_set1_epi32(_yTileRelative);

	__m256 zOverW = _mm256_fmadd_ps(_mm256_cvtepi32_ps(yTileSimd), _zOverW.dy, _zOverW.tileTopLeft);
	zOverW = _mm256_add_ps(zOverW, _mm256_mul_ps(_mm256_cvtepi32_ps(xTileSimd), _zOverW.dx));

	// offset from top left tile to this 8x8 block

	for (uint32_t i = 0; i < 3; ++i)
	{
		edges[i] = _mm256_add_epi32(_edges.tileTopLeftEdge[i], _mm256_add_epi32(_mm256_mullo_epi32(yTileSimd, _edges.dx[i]), _mm256_mullo_epi32(xTileSimd, _edges.dy[i])));
	}

	// Todo: look at other coverage mask generation approaches. Eg LUT based.

	// Manually unroll ?
	for (uint32_t i = 0; i < 8; ++i)
	{
		__m256i edgeMask = _mm256_or_si256(_mm256_or_si256(edges[0], edges[1]), edges[2]);
		edgeMask = _mm256_srai_epi32(edgeMask, 31); // or sign bits, arith shift right so sign bit = ~0 and no sign bit = 0
		edgeMask = _mm256_xor_si256(edgeMask, _mm256_cmpeq_epi32(edgeMask, edgeMask)); // flip bits so negative (~0) is 0.

		// test Z
		float* depthPtr = _depth->m_depth + (_xTileRelative + (_yTileRelative + i) * Config::c_binWidth);
		__m256 const depthGather = _mm256_loadu_ps(depthPtr);

		__m256 depthCmpMask = _mm256_and_ps(_mm256_cmp_ps(zOverW, _mm256_setzero_ps(), _CMP_GT_OQ), _mm256_cmp_ps(depthGather, zOverW, _CMP_GT_OQ));
		depthCmpMask = _mm256_and_ps(_mm256_castsi256_ps(edgeMask), depthCmpMask);

		uint64_t const laneMask = _mm256_movemask_ps(depthCmpMask);

		__m256 const newDepth = _mm256_blendv_ps(depthGather, zOverW, depthCmpMask);
		_mm256_storeu_ps(depthPtr, newDepth);

		mask8x8 |= (laneMask << (i * 8));

		edges[0] = _mm256_add_epi32(edges[0], _edges.dx[0]);
		edges[1] = _mm256_add_epi32(edges[1], _edges.dx[1]);
		edges[2] = _mm256_add_epi32(edges[2], _edges.dx[2]);

		zOverW = _mm256_add_ps(zOverW, _zOverW.dy);
	}

	return mask8x8;
}

static void RasterizeTrisInBin_OutputFragments(DrawCall const& _call, DepthTile* _depth, BinChunk const& _chunk, uint32_t _chunkIdx, FragmentBuffer& o_buffer)
{
	__m256i const rampi = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
	__m256 const rampf = _mm256_setr_ps(0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f);

	for (uint32_t triIdx = 0; triIdx < _chunk.m_numTris; ++triIdx)
	{
		BinChunk::EdgeEq const& edges = _chunk.m_edgeEq[triIdx];

		uint32_t const xBlockBegin = edges.blockMinX & ~7;
		uint32_t const yBlockBegin = edges.blockMinY & ~7;

		uint32_t const xBlockEnd = edges.blockMaxX;
		uint32_t const yBlockEnd = edges.blockMaxY;

		EdgeEquations8x8_New blockEdgesSimd;

		ZOverW8x8 zOverWPlaneSimd;

		// Todo: if we only ever do early Z should maybe pack zOverW plane and edge eqs for cache perf as they are always fetched together?
		// although probably want toggleable earlyZ for blending
		zOverWPlaneSimd.dx = _mm256_broadcast_ss(&_chunk.m_zOverW[triIdx].dx);
		zOverWPlaneSimd.dy = _mm256_broadcast_ss(&_chunk.m_zOverW[triIdx].dy);
		zOverWPlaneSimd.tileTopLeft = _mm256_fmadd_ps(rampf, zOverWPlaneSimd.dx, _mm256_broadcast_ss(&_chunk.m_zOverW[triIdx].c0));

		for (uint32_t i = 0; i < 3; ++i)
		{
			blockEdgesSimd.dx[i] = _mm256_set1_epi32(edges.dx[i]);
			blockEdgesSimd.dy[i] = _mm256_set1_epi32(edges.dy[i]);
			blockEdgesSimd.tileTopLeftEdge[i] = _mm256_add_epi32(_mm256_mullo_epi32(blockEdgesSimd.dy[i], rampi), _mm256_set1_epi32(edges.c[i]));
		}

		for (uint32_t yBlock = yBlockBegin; yBlock < yBlockEnd; yBlock += 8)
		{
			for (uint32_t xBlock = xBlockBegin; xBlock < xBlockEnd; xBlock += 8)
			{
				int32_t const binScreenX0 = xBlock;
				int32_t const binScreenX1 = xBlock + Config::c_binWidth;

				int32_t const binScreenY0 = yBlock;
				int32_t const binScreenY1 = yBlock + Config::c_binHeight;

				// Todo: slow, can offset edges and just test upper corner
				int32_t const e0_x0y0 = edges.c[0] + edges.dy[0] * binScreenX0 + edges.dx[0] * binScreenY0;
				int32_t const e0_x0y1 = edges.c[0] + edges.dy[0] * binScreenX0 + edges.dx[0] * binScreenY1;
				int32_t const e0_x1y0 = edges.c[0] + edges.dy[0] * binScreenX1 + edges.dx[0] * binScreenY0;
				int32_t const e0_x1y1 = edges.c[0] + edges.dy[0] * binScreenX1 + edges.dx[0] * binScreenY1;

				uint32_t const e0_allOut = (e0_x0y0 > 0) | ((e0_x0y1 > 0) << 1) | ((e0_x1y0 > 0) << 2) | ((e0_x1y1 > 0) << 3);

				int32_t const e1_x0y0 = edges.c[1] + edges.dy[1] * binScreenX0 + edges.dx[1] * binScreenY0;
				int32_t const e1_x0y1 = edges.c[1] + edges.dy[1] * binScreenX0 + edges.dx[1] * binScreenY1;
				int32_t const e1_x1y0 = edges.c[1] + edges.dy[1] * binScreenX1 + edges.dx[1] * binScreenY0;
				int32_t const e1_x1y1 = edges.c[1] + edges.dy[1] * binScreenX1 + edges.dx[1] * binScreenY1;

				uint32_t const e1_allOut = (e1_x0y0 > 0) | ((e1_x0y1 > 0) << 1) | ((e1_x1y0 > 0) << 2) | ((e1_x1y1 > 0) << 3);

				int32_t const e2_x0y0 = edges.c[2] + edges.dy[2] * binScreenX0 + edges.dx[2] * binScreenY0;
				int32_t const e2_x0y1 = edges.c[2] + edges.dy[2] * binScreenX0 + edges.dx[2] * binScreenY1;
				int32_t const e2_x1y0 = edges.c[2] + edges.dy[2] * binScreenX1 + edges.dx[2] * binScreenY0;
				int32_t const e2_x1y1 = edges.c[2] + edges.dy[2] * binScreenX1 + edges.dx[2] * binScreenY1;

				uint32_t const e2_allOut = (e2_x0y0 > 0) | ((e2_x0y1 > 0) << 1) | ((e2_x1y0 > 0) << 2) | ((e2_x1y1 > 0) << 3);

				if (!e0_allOut || !e1_allOut || !e2_allOut)
				{
					continue;
				}

				uint64_t mask8x8 = 0;

				//if (e0_allOut == 0xF && e1_allOut == 0xF && e2_allOut == 0xF)
				//{
				//	//// Todo: mark block instead of outputting 64 frags
				//	mask8x8 = ~0ull;
				//	//KT_ASSERT(o_buffer.m_numFragments < FragmentBuffer::c_maxFragments);
				//	//FragmentBuffer::Frag& frag = o_buffer.m_fragments[o_buffer.m_numFragments++];

				//	//frag.x = uint8_t(binScreenX0);
				//	//frag.y = uint8_t(binScreenY0);
				//	//frag.flags = FragmentBuffer::Flags::Full8x8Block;
				//	//frag.triIdx = triIdx;
				//	//frag.chunkIdx = _chunkIdx;
				//	//continue;

				//}
				//else
				//{
					mask8x8 = ComputeBlockMask8x8(blockEdgesSimd, zOverWPlaneSimd, _call, _depth, xBlock, yBlock);
				//}
				
				// KT_ASSERT(mask8x8 != ~0ull && "Rasterization routine found full 8x8 block, but above code did not detect it!");

				while (mask8x8)
				{
					// Todo: Maybe could do some fancy simd left packing?
					uint64_t const bitIdx = kt::Cnttz(mask8x8);
					uint8_t const bitY = uint8_t(bitIdx / 8) + binScreenY0;
					uint8_t const bitX = uint8_t(bitIdx & 7) + binScreenX0;

					KT_ASSERT(o_buffer.m_numFragments < FragmentBuffer::c_maxFragments);
					FragmentBuffer::Frag& frag = o_buffer.m_fragments[o_buffer.m_numFragments++];
					frag.x = bitX;
					frag.y = bitY;
					frag.flags = FragmentBuffer::Flags::None;
					frag.triIdx = triIdx;
					frag.chunkIdx = _chunkIdx;

					mask8x8 ^= (1ull << bitIdx);
				}
			}
		}
	}
}

static void RasterTrisInBin(DrawCall const& _call, ColourTile* _colour, DepthTile* _depth, BinChunk const& _chunk)
{
	BlockInterpolants8x8 blockinterpolant;
	EdgeEquations8x8 edges8x8;

	for (uint32_t triIdx = 0; triIdx < _chunk.m_numTris; ++triIdx)
	{
		BinChunk::EdgeEq const& edges = _chunk.m_edgeEq[triIdx];

		uint32_t const xBlockBegin = edges.blockMinX & ~7;
		uint32_t const yBlockBegin = edges.blockMinY & ~7;

		uint32_t const xBlockEnd = edges.blockMaxX;
		uint32_t const yBlockEnd = edges.blockMaxY;

		for (uint32_t yBlock = yBlockBegin; yBlock < yBlockEnd; yBlock += 8)
		{
			for (uint32_t xBlock = xBlockBegin; xBlock < xBlockEnd; xBlock += 8)
			{
				int32_t const binScreenX0 = xBlock;
				int32_t const binScreenX1 = xBlock + Config::c_binWidth;

				int32_t const binScreenY0 = yBlock;
				int32_t const binScreenY1 = yBlock + Config::c_binHeight;

				// Todo: slow, can offset edges and just test upper corner
				int32_t const e0_x0y0 = edges.c[0] + edges.dy[0] * binScreenX0 + edges.dx[0] * binScreenY0;
				int32_t const e0_x0y1 = edges.c[0] + edges.dy[0] * binScreenX0 + edges.dx[0] * binScreenY1;
				int32_t const e0_x1y0 = edges.c[0] + edges.dy[0] * binScreenX1 + edges.dx[0] * binScreenY0;
				int32_t const e0_x1y1 = edges.c[0] + edges.dy[0] * binScreenX1 + edges.dx[0] * binScreenY1;

				uint32_t const e0_allOut = (e0_x0y0 > 0) | ((e0_x0y1 > 0) << 1) | ((e0_x1y0 > 0) << 2) | ((e0_x1y1 > 0) << 3);

				int32_t const e1_x0y0 = edges.c[1] + edges.dy[1] * binScreenX0 + edges.dx[1] * binScreenY0;
				int32_t const e1_x0y1 = edges.c[1] + edges.dy[1] * binScreenX0 + edges.dx[1] * binScreenY1;
				int32_t const e1_x1y0 = edges.c[1] + edges.dy[1] * binScreenX1 + edges.dx[1] * binScreenY0;
				int32_t const e1_x1y1 = edges.c[1] + edges.dy[1] * binScreenX1 + edges.dx[1] * binScreenY1;

				uint32_t const e1_allOut = (e1_x0y0 > 0) | ((e1_x0y1 > 0) << 1) | ((e1_x1y0 > 0) << 2) | ((e1_x1y1 > 0) << 3);

				int32_t const e2_x0y0 = edges.c[2] + edges.dy[2] * binScreenX0 + edges.dx[2] * binScreenY0;
				int32_t const e2_x0y1 = edges.c[2] + edges.dy[2] * binScreenX0 + edges.dx[2] * binScreenY1;
				int32_t const e2_x1y0 = edges.c[2] + edges.dy[2] * binScreenX1 + edges.dx[2] * binScreenY0;
				int32_t const e2_x1y1 = edges.c[2] + edges.dy[2] * binScreenX1 + edges.dx[2] * binScreenY1;

				uint32_t const e2_allOut = (e2_x0y0 > 0) | ((e2_x0y1 > 0) << 1) | ((e2_x1y0 > 0) << 2) | ((e2_x1y1 > 0) << 3);

				if (!e0_allOut || !e1_allOut || !e2_allOut)
				{
					continue;
				}

				//if(e0_allOut == 0xF && e1_allOut == 0xF && e2_allOut == 0xF) // todo, full block

				GatherPartialBlockState(blockinterpolant, edges8x8, _chunk, xBlock , yBlock, triIdx);
				ShadePartialBlockAVX_8x8(_call, _colour, _depth, xBlock , yBlock, blockinterpolant, edges8x8);
			}
		}
	}
}

static void ShadeFragmentBuffer(ThreadRasterCtx const& _ctx, uint32_t const _tileIdx, BinChunk const* const* _chunks, FragmentBuffer& _buffer)
{
	// Fill interpolants.
	if (!_buffer.m_numFragments)
	{
		return;
	}

	uint32_t* fragsPerCall = (uint32_t*)_ctx.m_allocator->Alloc(sizeof(uint32_t) * _ctx.m_numDrawCalls);
	memset(fragsPerCall, 0, sizeof(uint32_t) * _ctx.m_numDrawCalls);

	{
		uint32_t fragIdx = 0;

		FragmentBuffer::Frag const* frag = &_buffer.m_fragments[fragIdx++];

		float* outAttribs = _buffer.m_interpolants;

		for (;;)
		{
			uint32_t packedTriChunkIdx = frag->packedChunkTriIdx;

			BinChunk const& chunk = *_chunks[frag->chunkIdx];
			uint32_t curDrawCallIdx = chunk.m_drawCallIdx;

			BinChunk::PlaneEq const& recipW = chunk.m_recipW[frag->triIdx];
			uint32_t const numAttribs = chunk.m_attribsPerTri;
			uint32_t const numAttribsAvx = (numAttribs + 7) / 8;

			__m256 attribPlaneDx[2];
			__m256 attribPlaneDy[2];
			__m256 attribPlaneC[2];

			for (uint32_t i = 0; i < numAttribsAvx; ++i)
			{
				attribPlaneDx[i] = _mm256_load_ps(&chunk.m_attribsDx[i * 8 + numAttribs * frag->triIdx]);
				attribPlaneDy[i] = _mm256_load_ps(&chunk.m_attribsDy[i * 8 + numAttribs * frag->triIdx]);
				attribPlaneC[i] = _mm256_load_ps(&chunk.m_attribsC[i * 8 + numAttribs * frag->triIdx]);
			}

			static_assert(Config::c_maxVaryings <= 16, "SIMD code assumed maxvaryings <= 16.");

			do
			{
				__m256 const fragX = _mm256_set1_ps(float(frag->x));
				__m256 const fragY = _mm256_set1_ps(float(frag->y));
				__m256 const recipWavx = _mm256_div_ps(_mm256_set1_ps(1.0f), _mm256_set1_ps(recipW.c0 + float(frag->x) * recipW.dx + float(frag->y) * recipW.dy));
				for (uint32_t i = 0; i < numAttribsAvx; ++i)
				{
					// Todo: derivatives
					__m256 const attribs = _mm256_fmadd_ps(attribPlaneDy[i], fragY, _mm256_fmadd_ps(attribPlaneDx[i], fragX, attribPlaneC[i]));
					_mm256_storeu_ps(outAttribs + i * 8, _mm256_mul_ps(recipWavx, attribs));
				}

				outAttribs += numAttribs;

				++fragsPerCall[curDrawCallIdx];
				frag = &_buffer.m_fragments[fragIdx++];

			} while (fragIdx < _buffer.m_numFragments
					 && frag->packedChunkTriIdx == packedTriChunkIdx);

			if (fragIdx >= _buffer.m_numFragments)
			{
				break;
			}
		}
	}


	// Now shade
	float const* interpolants = _buffer.m_interpolants;
	uint32_t globalFragIdx = 0;

	for (uint32_t drawCallIdx = 0; drawCallIdx < _ctx.m_numDrawCalls; ++drawCallIdx)
	{
		uint32_t const numFragsForCall = fragsPerCall[drawCallIdx];

		DrawCall const& call = _ctx.m_drawCalls[drawCallIdx];
		uint8_t* pixelWrite = call.m_frameBuffer->m_colourTiles[_tileIdx].m_colour;

		for (uint32_t drawCallFrag = 0; drawCallFrag < numFragsForCall; drawCallFrag += 8)
		{
			KT_ALIGNAS(32) float colourRGBA[4 * 8];
			call.m_pixelShader(call.m_pixelUniforms, interpolants, colourRGBA, _mm256_setzero_ps());

			uint32_t const writePixels = kt::Min(8u, numFragsForCall - drawCallFrag);

			interpolants += writePixels * (call.m_attributeBuffer.m_stride / sizeof(float));

			for (uint32_t i = 0; i < writePixels; ++i)
			{
				FragmentBuffer::Frag const& frag = _buffer.m_fragments[globalFragIdx++];
				uint32_t const pixIdx = frag.y * Config::c_binWidth + frag.x;
				pixelWrite[pixIdx * 4 + 0] = uint8_t(kt::Min(1.0f, colourRGBA[i * 4]) * 255.0f);
				pixelWrite[pixIdx * 4 + 1] = uint8_t(kt::Min(1.0f, colourRGBA[i * 4 + 1]) * 255.0f);
				pixelWrite[pixIdx * 4 + 2] = uint8_t(kt::Min(1.0f, colourRGBA[i * 4 + 2]) * 255.0f);
				pixelWrite[pixIdx * 4 + 3] = 255;

			}
		}
	}
}

void RasterAndShadeBin2(ThreadRasterCtx const& _ctx)
{
	// sort draw calls
	uint32_t numChunks = 0;

	// gather draw calls from all threads
	for (uint32_t threadBinIdx = 0; threadBinIdx < _ctx.m_binner->m_numThreads; ++threadBinIdx)
	{
		ThreadBin& bin = _ctx.m_binner->LookupThreadBin(threadBinIdx, _ctx.m_tileX, _ctx.m_tileY);

		numChunks += bin.m_numChunks;
	}

	if (!numChunks)
	{
		return;
	}

	BinChunk** sortedChunks = (BinChunk**)KT_ALLOCA(sizeof(BinChunk*) * numChunks);

	{
		uint32_t chunkIdx = 0;
		for (uint32_t threadBinIdx = 0; threadBinIdx < _ctx.m_binner->m_numThreads; ++threadBinIdx)
		{
			ThreadBin& bin = _ctx.m_binner->LookupThreadBin(threadBinIdx, _ctx.m_tileX, _ctx.m_tileY);
			memcpy(sortedChunks + chunkIdx, bin.m_binChunks, bin.m_numChunks * sizeof(BinChunk*));
			chunkIdx += bin.m_numChunks;
			KT_ASSERT(chunkIdx <= numChunks);
		}
	}

	BinChunk** radixTemp = (BinChunk**)KT_ALLOCA(sizeof(BinChunk**) * numChunks);

	kt::RadixSort(sortedChunks, sortedChunks + numChunks, radixTemp, [](BinChunk const* _c) { return _c->m_drawCallIdx; });

	uint32_t const tileIdx = _ctx.m_tileY * _ctx.m_binner->m_numBinsX + _ctx.m_tileX;

	FragmentBuffer* buffer = kt::New<FragmentBuffer>(_ctx.m_allocator);

	for (uint32_t chunkIdx = 0; chunkIdx < numChunks; ++chunkIdx)
	{
		buffer->Reset();
		BinChunk& curChunk = *sortedChunks[chunkIdx];
		DrawCall const& call = _ctx.m_drawCalls[curChunk.m_drawCallIdx];
		RasterizeTrisInBin_OutputFragments(call, &call.m_frameBuffer->m_depthTiles[tileIdx], curChunk, chunkIdx, *buffer);
		ShadeFragmentBuffer(_ctx, tileIdx, sortedChunks, *buffer); // Todo: shouldn't reset every chunk
	}

}

void RasterAndShadeBin(ThreadRasterCtx const& _ctx)
{
	// TOdo: sort draw calls
#define TEST_NEW 1

#if TEST_NEW
	RasterAndShadeBin2(_ctx);
					   
#else
	for (uint32_t threadBinIdx = 0; threadBinIdx < _ctx.m_binner->m_numThreads; ++threadBinIdx)
	{
		ThreadBin& bin = _ctx.m_binner->LookupThreadBin(threadBinIdx, _ctx.m_tileX, _ctx.m_tileY);

		for (uint32_t j = 0; j < bin.m_numChunks; ++j)
		{
			uint32_t tileIdx = _ctx.m_tileY * _ctx.m_binner->m_numBinsX + _ctx.m_tileX;
			DrawCall& call = _ctx.m_drawCalls[bin.m_binChunks[j]->m_drawCallIdx];

			RasterTrisInBin(call, &call.m_frameBuffer->m_colourTiles[tileIdx], &call.m_frameBuffer->m_depthTiles[tileIdx], *bin.m_binChunks[j]);
		}

	}
#endif

}

}
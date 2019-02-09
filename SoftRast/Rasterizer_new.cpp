#include "Rasterizer_new.h"
#include "Binning.h"
#include "Renderer.h"

namespace sr
{

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

	for (uint32_t i = 0; i < Config::c_maxVaryings; ++i)
	{
		o_block.attrib_eval[i] = _mm256_broadcast_ss(&_chunk.m_attribs[_triIdx].c0[i]);
		__m256 const attribdx = _mm256_broadcast_ss(&_chunk.m_attribs[_triIdx].dx[i]);
		o_block.attrib_dy[i] = _mm256_broadcast_ss(&_chunk.m_attribs[_triIdx].dy[i]);
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

				_call.m_pixelShader(_call.m_pixelUniforms, interpolants, colourRGBA, depthCmpMask);

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


void RasterAndShadeBin(ThreadRasterCtx const& _ctx)
{
	// TOdo: sort draw calls

	for (uint32_t threadBinIdx = 0; threadBinIdx < _ctx.m_binner->m_numThreads; ++threadBinIdx)
	{
		ThreadBin& bin = _ctx.m_binner->LookupThreadBin(threadBinIdx, _ctx.m_tileX, _ctx.m_tileY);

		for (uint32_t j = 0; j < bin.m_numChunks; ++j)
		{
			uint32_t tileIdx = _ctx.m_tileY * _ctx.m_binner->m_numBinsX + _ctx.m_tileX;
			DrawCall& call = _ctx.m_drawCalls[bin.m_drawCallIndicies[j]];

			RasterTrisInBin(call, &call.m_frameBuffer->m_colourTiles[tileIdx], &call.m_frameBuffer->m_depthTiles[tileIdx], *bin.m_binChunks[j]);
		}

	}


}

}
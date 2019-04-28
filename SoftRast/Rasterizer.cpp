#include "Rasterizer.h"
#include "Binning.h"
#include "Renderer.h"
#include "kt/Sort.h"

#include "microprofile.h"

MICROPROFILE_DEFINE(RasterAndShadeTile, "Backend", "RasterAndShadeTile", MP_GREEN);
MICROPROFILE_DEFINE(RasterFragments, "Backend", "RasterFragments", MP_GREEN1);
MICROPROFILE_DEFINE(ComputeInterpolants, "Backend", "ComputeInterpolants", MP_GREEN2);
MICROPROFILE_DEFINE(ShadeFragments, "Backend", "ShadeFragments", MP_GREEN3);

namespace sr
{


struct FragmentBuffer
{
	struct Frag
	{
		uint8_t x;
		uint8_t y;
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

	ThreadScratchAllocator* m_allocator = nullptr;

	void ReserveFragments(uint32_t _count)
	{
		KT_ASSERT(m_fragments);
		KT_ASSERT(!m_interpolants);
		void* p = m_allocator->Alloc(sizeof(Frag) * _count, 1);
		KT_ASSERT(p);
		KT_UNUSED(p);
	}

	Frag* m_fragments = nullptr;
	float* m_interpolants = nullptr;

	uint32_t m_numFragments = 0;
	uint32_t m_interpolantsAllocSize = 0;
};

// Todo: naming is really confusing as edge eq dx*y and dy*x are used but the opposite for planes. Should maybe clarify naming.
struct EdgeEquations8x8
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

KT_FORCEINLINE __m256 DepthCmpMask(__m256 _old, __m256 _new)
{
#if SR_USE_REVERSE_Z
	return _mm256_and_ps(_mm256_cmp_ps(_new, _mm256_setzero_ps(), _CMP_GT_OQ), _mm256_cmp_ps(_new, _old, _CMP_GT_OQ));
#else
	return _mm256_and_ps(_mm256_cmp_ps(_new, _mm256_setzero_ps(), _CMP_GT_OQ), _mm256_cmp_ps(_new, _old, _CMP_LT_OQ));
#endif
}

static uint64_t ComputeBlockMask8x8_DepthOnly
(
	ZOverW8x8 const& _zOverW,
	DepthTile* _depth,
	int32_t _xTileRelative,
	int32_t _yTileRelative
)
{
	__m256i const xTileSimd = _mm256_set1_epi32(_xTileRelative);
	__m256i const yTileSimd = _mm256_set1_epi32(_yTileRelative);

	__m256 zOverW = _mm256_fmadd_ps(_mm256_cvtepi32_ps(yTileSimd), _zOverW.dy, _zOverW.tileTopLeft);
	zOverW = _mm256_add_ps(zOverW, _mm256_mul_ps(_mm256_cvtepi32_ps(xTileSimd), _zOverW.dx));

	uint64_t mask8x8 = 0;

	for (uint32_t i = 0; i < 8; ++i)
	{
		float* depthPtr = _depth->m_depth + (_xTileRelative + (_yTileRelative + i) * Config::c_binWidth);
		__m256 const depthGather = _mm256_loadu_ps(depthPtr);

		__m256 depthCmpMask = DepthCmpMask(depthGather, zOverW);

		uint64_t const laneMask = _mm256_movemask_ps(depthCmpMask);

		__m256 const newDepth = _mm256_blendv_ps(depthGather, zOverW, depthCmpMask);
		_mm256_storeu_ps(depthPtr, newDepth);

		mask8x8 |= (laneMask << (i * 8ull));

		zOverW = _mm256_add_ps(zOverW, _zOverW.dy);
	}
	return mask8x8;
}

static uint64_t ComputeBlockMask8x8
(
	EdgeEquations8x8 const& _edges,
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

	__m256i const signBit = _mm256_set1_epi32(0x80000000);

	for (uint32_t i = 0; i < 8; ++i)
	{
		__m256i const edgeMask = _mm256_xor_si256(_mm256_or_si256(_mm256_or_si256(edges[0], edges[1]), edges[2]), signBit);

		// test Z
		float* depthPtr = _depth->m_depth + (_xTileRelative + (_yTileRelative + i) * Config::c_binWidth);
		__m256 const depthGather = _mm256_loadu_ps(depthPtr);

		__m256 depthCmpMask = DepthCmpMask(depthGather, zOverW);
		depthCmpMask = _mm256_and_ps(_mm256_castsi256_ps(edgeMask), depthCmpMask);

		uint64_t const laneMask = _mm256_movemask_ps(depthCmpMask);

		__m256 const newDepth = _mm256_blendv_ps(depthGather, zOverW, depthCmpMask);
		_mm256_storeu_ps(depthPtr, newDepth);

		mask8x8 |= (laneMask << (i * 8ull));

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

		EdgeEquations8x8 blockEdgesSimd;

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

				if (e0_allOut == 0xF && e1_allOut == 0xF && e2_allOut == 0xF)
				{
					mask8x8 = ComputeBlockMask8x8_DepthOnly(zOverWPlaneSimd, _depth, xBlock, yBlock);
				}
				else
				{
					mask8x8 = ComputeBlockMask8x8(blockEdgesSimd, zOverWPlaneSimd, _call, _depth, xBlock, yBlock);
				}

				if (!mask8x8)
				{
					continue;
				}

				uint32_t const numFragsToOutput = uint32_t(kt::Popcnt(mask8x8));
				o_buffer.ReserveFragments(numFragsToOutput);

				o_buffer.m_interpolantsAllocSize += numFragsToOutput * _call.m_attributeBuffer.m_stride + 4; // 4 for derivs, make constant?

				do
				{
					// Todo: Maybe could do some fancy simd left packing?
					// Todo: should maybe compress fragment stream?
					uint64_t const bitIdx = kt::Cnttz(mask8x8);
					uint8_t const bitY = uint8_t(bitIdx / 8) + binScreenY0;
					uint8_t const bitX = uint8_t(bitIdx & 7) + binScreenX0;

					KT_ASSERT(bitY < Config::c_binHeight);
					KT_ASSERT(bitX < Config::c_binWidth);

					FragmentBuffer::Frag& frag = o_buffer.m_fragments[o_buffer.m_numFragments++];
					frag.x = bitX;
					frag.y = bitY;
					frag.triIdx = triIdx;
					frag.chunkIdx = _chunkIdx;

					mask8x8 ^= (1ull << bitIdx);
				} while (mask8x8);
			}
		}
	}
}

template <uint32_t NumAttribsAVXNoDerivT>
static bool ComputeInterpolantsDrawCallImpl
(
	ThreadRasterCtx const& _ctx, 
	FragmentBuffer& _buffer, 
	float*& io_attribs,
	BinChunk const& _chunk,
	uint32_t* o_fragsPerDrawCall,
	uint32_t const _numAttribsNoDeriv,
	FragmentBuffer::Frag const*& io_frag,
	uint32_t& io_nextFragIdx
)
{
	uint32_t const drawCallIdx = _chunk.m_drawCallIdx;
	uint32_t const chunkIdx = io_frag->chunkIdx;
	uint32_t const packedTriChunkIdx = io_frag->packedChunkTriIdx;
	do
	{
		__m256 attribPlaneDx[NumAttribsAVXNoDerivT];
		__m256 attribPlaneDy[NumAttribsAVXNoDerivT];
		__m256 attribPlaneC[NumAttribsAVXNoDerivT];

		for (uint32_t i = 0; i < NumAttribsAVXNoDerivT; ++i)
		{
			attribPlaneDx[i] = _mm256_load_ps(&_chunk.m_attribsDx[i * 8 + _numAttribsNoDeriv * io_frag->triIdx]);
			attribPlaneDy[i] = _mm256_load_ps(&_chunk.m_attribsDy[i * 8 + _numAttribsNoDeriv * io_frag->triIdx]);
			attribPlaneC[i] = _mm256_load_ps(&_chunk.m_attribsC[i * 8 + _numAttribsNoDeriv * io_frag->triIdx]);
		}

		BinChunk::PlaneEq const& recipW = _chunk.m_recipW[io_frag->triIdx];

		__m256 recipW_dx = _mm256_broadcast_ss(&recipW.dx);
		__m256 recipW_dy = _mm256_broadcast_ss(&recipW.dy);
		__m256 recipW_c = _mm256_broadcast_ss(&recipW.c0);

		__m256 u4_v4_dx;
		__m256 u4_v4_dy;
		__m256 u4_v4_c;

		uint32_t const uOffset = _ctx.m_drawCalls[drawCallIdx].m_uvOffset;

		{
			__m256 const udx = _mm256_broadcast_ss(&_chunk.m_attribsDx[io_frag->triIdx * _numAttribsNoDeriv + uOffset]);
			__m256 const vdx = _mm256_broadcast_ss(&_chunk.m_attribsDx[io_frag->triIdx * _numAttribsNoDeriv + uOffset + 1]);

			u4_v4_dx = _mm256_permute2f128_ps(udx, vdx, (0 | (2 << 4)));
		}
		{
			__m256 const udy = _mm256_broadcast_ss(&_chunk.m_attribsDy[io_frag->triIdx * _numAttribsNoDeriv + uOffset]);
			__m256 const vdy = _mm256_broadcast_ss(&_chunk.m_attribsDy[io_frag->triIdx * _numAttribsNoDeriv + uOffset + 1]);
			u4_v4_dy = _mm256_permute2f128_ps(udy, vdy, (0 | (2 << 4)));
		}
		{
			__m256 const uc = _mm256_broadcast_ss(&_chunk.m_attribsC[io_frag->triIdx * _numAttribsNoDeriv + uOffset]);
			__m256 const vc = _mm256_broadcast_ss(&_chunk.m_attribsC[io_frag->triIdx * _numAttribsNoDeriv + uOffset + 1]);
			u4_v4_c = _mm256_permute2f128_ps(uc, vc, (0 | (2 << 4)));
		}

		static_assert(Config::c_maxVaryings <= 16, "SIMD code assumed maxvaryings <= 16.");

		do
		{
			static __m256 const interleaved_0_1_0_1 = _mm256_setr_ps(0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f);
			__m256 const interleaved_0_0_1_1 = _mm256_shuffle_ps(interleaved_0_1_0_1, interleaved_0_1_0_1, _MM_SHUFFLE(1, 1, 0, 0));

			__m256 const fragX = _mm256_set1_ps(float(io_frag->x));
			__m256 const fragY = _mm256_set1_ps(float(io_frag->y));

			__m256 const frag_x0_x1_x0_x1 = _mm256_add_ps(interleaved_0_1_0_1, fragX);
			__m256 const frag_y0_y0_y1_y1 = _mm256_add_ps(interleaved_0_0_1_1, fragY);

			// duvdx = uv_x1y0 - uv_x0y0 
			// duvdy = uv_x0y1 - uv_x0y0

			__m256 const uvOverW_x0y0_x1y0_x0y1_x1y1 = _mm256_fmadd_ps(u4_v4_dy, frag_y0_y0_y1_y1, _mm256_fmadd_ps(u4_v4_dx, frag_x0_x1_x0_x1, u4_v4_c));

			// evaluate recipW for quad (though we only need 3 vals denoted in above comment)
			__m256 const recipW_x0y0_x1y0_x0y1_x1y1 = _mm256_div_ps(_mm256_set1_ps(1.0f), _mm256_fmadd_ps(recipW_dx, frag_x0_x1_x0_x1, _mm256_fmadd_ps(recipW_dy, frag_y0_y0_y1_y1, recipW_c)));

			__m256 const uv_eval_x0y0_x1y0_x0y1_x1y1 = _mm256_mul_ps(recipW_x0y0_x1y0_x0y1_x1y1, uvOverW_x0y0_x1y0_x0y1_x1y1);

			// swizzle and calculate derivatives.
			__m256 const uv_eval_x0y0 = _mm256_shuffle_ps(uv_eval_x0y0_x1y0_x0y1_x1y1, uv_eval_x0y0_x1y0_x0y1_x1y1, _MM_SHUFFLE(0, 0, 0, 0));
			__m256 const derivs = _mm256_sub_ps(uv_eval_x0y0_x1y0_x0y1_x1y1, uv_eval_x0y0);
			// 4 du then 4 dv -> [0, dx, dy]
			static const __m256i deriv_permute_mask = _mm256_setr_epi32(1, 2, 1 + 4, 2 + 4, 0, 0, 0, 0);
			__m128 const finalDerivs = _mm256_castps256_ps128(_mm256_permutevar8x32_ps(derivs, deriv_permute_mask));

			// Storing: dudx, dudy, dvdx, dvdy
			_mm_storeu_ps(io_attribs, finalDerivs);
			io_attribs += 4;

			// Broadcast original pixel 1/w for rest of attributes.
			__m256 const recipW_x0y0 = _mm256_shuffle_ps(recipW_x0y0_x1y0_x0y1_x1y1, recipW_x0y0_x1y0_x0y1_x1y1, _MM_SHUFFLE(0, 0, 0, 0));

			for (uint32_t i = 0; i < NumAttribsAVXNoDerivT; ++i)
			{
				__m256 const attribs = _mm256_fmadd_ps(attribPlaneDy[i], fragY, _mm256_fmadd_ps(attribPlaneDx[i], fragX, attribPlaneC[i]));
				_mm256_storeu_ps(io_attribs + i * 8, _mm256_mul_ps(recipW_x0y0, attribs));
			}

			io_attribs += _numAttribsNoDeriv;

			++o_fragsPerDrawCall[drawCallIdx];

			if (io_nextFragIdx == _buffer.m_numFragments)
			{
				return true;
			}

			io_frag = &_buffer.m_fragments[io_nextFragIdx++];

		} while (io_frag->packedChunkTriIdx == packedTriChunkIdx);
	} while (io_frag->chunkIdx == chunkIdx);

	return false;
}


static void ComputeInterpolants(ThreadRasterCtx const& _ctx, BinChunk const* const* _chunks, FragmentBuffer& _buffer, uint32_t* o_fragsPerDrawCall)
{
	MICROPROFILE_SCOPE(ComputeInterpolants);
	uint32_t nextFragIdx = 0;
	FragmentBuffer::Frag const* frag = &_buffer.m_fragments[nextFragIdx++];

	float* outAttribs = _buffer.m_interpolants;

	for (;;)
	{
		BinChunk const& chunk = *_chunks[frag->chunkIdx];
		KT_ASSERT(chunk.m_drawCallIdx < _ctx.m_numDrawCalls);

		uint32_t const numAttribsNoDeriv = chunk.m_attribsPerTri;
		uint32_t const numAttribsAvxNoDeriv = (numAttribsNoDeriv + 7) / 8;

		if (numAttribsAvxNoDeriv == 1)
		{
			if (ComputeInterpolantsDrawCallImpl<1>(_ctx, _buffer, outAttribs, chunk, o_fragsPerDrawCall, numAttribsNoDeriv, frag, nextFragIdx))
			{
				break;
			}
		}
		else if (numAttribsAvxNoDeriv == 2)
		{
			if (ComputeInterpolantsDrawCallImpl<2>(_ctx, _buffer, outAttribs, chunk, o_fragsPerDrawCall, numAttribsNoDeriv, frag, nextFragIdx))
			{
				break;
			}
		}
		else
		{
			KT_UNREACHABLE;
		}
	}
	KT_ASSERT(nextFragIdx == _buffer.m_numFragments);
}

static void ShadeFragmentBuffer(ThreadRasterCtx const& _ctx, uint32_t const _tileIdx, BinChunk const* const* _chunks, FragmentBuffer& _buffer)
{
	// Fill interpolants.
	if (!_buffer.m_numFragments)
	{
		return;
	}

	uint32_t* fragsPerCall = (uint32_t*)KT_ALLOCA(sizeof(uint32_t) * _ctx.m_numDrawCalls);
	memset(fragsPerCall, 0, sizeof(uint32_t) * _ctx.m_numDrawCalls);

	ComputeInterpolants(_ctx, _chunks, _buffer, fragsPerCall);

#if KT_DEBUG
	{
		uint32_t totalFrags = 0;
		for (uint32_t i = 0; i < _ctx.m_numDrawCalls; ++i)
		{
			totalFrags += fragsPerCall[i];
		}
		KT_ASSERT(totalFrags == _buffer.m_numFragments);
	}
#endif
	MICROPROFILE_SCOPE(ShadeFragments);

	float const* interpolants = _buffer.m_interpolants;
	uint32_t globalFragIdx = 0;

	for (uint32_t drawCallIdx = 0; drawCallIdx < _ctx.m_numDrawCalls; ++drawCallIdx)
	{
		uint32_t const numFragsForCall = fragsPerCall[drawCallIdx];

		DrawCall const& call = _ctx.m_drawCalls[drawCallIdx];
		uint32_t* pixelWrite = (uint32_t*)call.m_frameBuffer->m_colourTiles[_tileIdx].m_colour;

		for (uint32_t drawCallFrag = 0; drawCallFrag < numFragsForCall; drawCallFrag += 8)
		{
			KT_ALIGNAS(32) uint32_t colourRGBA[8];

			uint32_t const writePixels = kt::Min(8u, numFragsForCall - drawCallFrag);
			uint32_t const laneMask = (1 << writePixels) - 1;

			call.m_pixelShader(call.m_pixelUniforms, interpolants, colourRGBA, laneMask);

			// Todo: 4 constant for derivatives
			interpolants += writePixels * (4 + call.m_attributeBuffer.m_stride / sizeof(float));

			for (uint32_t i = 0; i < writePixels; ++i)
			{
				FragmentBuffer::Frag const& frag = _buffer.m_fragments[globalFragIdx++];
				uint32_t const pixIdx = frag.y * Config::c_binWidth + frag.x;
				pixelWrite[pixIdx] = colourRGBA[i];
			}
		}
	}
	KT_ASSERT(globalFragIdx == _buffer.m_numFragments);
}

void RasterAndShadeBin(ThreadRasterCtx const& _ctx)
{
	MICROPROFILE_SCOPE(RasterAndShadeTile);

	ThreadScratchAllocator& threadAllocator = _ctx.m_ctx->ThreadAllocator();
	ThreadScratchAllocator::AllocScope const allocScope(threadAllocator);

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

	FragmentBuffer buffer;
	buffer.m_fragments = (FragmentBuffer::Frag*)threadAllocator.Align(KT_ALIGNOF(FragmentBuffer::Frag));
	buffer.m_allocator = &threadAllocator;
	KT_ASSERT(buffer.m_fragments);

	{
		//MICROPROFILE_SCOPE(RasterFragments);
		MICROPROFILE_SCOPEI("Backend", "RasterFragments", MP_BLUE1);
		for (uint32_t chunkIdx = 0; chunkIdx < numChunks; ++chunkIdx)
		{
			BinChunk& curChunk = *sortedChunks[chunkIdx];
			DrawCall const& call = _ctx.m_drawCalls[curChunk.m_drawCallIdx];
			RasterizeTrisInBin_OutputFragments(call, &call.m_frameBuffer->m_depthTiles[tileIdx], curChunk, chunkIdx, buffer);
		}
	}

	{
		buffer.m_interpolants = (float*)threadAllocator.Alloc(buffer.m_interpolantsAllocSize, 32);
		KT_ASSERT(buffer.m_interpolants);
		ShadeFragmentBuffer(_ctx, tileIdx, sortedChunks, buffer);
	}
}

}
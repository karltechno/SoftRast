#include "Rasterizer.h"
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

	static_assert((1 << 7) >= Config::c_binHeight, "Cant fit height in 7 bits");
	static_assert((1 << 7) >= Config::c_binWidth, "Cant fit width in 7 bits");

	// Todo: chunk the buffer
	static uint32_t constexpr c_maxFragments = 2048 * 16;

	KT_ALIGNAS(32) float m_interpolants[Config::c_maxVaryings * c_maxFragments];

	KT_ALIGNAS(32) Frag m_fragments[c_maxFragments];
	uint32_t m_numFragments = 0;

	void Reset()
	{
		m_numFragments = 0;
	}
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

		__m256 depthCmpMask = _mm256_and_ps(_mm256_cmp_ps(zOverW, _mm256_setzero_ps(), _CMP_GT_OQ), _mm256_cmp_ps(depthGather, zOverW, _CMP_GT_OQ));

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

	// Todo: look at other coverage mask generation approaches. Eg LUT based.

	// Manually unroll ?

	__m256i const signBit = _mm256_set1_epi32(0x80000000);

	for (uint32_t i = 0; i < 8; ++i)
	{
		__m256i const edgeMask = _mm256_xor_si256(_mm256_or_si256(_mm256_or_si256(edges[0], edges[1]), edges[2]), signBit);

		// test Z
		float* depthPtr = _depth->m_depth + (_xTileRelative + (_yTileRelative + i) * Config::c_binWidth);
		__m256 const depthGather = _mm256_loadu_ps(depthPtr);

		__m256 depthCmpMask = _mm256_and_ps(_mm256_cmp_ps(zOverW, _mm256_setzero_ps(), _CMP_GT_OQ), _mm256_cmp_ps(depthGather, zOverW, _CMP_GT_OQ));
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
				
				while (mask8x8)
				{
					// Todo: Maybe could do some fancy simd left packing?
					// Todo: should maybe compress fragment stream?
					uint64_t const bitIdx = kt::Cnttz(mask8x8);
					uint8_t const bitY = uint8_t(bitIdx / 8) + binScreenY0;
					uint8_t const bitX = uint8_t(bitIdx & 7) + binScreenX0;

					KT_ASSERT(bitY < Config::c_binHeight);
					KT_ASSERT(bitX < Config::c_binWidth);

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

static void ShadeFragmentBuffer(ThreadRasterCtx const& _ctx, uint32_t const _tileIdx, BinChunk const* const* _chunks, FragmentBuffer& _buffer)
{
	// Fill interpolants.
	if (!_buffer.m_numFragments)
	{
		return;
	}

	uint32_t* fragsPerCall = (uint32_t*)KT_ALLOCA(sizeof(uint32_t) * _ctx.m_numDrawCalls);
	memset(fragsPerCall, 0, sizeof(uint32_t) * _ctx.m_numDrawCalls);

	{
		uint32_t fragIdx = 0;
		FragmentBuffer::Frag const* frag = &_buffer.m_fragments[fragIdx++];

		float* outAttribs = _buffer.m_interpolants;

		for (;;)
		{
			uint32_t packedTriChunkIdx = frag->packedChunkTriIdx;

			BinChunk const& chunk = *_chunks[frag->chunkIdx];
			uint32_t const curDrawCallIdx = chunk.m_drawCallIdx;
			KT_ASSERT(curDrawCallIdx < _ctx.m_numDrawCalls);

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
				KT_ASSERT(outAttribs < _buffer.m_interpolants + KT_ARRAY_COUNT(_buffer.m_interpolants));

				++fragsPerCall[curDrawCallIdx];
				
				if (fragIdx == _buffer.m_numFragments)
				{
					goto do_shade;
				}

				frag = &_buffer.m_fragments[fragIdx++];

			} while (frag->packedChunkTriIdx == packedTriChunkIdx);
		}
	}

do_shade:
	   
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
	KT_ASSERT(globalFragIdx == _buffer.m_numFragments);
}

void RasterAndShadeBin(ThreadRasterCtx const& _ctx)
{
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

	FragmentBuffer* buffer = kt::New<FragmentBuffer>(&threadAllocator);

	for (uint32_t chunkIdx = 0; chunkIdx < numChunks; ++chunkIdx)
	{
		buffer->Reset();
		BinChunk& curChunk = *sortedChunks[chunkIdx];
		DrawCall const& call = _ctx.m_drawCalls[curChunk.m_drawCallIdx];
		RasterizeTrisInBin_OutputFragments(call, &call.m_frameBuffer->m_depthTiles[tileIdx], curChunk, chunkIdx, *buffer);
		ShadeFragmentBuffer(_ctx, tileIdx, sortedChunks, *buffer); // Todo: shouldn't reset every chunk
	}

}

}
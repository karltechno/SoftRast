#pragma once
#include <stdint.h>
#include "SoftRastTypes.h"
#include "kt/Macros.h"
#include "TaskSystem.h"

namespace sr
{

struct DrawCall;

static uint32_t const c_trisPerBinChunk = 512;
static uint32_t const c_maxThreadBinChunks = 256;

struct BinChunk
{
	struct EdgeEq
	{
		int32_t c[3];
		int32_t dx[3];
		int32_t dy[3];
	};

	struct PlaneEq
	{
		float c0;
		float dx;
		float dy;
	};

	struct AttribPlaneChunk
	{
		KT_ALIGNAS(32) float c0[Config::c_maxVaryings];
		KT_ALIGNAS(32) float dx[Config::c_maxVaryings];
		KT_ALIGNAS(32) float dy[Config::c_maxVaryings];
	};

	EdgeEq m_edgeEq[c_trisPerBinChunk];

	PlaneEq m_recipW[c_trisPerBinChunk];

	PlaneEq m_zOverW[c_trisPerBinChunk];

	AttribPlaneChunk m_attribs[c_trisPerBinChunk];

	uint32_t m_numTris = 0;
};

struct ThreadBin
{
	uint32_t m_drawCallIndicies[c_maxThreadBinChunks];
	BinChunk* m_binChunks[c_maxThreadBinChunks];

	uint32_t m_numChunks = 0;
};

struct BinContext
{
	~BinContext();

	void Init(uint32_t _numThreads, uint32_t _binsX, uint32_t _binsY);

	ThreadBin& LookupThreadBin(uint32_t _tileX, uint32_t _tileY, uint32_t _threadIdx);

	ThreadBin* m_bins = nullptr;

	uint32_t m_numBinsX = 0;
	uint32_t m_numBinsY = 0;
	uint32_t m_numThreads = 0;
};

void BinTrisEntry(BinContext& _ctx, ThreadScratchAllocator& _alloc, uint32_t _threadIdx, uint32_t _triIdxBegin, uint32_t _triIdxEnd, DrawCall const& _drawCall);

}
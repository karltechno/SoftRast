#pragma once
#include <stdint.h>
#include "TaskSystem.h"

namespace sr
{

struct BinChunk;
struct DepthTile;
struct ColourTile;
struct DrawCall;
struct BinContext;

struct ThreadRasterCtx
{
	ThreadScratchAllocator* m_allocator = nullptr;

	DrawCall* m_drawCalls = nullptr;

	BinContext* m_binner = nullptr;

	uint32_t m_numDrawCalls = 0;
	uint32_t m_tileX = 0;
	uint32_t m_tileY = 0;
};

void RasterAndShadeBin(ThreadRasterCtx const& _ctx);

}
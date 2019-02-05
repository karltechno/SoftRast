#pragma once

namespace sr
{

struct BinChunk;
struct DepthTile;
struct ColourTile;
struct DrawCall;

void RasterTrisInBin(DrawCall const& _call, BinChunk const& _chunk, DepthTile* _depth, ColourTile* _colour);

}
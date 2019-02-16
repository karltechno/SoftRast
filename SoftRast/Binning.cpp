#include "Binning.h"
#include "kt/Macros.h"
#include "kt/Memory.h"
#include "Renderer.h"
#include "TaskSystem.h"

namespace sr
{

BinContext::~BinContext()
{
	kt::Free(m_bins);
}

void BinContext::Init(uint32_t _numThreads, uint32_t _binsX, uint32_t _binsY)
{
	KT_ASSERT(!m_bins);

	m_numBinsX = _binsX;
	m_numBinsY = _binsY;
	m_numThreads = _numThreads;

	m_bins = (ThreadBin*)kt::Malloc(sizeof(ThreadBin) * _numThreads * _binsX * _binsY);
}

ThreadBin& BinContext::LookupThreadBin(uint32_t _threadIdx, uint32_t _binX, uint32_t _binY)
{
	KT_ASSERT(_binX < m_numBinsX);
	KT_ASSERT(_binY < m_numBinsY);
	KT_ASSERT(_threadIdx < m_numThreads);
	uint32_t const linearIdx = m_numThreads * (_binY * m_numBinsX + _binX) + _threadIdx;
	return m_bins[linearIdx];
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

// Each frustum plane can turn a point into an edge (1 vert -> 2 verts). Therefore each clip plane can add 1 vertex. 9 total (including 3 from initial tri)
constexpr uint32_t CLIP_VERT_BUFFER_SIZE = 3 + 6;

constexpr uint32_t CLIP_OUT_MAX_TRIS = (CLIP_VERT_BUFFER_SIZE - 2);

struct ClipBuffer
{
	float attribs[2][CLIP_VERT_BUFFER_SIZE][Config::c_maxVaryings];
	kt::Vec4 verts[2][CLIP_VERT_BUFFER_SIZE];

	uint32_t numInputVerts = 0;
	uint32_t numOutputVerts = 0;

	// 0 or 1
	uint32_t inputIdx = 0;
};

static void ClipPlane(ClipBuffer& _buffer, uint32_t _clipPlaneIdx)
{
	static kt::Vec4 const c_clipPlanes[6] =
	{
		{ 1.0f, 0.0f, 0.0f, 1.0f }, // X_Neg
		{ -1.0f, 0.0f, 0.0f, 1.0f  }, // X_Pos

		{ 0.0f, 1.0f, 0.0f, 1.0f }, // Y_Neg
		{ 0.0f, -1.0f, 0.0f, 1.0f  }, // Y_Pos

		{ 0.0f, 0.0f, 1.0f, 1.0f }, // Z_Near
		{ 0.0f, 0.0f, -1.0f, 1.0f }, // Z_Far
	};

	KT_ASSERT(_buffer.numInputVerts);

	kt::Vec4 const& clipPlane = c_clipPlanes[_clipPlaneIdx];

	uint32_t v0_input_idx = _buffer.numInputVerts - 1;

	kt::Vec4(&input_vec)[CLIP_VERT_BUFFER_SIZE] = _buffer.verts[_buffer.inputIdx];
	float(&input_attrib)[CLIP_VERT_BUFFER_SIZE][Config::c_maxVaryings] = _buffer.attribs[_buffer.inputIdx];

	kt::Vec4(&output_vec)[CLIP_VERT_BUFFER_SIZE] = _buffer.verts[_buffer.inputIdx ^ 1];
	float(&output_attrib)[CLIP_VERT_BUFFER_SIZE][Config::c_maxVaryings] = _buffer.attribs[_buffer.inputIdx ^ 1];

	float v0_dot = kt::Dot(clipPlane, input_vec[v0_input_idx]);

	for (uint32_t nextClipVertIdx = 0; nextClipVertIdx < _buffer.numInputVerts; ++nextClipVertIdx)
	{
		uint32_t const v1_input_idx = nextClipVertIdx;
		float const v1_dot = kt::Dot(clipPlane, input_vec[v1_input_idx]);
		bool const v0_inside = v0_dot >= 0.0f;
		bool const v1_inside = v1_dot >= 0.0f;

		if (v0_inside)
		{
			KT_ASSERT(_buffer.numOutputVerts < CLIP_VERT_BUFFER_SIZE);
			uint32_t const outIdx = _buffer.numOutputVerts++;

			output_vec[outIdx] = input_vec[v0_input_idx];
			memcpy(output_attrib[outIdx], input_attrib[v0_input_idx], sizeof(float) * Config::c_maxVaryings);
		}

		if (v0_inside ^ v1_inside)
		{
			if (v1_inside)
			{
				KT_ASSERT(_buffer.numOutputVerts < CLIP_VERT_BUFFER_SIZE);
				float const tInterp = v1_dot / (v1_dot - v0_dot);
				uint32_t const outputVertIdx = _buffer.numOutputVerts++;
				output_vec[outputVertIdx] = kt::Lerp(input_vec[v1_input_idx], input_vec[v0_input_idx], tInterp);

				for (uint32_t varI = 0; varI < Config::c_maxVaryings; ++varI)
				{
					output_attrib[outputVertIdx][varI] = kt::Lerp(input_attrib[v1_input_idx][varI], input_attrib[v0_input_idx][varI], tInterp);
				}
			}
			else
			{
				KT_ASSERT(_buffer.numOutputVerts < CLIP_VERT_BUFFER_SIZE);
				float const tInterp = v0_dot / (v0_dot - v1_dot);

				uint32_t const outputVertIdx = _buffer.numOutputVerts++;
				output_vec[outputVertIdx] = kt::Lerp(input_vec[v0_input_idx], input_vec[v1_input_idx], tInterp);

				for (uint32_t varI = 0; varI < Config::c_maxVaryings; ++varI)
				{
					output_attrib[outputVertIdx][varI] = kt::Lerp(input_attrib[v0_input_idx][varI], input_attrib[v1_input_idx][varI], tInterp);
				}
			}
		}

		v0_dot = v1_dot;
		v0_input_idx = v1_input_idx;
	}

	_buffer.numInputVerts = _buffer.numOutputVerts;
	_buffer.numOutputVerts = 0;
	_buffer.inputIdx ^= 1;
}

KT_FORCEINLINE static void FetchIndicies(DrawCall const& _call, uint32_t _triIdx, uint32_t o_indicies[3])
{
	switch (_call.m_indexBuffer.m_stride)
	{
		case 1:
		{
			uint8_t* ptr = (uint8_t*)_call.m_indexBuffer.m_ptr;
			o_indicies[0] = ptr[_triIdx * 3];
			o_indicies[1] = ptr[_triIdx * 3 + 1];
			o_indicies[2] = ptr[_triIdx * 3 + 2];
		} break;

		case 2:
		{
			uint16_t* ptr = (uint16_t*)_call.m_indexBuffer.m_ptr;
			o_indicies[0] = ptr[_triIdx * 3];
			o_indicies[1] = ptr[_triIdx * 3 + 1];
			o_indicies[2] = ptr[_triIdx * 3 + 2];
		} break;

		case 4:
		{
			uint32_t* ptr = (uint32_t*)_call.m_indexBuffer.m_ptr;
			o_indicies[0] = ptr[_triIdx * 3];
			o_indicies[1] = ptr[_triIdx * 3 + 1];
			o_indicies[2] = ptr[_triIdx * 3 + 2];
		} break;

		default:
		{
			KT_ASSERT(false);
		} break;
	}

	KT_ASSERT(o_indicies[0] < _call.m_positionBuffer.m_num);
	KT_ASSERT(o_indicies[1] < _call.m_positionBuffer.m_num);
	KT_ASSERT(o_indicies[2] < _call.m_positionBuffer.m_num);
}

KT_FORCEINLINE static void FetchPositions(DrawCall const& _call, uint32_t const i_indicies[3], kt::Vec4 o_positions[3])
{
	uint8_t* buff = (uint8_t*)_call.m_positionBuffer.m_ptr;
	o_positions[0] = kt::Vec4(*(kt::Vec3*)(buff + i_indicies[0] * _call.m_positionBuffer.m_stride), 1.0f);
	o_positions[1] = kt::Vec4(*(kt::Vec3*)(buff + i_indicies[1] * _call.m_positionBuffer.m_stride), 1.0f);
	o_positions[2] = kt::Vec4(*(kt::Vec3*)(buff + i_indicies[2] * _call.m_positionBuffer.m_stride), 1.0f);
}

KT_FORCEINLINE static void FetchAttribPointers(DrawCall const& _call, uint32_t const i_indicies[3], float const* (&o_attribs)[3])
{
	uint8_t* buff = (uint8_t*)_call.m_attributeBuffer.m_ptr;

	o_attribs[0] = (float*)(buff + i_indicies[0] * _call.m_attributeBuffer.m_stride);
	o_attribs[1] = (float*)(buff + i_indicies[1] * _call.m_attributeBuffer.m_stride);
	o_attribs[2] = (float*)(buff + i_indicies[2] * _call.m_attributeBuffer.m_stride);
}

static BinChunk& GetOrCreateBinForDrawCall(ThreadScratchAllocator& _alloc, BinContext& _ctx, ThreadBin& _bin, DrawCall const& _call)
{
	if (_bin.m_numChunks
		&& _bin.m_drawCallIndicies[_bin.m_numChunks - 1] == _call.m_drawCallIdx
		&& _bin.m_binChunks[_bin.m_numChunks - 1]->m_numTris < c_trisPerBinChunk)
	{
		return *_bin.m_binChunks[_bin.m_numChunks - 1];
	}

	KT_ASSERT(_bin.m_numChunks < c_maxThreadBinChunks);
	BinChunk* newChunk = (BinChunk*)_alloc.Alloc(sizeof(BinChunk), KT_ALIGNOF(BinChunk));
	uint32_t const chunkIdx = _bin.m_numChunks++;
	newChunk->m_numTris = 0;
	newChunk->m_attribStride = _call.m_attributeBuffer.m_stride / sizeof(float);
	_bin.m_drawCallIndicies[chunkIdx] = _call.m_drawCallIdx;
	_bin.m_binChunks[chunkIdx] = newChunk;
	return *newChunk;
}

static void SetupEdge(BinChunk::EdgeEq& _e, uint32_t const _idx, int32_t const (&_v0)[2], int32_t const (&_v1)[2])
{
	int32_t const dy = (_v1[1] - _v0[1]);
	int32_t const dx = (_v0[0] - _v1[0]);

	int64_t c = _v0[1] * (_v1[0] - _v0[0]) - _v0[0] * (_v1[1] - _v0[1]);

	// Left/horizontal fill rule
	if (dy < 0 || dx == 0 && dy > 0)
	{
		c += 1;
	}

	_e.dy[_idx] = (_v1[1] - _v0[1]);
	_e.dx[_idx] = (_v0[0] - _v1[0]);
	_e.c[_idx] = uint32_t(c >> Config::c_subPixelBits);
}

static void SetupPlane
(
	float const _constantTerm,
	kt::Vec2 const& _d10,
	kt::Vec2 const& _d20,
	float const _attrib_d10,
	float const _attrib_d20,
	float& o_dx,
	float& o_dy
)
{
	float const A = _d10.y * _attrib_d20 - _attrib_d10 * _d20.y;
	float const B = _d20.x * _attrib_d10 - _d10.x * _attrib_d20;

	o_dx = -A / _constantTerm;
	o_dy = -B / _constantTerm;
}

static void BinTransformedAndClippedTri
(
	BinContext& _ctx, 
	ThreadScratchAllocator& _alloc, 
	uint32_t _threadIdx, 
	kt::Vec4 const& _v0, 
	kt::Vec4 const& _v1, 
	kt::Vec4 const& _v2, 
	float const* (&_attribPtrs)[3], 
	DrawCall const& _call
)
{
	uint32_t const height = _call.m_frameBuffer->m_height;
	uint32_t const width = _call.m_frameBuffer->m_width;
	kt::Vec2 const halfScreenCoords = kt::Vec2((float)width, (float)height) * 0.5f;

	float const invW[3] = { 1.0f / _v0.w, 1.0f / _v1.w, 1.0f / _v2.w };

	kt::Vec2 const v0raster = kt::Vec2(invW[0] * _v0.x * halfScreenCoords.x + halfScreenCoords.x, invW[0] * _v0.y * -halfScreenCoords.y + halfScreenCoords.y);
	kt::Vec2 const v1raster = kt::Vec2(invW[1] * _v1.x * halfScreenCoords.x + halfScreenCoords.x, invW[1] * _v1.y * -halfScreenCoords.y + halfScreenCoords.y);
	kt::Vec2 const v2raster = kt::Vec2(invW[2] * _v2.x * halfScreenCoords.x + halfScreenCoords.x, invW[2] * _v2.y * -halfScreenCoords.y + halfScreenCoords.y);

	int32_t const v0_fp[2] = { int32_t(v0raster.x * Config::c_subPixelStep + 0.5f), int32_t(v0raster.y * Config::c_subPixelStep + 0.5f) };
	int32_t const v1_fp[2] = { int32_t(v1raster.x * Config::c_subPixelStep + 0.5f), int32_t(v1raster.y * Config::c_subPixelStep + 0.5f) };
	int32_t const v2_fp[2] = { int32_t(v2raster.x * Config::c_subPixelStep + 0.5f), int32_t(v2raster.y * Config::c_subPixelStep + 0.5f) };

	int64_t triArea2_fp = (v2_fp[0] - v0_fp[0]) * (v1_fp[1] - v0_fp[1]) - (v2_fp[1] - v0_fp[1]) * (v1_fp[0] - v0_fp[0]);

	// Todo: allow switching winding order
	if (triArea2_fp <= 0)
	{
		return;
	}

	uint32_t xmin, ymin;
	uint32_t xmax, ymax;

	BinChunk::EdgeEq edges;

	xmin = (uint16_t)kt::Clamp(((kt::Min(kt::Min(v0_fp[0], v1_fp[0]), v2_fp[0]) + Config::c_subPixelMask) >> Config::c_subPixelBits), 0, (int32_t)width - 1);
	ymin = (uint16_t)kt::Clamp(((kt::Min(kt::Min(v0_fp[1], v1_fp[1]), v2_fp[1]) + Config::c_subPixelMask) >> Config::c_subPixelBits), 0, (int32_t)height - 1);

	xmax = (uint16_t)kt::Clamp(((kt::Max(kt::Max(v0_fp[0], v1_fp[0]), v2_fp[0]) + Config::c_subPixelMask) >> Config::c_subPixelBits), 0, (int32_t)width - 1);
	ymax = (uint16_t)kt::Clamp(((kt::Max(kt::Max(v0_fp[1], v1_fp[1]), v2_fp[1]) + Config::c_subPixelMask) >> Config::c_subPixelBits), 0, (int32_t)height - 1);

	SetupEdge(edges, 0, v0_fp, v1_fp);
	SetupEdge(edges, 1, v1_fp, v2_fp);
	SetupEdge(edges, 2, v2_fp, v0_fp);

	kt::Vec2 const d10_raster = v1raster - v0raster;
	kt::Vec2 const d20_raster = v2raster - v0raster;

	float const plane_eq_c = d10_raster.x * d20_raster.y - d10_raster.y * d20_raster.x;

	BinChunk::PlaneEq zOverW_plane;

	SetupPlane(plane_eq_c, d10_raster, d20_raster, _v1.z * invW[1] - _v0.z * invW[0], _v2.z * invW[2] - _v0.z * invW[0], zOverW_plane.dx, zOverW_plane.dy);

	BinChunk::PlaneEq recipW_plane;

	SetupPlane(plane_eq_c, d10_raster, d20_raster, invW[1] - invW[0], invW[2] - invW[0], recipW_plane.dx, recipW_plane.dy);

	BinChunk::PlaneEq attribPlanes[Config::c_maxVaryings];

	for (uint32_t i = 0; i < _call.m_attributeBuffer.m_stride / sizeof(float); ++i)
	{
		float const attrib_d10 = _attribPtrs[1][i] * invW[1] - _attribPtrs[0][i] * invW[0];
		float const attrib_d20 = _attribPtrs[2][i] * invW[2] - _attribPtrs[0][i] * invW[0];
		SetupPlane(plane_eq_c, d10_raster, d20_raster, attrib_d10, attrib_d20, attribPlanes[i].dx, attribPlanes[i].dy);
	}

	uint32_t const binYmin = ymin >> Config::c_binHeightLog2;
	uint32_t const binYmax = ymax >> Config::c_binHeightLog2;

	uint32_t const binXmax = xmax >> Config::c_binWidthLog2;
	uint32_t const binXmin = xmin >> Config::c_binWidthLog2;

	uint32_t const numYbins = binYmax - binYmin + 1;
	uint32_t const numXbins = binYmax - binYmin + 1;

	bool doTileCoverageCheck = true;

	if (numYbins <= 2 && numXbins <= 2)
	{
		// 2x1 or 1x2 or 2x2.
		// 2x1 and 1x2 must cover all bins.
		// 2x2 common case straddling corner of 4 bins.
		// In all these cases we just skip check and directly bin.
		doTileCoverageCheck = false;
	}

	for (uint32_t binY = binYmin; binY <= binYmax; ++binY)
	{
		for (uint32_t binX = binXmin; binX <= binXmax; ++binX)
		{
			int32_t const binScreenX0 = binX * Config::c_binWidth;
			int32_t const binScreenX1 = binScreenX0 + Config::c_binWidth;

			int32_t const binScreenY0 = binY * Config::c_binHeight;
			int32_t const binScreenY1 = binScreenY0 + Config::c_binHeight;

			if (doTileCoverageCheck)
			{
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
			}

			// todo full block
			ThreadBin& bin = _ctx.LookupThreadBin(_threadIdx, binX, binY);

			BinChunk& chunk = GetOrCreateBinForDrawCall(_alloc, _ctx, bin, _call);

			KT_ASSERT(chunk.m_numTris < c_trisPerBinChunk);
			uint32_t const chunkTriIdx = chunk.m_numTris++;

			BinChunk::EdgeEq& outEdge = chunk.m_edgeEq[chunkTriIdx];
			outEdge = edges;

			for (uint32_t i = 0; i < 3; ++i)
			{
				// pre shift edge constant terms to top of tile
				outEdge.c[i] += (outEdge.dx[i] * binScreenY0) + (outEdge.dy[i] * binScreenX0);
			}

			outEdge.blockMinX = uint8_t(kt::Clamp<int32_t>(int32_t(xmin) - int32_t(binScreenX0), 0, Config::c_binWidth));
			outEdge.blockMaxX = uint8_t(kt::Clamp<int32_t>(int32_t(xmax) - int32_t(binScreenX0), 0, Config::c_binWidth));

			outEdge.blockMinY = uint8_t(kt::Clamp<int32_t>(int32_t(ymin) - int32_t(binScreenY0), 0, Config::c_binHeight));
			outEdge.blockMaxY = uint8_t(kt::Clamp<int32_t>(int32_t(ymax) - int32_t(binScreenY0), 0, Config::c_binHeight));

			// pre shift plane constants to tile top left
			float const screenV0dx = (float)binScreenX0 - v0raster.x;
			float const screenV0dy = (float)binScreenY0 - v0raster.y;

			chunk.m_recipW[chunkTriIdx] = recipW_plane;
			chunk.m_recipW[chunkTriIdx].c0 = recipW_plane.dx * screenV0dx + recipW_plane.dy * screenV0dy + invW[0];

			chunk.m_zOverW[chunkTriIdx] = zOverW_plane;
			chunk.m_zOverW[chunkTriIdx].c0 = zOverW_plane.dx * screenV0dx + zOverW_plane.dy * screenV0dy + _v0.z * invW[0];

			uint32_t const attribElementStride = chunk.m_attribStride;

			for (uint32_t i = 0; i < attribElementStride; ++i)
			{
				chunk.m_attribPlanes[chunkTriIdx * attribElementStride + i].dx = attribPlanes[i].dx;
				chunk.m_attribPlanes[chunkTriIdx * attribElementStride + i].dy = attribPlanes[i].dy;
				chunk.m_attribPlanes[chunkTriIdx * attribElementStride + i].c0 = attribPlanes[i].dx * screenV0dx + attribPlanes[i].dy * screenV0dy + _attribPtrs[0][i] * invW[0];
			}
		}
	}
}

void BinTrisEntry(BinContext& _ctx, ThreadScratchAllocator& _alloc, uint32_t _threadIdx, uint32_t _triIdxBegin, uint32_t _triIdxEnd, DrawCall const& _drawCall)
{
	for (uint32_t triIdx = _triIdxBegin; triIdx < _triIdxEnd; ++triIdx)
	{
		KT_ASSERT(triIdx < _drawCall.m_indexBuffer.m_num);

		uint32_t indicies[3];
		kt::Vec4 vtx[3];

		FetchIndicies(_drawCall, triIdx, indicies);
		FetchPositions(_drawCall, indicies, vtx);

		vtx[0] = _drawCall.m_mvp * vtx[0];
		vtx[1] = _drawCall.m_mvp * vtx[1];
		vtx[2] = _drawCall.m_mvp * vtx[2];

		uint32_t maskOr = 0;

		uint8_t clipv0 = ComputeClipMask(vtx[0]);
		uint8_t clipv1 = ComputeClipMask(vtx[1]);
		uint8_t clipv2 = ComputeClipMask(vtx[2]);

		maskOr = clipv0 | clipv1 | clipv2;

		if (maskOr == 0)
		{
			float const* originalAttribs[3];
			FetchAttribPointers(_drawCall, indicies, originalAttribs);
			BinTransformedAndClippedTri(_ctx, _alloc, _threadIdx, vtx[0], vtx[1], vtx[2], originalAttribs, _drawCall);
			continue;
		}
		else if (clipv0 & clipv1 & clipv2)
		{
			// If clip AND mask has any bits set, all verts are the wrong side of a clip plane, so the whole triangle can be culled.
			continue;
		}

		float const* originalAttribs[3];
		FetchAttribPointers(_drawCall, indicies, originalAttribs);

		ClipBuffer buf;

		buf.numInputVerts = 3;
		memcpy(buf.attribs[buf.inputIdx][0], originalAttribs[0], _drawCall.m_attributeBuffer.m_stride);
		memcpy(buf.attribs[buf.inputIdx][1], originalAttribs[1], _drawCall.m_attributeBuffer.m_stride);
		memcpy(buf.attribs[buf.inputIdx][2], originalAttribs[2], _drawCall.m_attributeBuffer.m_stride);
		buf.verts[buf.inputIdx][0] = vtx[0];
		buf.verts[buf.inputIdx][1] = vtx[1];
		buf.verts[buf.inputIdx][2] = vtx[2];

		do
		{
			uint32_t clipIdx = kt::Cnttz(maskOr);
			maskOr ^= (1 << clipIdx);
			ClipPlane(buf, clipIdx);
		} while (maskOr && buf.numInputVerts);

		// Fan triangulation
		kt::Vec4(&input_vec)[CLIP_VERT_BUFFER_SIZE] = buf.verts[buf.inputIdx];
		float(&input_attribs)[CLIP_VERT_BUFFER_SIZE][Config::c_maxVaryings] = buf.attribs[buf.inputIdx];
		for (uint32_t i = 2; i < buf.numInputVerts; ++i)
		{
			float const* attribPtrs[3] = { input_attribs[0], input_attribs[i - 1], input_attribs[i] };
			BinTransformedAndClippedTri(_ctx, _alloc, _threadIdx, input_vec[0], input_vec[i - 1], input_vec[i], attribPtrs, _drawCall);
		}
	}
}

}
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


void FillScreenTest(Framebuffer& _buffer, uint8_t const color[3])
{
	uint8_t* p = _buffer.ptr;
	uint8_t* pEnd = _buffer.ptr + _buffer.width * _buffer.height * 3;

	while (p != pEnd)
	{
		*p++ = color[0];
		*p++ = color[1];
		*p++ = color[2];
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
		x = (_v1[1] - _v0[1]);
		y = (_v0[0] - _v1[0]);

		// Left/horizontal fill rule
		if (x < 0 || y == 0 && x > 0)
		{
			c += c_subPixelStep;
		}

		cur_base = c + x * (_xmin << c_subPixelBits) + y * (_ymin << c_subPixelBits);
		eval = cur_base;
		x <<= c_subPixelBits; // so they become subpixel deltas
		y <<= c_subPixelBits;
	}

	void IncX()
	{
		eval += x;
	}

	void NextLine()
	{
		cur_base += y;
		eval = cur_base;
	}

	int32_t c;
	int32_t x;
	int32_t y;

	int32_t cur_base;

	int32_t eval;
};

static void RasterTransformedTri(Framebuffer& _buffer, DepthBuffer& _depthBuffer, Transformed_Tri const& _tri)
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
			uint8_t* pixel = _buffer.ptr + y * _buffer.width * 3 + x * 3;

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
#elif 1
					//pixel[0] = (uint8_t)(v0_persp * 255);
					//pixel[1] = (uint8_t)(v1_persp * 255);
					//pixel[2] = (uint8_t)(v2_persp * 255);
					pixel[0] = uint8_t((vary0[0] * v0_persp + vary1[0] * v1_persp + vary2[0] * v2_persp) * 255.0f);
					pixel[1] = uint8_t((vary0[1] * v0_persp + vary1[1] * v1_persp + vary2[1] * v2_persp) * 255.0f);
					pixel[2] = uint8_t((vary0[2] * v0_persp + vary1[2] * v1_persp + vary2[2] * v2_persp) * 255.0f);
#else
					pixel[0] = thecol;
					pixel[1] = thecol;
					pixel[2] = thecol;
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


static void TransformAndClip(kt::Vec4 const& _v0, kt::Vec4 const& _v1, kt::Vec4 const& _v2, kt::Mat4 const& _mtx, Transformed_Tri o_tris[CLIP_OUT_TRI_BUFF_SIZE], uint32_t& o_numOutTris)
{

}

static bool SetupTri(Transformed_Vert const& _v0, Transformed_Vert const& _v1, Transformed_Vert const& _v2, Transformed_Tri& o_tri, uint32_t const _screenWidth, uint32_t const _screenHeight)
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

	if (triArea2_fp <= 0.0f)
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

void ClipAndSetup(Transformed_Vert clipSpaceVerts[3], Transformed_Tri o_tris[CLIP_OUT_TRI_BUFF_SIZE], uint32_t& o_numOutTris)
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

	auto copyVert = [](Transformed_Vert& out, Transformed_Vert const& _in) 
	{
		memcpy(&out, &_in, sizeof(Transformed_Vert));
	};

	Transformed_Vert input[CLIP_BUFFER_SIZE];
	copyVert(input[0], clipSpaceVerts[0]);
	copyVert(input[1], clipSpaceVerts[1]);
	copyVert(input[2], clipSpaceVerts[2]);

	Transformed_Vert output[CLIP_BUFFER_SIZE];

	uint32_t numInputVerts = 3;
	uint32_t numOutputVerts = 0;

	do
	{
		uint32_t clipIdx = kt::Cnttz(maskOr);
		maskOr ^= (1 << clipIdx);

		kt::Vec4 const& clipPlane = c_clipPlanes[clipIdx];

		Transformed_Vert const* v0 = &input[numInputVerts - 1];
		float v0_dot = kt::Dot(clipPlane, v0->transformedPos);

		for (uint32_t nextClipVertIdx = 0; nextClipVertIdx < numInputVerts; ++nextClipVertIdx)
		{
			Transformed_Vert const*const v1 = &input[nextClipVertIdx];
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
					Transformed_Vert* out = &output[numOutputVerts++];
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
					Transformed_Vert* out = &output[numOutputVerts++];
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

		memcpy(input, output, sizeof(Transformed_Vert) * numOutputVerts); // this is slow - should have pointer indirection into buffers so no copy is needed.
		numInputVerts = numOutputVerts;
		numOutputVerts = 0;
	} while (maskOr);

	// Fan triangulation
	o_numOutTris = 0;
	for (uint32_t i = 2; i < numInputVerts; ++i)
	{
		KT_ASSERT(o_numOutTris <= CLIP_OUT_TRI_BUFF_SIZE);
		Transformed_Tri& tri = o_tris[o_numOutTris];
		if (SetupTri(input[0], input[i - 1], input[i], o_tris[o_numOutTris], 1280, 720)) // ahhh hardcoded
		{
			++o_numOutTris;
		}
	}
}


void RasterClippedTri(Framebuffer& _buffer, DepthBuffer& _depthBuffer, kt::Vec4 const& _v0, kt::Vec4 const& _v1, kt::Vec4 const& _v2)
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
}

void SetupAndRasterTriTest(Framebuffer& _buffer, DepthBuffer& _depthBuffer, kt::Mat4 const& _mtx, kt::Vec3 const& _v0, kt::Vec3 const& _v1, kt::Vec3 const& _v2)
{
	Transformed_Tri outTris[CLIP_OUT_TRI_BUFF_SIZE];
	uint32_t numClipTris;
	
	Transformed_Vert verts[3];
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
		RasterTransformedTri(_buffer, _depthBuffer, outTris[i]);
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
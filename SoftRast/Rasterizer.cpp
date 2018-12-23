#include "Rasterizer.h"
#include <float.h>

namespace Raster
{


void FillScreenTest(Framebuffer const& _buffer, uint8_t const color[3])
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

struct EdgeEq
{
	void Set(kt::Vec2 const& _v0, kt::Vec2 const& _v1)
	{
		ec = (-_v0.x * (_v1.y - _v0.y) + _v0.y * (_v1.x - _v0.x));
		ex = (_v1.y - _v0.y);
		ey = (_v0.x - _v1.x);
	}

	float Eval(int32_t _x, int32_t _y)
	{
		return ec + ex * _x + ey * _y;
	}

	float ec;
	float ex;
	float ey;
};

float EvalEdge(kt::Vec2 const& _v0, kt::Vec2 const& _v1, int32_t _x, int32_t _y)
{
	float a = (_v1.x - _v0.x) * (_y - _v0.y);
	float b = (_v1.y - _v0.y) * (_x - _v0.x);
	return b - a;
}

void RasterTriTest(Framebuffer const& _buffer, kt::Vec3 const& _v0, kt::Vec3 const& _v1, kt::Vec3 const& _v2)
{
	kt::Vec2 bboxMin(FLT_MAX);
	kt::Vec2 bboxMax(-FLT_MAX);

	kt::Vec2 const halfScreenCoords(_buffer.width*0.5f, _buffer.height*0.5f);

	// v2 and v1 are swapped to preserve winding order when flipping y (from NDC -> screen)
	kt::Vec2 const v0raster = kt::Vec2(_v0.x * halfScreenCoords.x + halfScreenCoords.x, _v0.y * -halfScreenCoords.y + halfScreenCoords.y);
	kt::Vec2 const v1raster = kt::Vec2(_v2.x * halfScreenCoords.x + halfScreenCoords.x, _v2.y * -halfScreenCoords.y + halfScreenCoords.y);
	kt::Vec2 const v2raster = kt::Vec2(_v1.x * halfScreenCoords.x + halfScreenCoords.x, _v1.y * -halfScreenCoords.y + halfScreenCoords.y);

	kt::Vec2 const edge0 = v1raster - v0raster;
	kt::Vec2 const edge1 = v2raster - v1raster;
	kt::Vec2 const edge2 = v0raster - v2raster;

	EdgeEq e0, e1, e2;
	e0.Set(v0raster, v1raster);
	e1.Set(v1raster, v2raster);
	e2.Set(v2raster, v0raster);

	bboxMin.x = kt::Min(kt::Min(v0raster.x, v1raster.x), v2raster.x);
	bboxMax.x = kt::Max(kt::Max(v0raster.x, v1raster.x), v2raster.x);

	bboxMin.y = kt::Min(kt::Min(v0raster.y, v1raster.y), v2raster.y);
	bboxMax.y = kt::Max(kt::Max(v0raster.y, v1raster.y), v2raster.y);

	bboxMin = kt::Clamp(bboxMin, kt::Vec2(0.0f), kt::Vec2((float)_buffer.width - 1, (float)_buffer.height - 1));
	bboxMax = kt::Clamp(bboxMax, kt::Vec2(0.0f), kt::Vec2((float)_buffer.width - 1, (float)_buffer.height - 1));

	int32_t xmin, ymin;
	int32_t xmax, ymax;

	xmin = (int32_t)bboxMin.x;
	ymin = (int32_t)bboxMin.y;
	xmax = (int32_t)bboxMax.x;
	ymax = (int32_t)bboxMax.y;

	float const recipTriArea2 = 1.0f / kt::Cross(v2raster - v0raster, v1raster - v0raster);

	for (int32_t y = ymin; y <= ymax; ++y)
	{
		for (int32_t x = xmin; x <= xmax; ++x)
		{
			bool inTri = true;
#if 1
			float const e0Eval = e0.Eval(x, y);
			float const e1Eval = e1.Eval(x, y);
			float const e2Eval = e2.Eval(x, y);

			inTri &= (e0Eval == 0.0f) ? (edge0.y == 0.0f && edge0.x > 0.0f) || edge0.y > 0.0f : (e0Eval > 0.0f);
			inTri &= (e1Eval == 0.0f) ? (edge1.y == 0.0f && edge1.x > 0.0f) || edge1.y > 0.0f : (e1Eval > 0.0f);
			inTri &= (e2Eval == 0.0f) ? (edge2.y == 0.0f && edge2.x > 0.0f) || edge2.y > 0.0f : (e2Eval > 0.0f);

			//inTri &= e0Eval >= 0.0f;
			//inTri &= e1Eval >= 0.0f;
			//inTri &= e2Eval >= 0.0f;
#else
			float const eval0 = EvalEdge(v0raster, v1raster, x, y);
			float const eval1 = EvalEdge(v1raster, v2raster, x, y);
			float const eval2 = EvalEdge(v2raster, v0raster, x, y);
			inTri = eval0 >= 0.0f && eval1 >= 0.0f && eval2 >= 0.0f;
#endif
			uint8_t* pixel = _buffer.ptr + y * _buffer.width * 3 + x * 3;
			
			if (inTri)
			{
				float const baryV0 = e1Eval * recipTriArea2;
				float const baryV1 = e2Eval * recipTriArea2;
				float const baryV2 = e0Eval * recipTriArea2;

				pixel[0] = (uint8_t)(baryV0 * 255);
				pixel[1] = (uint8_t)(baryV1 * 255);
				pixel[2] = (uint8_t)(baryV2 * 255);
			}
		}
	}
}

}
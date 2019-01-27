#pragma once
#include <kt/kt.h>
#include <kt/Vec4.h>
#include "Renderer.h"

namespace kt
{
struct Vec3;
struct Mat4;
struct IAllocator;
}

namespace sr
{

namespace Tex
{
struct TextureData;
}

namespace Raster
{
struct FrameBuffer;
struct DepthBuffer;

struct RasterContext
{
	FrameBuffer* frameBuffer;
	DepthBuffer* depthBuffer;

	uint32_t vpHeight;
	uint32_t vpWidth;
};


struct PipelineVert
{
	PipelineVert(){}

	union
	{
		kt::Vec4 transformedPos;
		
		struct  
		{
			float raster_x;
			float raster_y;
			float zOverW;
			float recipW;
		};
	};

	KT_ALIGNAS(16) float varyings[c_maxVaryings];
};

struct PipelineTri
{
	PipelineVert verts[3];

	// Raster space recip determinant.
	float halfRecipTriArea_fp;

	// Fixed point verticies in raster space.
	int32_t v0_fp[2];
	int32_t v1_fp[2];
	int32_t v2_fp[2];

	// Raster space bbox.
	int32_t xmin, xmax;
	int32_t ymin, ymax;
};

struct FrameBuffer
{
	// R8G8B8A8 framebuffer
	uint8_t* ptr;

	// Framebuffer width
	uint32_t width;

	// Framebuffer height
	uint32_t height;
};

struct DepthBuffer
{
	void Init(kt::IAllocator* _allocator, uint32_t _width, uint32_t _height);

	~DepthBuffer();

	float* At(uint32_t _x, uint32_t _y)
	{
		KT_ASSERT(_x < width && _y < height);
		return ptr + _y * width + _x;
	}

	float* ptr;

	uint32_t width;

	uint32_t height;

	kt::IAllocator* allocator = nullptr;
};

enum class WindingOrder
{
	CW,
	CCW,
	Default = CCW
};

void ClearDepthBufferTest(DepthBuffer& _buffer);

void FillScreenTest(FrameBuffer& _buffer, uint8_t const color[3]);

// Raster a tri (CCW winding)
void SetupAndRasterTriTest(FrameBuffer& _buffer, DepthBuffer& _depthBuffer, kt::Mat4 const& _mtx, kt::Vec3 const& _v0, kt::Vec3 const& _v1, kt::Vec3 const& _v2);

void DrawSerial_Test(FrameBuffer& _buffer, DepthBuffer& _depthBuffer, kt::Mat4 const& _mtx, Renderer::DrawCall const& _call);


}

}
#pragma once
#include <kt/Vec3.h>

struct Framebuffer
{
	// R8G8B8 framebuffer
	uint8_t* ptr;

	// Framebuffer width
	uint32_t width;

	// Framebuffer height
	uint32_t height;
};

namespace Raster
{


void FillScreenTest(Framebuffer const& _buffer, uint8_t const color[3]);
void RasterTriTest(Framebuffer const& _buffer, kt::Vec3 const& _v0, kt::Vec3 const& _v1, kt::Vec3 const& _v2);

}
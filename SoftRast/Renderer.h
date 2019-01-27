#pragma once
#include <immintrin.h>
#include <stdint.h>

#include "SoftRastTypes.h"
#include "kt/Array.h"


namespace sr
{

namespace Renderer
{

// Always RGBA8
struct ColourTile
{
	static uint32_t const c_bytesPerPixel = 4;

	KT_ALIGNAS(32) uint8_t m_colour[Config::c_tileHeight * Config::c_tileWidth * c_bytesPerPixel];
	char _pad_[64];
};

struct DepthTile
{
	KT_ALIGNAS(32) float m_depth[Config::c_tileWidth * Config::c_tileWidth];
	char _pad_[64];
};

struct FrameBuffer
{
	KT_NO_COPY(FrameBuffer);
	
	FrameBuffer(uint32_t _width, uint32_t _height, bool _colour = true, bool _depth = true);
	~FrameBuffer();

	void Init(uint32_t _width, uint32_t _height, bool _colour = true, bool _depth = true);

	ColourTile* m_colourTiles = nullptr;
	DepthTile* m_depthTiles = nullptr;

	uint32_t m_height = 0;
	uint32_t m_width = 0;

	uint32_t m_tilesY = 0;
	uint32_t m_tilesX = 0;
};

using PixelShaderFn = void(*)(void const* _uniforms, __m256 const _varyings[sr::Config::c_maxVaryings], float o_colour[4 * 8], __m256 const& _execMask);

struct GenericDrawBuffer
{
	void const* m_ptr = nullptr;
	uint32_t m_num = 0;
	uint32_t m_stride = 0;
};

struct DrawCall
{
	DrawCall();

	DrawCall& SetPixelShader(PixelShaderFn _fn, void const* _uniform);
	DrawCall& SetIndexBuffer(void const* _buffer, uint32_t const _stride, uint32_t const _num);
	DrawCall& SetPositionBuffer(void const* _buffer, uint32_t const _stride, uint32_t const _num);
	DrawCall& SetAttributeBuffer(void const* _buffer, uint32_t const _stride, uint32_t const _num);
	DrawCall& SetFrameBuffer(FrameBuffer const* _buffer);

	PixelShaderFn m_pixelShader = nullptr;
	void const* m_pixelUniforms = nullptr;

	GenericDrawBuffer m_indexBuffer;
	GenericDrawBuffer m_positionBuffer;
	GenericDrawBuffer m_attributeBuffer;
	
	FrameBuffer const* m_frameBuffer = nullptr;

	uint32_t m_colourWrite		: 1;
	uint32_t m_depthWrite		: 1;
	uint32_t m_depthRead		: 1;
};


class Context
{
public:

	Context();
	~Context();

	void DrawIndexed(DrawCall const& _call);

	void BeginFrame();
	void EndFrame();

private:
	kt::Array<DrawCall> m_drawCalls;
};




}
}
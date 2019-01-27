#include "Renderer.h"
#include "kt/Memory.h"

namespace sr
{

namespace Renderer
{


FrameBuffer::FrameBuffer(uint32_t _width, uint32_t _height, bool _colour /*= true*/, bool _depth /*= true*/)
{
	
}

FrameBuffer::~FrameBuffer()
{
	kt::Free(m_depthTiles);
	kt::Free(m_colourTiles);
}

void FrameBuffer::Init(uint32_t _width, uint32_t _height, bool _colour /*= true*/, bool _depth /*= true*/)
{
	m_height = _height;
	m_width = _width;

	m_tilesX = (uint32_t)kt::AlignValue(m_width, Config::c_tileWidth) >> Config::c_tileWidthLog2;
	m_tilesY = (uint32_t)kt::AlignValue(m_height, Config::c_tileHeight) >> Config::c_tileHeightLog2;

	if (_colour)
	{
		m_colourTiles = (ColourTile*)kt::Malloc(sizeof(ColourTile) * m_tilesX * m_tilesY, KT_ALIGNOF(ColourTile));
	}

	if (_depth)
	{
		m_depthTiles = (DepthTile*)kt::Malloc(sizeof(DepthTile) * m_tilesX * m_tilesY, KT_ALIGNOF(DepthTile));
	}
}

DrawCall::DrawCall()
	: m_colourWrite(1)
	, m_depthRead(1)
	, m_depthWrite(1)
{

}

DrawCall& DrawCall::SetPixelShader(PixelShaderFn _fn, void const* _uniform)
{
	m_pixelShader = _fn;
	m_pixelUniforms = _uniform;
	return *this;
}

DrawCall& DrawCall::SetIndexBuffer(void const* _buffer, uint32_t const _stride, uint32_t const _num)
{
	m_indexBuffer.m_num = _num;
	m_indexBuffer.m_ptr = _buffer;
	m_indexBuffer.m_stride = _stride;
	return *this;
}

DrawCall& DrawCall::SetPositionBuffer(void const* _buffer, uint32_t const _stride, uint32_t const _num)
{
	m_positionBuffer.m_num = _num;
	m_positionBuffer.m_ptr = _buffer;
	m_positionBuffer.m_stride = _stride;
	return *this;
}

DrawCall& DrawCall::SetAttributeBuffer(void const* _buffer, uint32_t const _stride, uint32_t const _num)
{
	m_attributeBuffer.m_num = _num;
	m_attributeBuffer.m_ptr = _buffer;
	m_attributeBuffer.m_stride = _stride;
	return *this;
}

DrawCall& DrawCall::SetFrameBuffer(FrameBuffer const* _buffer)
{
	m_frameBuffer = _buffer;
	return *this;
}


void Context::DrawIndexed(DrawCall const& _call)
{
	KT_ASSERT(_call.m_indexBuffer.m_ptr && "No index buffer bound.");
	m_drawCalls.PushBack(_call);
}

}

}
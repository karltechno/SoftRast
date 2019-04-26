#pragma once
#include <immintrin.h>
#include <stdint.h>

#include <kt/Array.h>
#include <kt/Mat4.h>

#include <atomic>

#include "SoftRastTypes.h"
#include "TaskSystem.h"
#include "Binning.h"
#include "Config.h"


namespace sr
{


// Always RGBA8
struct ColourTile
{
	static uint32_t const c_bytesPerPixel = 4;

	void Clear(uint32_t _col);

	KT_ALIGNAS(32) uint8_t m_colour[Config::c_binHeight * Config::c_binWidth * c_bytesPerPixel];
};

struct DepthTile
{
	void Clear();

	KT_ALIGNAS(32) float m_depth[Config::c_binWidth * Config::c_binWidth];
	float m_hiZmin;
	float m_hiZmax;
};

struct FrameBufferPlane
{
	FrameBufferPlane() = default;
	KT_NO_COPY(FrameBufferPlane);

	void Init(uint32_t _width, uint32_t _height, bool _colour = true, bool _depth = true);

	ColourTile* m_colourTiles = nullptr;
	DepthTile* m_depthTiles = nullptr;

	uint32_t m_height = 0;
	uint32_t m_width = 0;

	uint32_t m_tilesY = 0;
	uint32_t m_tilesX = 0;
};

struct FrameBuffer
{
	KT_NO_COPY(FrameBuffer);
	FrameBuffer(uint32_t _width, uint32_t _height, bool _colour = true, bool _depth = true);
	~FrameBuffer();

	FrameBufferPlane* WritePlane()
	{
		return &m_bufferedPlanes[m_writePlane];
	}

	FrameBufferPlane const* ReadPlane() const
	{
		return &m_bufferedPlanes[m_writePlane ^ 1];
	}

	void SwapPlanes()
	{
		m_writePlane ^= 1;
	}

	FrameBufferPlane m_bufferedPlanes[2];
	uint32_t m_writePlane = 0;

	struct JobData
	{
		Task m_task;
		FrameBufferPlane* m_plane = nullptr;
		uint8_t* m_linearPixels = nullptr;
		std::atomic<uint32_t> m_counter{ 0 };
		void(*m_onFinishBlit)(void*) = nullptr;
		void* m_onFinishBlitUser = nullptr;
	} m_jobs[2];
};

using PixelShaderFn = void(void const* _uniforms, float const* _varyings, uint32_t o_texels[8], uint32_t _execMask);

using VertexShaderFn = void(kt::Vec3 const& _vtx, void const* i_uniforms, void const* i_attribs, float* o_attribs);

struct GenericDrawBuffer
{
	void const* m_ptr = nullptr;
	uint32_t m_num = 0;
	uint32_t m_stride = 0;
};

struct DrawCall
{
	static const uint32_t UV_OFFSET_INVALID = 0xFFFFFFFF;

	DrawCall();

	DrawCall& SetVertexShader(VertexShaderFn* _fn, void const* _uniforms, uint32_t const _outAttributeStrideBytes);
	DrawCall& SetPixelShader(PixelShaderFn* _fn, void const* _uniforms);
	DrawCall& SetIndexBuffer(void const* _buffer, uint32_t const _stride, uint32_t const _num);
	DrawCall& SetPositionBuffer(void const* _buffer, uint32_t const _stride, uint32_t const _num);
	DrawCall& SetAttributeBuffer(void const* _buffer, uint32_t const _stride, uint32_t const _num, uint32_t const _uvOffset = 0);
	DrawCall& SetFrameBuffer(FrameBuffer* _buffer);
	DrawCall& SetMVP(kt::Mat4 const& _mvp);

	VertexShaderFn* m_vertexShader = nullptr;
	void const* m_vertexUniforms = nullptr;
	uint32_t m_outAttributeStrideBytes = 0;

	PixelShaderFn* m_pixelShader = nullptr;
	void const* m_pixelUniforms = nullptr;

	GenericDrawBuffer m_indexBuffer;
	GenericDrawBuffer m_positionBuffer;
	GenericDrawBuffer m_attributeBuffer;

	uint32_t m_uvOffset = 0;
	
	FrameBufferPlane const* m_frameBuffer = nullptr;

	kt::Mat4 m_mvp = kt::Mat4::Identity();

	uint32_t m_drawCallIdx = 0;

	uint32_t m_colourWrite		: 1;
	uint32_t m_depthWrite		: 1;
	uint32_t m_depthRead		: 1;
};


class RenderContext
{
public:
	RenderContext();
	~RenderContext();

	void DrawIndexed(DrawCall const& _call);

	void ClearFrameBuffer(FrameBuffer& _buffer, uint32_t _color = 0x00000000, bool _clearColour = true, bool _clearDepth = true);

	ThreadScratchAllocator& ThreadAllocator();

	void BeginFrame();
	void EndFrame();

	void Blit(FrameBuffer& _fb, uint8_t* _linearPixels, void(*_onFinishBlit)(void*) = nullptr, void* _onFinishUser = nullptr);

private:
	TaskSystem m_taskSystem;

	BinContext m_binner;
	kt::Array<DrawCall> m_drawCalls;
};




}

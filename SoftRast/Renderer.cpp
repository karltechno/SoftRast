#include "Renderer.h"
#include "Rasterizer.h"
#include "kt/Memory.h"
#include "kt/Logging.h"

namespace sr
{


FrameBuffer::FrameBuffer(uint32_t _width, uint32_t _height, bool _colour /*= true*/, bool _depth /*= true*/)
{
	for (uint32_t i = 0; i < KT_ARRAY_COUNT(m_bufferedPlanes); ++i)
	{
		m_bufferedPlanes[i].Init(_width, _height, _colour, _depth);
	}
}

FrameBuffer::~FrameBuffer()
{
	for (uint32_t i = 0; i < KT_ARRAY_COUNT(m_bufferedPlanes); ++i)
	{
		kt::Free(m_bufferedPlanes[i].m_depthTiles);
		kt::Free(m_bufferedPlanes[i].m_colourTiles);
	}
}

void FrameBufferPlane::Init(uint32_t _width, uint32_t _height, bool _colour /*= true*/, bool _depth /*= true*/)
{
	m_height = _height;
	m_width = _width;

	m_tilesX = (uint32_t)kt::AlignUp(m_width, Config::c_binWidth) >> Config::c_binWidthLog2;
	m_tilesY = (uint32_t)kt::AlignUp(m_height, Config::c_binHeight) >> Config::c_binHeightLog2;

	if (_colour)
	{
		m_colourTiles = (ColourTile*)kt::Malloc(sizeof(ColourTile) * m_tilesX * m_tilesY, KT_ALIGNOF(ColourTile));
	}

	if (_depth)
	{
		m_depthTiles = (DepthTile*)kt::Malloc(sizeof(DepthTile) * m_tilesX * m_tilesY, KT_ALIGNOF(DepthTile));
	}
}
//
//void FrameBufferPlane::BlitDepth(uint8_t* _linearFramebuffer)
//{
//	// hack/slow
//	uint32_t* fb32 = (uint32_t*)_linearFramebuffer;
//	for (uint32_t tileY = 0; tileY < m_tilesY; ++tileY)
//	{
//		for (uint32_t tileX = 0; tileX < m_tilesX; ++tileX)
//		{
//			DepthTile const& tile = m_depthTiles[tileY * m_tilesX + tileX];
//
//			uint32_t const yEnd = kt::Min(Config::c_binHeight, m_height - tileY * Config::c_binHeight);
//			uint32_t const widthCopySize = kt::Min(Config::c_binWidth, m_width - tileX * Config::c_binWidth);
//
//			for (uint32_t y = 0; y < yEnd; ++y)
//			{
//				float const* src = &tile.m_depth[y * Config::c_binWidth];
//				uint32_t* dest = fb32 + tileY * Config::c_binHeight * m_width + tileX * Config::c_binWidth + y * m_width;
//				for (uint32_t x = 0; x < widthCopySize; ++x)
//				{
//					*dest++ = uint8_t(*src++ * 255.0f);
//				}
//			}
//		}
//	}
//}

DrawCall::DrawCall()
	: m_colourWrite(1)
	, m_depthWrite(1)
	, m_depthRead(1)
{
}


DrawCall& DrawCall::SetVertexShader(VertexShaderFn* _fn, void const* _uniforms, uint32_t const _outAttributeStrideBytes)
{
	KT_ASSERT(_outAttributeStrideBytes <= Config::c_maxVaryings * sizeof(float));
	m_vertexShader = _fn;
	m_vertexUniforms = _uniforms;
	m_outAttributeStrideBytes = _outAttributeStrideBytes;
	return *this;
}

DrawCall& DrawCall::SetPixelShader(PixelShaderFn _fn, void const* _uniforms)
{
	m_pixelShader = _fn;
	m_pixelUniforms = _uniforms;
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

DrawCall& DrawCall::SetAttributeBuffer(void const* _buffer, uint32_t const _stride, uint32_t const _num, uint32_t const _uvOffset)
{
	m_uvOffset = _uvOffset;
	m_attributeBuffer.m_num = _num;
	m_attributeBuffer.m_ptr = _buffer;
	m_attributeBuffer.m_stride = _stride;
	return *this;
}

DrawCall& DrawCall::SetFrameBuffer(FrameBuffer* _buffer)
{
	m_frameBuffer = _buffer->WritePlane();
	return *this;
}


DrawCall& DrawCall::SetMVP(kt::Mat4 const& _mvp)
{
	m_mvp = _mvp;
	return *this;
}

RenderContext::RenderContext()
{
#if !SR_DEBUG_SINGLE_THREADED
	m_taskSystem.InitFromMainThread(kt::LogicalCoreCount() - 1);
#else
	m_taskSystem.InitFromMainThread(0);
#endif

	// Todo: frame buffer size hardcoded!!
	m_binner.Init(m_taskSystem.TotalThreadsIncludingMainThread(), uint32_t(kt::AlignUp(Config::c_screenWidth, Config::c_binWidth)) / Config::c_binWidth, 
				  uint32_t(kt::AlignUp(Config::c_screenHeight, Config::c_binHeight)) / Config::c_binHeight);

}

RenderContext::~RenderContext()
{
	m_taskSystem.WaitAndShutdown();
}

void RenderContext::DrawIndexed(DrawCall const& _call)
{
	KT_ASSERT(_call.m_indexBuffer.m_ptr && "No index buffer bound.");
	m_drawCalls.PushBack(_call);
	m_drawCalls.Back().m_drawCallIdx = m_drawCalls.Size() - 1;
}

void RenderContext::ClearFrameBuffer(FrameBuffer& _buffer, uint32_t _color, bool _clearColour /*= true*/, bool _clearDepth /*= true*/)
{
	// TODO: Fast clear
	FrameBufferPlane const& plane = *_buffer.WritePlane();
	if (_clearDepth)
	{
		for (uint32_t i = 0; i < (plane.m_tilesX * plane.m_tilesY); ++i)
		{
			for (uint32_t j = 0; j < (Config::c_binHeight * Config::c_binWidth); ++j)
			{
				plane.m_depthTiles[i].m_depth[j] = Config::c_depthMax;
			}
		}
	}

	if (_clearColour)
	{
		for (uint32_t i = 0; i < (plane.m_tilesX * plane.m_tilesY); ++i)
		{	
			memset(plane.m_colourTiles[i].m_colour, _color, sizeof(plane.m_colourTiles[i].m_colour));
		}
	}
}

ThreadScratchAllocator& RenderContext::ThreadAllocator()
{
	return m_taskSystem.ThreadAllocator();
}

void RenderContext::BeginFrame()
{
	m_drawCalls.Clear();
	m_taskSystem.ResetAllocators();
}

void RenderContext::EndFrame()
{
	for (uint32_t i = 0; i < m_binner.m_numBinsX * m_binner.m_numBinsY * m_binner.m_numThreads; ++i)
	{
		// todo frame number dirty
		m_binner.m_bins[i].Reset();
	}

	struct BinTrisTaskData
	{
		DrawCall const* call;
		RenderContext* ctx;
	};

	Task* drawCallTasks = (Task*)KT_ALLOCA(sizeof(Task) * m_drawCalls.Size());
	BinTrisTaskData* drawCallTasksData = (BinTrisTaskData*)KT_ALLOCA(sizeof(BinTrisTaskData) * m_drawCalls.Size());

	std::atomic<uint32_t> frontEndCounter(0);

	for (uint32_t i = 0; i < m_drawCalls.Size(); ++i)
	{
		DrawCall const& draw = m_drawCalls[i];

		auto drawCallTaskFn = [](Task const* _task, uint32_t _threadIdx, uint32_t _start, uint32_t _end)
		{
			BinTrisTaskData* data = (BinTrisTaskData*)_task->m_userData;
			BinTrisEntry(data->ctx->m_binner, data->ctx->ThreadAllocator(), _threadIdx, _start, _end, *data->call);
		};

		BinTrisTaskData* taskData = drawCallTasksData + i;

		Task* task = drawCallTasks + i;

		kt::PlacementNew(task, drawCallTaskFn, draw.m_indexBuffer.m_num / 3, 512, taskData);
		taskData->call = &draw;
		taskData->ctx = this;
		task->m_taskCounter = &frontEndCounter;

		m_taskSystem.PushTask(task);
	}

	m_taskSystem.WaitForCounter(&frontEndCounter);

	std::atomic<uint32_t> tileRasterCounter{ 0 };

	for (uint32_t binY = 0; binY < m_binner.m_numBinsY; ++binY)
	{
		for (uint32_t binX = 0; binX < m_binner.m_numBinsX; ++binX)
		{
			bool anyTris = false;
			for (uint32_t threadIdx = 0; threadIdx < m_binner.m_numThreads; ++threadIdx)
			{
				ThreadBin& bin = m_binner.LookupThreadBin(threadIdx, binX, binY);
				anyTris |= bin.m_numChunks != 0;
			}

			if (anyTris)
			{
				struct TileTaskData
				{
					Task t;
					ThreadRasterCtx rasterCtx;
				};

				auto tileRasterFn = [](Task const* _task, uint32_t _threadIdx, uint32_t _start, uint32_t _end)
				{
					TileTaskData* data = (TileTaskData*)_task->m_userData;
					RasterAndShadeBin(data->rasterCtx);
				};

				TileTaskData* t = (TileTaskData*)KT_ALLOCA(sizeof(TileTaskData));
				t->rasterCtx.m_binner = &m_binner;
				t->rasterCtx.m_tileX = binX;
				t->rasterCtx.m_tileY = binY;
				t->rasterCtx.m_drawCalls = m_drawCalls.Data();
				t->rasterCtx.m_numDrawCalls = m_drawCalls.Size();
				t->rasterCtx.m_ctx = this;
				kt::PlacementNew(&t->t, tileRasterFn, 1, 1, t);

				t->t.m_taskCounter = &tileRasterCounter;
				m_taskSystem.PushTask(&t->t);
			}

		}
	}

	m_taskSystem.WaitForCounter(&tileRasterCounter);
}

static void BlitJobFn(FrameBuffer::JobData const& _job)
{
	FrameBufferPlane const& plane = *_job.m_plane;
	uint32_t* fb32 = (uint32_t*)_job.m_linearPixels;
	for (uint32_t tileY = 0; tileY < plane.m_tilesY; ++tileY)
	{
		for (uint32_t tileX = 0; tileX < plane.m_tilesX; ++tileX)
		{
			ColourTile const& tile = plane.m_colourTiles[tileY * plane.m_tilesX + tileX];

			uint32_t const yEnd = kt::Min(Config::c_binHeight, plane.m_height - tileY * Config::c_binHeight);
			uint32_t const widthCopySize = kt::Min(Config::c_binWidth, plane.m_width - tileX * Config::c_binWidth);

			for (uint32_t y = 0; y < yEnd; ++y)
			{
				uint8_t const* src = &tile.m_colour[y * 4 * Config::c_binWidth];

				uint32_t* dest = fb32 + tileY * Config::c_binHeight * plane.m_width + tileX * Config::c_binWidth + y * plane.m_width;
				memcpy(dest, src, 4 * widthCopySize);
			}
		}
	}

	if (_job.m_onFinishBlit)
	{
		_job.m_onFinishBlit(_job.m_onFinishBlitUser);
	}
}


void RenderContext::Blit(FrameBuffer& _fb, uint8_t* _linearPixels, void(*_onFinishBlit)(void*), void* _onFinishUser)
{
	uint32_t const idx = _fb.m_writePlane;
	std::atomic_store_explicit(&_fb.m_jobs[idx].m_counter, 0, std::memory_order_relaxed);
	_fb.m_jobs[idx].m_linearPixels = _linearPixels;
	_fb.m_jobs[idx].m_plane = _fb.WritePlane();
	_fb.m_jobs[idx].m_onFinishBlit = _onFinishBlit;
	_fb.m_jobs[idx].m_onFinishBlitUser = _onFinishUser;

	auto taskFN = [](Task const* _task, uint32_t _threadIdx, uint32_t _start, uint32_t _end)
	{
		BlitJobFn(*(FrameBuffer::JobData*)_task->m_userData);
	};

	_fb.m_jobs[idx].m_task.Set(taskFN, 1, 1, &_fb.m_jobs[idx], &_fb.m_jobs[idx].m_counter);

	_fb.SwapPlanes();


	m_taskSystem.PushTask(&_fb.m_jobs[idx].m_task);

	// ensure we can start writing to the next frame buffer.
	m_taskSystem.WaitForCounter(&_fb.m_jobs[_fb.m_writePlane].m_counter);
}

}
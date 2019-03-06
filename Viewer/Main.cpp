#include "Platform/Window_Win32.h"
#include "Rasterizer.h"
#include <kt/Timer.h>
#include <kt/Logging.h>
#include <kt/Vec3.h>
#include <kt/Mat4.h>
#include <vector>

#include "Camera.h"
#include "Obj.h"

#include "Input.h"
#include "kt/FilePath.h"
#include "Renderer.h"


template <typename T>
static void DrawT(T const* _indexBuffer, uint32_t const _numIdx, sr::Obj::Vertex const* _vtxBuffer, sr::Raster::FrameBuffer& _fb, sr::Raster::DepthBuffer& _db, kt::Mat4 const& _mtx)
{
	for (uint32_t i = 0; i < _numIdx; i += 3)
	{
		T const i0 = _indexBuffer[i];
		T const i1 = _indexBuffer[i + 1];
		T const i2 = _indexBuffer[i + 2];

		sr::Raster::SetupAndRasterTriTest(_fb, _db, _mtx, _vtxBuffer[i0].pos, _vtxBuffer[i1].pos, _vtxBuffer[i2].pos);
	}
}

struct UniformTest
{
	sr::Tex::TextureData const* diffuse;
};

void DiffuseTest(void const* _uniforms, __m256 const _varyings[sr::Config::c_maxVaryings], float o_colour[4 * 8], __m256 const& _execMask)
{
	sr::Tex::TextureData* tex = (sr::Tex::TextureData*)_uniforms;
	
	if (!tex || !tex->m_data)
	{
		memset(o_colour, 0, 4 * 8 * sizeof(float));
		return;
	}

	// todo should pass in float ptr?
	KT_ALIGNAS(32) float u[8];
	KT_ALIGNAS(32) float v[8];
	_mm256_store_ps(u, _varyings[offsetof(sr::Obj::Vertex, uv) / sizeof(float)]);
	_mm256_store_ps(v, _varyings[offsetof(sr::Obj::Vertex, uv) / sizeof(float) + 1]);

	for (uint32_t i = 0; i < 8; ++i)
	{
		sr::Tex::SampleWrap_Slow(*tex, u[i], v[i], &o_colour[i * 4]);
	}
}

void NormalShaderTest(void const* _uniforms, __m256 const _varyings[sr::Config::c_maxVaryings], float o_colour[4 * 8], __m256 const& _execMask)
{

	// todo should pass in float ptr?

	uint32_t const redOffset = offsetof(sr::Obj::Vertex, norm) / sizeof(float);
	uint32_t const greenOffset = redOffset + 1;
	uint32_t const blueOffset = greenOffset + 1;

	__m256 const half = _mm256_set1_ps(0.5f);

	__m256 const red = _mm256_fmadd_ps(_varyings[redOffset], half, half);
	__m256 const green = _mm256_fmadd_ps(_varyings[greenOffset], half, half);
	__m256 const blue = _mm256_fmadd_ps(_varyings[blueOffset], half, half);

	KT_ALIGNAS(32) float red_store[8];
	KT_ALIGNAS(32) float green_store[8];
	KT_ALIGNAS(32) float blue_store[8];
	_mm256_store_ps(red_store, red);
	_mm256_store_ps(green_store, green);
	_mm256_store_ps(blue_store, blue);

	for (uint32_t i = 0; i < 8; ++i)
	{
		o_colour[i * 4 + 0] = red_store[i];
		o_colour[i * 4 + 1] = green_store[i];
		o_colour[i * 4 + 2] = blue_store[i];
		o_colour[i * 4 + 3] = 0xFF;
	}
}

int main(int argc, char** argv)
{
	sr::input::Init();
	sr::Window_Win32 window("SoftRast", 1280, 720);

	kt::TimePoint prevFrameTime = kt::TimePoint::Now();
	kt::Duration totalTime = kt::Duration::Zero();

	sr::Raster::FrameBuffer fb;
	fb.height = window.Height();
	fb.width = window.Width();
	fb.ptr = window.BackBufferData();

	sr::Raster::DepthBuffer depthBuff;
	depthBuff.Init(kt::GetDefaultAllocator(), 1280, 720);

	kt::Mat4 perspMtx = kt::Mat4::PerspectiveLH_ZO(kt::ToRadians(85.0f), 1280.0f / 720.0f, 0.01f, 1000.0f);

	sr::FreeCamController controller;
	
	sr::FreeCamController::ProjectionParams proj;
	proj.m_aspect = 1280.0f / 720.0f;
	proj.m_fov = kt::ToRadians(85.0f);
	proj.m_nearPlane = 0.1f;
	proj.m_farPlane = 1000.0f;

	controller.SetProjectionParams(proj);
	controller.SetPos({ 0.0f, 0.0f, -2.0f });

	uint32_t logDtCounter = 0;
	
	sr::Obj::Model model;
	//model.Load("Models/dragon.obj", kt::GetDefaultAllocator(), sr::Obj::LoadFlags::FlipWinding);
	//model.Load("Models/bunny.obj", kt::GetDefaultAllocator(), sr::Obj::LoadFlags::FlipWinding);
	model.Load("Models/sponza/sponza.obj", kt::GetDefaultAllocator(), sr::Obj::LoadFlags::FlipWinding | sr::Obj::LoadFlags::FlipUVs);
	//model.Load("Models/teapot/teapot.obj", kt::GetDefaultAllocator(), sr::Obj::LoadFlags::FlipWinding);

	kt::FilePath const f = kt::FilePath::WorkingDirectory();

	kt::Duration frameTime = kt::Duration::FromMicroseconds(16.0);

	sr::RenderContext renderCtx;
	
	sr::FrameBuffer framebuffer;
	framebuffer.Init(1280, 720);

	while (!window.WantsQuit())
	{
		window.PumpMessageLoop();
		sr::input::Tick((float)frameTime.Seconds());
		controller.UpdateViewGamepad((float)frameTime.Seconds());
#if 0
		uint8_t const c =  (uint8_t)(255 * fmod(totalTime.Seconds(), 1.0));
#else
		uint8_t const c = 0;
#endif
		uint8_t const col[3] = { c, c, c };
		sr::Raster::FillScreenTest(fb, col);
		sr::Raster::ClearDepthBufferTest(depthBuff);

		renderCtx.BeginFrame();
		renderCtx.ClearFrameBuffer(framebuffer, 0);

		for (sr::Obj::Mesh const& mesh : model.m_meshes)
		{
			sr::DrawCall call;
			call.m_frameBuffer = &framebuffer;
			call.m_mvp = controller.GetCam().GetCachedViewProj();
			call.m_attributeBuffer.m_ptr = (uint8_t*)mesh.m_vertexData;
			call.m_attributeBuffer.m_stride = sizeof(sr::Obj::Vertex);
			call.m_attributeBuffer.m_num = mesh.m_numVertices;

			call.m_positionBuffer.m_ptr = (uint8_t*)mesh.m_vertexData;
			call.m_positionBuffer.m_stride = sizeof(sr::Obj::Vertex);
			call.m_positionBuffer.m_num = mesh.m_numVertices;

			call.m_indexBuffer.m_ptr = mesh.m_indexData.index16;
			//call.m_indexBuffer.m_ptr = mesh.m_indexData.index16 + 18;
			//call.m_indexBuffer.m_num = 3;
			call.m_indexBuffer.m_num = mesh.m_numIndicies;
			call.m_indexBuffer.m_stride = mesh.m_indexType == sr::IndexType::u16 ? sizeof(uint16_t) : sizeof(uint32_t);
			
		
			if (mesh.m_matIdx < model.m_materials.Size())
			{
				call.m_pixelUniforms = &model.m_materials[mesh.m_matIdx].m_diffuse;
				call.m_pixelShader = DiffuseTest;
			}
			else
			{
				call.m_pixelShader = NormalShaderTest;
			}



			//sr::Raster::DrawSerial_Test(fb, depthBuff, controller.GetCam().GetCachedViewProj(), call);

			renderCtx.DrawIndexed(call);


		}

		renderCtx.EndFrame();

		static bool blitDepth = false;

		if (blitDepth)
		{
			framebuffer.BlitDepth(window.BackBufferData());
		}
		else
		{
			framebuffer.Blit(window.BackBufferData());
		}


		window.Flip();

		kt::TimePoint const timeNow = kt::TimePoint::Now();
		frameTime = timeNow - prevFrameTime;
		prevFrameTime = timeNow;
		totalTime += frameTime;

		if (++logDtCounter % 10 == 0)
		{
			KT_LOG_INFO("Frame took: %.3fms", frameTime.Milliseconds());
		}
	}

	sr::input::Shutdown();
	return 1;
}
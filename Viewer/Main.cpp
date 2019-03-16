#include <kt/Timer.h>
#include <kt/Logging.h>
#include <kt/Vec3.h>
#include <kt/Mat4.h>

#include "Camera.h"
#include "Obj.h"
#include "Input.h"
#include "Renderer.h"
#include "Platform/Window_Win32.h"
#include "Rasterizer.h"
#include "Config.h"

void DiffuseTest(void const* _uniforms, float const* _varyings, float o_colour[4 * 8], __m256 const& _execMask)
{
	sr::Tex::TextureData* tex = (sr::Tex::TextureData*)_uniforms;
	
	if (!tex || tex->m_texels.Size() == 0)
	{
		memset(o_colour, 0, 4 * 8 * sizeof(float));
		return;
	}

	uint32_t const stride = (sizeof(sr::Obj::Vertex) / sizeof(float)) + 4;

	KT_ALIGNAS(32) float u[8];
	KT_ALIGNAS(32) float v[8];

	KT_ALIGNAS(32) float dudx[8];
	KT_ALIGNAS(32) float dudy[8];
	KT_ALIGNAS(32) float dvdx[8];
	KT_ALIGNAS(32) float dvdy[8];

	uint32_t const uOffs = offsetof(sr::Obj::Vertex, uv) / sizeof(float);
	uint32_t const vOffs = 1 + offsetof(sr::Obj::Vertex, uv) / sizeof(float);

	for (uint32_t i = 0; i < 8; ++i)
	{
		u[i] = _varyings[4 + i * stride + uOffs];
		v[i] = _varyings[4 + i * stride + vOffs];

		dudx[i] = _varyings[i * stride];
		dudy[i] = _varyings[i * stride + 1];
		dvdx[i] = _varyings[i * stride + 2];
		dvdy[i] = _varyings[i * stride + 3];
	}



#if 0
	for (uint32_t i = 0; i < 8; ++i)
	{
		sr::Tex::SampleWrap_Slow(*tex, 0, u[i], v[i], &o_colour[i * 4]);
	}
#else
	sr::Tex::SampleWrap(*tex, _mm256_load_ps(u), _mm256_load_ps(v), _mm256_load_ps(dudx), _mm256_load_ps(dudy), _mm256_load_ps(dvdx), _mm256_load_ps(dvdy), o_colour);
#endif
}

void NormalShaderTest(void const* _uniforms, float const* _varyings, float o_colour[4 * 8], __m256 const& _execMask)
{
	uint32_t const stride = (sizeof(sr::Obj::Vertex) / sizeof(float)) + 4;
	uint32_t const redOffset = offsetof(sr::Obj::Vertex, norm) / sizeof(float) + 4;
	uint32_t const greenOffset = redOffset + 1;
	uint32_t const blueOffset = greenOffset + 1;

#define WHITE_OUT (0)

	for (uint32_t i = 0; i < 8; ++i)
	{
#if WHITE_OUT
		o_colour[i * 4 + 0] = 1.0f;
		o_colour[i * 4 + 1] = 1.0f;
		o_colour[i * 4 + 2] = 1.0f;
#else
		float const r = _varyings[i * stride + redOffset];
		float const g = _varyings[i * stride + greenOffset];
		float const b = _varyings[i * stride + blueOffset];

		float const mul = 1.0f / sqrtf(r * r + g * g + b * b);

		o_colour[i * 4 + 0] = r * mul * 0.5f + 0.5f;
		o_colour[i * 4 + 1] = g * mul * 0.5f + 0.5f;
		o_colour[i * 4 + 2] = b * mul * 0.5f + 0.5f;
#endif
		o_colour[i * 4 + 3] = 1.0f;
	}
}

int main(int argc, char** argv)
{
	sr::input::Init();
	sr::Window_Win32 window("SoftRast", sr::Config::c_screenWidth, sr::Config::c_screenHeight);

	kt::TimePoint prevFrameTime = kt::TimePoint::Now();
	kt::Duration totalTime = kt::Duration::Zero();

	sr::FreeCamController controller;
	
	sr::FreeCamController::ProjectionParams proj;
	proj.m_aspect = float(sr::Config::c_screenWidth) / float(sr::Config::c_screenHeight);
	proj.m_fov = kt::ToRadians(85.0f);
	proj.m_nearPlane = 0.1f;
	proj.m_farPlane = 10000.0f;

	controller.SetProjectionParams(proj);
	controller.SetPos({ 0.0f, 0.0f, -2.0f });

	uint32_t logDtCounter = 0;
	
	sr::Obj::Model model;
	//model.Load("Models/dragon.obj", kt::GetDefaultAllocator(), sr::Obj::LoadFlags::FlipWinding);
	//model.Load("Models/bunny.obj", kt::GetDefaultAllocator(), sr::Obj::LoadFlags::FlipWinding);
	//model.Load("Models/cube/cube.obj", kt::GetDefaultAllocator(), sr::Obj::LoadFlags::FlipWinding);
	//model.Load("Models/sponza/sponza.obj", kt::GetDefaultAllocator(), sr::Obj::LoadFlags::FlipWinding | sr::Obj::LoadFlags::FlipUVs);
	model.Load("Models/sponza-crytek/sponza.obj", kt::GetDefaultAllocator(), sr::Obj::LoadFlags::FlipWinding | sr::Obj::LoadFlags::FlipUVs);
	//model.Load("Models/teapot/teapot.obj", kt::GetDefaultAllocator(), sr::Obj::LoadFlags::FlipWinding);

	kt::Duration frameTime = kt::Duration::FromMicroseconds(16.0);

	sr::RenderContext renderCtx;
	
	sr::FrameBuffer framebuffer;
	framebuffer.Init(sr::Config::c_screenWidth, sr::Config::c_screenHeight);

	while (!window.WantsQuit())
	{
		window.PumpMessageLoop();
		sr::input::Tick((float)frameTime.Seconds());
		controller.UpdateViewGamepad((float)frameTime.Seconds());

		renderCtx.BeginFrame();
		renderCtx.ClearFrameBuffer(framebuffer, 0);

		for (sr::Obj::Mesh const& mesh : model.m_meshes)
		{
			sr::DrawCall call;
			call.m_frameBuffer = &framebuffer;
			call.m_mvp = controller.GetCam().GetCachedViewProj();

			call.SetAttributeBuffer(mesh.m_vertexData.Data(), sizeof(sr::Obj::Vertex), mesh.m_vertexData.Size(), offsetof(sr::Obj::Vertex, uv) / sizeof(float));

			call.m_positionBuffer.m_ptr = (uint8_t*)mesh.m_vertexData.Data();
			call.m_positionBuffer.m_stride = sizeof(sr::Obj::Vertex);
			call.m_positionBuffer.m_num = mesh.m_vertexData.Size();

			call.m_indexBuffer.m_ptr = mesh.m_indexData.Data();
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
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

void DiffuseTest(void const* _uniforms, float const* _varyings, uint32_t o_texels[8], uint32_t _execMask)
{
	sr::Tex::TextureData* tex = (sr::Tex::TextureData*)_uniforms;
	
	if (!tex || tex->m_texels.Size() == 0)
	{
		memset(o_texels, 0xFFFFFFFF, 8 * sizeof(uint32_t));
		return;
	}

	uint32_t const stride = (sizeof(sr::Obj::Vertex) / sizeof(float)) + 4;

	__m256 dudx;
	__m256 dudy;
	__m256 dvdx;
	__m256 dvdy;

	// load as rows and transpose sub 4x4 matricies.
	// [dudx0 dudy0 dvdx0 dvdy0 | dudx4 dudy4 dvdx4 dvdy4] 
	// [dudx1 dudy1 dvdx1 dvdy1 | dudx5 dudy5 dvdx5 dvdy5] 
	// [dudx2 dudy2 dvdx2 dvdy2 | dudx6 dudy6 dvdx6 dvdy6] 
	// [dudx3 dudy3 dvdx3 dvdy3 | dudx7 dudy7 dvdx7 dvdy7] 

	dudx = _mm256_loadu2_m128(_varyings + stride * 4, _varyings + stride * 0);
	dudy = _mm256_loadu2_m128(_varyings + stride * 5, _varyings + stride * 1);
	dvdx = _mm256_loadu2_m128(_varyings + stride * 6, _varyings + stride * 2);
	dvdy = _mm256_loadu2_m128(_varyings + stride * 7, _varyings + stride * 3);

	sr::simdutil::Transpose4x4SubMatricies(dudx, dudy, dvdx, dvdy);

	/*
	struct Vertex
	{
		kt::Vec3 pos;
		kt::Vec3 norm;
		kt::Vec2 uv;
	};
	*/

	__m256 pos_x, pos_y, pos_z;
	__m256 norm_x, norm_y, norm_z;
	__m256 u, v;

	float const* vertexVaryingBegin = _varyings + 4;
	pos_x	= _mm256_loadu_ps(vertexVaryingBegin);
	pos_y	= _mm256_loadu_ps(vertexVaryingBegin + stride);
	pos_z	= _mm256_loadu_ps(vertexVaryingBegin + stride * 2);
	norm_x	= _mm256_loadu_ps(vertexVaryingBegin + stride * 3);
	norm_y	= _mm256_loadu_ps(vertexVaryingBegin + stride * 4);
	norm_z	= _mm256_loadu_ps(vertexVaryingBegin + stride * 5);
	u		= _mm256_loadu_ps(vertexVaryingBegin + stride * 6);
	v		= _mm256_loadu_ps(vertexVaryingBegin + stride * 7);

	sr::simdutil::Transpose8x8(pos_x, pos_y, pos_z, norm_x, norm_y, norm_z, u, v);

	__m256 r;
	__m256 g;
	__m256 b;
	__m256 a;

	sr::Tex::SampleWrap(*tex, u, v, dudx, dudy, dvdx, dvdy, r, g, b, a, _execMask);
	sr::simdutil::RGBA32SoA_To_RGBA8AoS(r, g, b, a, o_texels);
}

void NormalShaderTest(void const* _uniforms, float const* _varyings, uint32_t o_texels[8], uint32_t _execMask)
{
	//uint32_t const stride = (sizeof(sr::Obj::Vertex) / sizeof(float)) + 4;
	//uint32_t const redOffset = offsetof(sr::Obj::Vertex, norm) / sizeof(float) + 4;
	//uint32_t const greenOffset = redOffset + 1;
	//uint32_t const blueOffset = greenOffset + 1;

#define WHITE_OUT (0)

	for (uint32_t i = 0; i < 8; ++i)
	{
//#if WHITE_OUT
//		o_colour[i * 4 + 0] = 1.0f;
//		o_colour[i * 4 + 1] = 1.0f;
//		o_colour[i * 4 + 2] = 1.0f;
//#else
//		float const r = _varyings[i * stride + redOffset];
//		float const g = _varyings[i * stride + greenOffset];
//		float const b = _varyings[i * stride + blueOffset];
//
//		float const mul = 1.0f / sqrtf(r * r + g * g + b * b);
//
//		o_colour[i * 4 + 0] = r * mul * 0.5f + 0.5f;
//		o_colour[i * 4 + 1] = g * mul * 0.5f + 0.5f;
//		o_colour[i * 4 + 2] = b * mul * 0.5f + 0.5f;
//#endif
//		o_colour[i * 4 + 3] = 1.0f;
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

#if SR_USE_REVERSE_Z
	proj.m_nearPlane = 10000.0f;
	proj.m_farPlane = 0.1f;
#else	
	proj.m_nearPlane = 0.1f;
	proj.m_farPlane = 10000.0f;
#endif

	controller.SetProjectionParams(proj);
	controller.SetPos({ 0.0f, 0.0f, -2.0f });

	uint32_t logDtCounter = 0;
	
	sr::Obj::Model model;
	//model.Load("Models/dragon.obj", kt::GetDefaultAllocator(), sr::Obj::LoadFlags::FlipWinding);
	//model.Load("Models/bunny.obj", kt::GetDefaultAllocator(), sr::Obj::LoadFlags::FlipWinding);
	//model.Load("Models/cube/cube.obj", kt::GetDefaultAllocator(), sr::Obj::LoadFlags::FlipWinding);
	model.Load("Models/sponza-crytek/sponza.obj", kt::GetDefaultAllocator(), sr::Obj::LoadFlags::FlipWinding | sr::Obj::LoadFlags::FlipUVs);
	//model.Load("Models/teapot/teapot.obj", kt::GetDefaultAllocator(), sr::Obj::LoadFlags::FlipWinding);

	kt::Duration frameTime = kt::Duration::FromMilliseconds(16.0);

	sr::RenderContext renderCtx;
	
	sr::FrameBuffer framebuffer(sr::Config::c_screenWidth, sr::Config::c_screenHeight);

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
			call.SetFrameBuffer(&framebuffer);
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

		//static bool blitDepth = false;

		//if (blitDepth)
		//{
		//	framebuffer.BlitDepth(window.BackBufferData());
		//}
		//else
		{
			renderCtx.Blit(framebuffer, window.BackBufferData(), [](void* _ptr) { ((sr::Window_Win32*)_ptr)->Flip(); }, &window);
		}


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
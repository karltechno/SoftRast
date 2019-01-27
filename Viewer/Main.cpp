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

void PixelTest(void const* _uniforms, __m256 const _varyings[sr::Config::c_maxVaryings], float o_colour[4 * 8], __m256 const& _execMask)
{
	UniformTest* uniform = (UniformTest*)_uniforms;
	memset(o_colour, 0, 4 * 8 * sizeof(float));

	// todo should pass in float ptr?
	float u[8];
	float v[8];
	_mm256_store_ps(u, _varyings[0]);
	_mm256_store_ps(v, _varyings[1]);

	for (uint32_t i = 0; i < 8; ++i)
	{
		sr::Tex::SampleClamp_Slow(*uniform->diffuse, u[i], v[i], &o_colour[i * 4]);
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
	proj.m_farPlane = 10000.0f;

	controller.SetProjectionParams(proj);
	controller.SetPos({ 0.0f, 0.0f, -2.0f });


	kt::Vec3 eye = { 75.0f, 75.0f, 75.0f };
	kt::Mat4 fullMtx = perspMtx * kt::Mat4::LookAtLH(eye * 10.0f, -kt::Normalize(eye));
	uint32_t logDtCounter = 0;
	
	sr::Obj::Model model;
	model.Load("Models/cube/cube.obj", kt::GetDefaultAllocator(), sr::Obj::LoadFlags::FlipWinding);

	kt::FilePath const f = kt::FilePath::WorkingDirectory();

	kt::Duration frameTime = kt::Duration::FromMicroseconds(16.0);

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


		for (sr::Obj::Mesh const& mesh : model.m_meshes)
		{
			sr::Renderer::DrawCall call;
			call.m_attributeBuffer.m_ptr = (uint8_t*)mesh.m_vertexData + offsetof(sr::Obj::Vertex, uv);
			call.m_attributeBuffer.m_stride = sizeof(sr::Obj::Vertex);
			call.m_attributeBuffer.m_num = mesh.m_numVertices;

			call.m_positionBuffer.m_ptr = (uint8_t*)mesh.m_vertexData;
			call.m_positionBuffer.m_stride = sizeof(sr::Obj::Vertex);
			call.m_positionBuffer.m_num = mesh.m_numVertices;

			call.m_indexBuffer.m_ptr = mesh.m_indexData.index16;
			call.m_indexBuffer.m_num = mesh.m_numIndicies;
			call.m_indexBuffer.m_stride = mesh.m_indexType == sr::IndexType::u16 ? sizeof(uint16_t) : sizeof(uint32_t);
			
			call.m_pixelShader = PixelTest;
		
			UniformTest uniforms;
			if (model.m_materials.Size())
			{
				uniforms.diffuse = &model.m_materials[0].m_diffuse;
			}
			call.m_pixelUniforms = &uniforms;
			sr::Raster::DrawSerial_Test(fb, depthBuff, controller.GetCam().GetCachedViewProj(), call);

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
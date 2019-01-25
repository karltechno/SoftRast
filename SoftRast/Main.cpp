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

	kt::Duration frameTime = kt::Duration::FromMicroseconds(16);

	while (!window.WantsQuit())
	{
		window.PumpMessageLoop();
		sr::input::Tick(frameTime.Seconds());
		controller.UpdateViewGamepad(frameTime.Seconds());
#if 0
		uint8_t const c =  (uint8_t)(255 * fmod(totalTime.Seconds(), 1.0));
#else
		uint8_t const c = 0;
#endif
		uint8_t const col[3] = { c, c, c };
		sr::Raster::FillScreenTest(fb, col);
		sr::Raster::ClearDepthBufferTest(depthBuff);

#if 0
		float myZ = (sinf(totalTime.Seconds() / 2.0f) + 1.0f) * 10.0f + 1.0f;
		//myZ = 1.0f;
		//Raster::RasterTriTest(fb, depthBuff, perspMtx, kt::Vec3(-0.5f, -50.0f, myZ), kt::Vec3(10.0f, 0.0f, 5.0f), kt::Vec3(0.0f, 0.5f, 5.0f));
		kt::Mat3 rotMtx = kt::Mat3::RotX(totalTime.Seconds());

		//Raster::RasterTriTest(fb, depthBuff, controller.GetCam().GetCachedViewProj(), kt::Vec3(-0.5f, -0.5f, myZ), kt::Vec3(0.5f, 0.0f, 1.0f), kt::Vec3(0.0f, 0.5f, myZ));


		//Raster::RasterTriTest(fb, depthBuff, perspMtx, kt::Vec3(-0.75f, 0.0f, 10.0f), kt::Vec3(1.0f, 0.0f, 10.0f), kt::Vec3(0.0f, 0.25f, 10.0f));
		//Raster::RasterTriTest(fb, depthBuff, perspMtx, kt::Vec3(-0.5f, 0.0f, 15.0f), kt::Vec3(10.0f, 0.0f, 15.0f), kt::Vec3(0.0f, 0.5f, 15.0f));
#else
		for (sr::Obj::Mesh const& mesh : model.m_meshes)
		{
			if (mesh.m_indexType == sr::IndexType::u16)
			{
				DrawT(mesh.m_indexData.index16, mesh.m_numIndicies, mesh.m_vertexData, fb, depthBuff, controller.GetCam().GetCachedViewProj());
			}
			else
			{
				DrawT(mesh.m_indexData.index32, mesh.m_numIndicies, mesh.m_vertexData, fb, depthBuff, controller.GetCam().GetCachedViewProj());
			}
		}

		//sr::Raster::SetupAndRasterTriTest(fb, depthBuff, controller.GetCam().GetCachedViewProj(), kt::Vec3(-1.0f, -1.0f, 0.0f), kt::Vec3(1.0f, -1.0f, 0.0f), kt::Vec3(0.0f, 1.0f, 0.0f));
#endif
		
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
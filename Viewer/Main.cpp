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
#include "Shaders.h"
#include "Scene.h"
#include "SponzaScene.h"

#include "microprofile.h"

int main(int argc, char** argv)
{
	sr::input::Init();
	sr::Window_Win32 window("SoftRast", sr::Config::c_screenWidth, sr::Config::c_screenHeight);

	kt::TimePoint prevFrameTime = kt::TimePoint::Now();
	kt::Duration totalTime = kt::Duration::Zero();

	uint32_t logDtCounter = 0;
	sr::Scene* scene = nullptr;
	sr::Obj::Model model;
	//scene = new sr::SimpleModelScene("Models/dragon.obj", sr::Obj::LoadFlags::FlipWinding);
	//scene = new sr::SimpleModelScene("Models/bunny.obj", sr::Obj::LoadFlags::FlipWinding);
	//scene = new sr::SimpleModelScene("Models/cube/cube.obj", sr::Obj::LoadFlags::FlipWinding);
	scene = new sr::SponzaScene("Models/sponza-crytek/sponza.obj", sr::Obj::LoadFlags::FlipWinding | sr::Obj::LoadFlags::FlipUVs);
	//scene = new sr::SimpleModelScene("Models/teapot/teapot.obj", sr::Obj::LoadFlags::FlipWinding);
	//scene = new sr::SimpleModelScene("Models/lost_empire/lost_empire.obj", sr::Obj::LoadFlags::FlipWinding | sr::Obj::LoadFlags::FlipUVs);

	scene->Init(sr::Config::c_screenWidth, sr::Config::c_screenHeight);

	kt::Duration frameTime = kt::Duration::FromMilliseconds(16.0);

	sr::RenderContext renderCtx;
	sr::FrameBuffer framebuffer(sr::Config::c_screenWidth, sr::Config::c_screenHeight);

	MicroProfileOnThreadCreate("Main");
	MicroProfileSetEnableAllGroups(true);
	MicroProfileSetForceMetaCounters(true);
	MicroProfileStartContextSwitchTrace();
	

	while (!window.WantsQuit())
	{
		window.PumpMessageLoop();
		sr::input::Tick((float)frameTime.Seconds());

		renderCtx.BeginFrame();
	
		scene->Update(renderCtx, framebuffer, frameTime.Seconds());

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

		MicroProfileFlip(0);
	}

	delete scene;
	sr::input::Shutdown();

	MicroProfileShutdown();

	
}
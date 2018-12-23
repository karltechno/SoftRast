#include "Platform/Window_Win32.h"
#include "Rasterizer.h"
#include <kt/Timer.h>
#include <kt/Logging.h>

int main(int argc, char** argv)
{
	Window_Win32 window("SoftRast", 1280, 720);

	kt::TimePoint prevFrameTime = kt::TimePoint::Now();
	kt::Duration totalTime = kt::Duration::Zero();

	Framebuffer fb;
	fb.height = window.Height();
	fb.width = window.Width();
	fb.ptr = window.BackBufferData();

	uint32_t logDtCounter = 0;

	while (!window.WantsQuit())
	{
		window.PumpMessageLoop();
#if 0
		uint8_t const c =  (uint8_t)(255 * fmod(totalTime.Seconds(), 1.0));
#else
		uint8_t const c = 0;
#endif
		uint8_t const col[3] = { c, c, c };
		Raster::FillScreenTest(fb, col);

#if 1
		Raster::RasterTriTest(fb, kt::Vec3(-0.5f, 0.0f, 0.0f), kt::Vec3(0.0f, 0.5f, 0.0f), kt::Vec3(0.5f, 0.0f, 0.0f));
		Raster::RasterTriTest(fb, kt::Vec3(-0.75f, 0.0f, 0.0f), kt::Vec3(0.0f, 0.25f, 0.0f), kt::Vec3(1.0f, 0.0f, 0.0f));
#else
		Raster::RasterTriTest(fb, kt::Vec3(-1.0f, -1.0f, 0.0f), kt::Vec3(0.0f, 1.0f, 0.0f), kt::Vec3(1.0f, -1.0f, 0.0f));
#endif
		
		window.Flip();

		kt::TimePoint const timeNow = kt::TimePoint::Now();
		kt::Duration const frameTime = timeNow - prevFrameTime;
		prevFrameTime = timeNow;
		totalTime += frameTime;

		if (++logDtCounter % 10 == 0)
		{
			KT_LOG_INFO("Frame took: %.3fms", frameTime.Milliseconds());
		}
	}

	return 1;
}
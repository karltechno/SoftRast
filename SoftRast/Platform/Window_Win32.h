#pragma once
#include <kt/kt.h>
#include <kt/Memory.h>

#define WIN32_LEAN_AND_MEAN
#include <Windows.h>

namespace sr
{

class Window_Win32
{
public:
	Window_Win32(char const* _windowName, uint32_t _width, uint32_t _height);
	~Window_Win32();

	void PumpMessageLoop();
	bool WantsQuit() const;

	void Shutdown();

	void Flip();
	uint8_t* BackBufferData();

	uint32_t Height() const;
	uint32_t Width() const;

	HWND Nwh() { return m_hwnd; }

private:
	bool InitWindow(char const* _name, uint32_t _height, uint32_t _width);
	bool InitBackBuffer(uint32_t _height, uint32_t _width);

	static LRESULT CALLBACK WndProcHook(HWND _hwnd, UINT _msg, WPARAM _wparam, LPARAM _lparam);

	HWND m_hwnd = 0;
	HDC m_windowHDC = 0;
	HDC m_backBufferHDC = 0;
	HBITMAP m_drawDibSection = 0;
	void* m_backBufferBitmapPtr = nullptr;

	uint32_t m_width;
	uint32_t m_height;

	ATOM m_windowAtom = 0;

	bool m_wantsQuit = false;
};

}
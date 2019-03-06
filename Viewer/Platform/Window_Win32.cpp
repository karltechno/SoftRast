#include "Window_Win32.h"
#include <kt/Logging.h>

#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <immintrin.h>

namespace sr
{


LRESULT CALLBACK Window_Win32::WndProcHook(HWND _hwnd, UINT _msg, WPARAM _wparam, LPARAM _lparam)
{
	switch (_msg)
	{
		case WM_NCCREATE:
		{
			CREATESTRUCT* create = (CREATESTRUCT*)_lparam;
			::SetWindowLongPtr(_hwnd, GWLP_USERDATA, (LPARAM)create->lpCreateParams);
			return TRUE;
		} break;

		case WM_CLOSE:
		{
			Window_Win32* pState = (Window_Win32*)::GetWindowLongPtr(_hwnd, GWLP_USERDATA);
			pState->m_wantsQuit = true;
			return TRUE;
		} break;

		//case WM_DESTROY:
		//{

		//	return FALSE;
		//} break;
	}
	return ::DefWindowProcA(_hwnd, _msg, _wparam, _lparam);
}

Window_Win32::Window_Win32(char const* _windowName, uint32_t _width, uint32_t _height)
	: m_width(_width)
	, m_height(_height)
{
	if (!InitWindow(_windowName, _height, _width))
	{
		return;
	}

	InitBackBuffer(_height, _width);
}

Window_Win32::~Window_Win32()
{
	Shutdown();
}

void Window_Win32::PumpMessageLoop()
{
	MSG m;
	while (::PeekMessage(&m, nullptr, 0, 0, PM_REMOVE) != 0)
	{
		::TranslateMessage(&m);
		::DispatchMessage(&m);
	}
}

bool Window_Win32::WantsQuit() const
{
	return m_wantsQuit;
}

void Window_Win32::Shutdown()
{
	if (m_hwnd)
	{
		::ReleaseDC(m_hwnd, m_windowHDC);
		::DestroyWindow(m_hwnd);
	}

	if (m_backBufferHDC)
	{
		::DeleteDC(m_backBufferHDC);
	}
}

void Window_Win32::Flip()
{
	// Swizzle back buffer rgba -> bgra
	uint8_t* ptr = (uint8_t*)m_backBufferBitmapPtr;
	uint8_t* pEnd = ptr + m_width * m_height * 4;

	__m256i const shufmask = _mm256_setr_epi8
	(
		2, 1, 0, 3, 6, 5, 4, 7, 10, 9, 8, 11, 14, 13, 12, 15,
		2, 1, 0, 3, 6, 5, 4, 7, 10, 9, 8, 11, 14, 13, 12, 15
	);

	while (pEnd - ptr >= 32)
	{
		__m256i const bgra_x8 = _mm256_loadu_si256((__m256i*)ptr);
		__m256i const rgba_x8 = _mm256_shuffle_epi8(bgra_x8, shufmask);
		_mm256_storeu_si256((__m256i*)ptr, rgba_x8);
		ptr += 32;
	}

	while (ptr != pEnd)
	{
		uint8_t p0 = ptr[0];
		ptr[0] = ptr[2];
		ptr[2] = p0;
		ptr += 4;
	}


	BOOL ok = ::BitBlt(m_windowHDC, 0, 0, m_width, m_height, m_backBufferHDC, 0, 0, SRCCOPY);
	KT_UNUSED(ok);
	KT_ASSERT(ok == TRUE);
}

uint8_t* Window_Win32::BackBufferData()
{
	return (uint8_t*)m_backBufferBitmapPtr;
}

uint32_t Window_Win32::Height() const
{
	return m_height;
}

uint32_t Window_Win32::Width() const
{
	return m_width;
}
	
bool Window_Win32::InitWindow(char const* _name, uint32_t _height, uint32_t _width)
{
	static char const* c_className = "SOFT_RAST";

	WNDCLASSEXA wndClass{};
	wndClass.cbSize = sizeof(WNDCLASSEXA);
	wndClass.style = CS_HREDRAW | CS_VREDRAW;
	wndClass.lpfnWndProc = &WndProcHook;
	wndClass.hInstance = ::GetModuleHandle(nullptr);
	wndClass.lpszClassName = c_className;
	m_windowAtom = RegisterClassExA(&wndClass);
	if (!m_windowAtom)
	{
		KT_LOG_ERROR("Failed to register window class.");
		return false;
	}

	RECT windowRectAdjusted;
	windowRectAdjusted.top = windowRectAdjusted.left = 0;
	windowRectAdjusted.bottom = _height;
	windowRectAdjusted.right = _width;
	::AdjustWindowRect(&windowRectAdjusted, WS_OVERLAPPEDWINDOW, false);
	int32_t const adjustedWidth = windowRectAdjusted.right - windowRectAdjusted.left;
	int32_t const adjustedHeight = windowRectAdjusted.bottom - windowRectAdjusted.top;
	
	m_hwnd = ::CreateWindowExA(0, c_className, _name, WS_OVERLAPPEDWINDOW, 0, 0, adjustedWidth, adjustedHeight, 0, 0, ::GetModuleHandle(nullptr), this);
	KT_ASSERT(m_hwnd);
	if (!m_hwnd)
	{
		KT_LOG_ERROR("Failed to create window, error: %#08x", GetLastError());
		return false;
	}

	::ShowWindow((HWND)m_hwnd, SW_SHOW);
	return true;
}

bool Window_Win32::InitBackBuffer(uint32_t _height, uint32_t _width)
{
	// Setup bitmap info structure
	BITMAPINFO bmi = {};

	bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	bmi.bmiHeader.biWidth = (LONG)_width;
	bmi.bmiHeader.biHeight = -(LONG)_height; // flip so 0,0 is top
	bmi.bmiHeader.biPlanes = 1;
	bmi.bmiHeader.biBitCount = 32;
	bmi.bmiHeader.biCompression = BI_RGB;
	
	m_windowHDC = ::GetDC(m_hwnd);
	m_backBufferHDC = ::CreateCompatibleDC(m_windowHDC);
	KT_ASSERT(m_backBufferHDC);

	m_drawDibSection = ::CreateDIBSection(m_windowHDC, &bmi, DIB_RGB_COLORS, &m_backBufferBitmapPtr, 0, 0);
	KT_ASSERT(m_drawDibSection);
	::SelectObject(m_backBufferHDC, m_drawDibSection);

	Flip();
	return true;
}

}
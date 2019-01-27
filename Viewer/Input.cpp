#include <kt/kt.h>
#include <kt/Memory.h>

#include "Input.h"
#include "Platform/Window_Win32.h"

#if KT_PLATFORM_WINDOWS
#include "Platform/Input_Win32.h"
namespace sr { namespace input { using Impl = Impl_Windows; } }
#else
#	error Platform not implemented
#endif

namespace sr
{


namespace input
{

KT_ALIGNAS(KT_ALIGNOF(Impl)) static char s_implStorage[sizeof(Impl)];

static Impl* s_impl = nullptr;

void Init()
{
	if (s_impl)
	{
		return;
	}

	s_impl = kt::PlacementNew<Impl>((Impl*)s_implStorage);
}

void Shutdown()
{
	if (s_impl)
	{
		s_impl->~Impl();
		s_impl = nullptr;
	}
}

void Tick(float const _dt)
{
	KT_ASSERT(s_impl);
	s_impl->Tick(_dt);
}

bool IsDown(GamePadButton const _button, uint8_t const _controllerIdx)
{
	KT_ASSERT(_controllerIdx < s_maxGamepads);
	return s_impl->IsDown(_button, _controllerIdx);
}

bool WasPressed(GamePadButton const _button, uint8_t const _controllerIdx)
{
	KT_ASSERT(_controllerIdx < s_maxGamepads);
	return s_impl->WasPressed(_button, _controllerIdx);
}

bool WasReleased(GamePadButton const _button, uint8_t const _controllerIdx)
{
	KT_ASSERT(_controllerIdx < s_maxGamepads);
	return s_impl->WasReleased(_button, _controllerIdx);
}

float GetGamepadAxis(GamepadAxis const _axis, uint8_t const _controllerIdx)
{
	return s_impl->GetGamepadAxis(_axis, _controllerIdx);
}

bool IsGamepadConnected(uint8_t const _idx)
{
	KT_ASSERT(_idx < s_maxGamepads);
	return s_impl->IsGamepadConnected(_idx);
}


}

}
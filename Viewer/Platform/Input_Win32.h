#pragma once
#include "../Input.h"
#define WIN32_LEAN_AND_MEAN 
#include <Windows.h>
#include <kt/Array.h>

class App;

namespace sr
{


namespace input
{

struct Impl_Windows
{
	Impl_Windows();

	void Tick(float const _dt);

	bool IsDown(GamePadButton const _button, uint8_t const _controllerIdx);
	bool WasPressed(GamePadButton const _button, uint8_t const _controllerIdx);
	bool WasReleased(GamePadButton const _button, uint8_t const _controllerIdx);
	float GetGamepadAxis(GamepadAxis const _axis, uint8_t const _controllerIdx);
	bool IsGamepadConnected(uint8_t const _idx);

private:
	struct InternalXinputState
	{
		DWORD m_lastPacket = 0;
		float m_axisState[(uint32_t)GamepadAxis::MaxAxis] = {};
		uint8_t m_connected;

		uint8_t m_buttonState[(uint32_t)GamePadButton::NumButtons] = {};
	};

	void UpdateXInput();
	void HandlePadButtonState(InternalXinputState& _state, GamePadButton const _button, bool const _isDown);

	enum : uint8_t
	{
		c_isDownBit = 1 << 0,
		c_wasPressedBit = 1 << 1,
		c_wasReleasedBit = 1 << 2,
	};

	InternalXinputState m_gamepads[s_maxGamepads];
};

}

}
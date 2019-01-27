#include "Input_Win32.h"

#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <Xinput.h>
#include <WinUser.h>
#include <hidusage.h>
#include <Dbt.h>

#include <string.h>
#include <limits.h>
#include "../Input.h"


#pragma comment(lib, "Xinput.lib") // dll ?

namespace sr
{

namespace input
{

Impl_Windows::Impl_Windows()
{
	for (uint32_t i = 0; i < s_maxGamepads; ++i)
	{
		XINPUT_STATE state;
		m_gamepads[i].m_connected = (XInputGetState(i, &state) == ERROR_SUCCESS);
	}
}



void Impl_Windows::Tick(float const _dt)
{
	KT_UNUSED(_dt);

	UpdateXInput();
}

bool Impl_Windows::IsDown(GamePadButton const _button, uint8_t const _controllerIdx)
{
	return m_gamepads[_controllerIdx].m_buttonState[(uint32_t)_button] & c_isDownBit;
}

bool Impl_Windows::WasPressed(GamePadButton const _button, uint8_t const _controllerIdx)
{
	return m_gamepads[_controllerIdx].m_buttonState[(uint32_t)_button] & c_wasPressedBit;
}

bool Impl_Windows::WasReleased(GamePadButton const _button, uint8_t const _controllerIdx)
{
	return m_gamepads[_controllerIdx].m_buttonState[(uint32_t)_button] & c_wasReleasedBit;
}

float Impl_Windows::GetGamepadAxis(GamepadAxis const _axis, uint8_t const _controllerIdx)
{
	return m_gamepads[_controllerIdx].m_axisState[(uint32_t)_axis];
}

bool Impl_Windows::IsGamepadConnected(uint8_t const _idx)
{
	return m_gamepads[_idx].m_connected;
}


struct XinputBitmaskToGamepadButton
{
	uint32_t m_mask;
	GamePadButton m_button;
};

static XinputBitmaskToGamepadButton const s_xinputBitmaskToButton[] =
{
	{ XINPUT_GAMEPAD_DPAD_UP,			GamePadButton::DpadUp		},
	{ XINPUT_GAMEPAD_DPAD_DOWN,			GamePadButton::DpadDown		},
	{ XINPUT_GAMEPAD_DPAD_LEFT,			GamePadButton::DpadLeft		},
	{ XINPUT_GAMEPAD_DPAD_RIGHT,		GamePadButton::DpadRight	},
	{ XINPUT_GAMEPAD_START,				GamePadButton::Start		},
	{ XINPUT_GAMEPAD_BACK,				GamePadButton::Back			},
	{ XINPUT_GAMEPAD_LEFT_THUMB,		GamePadButton::LeftThumb	},
	{ XINPUT_GAMEPAD_RIGHT_THUMB,		GamePadButton::RightThumb	},
	{ XINPUT_GAMEPAD_LEFT_SHOULDER,		GamePadButton::LeftBumper	},
	{ XINPUT_GAMEPAD_RIGHT_SHOULDER,	GamePadButton::RightBumper	},
	{ XINPUT_GAMEPAD_A,					GamePadButton::Pad_A		},
	{ XINPUT_GAMEPAD_B,					GamePadButton::Pad_B		},
	{ XINPUT_GAMEPAD_X,					GamePadButton::Pad_X		},
	{ XINPUT_GAMEPAD_Y,					GamePadButton::Pad_Y		}
};

static float XinputTriggerToFloat(BYTE const _in, BYTE const _deadZone)
{
	BYTE const range = UINT8_MAX - _deadZone;
	BYTE const val = _in <= _deadZone ? 0 : (_in - _deadZone);
	return kt::Clamp(val / (float)range, 0.0f, 1.0f);
}

static float XinputThumbAxisToFloat(SHORT const _in, SHORT const _deadZone)
{
	if (_in <= _deadZone && _in >= -_deadZone)
	{
		return 0.0f;
	}

	SHORT const deadZoneApplied = _in < 0.0f ? _in + _deadZone : _in - _deadZone;
	return kt::Clamp(deadZoneApplied / (float)(SHRT_MAX - _deadZone), -1.0f, 1.0f);
}


void Impl_Windows::UpdateXInput()
{
	for (uint32_t padIdx = 0; padIdx < s_maxGamepads; ++padIdx)
	{
		InternalXinputState& pad = m_gamepads[padIdx];

		XINPUT_STATE xinputState;
		if (XInputGetState(padIdx, &xinputState) != ERROR_SUCCESS)
		{
			pad.m_connected = false;
			continue;
		}

		if (xinputState.dwPacketNumber == pad.m_lastPacket)
		{
			continue;
		}

		pad.m_lastPacket = xinputState.dwPacketNumber;

		for (uint32_t i = 0; i < KT_ARRAY_COUNT(s_xinputBitmaskToButton); ++i)
		{
			XinputBitmaskToGamepadButton const& mask = s_xinputBitmaskToButton[i];
			HandlePadButtonState(pad, mask.m_button, mask.m_mask & xinputState.Gamepad.wButtons);
		}

		pad.m_axisState[(uint32_t)GamepadAxis::LeftTrigger] = XinputTriggerToFloat(xinputState.Gamepad.bLeftTrigger, XINPUT_GAMEPAD_TRIGGER_THRESHOLD);
		pad.m_axisState[(uint32_t)GamepadAxis::RightTrigger] = XinputTriggerToFloat(xinputState.Gamepad.bRightTrigger, XINPUT_GAMEPAD_TRIGGER_THRESHOLD);
		pad.m_axisState[(uint32_t)GamepadAxis::LeftStick_X] = XinputThumbAxisToFloat(xinputState.Gamepad.sThumbLX, XINPUT_GAMEPAD_LEFT_THUMB_DEADZONE);
		pad.m_axisState[(uint32_t)GamepadAxis::LeftStick_Y] = XinputThumbAxisToFloat(xinputState.Gamepad.sThumbLY, XINPUT_GAMEPAD_LEFT_THUMB_DEADZONE);
		pad.m_axisState[(uint32_t)GamepadAxis::RightStick_X] = XinputThumbAxisToFloat(xinputState.Gamepad.sThumbRX, XINPUT_GAMEPAD_RIGHT_THUMB_DEADZONE);
		pad.m_axisState[(uint32_t)GamepadAxis::RightStick_Y] = XinputThumbAxisToFloat(xinputState.Gamepad.sThumbRY, XINPUT_GAMEPAD_RIGHT_THUMB_DEADZONE);
	}
}

void Impl_Windows::HandlePadButtonState(InternalXinputState& _state, GamePadButton const _button, bool const _isDown)
{
	uint8_t& internalState = _state.m_buttonState[(uint32_t)_button];
	if (_isDown)
	{
		if ((internalState & c_isDownBit) == 0)
		{
			internalState |= c_isDownBit;
			internalState |= c_wasPressedBit;
		}
	}
	else
	{
		if (internalState & c_isDownBit)
		{
			internalState ^= c_isDownBit;
			internalState |= c_wasReleasedBit;
		}
	}
}

}

}
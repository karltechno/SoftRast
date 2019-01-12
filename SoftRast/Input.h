#pragma once
#include <kt/kt.h>
#include <kt/Vec2.h>

class Window_Win32;

namespace sr
{


namespace input
{

static constexpr uint32_t s_maxGamepads = 4u;
enum class Key : uint8_t;
enum class GamepadAxis : uint8_t;
enum class GamePadButton : uint8_t;
enum class MouseButton : uint8_t;

void Init();

void Shutdown();

void Tick(float const _dt);

bool IsDown(GamePadButton const _key, uint8_t const _controllerIdx = 0);
bool WasPressed(GamePadButton const _key, uint8_t const _controllerIdx = 0);
bool WasReleased(GamePadButton const _key, uint8_t const _controllerIdx = 0);
float GetGamepadAxis(GamepadAxis const _axis, uint8_t const _controllerIdx = 0);
bool IsGamepadConnected(uint8_t const _idx);


enum class GamepadAxis : uint8_t
{
	LeftStick_X,
	LeftStick_Y,
	RightStick_X,
	RightStick_Y,

	LeftTrigger,
	RightTrigger,

	MaxAxis
};

enum class GamePadButton : uint8_t
{
	Pad_A,
	Pad_B,
	Pad_X,
	Pad_Y,
	LeftBumper,
	RightBumper,
	Back,
	Start,
	Guide,
	LeftThumb,
	RightThumb,
	DpadUp,
	DpadRight,
	DpadDown,
	DpadLeft,

	NumButtons
};

}

}
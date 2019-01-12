#include "Camera.h"
#include "Input.h"

namespace sr
{

void Camera::SetProjectionMatrix(kt::Mat4 const& _mat)
{
	m_viewToClip = _mat;
	m_cachedWorldToClipDirty = true;
}

void Camera::SetCameraPos(kt::Vec3 const& _pos)
{
	m_viewToWorld.SetPos(_pos);
	m_worldToView = kt::InverseOrthoAffine(m_viewToWorld);
	m_cachedWorldToClipDirty = true;
}

void Camera::SetCameraMatrix(kt::Mat4 const& _viewToWorld)
{
	m_viewToWorld = _viewToWorld;
	m_worldToView = kt::InverseOrthoAffine(_viewToWorld);
	m_cachedWorldToClipDirty = true;
}

kt::Mat4 const& Camera::GetCachedViewProj() const
{
	if (m_cachedWorldToClipDirty)
	{
		m_cachedWorldToClip = m_viewToClip * m_worldToView;
	}
	return m_cachedWorldToClip;
}

kt::Mat4 const& Camera::GetCameraMatrix() const
{
	return m_viewToWorld;
}

kt::Mat4 const& Camera::GetInverseCameraMatrix() const
{
	return m_worldToView;
}

kt::Mat4 const& Camera::GetProjectionMatrix() const
{
	return m_viewToClip;
}

void FreeCamController::SetPos(kt::Vec3 const& _pos)
{
	m_camera.SetCameraPos(_pos);
}

void FreeCamController::Move(kt::Vec3 const& _movement)
{
	m_frameMovement += _movement;
}

void FreeCamController::RotateXY(kt::Vec2 const& _xy)
{
	kt::Quat rotx;
	rotx.FromNormalizedAxisAngle(kt::Vec3(0.0f, 1.0f, 0.0f), _xy.x);
	kt::Quat roty;
	roty.FromNormalizedAxisAngle(kt::Vec3(1.0f, 0.0f, 0.0f), _xy.y);

	m_camQuat = m_camQuat * roty * rotx;
}

void FreeCamController::RotateByMatrix(kt::Mat3 const& _mtx)
{
	m_camQuat = kt::ToQuat(_mtx) * m_camQuat;
}

void FreeCamController::SetRotation(kt::Mat3 const& _rot)
{
	m_camQuat = kt::ToQuat(_rot);
}

void FreeCamController::UpdateView()
{
	m_camQuat = kt::Normalize(m_camQuat);
	kt::Mat4 camMtx = kt::ToMat4(m_camQuat);
	kt::Vec4 move = camMtx.m_cols[0] * m_frameMovement.x;
	move += camMtx.m_cols[1] * m_frameMovement.y;
	move += camMtx.m_cols[2] * m_frameMovement.z;

	move += kt::Vec4(m_camPos, 0.0f);
	move.w = 1.0f;

	m_camPos = kt::Vec3(move.x, move.y, move.z);

	camMtx.m_cols[3] = move;

	m_camera.SetCameraMatrix(camMtx);

	m_frameMovement = kt::Vec3(0.0f);
}

void FreeCamController::UpdateViewGamepad(float const _dt)
{
	kt::Vec3 cameraMove(0.0f);

	//if (input::IsDown(input::Key::KeyW))
	//{
	//	cameraMove.z += 1.0f * _dt;
	//}

	//if (input::IsDown(input::Key::KeyS))
	//{
	//	cameraMove.z -= 1.0f * _dt;
	//}

	//if (input::IsDown(input::Key::KeyA))
	//{
	//	cameraMove.x -= 1.0f * _dt;
	//}

	//if (input::IsDown(input::Key::KeyD))
	//{
	//	cameraMove.x += 1.0f * _dt;
	//}

	kt::Vec2 gamepadMove(0.0f);

	gamepadMove.x = input::GetGamepadAxis(input::GamepadAxis::LeftStick_X);
	gamepadMove.y = input::GetGamepadAxis(input::GamepadAxis::LeftStick_Y);

	float const gamepadMoveY = input::GetGamepadAxis(input::GamepadAxis::RightTrigger) * _dt + input::GetGamepadAxis(input::GamepadAxis::LeftTrigger) * -_dt;

	float const gamepadMoveLength = kt::Length(gamepadMove);
	if (gamepadMoveLength > 1.0f)
	{
		gamepadMove /= gamepadMoveLength;
	}

	gamepadMove *= _dt;

	Move(cameraMove + kt::Vec3(gamepadMove.x, gamepadMoveY, gamepadMove.y));

	kt::Vec2 gamepadRot(0.0f);
	gamepadRot.y = -input::GetGamepadAxis(input::GamepadAxis::RightStick_Y);
	gamepadRot.x = input::GetGamepadAxis(input::GamepadAxis::RightStick_X);

	float const gamepadRotLength = kt::Length(gamepadRot);
	if (gamepadRotLength > 1.0f)
	{
		gamepadRot /= gamepadRotLength;
	}

	RotateXY(gamepadRot * _dt);
	UpdateView();
}

void FreeCamController::SetProjectionParams(ProjectionParams const& _params)
{
	m_projectionParams = _params;
	m_camera.SetProjectionMatrix(kt::Mat4::PerspectiveLH_ZO(_params.m_fov, _params.m_aspect, _params.m_nearPlane, _params.m_farPlane));
}

FreeCamController::ProjectionParams const& FreeCamController::GetProjectionParams() const
{
	return m_projectionParams;
}


Camera const& FreeCamController::GetCam() const
{
	return m_camera;
}

Camera& FreeCamController::GetCam()
{
	return m_camera;
}

}
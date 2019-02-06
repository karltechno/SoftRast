#pragma once 
#include <kt/Vec2.h>
#include <kt/Vec3.h>
#include <kt/Mat4.h>
#include <kt/Quat.h>

namespace sr
{


struct Camera
{
	void SetProjectionMatrix(kt::Mat4 const& _mat);
	void SetCameraPos(kt::Vec3 const& _pos);
	void SetCameraMatrix(kt::Mat4 const& _viewToWorld);

	kt::Mat4 const& GetCachedViewProj() const;
	kt::Mat4 const& GetCameraMatrix() const;
	kt::Mat4 const& GetInverseCameraMatrix() const;
	kt::Mat4 const& GetProjectionMatrix() const;

private:
	kt::Mat4 m_viewToWorld = kt::Mat4::Identity(); // Camera in world space
	kt::Mat4 m_worldToView = kt::Mat4::Identity(); // Inverse camera
	kt::Mat4 m_viewToClip = kt::Mat4::Identity(); // projection

	mutable kt::Mat4 m_cachedWorldToClip; // world -> view -> proj

	bool m_cachedWorldToClipDirty = true;
};

struct FreeCamController
{
	struct ProjectionParams
	{
		float m_nearPlane;
		float m_farPlane;
		float m_fov;
		float m_aspect;
	};

	void SetPos(kt::Vec3 const& _pos);
	void Move(kt::Vec3 const& _movement);
	void RotateXY(kt::Vec2 const& _xy);
	void RotateByMatrix(kt::Mat3 const& _mtx);
	void SetRotation(kt::Mat3 const& _rot);

	void UpdateView();

	void UpdateViewGamepad(float const _dt);

	void SetProjectionParams(ProjectionParams const& _params);
	ProjectionParams const& GetProjectionParams() const;

	Camera& GetCam();
	Camera const& GetCam() const;

private:
	float m_speedMult = 1.0f;

	kt::Quat m_camQuat = kt::Quat::Identity();
	kt::Vec3 m_camPos = kt::Vec3(0.0f);

	kt::Vec3 m_frameMovement = kt::Vec3(0.0f);

	ProjectionParams m_projectionParams;
	Camera m_camera;
};

}
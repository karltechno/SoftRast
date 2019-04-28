#pragma once
#include "Scene.h"
#include "Camera.h"
#include "Obj.h"

namespace sr
{

struct SponzaScene : Scene
{
	SponzaScene(char const* _modelPath, uint32_t _loadFlags);

	void Init(uint32_t _screenHeight, uint32_t _screenWidth) override;
	void Update(RenderContext& _ctx, FrameBuffer& _fb, float _dt) override;

	struct PointLight
	{
		kt::Vec3 m_pos;
		kt::Vec3 m_colour;
		float m_intensity;
		float m_falloff;

		// Anim:
		kt::Vec3 m_basePos;
		kt::Vec3 m_rotOffset;
		kt::Vec3 m_rotAxis;
		float m_angle = 0.0f;

		kt::Vec3 m_colourA;
		kt::Vec3 m_colourB;
	};

	struct Constants
	{
		static uint32_t constexpr c_numPointLights = 20;

		__m256 m_sunDir[3];
		__m256 m_ambCol[3];


		PointLight m_pointLights[c_numPointLights];
	};

	sr::Obj::Model m_model;
	FreeCamController m_camController;

	float m_animPhase = 0.0f;
};

}
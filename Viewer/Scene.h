#pragma once
#include "Obj.h"
#include "Camera.h"

namespace sr
{
class RenderContext;
struct FrameBuffer;

struct Scene
{
	virtual ~Scene() {}

	virtual void Init(uint32_t _screenHeight, uint32_t _screenWidth) {}
	virtual void Update(RenderContext& _ctx, FrameBuffer& _fb, float _dt) = 0;
};

struct SimpleModelScene : Scene
{
	SimpleModelScene(char const* _modelPath, uint32_t _loadFlags);

	void Init(uint32_t _screenHeight, uint32_t _screenWidth) override;
	void Update(RenderContext& _ctx, FrameBuffer& _fb, float _dt) override;

	sr::Obj::Model m_model;
	FreeCamController m_camController;
};

}
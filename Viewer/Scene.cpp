#include "Scene.h"
#include "Renderer.h"
#include "Shaders.h"


namespace sr
{

SimpleModelScene::SimpleModelScene(char const* _modelPath, uint32_t _loadFlags)
{
	m_model.Load(_modelPath, kt::GetDefaultAllocator(), _loadFlags);
}

void SimpleModelScene::Init(uint32_t _screenHeight, uint32_t _screenWidth)
{
	sr::FreeCamController::ProjectionParams proj;
	proj.m_aspect = float(sr::Config::c_screenWidth) / float(sr::Config::c_screenHeight);
	proj.m_fov = kt::ToRadians(85.0f);

#if SR_USE_REVERSE_Z
	proj.m_nearPlane = 10000.0f;
	proj.m_farPlane = 0.1f;
#else	
	proj.m_nearPlane = 0.1f;
	proj.m_farPlane = 10000.0f;
#endif

	m_camController.SetProjectionParams(proj);
	m_camController.SetPos({ 0.0f, 0.0f, -2.0f });
}

void SimpleModelScene::Update(RenderContext& _ctx, FrameBuffer& _fb, float _dt)
{
	m_camController.UpdateViewGamepad(_dt);
	_ctx.ClearFrameBuffer(_fb, 0);

	for (sr::Obj::Mesh const& mesh : m_model.m_meshes)
	{
		sr::DrawCall call;
		call.SetFrameBuffer(&_fb);
		call.m_mvp = m_camController.GetCam().GetCachedViewProj();

		call.SetAttributeBuffer(mesh.m_vertexData.Data(), sizeof(sr::Obj::Vertex), mesh.m_vertexData.Size(), offsetof(sr::Obj::Vertex, uv) / sizeof(float));

		call.m_positionBuffer.m_ptr = (uint8_t*)mesh.m_vertexData.Data();
		call.m_positionBuffer.m_stride = sizeof(sr::Obj::Vertex);
		call.m_positionBuffer.m_num = mesh.m_vertexData.Size();

		call.m_indexBuffer.m_ptr = mesh.m_indexData.Data();
		call.m_indexBuffer.m_num = mesh.m_numIndices;
		call.m_indexBuffer.m_stride = mesh.m_indexType == sr::IndexType::u16 ? sizeof(uint16_t) : sizeof(uint32_t);

		if (mesh.m_matIdx < m_model.m_materials.Size())
		{
			call.m_pixelUniforms = &m_model.m_materials[mesh.m_matIdx].m_diffuse;
			call.m_pixelShader = sr::shader::UnlitDiffuseShader;
		}
		else
		{
			call.m_pixelShader = sr::shader::VisualizeNormalsShader;
		}

		_ctx.DrawIndexed(call);
	}
}

}
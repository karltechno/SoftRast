#include "SponzaScene.h"
#include "Renderer.h"
#include "Shaders.h"
#include "kt/Random.h"
#include "kt/inl/Quat.inl"

namespace sr
{

// hack as a global for now.
static SponzaScene::Constants g_constants;

static void SponzaShader(void const* _uniforms, float const* _varyings, uint32_t o_texels[8], uint32_t _execMask)
{
	sr::Tex::TextureData* tex = (sr::Tex::TextureData*)_uniforms;

	if (!tex || tex->m_texels.Size() == 0)
	{
		memset(o_texels, 0xFFFFFFFF, 8 * sizeof(uint32_t));
		return;
	}

	sr::shader::OBJVaryings const objVaryings = sr::shader::UnpackOBJVaryings(_varyings);
	sr::shader::Derivatives const derivs = sr::shader::UnpackDerivatives(_varyings);

	__m256 radiance[3];

	// sun
	{
		__m256 nDotL = simdutil::Dot3SoA(objVaryings.norm_x, objVaryings.norm_y, objVaryings.norm_z, g_constants.m_sunDir[0], g_constants.m_sunDir[1], g_constants.m_sunDir[2]);

		// magic bias
		nDotL = _mm256_max_ps(_mm256_set1_ps(0.1f), nDotL);
		radiance[0] = nDotL;
		radiance[1] = nDotL;
		radiance[2] = nDotL;
	}


	for (uint32_t i = 0; i < SponzaScene::Constants::c_numPointLights; ++i)
	{
		SponzaScene::PointLight const& light = g_constants.m_pointLights[i];
	
		__m256 const lightPos_x = _mm256_broadcast_ss(&light.m_pos.x);
		__m256 const lightPos_y = _mm256_broadcast_ss(&light.m_pos.y);
		__m256 const lightPos_z = _mm256_broadcast_ss(&light.m_pos.z);
	
		__m256 const pToL_x = _mm256_sub_ps(lightPos_x, objVaryings.pos_x);
		__m256 const pToL_y = _mm256_sub_ps(lightPos_y, objVaryings.pos_y);
		__m256 const pToL_z = _mm256_sub_ps(lightPos_z, objVaryings.pos_z);

		__m256 const distSq = simdutil::Dot3SoA(pToL_x, pToL_y, pToL_z, pToL_x, pToL_y, pToL_z);
		__m256 const dist = _mm256_sqrt_ps(distSq);
		__m256 const recipDist = _mm256_rcp_ps(dist);
		
		__m256 const l_x = _mm256_mul_ps(pToL_x, recipDist);
		__m256 const l_y = _mm256_mul_ps(pToL_y, recipDist);
		__m256 const l_z = _mm256_mul_ps(pToL_z, recipDist);

		__m256 const nDotL = _mm256_max_ps(_mm256_setzero_ps(), simdutil::Dot3SoA(l_x, l_y, l_z, objVaryings.norm_x, objVaryings.norm_y, objVaryings.norm_z));

		__m256 const one = _mm256_set1_ps(1.0f);

		// Inverse square law
		
		__m256 const atten = _mm256_rcp_ps(_mm256_add_ps(one, _mm256_fmadd_ps(_mm256_set1_ps(0.1f), dist, _mm256_mul_ps(distSq, _mm256_set1_ps(0.01f)))));
		__m256 const lightRadiance = _mm256_mul_ps(nDotL, _mm256_mul_ps(_mm256_broadcast_ss(&light.m_intensity), atten));

		radiance[0] = _mm256_add_ps(radiance[0], _mm256_mul_ps(lightRadiance, _mm256_broadcast_ss(&light.m_colour.x)));
		radiance[1] = _mm256_add_ps(radiance[1], _mm256_mul_ps(lightRadiance, _mm256_broadcast_ss(&light.m_colour.y)));
		radiance[2] = _mm256_add_ps(radiance[2], _mm256_mul_ps(lightRadiance, _mm256_broadcast_ss(&light.m_colour.z)));
	}


	// Apply ambient term.
	radiance[0] = _mm256_add_ps(radiance[0], g_constants.m_ambCol[0]);
	radiance[1] = _mm256_add_ps(radiance[1], g_constants.m_ambCol[1]);
	radiance[2] = _mm256_add_ps(radiance[2], g_constants.m_ambCol[2]);

	__m256 r;
	__m256 g;
	__m256 b;
	__m256 a;

	sr::Tex::SampleWrap(*tex, objVaryings.u, objVaryings.v, derivs.dudx, derivs.dudy, derivs.dvdx, derivs.dvdy, r, g, b, a, _execMask);

	r = _mm256_mul_ps(radiance[0], r);
	g = _mm256_mul_ps(radiance[1], g);
	b = _mm256_mul_ps(radiance[2], b);

	sr::simdutil::RGBA32SoA_To_RGBA8AoS(r, g, b, a, o_texels);
}
SponzaScene::SponzaScene(char const* _modelPath, uint32_t _loadFlags)
{
	m_model.Load(_modelPath, kt::GetDefaultAllocator(), _loadFlags);
}

void SponzaScene::Init(uint32_t _screenHeight, uint32_t _screenWidth)
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
	m_camController.SetPos({ 0.0f, 0.0f, 2.0f });

	{
		g_constants.m_ambCol[0] = _mm256_set1_ps(0.1f);
		g_constants.m_ambCol[1] = _mm256_set1_ps(0.1f);
		g_constants.m_ambCol[2] = _mm256_set1_ps(0.1f);
	}

	{
		kt::Vec3 const sunDir = kt::Normalize(kt::Vec3(0.4f, 0.7f, 0.1f));
		g_constants.m_sunDir[0] = _mm256_broadcast_ss(&sunDir.x);
		g_constants.m_sunDir[1] = _mm256_broadcast_ss(&sunDir.x);
		g_constants.m_sunDir[2] = _mm256_broadcast_ss(&sunDir.x);
	}

	kt::XorShift32 rng;

	for (PointLight& light : g_constants.m_pointLights)
	{
		light.m_colourA = kt::Vec3(kt::RandomUnitFloat(rng), kt::RandomUnitFloat(rng), kt::RandomUnitFloat(rng));
		light.m_colourB = kt::Vec3(kt::RandomUnitFloat(rng), kt::RandomUnitFloat(rng), kt::RandomUnitFloat(rng));
	
		light.m_basePos.x = kt::Lerp(-1000.0f, 1000.0f, kt::RandomUnitFloat(rng));
		light.m_basePos.y = kt::Lerp(50.0f, 250.0f, kt::RandomUnitFloat(rng));
		light.m_basePos.z = kt::Lerp(-150.0f, 150.0f, kt::RandomUnitFloat(rng));

		light.m_falloff = kt::Lerp(500.0f, 2500.0f, kt::RandomUnitFloat(rng));
		light.m_intensity = kt::Lerp(150.0f, 350.0f, kt::RandomUnitFloat(rng));

		light.m_rotAxis = kt::Vec3(kt::RandomUnitFloat(rng), kt::RandomUnitFloat(rng), kt::RandomUnitFloat(rng));
		light.m_rotAxis = kt::Normalize(light.m_rotAxis * 2.0f - kt::Vec3(1.0f));
		light.m_rotOffset = 25.0f * (kt::Vec3(kt::RandomUnitFloat(rng), kt::RandomUnitFloat(rng), kt::RandomUnitFloat(rng)) * 2.0f - kt::Vec3(1.0f));
	}
}

void SponzaScene::Update(RenderContext& _ctx, FrameBuffer& _fb, float _dt)
{
	m_camController.UpdateViewGamepad(_dt);
	_ctx.ClearFrameBuffer(_fb, 0);

	float const sinT = (sinf(m_animPhase) * 0.5f + 0.5f);
	
	m_animPhase += _dt;

	for (PointLight& pointLight : g_constants.m_pointLights)
	{
		pointLight.m_colour = kt::Lerp(pointLight.m_colourA, pointLight.m_colourB, sinT);

		kt::Quat rot; 
		rot.FromNormalizedAxisAngle(pointLight.m_rotAxis, pointLight.m_angle);

		kt::Vec3 pos = pointLight.m_basePos - pointLight.m_rotOffset;
		pos = kt::Mul(rot, pos);
		pointLight.m_pos = pos + pointLight.m_rotOffset;

		pointLight.m_angle += _dt;
	}

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
		call.m_indexBuffer.m_num = mesh.m_numIndicies;
		call.m_indexBuffer.m_stride = mesh.m_indexType == sr::IndexType::u16 ? sizeof(uint16_t) : sizeof(uint32_t);
		call.m_pixelShader = SponzaShader;

		if (mesh.m_matIdx < m_model.m_materials.Size())
		{
			call.m_pixelUniforms = &m_model.m_materials[mesh.m_matIdx].m_diffuse;
		}
		else
		{
			call.m_pixelUniforms = nullptr;
		}

		_ctx.DrawIndexed(call);
	}
}

}
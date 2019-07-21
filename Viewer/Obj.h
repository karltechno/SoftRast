#pragma once
#include <stdint.h>
#include <kt/Array.h>
#include <kt/Strings.h>

#include "SoftRastTypes.h"
#include "Texture.h"

namespace sr
{


namespace Obj
{

struct Vertex
{
	kt::Vec3 pos;
	kt::Vec3 norm;
	kt::Vec2 uv;
};


struct Mesh
{
	KT_NO_COPY(Mesh);

	Mesh() = default;
	Mesh(Mesh&&) = default;
	Mesh& operator=(Mesh&&) = default;

	void Clear();

	kt::Array<uint8_t> m_indexData;

	IndexType m_indexType = IndexType::u16;
	uint32_t m_numIndices = 0;

	kt::Array<Vertex> m_vertexData;

	uint32_t m_matIdx = 0;
};

struct Material
{
	Material() = default;
	Material(Material&&) = default;
	Material& operator=(Material&&) = default;

	kt::String128 m_name;
	Tex::TextureData m_diffuse;
};

enum LoadFlags : uint32_t
{
	None = 0x0,
	FlipWinding = 0x1,
	GenNormals = 0x2, // todo
	FlipUVs = 0x4
};

struct Model
{
	~Model();

	bool Load(char const* _path, kt::IAllocator* _tempAllocator, uint32_t const _flags);
	void Clear();

	kt::Array<Mesh> m_meshes;
	kt::Array<Material> m_materials;
};

}


}

namespace kt
{

template<>
void Serialize(ISerializer* _s, sr::Obj::Model& _model);

template<>
void Serialize(ISerializer* _s, sr::Obj::Mesh& _mesh);

template<>
void Serialize(ISerializer* _s, sr::Obj::Material& _mat);
}
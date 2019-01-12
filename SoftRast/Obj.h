#pragma once
#include <stdint.h>
#include <kt/Array.h>

#include "SoftRastTypes.h"

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
	void Clear();

	union
	{
		uint16_t* index16 = nullptr;
		uint32_t* index32;
	} m_indexData;

	IndexType m_indexType = IndexType::u16;

	uint32_t m_numIndicies = 0;

	Vertex* m_vertexData = nullptr;
	uint32_t m_numVertices = 0;
};

enum LoadFlags : uint32_t
{
	FlipWinding = 0x1,
	GenNormals = 0x2
};

struct Model
{
	bool Load(char const* _path, kt::IAllocator* _tempAllocator, LoadFlags const _flags = LoadFlags::FlipWinding);
	void Clear();

	kt::Array<Mesh> m_meshes;
};


}


}
#include "Obj.h"

#include <stdio.h>

#include <kt/Memory.h>
#include <kt/HashMap.h>
#include <kt/Logging.h>
#include <kt/FilePath.h>
#include <kt/File.h>
#include <kt/Serialization.h>

namespace kt
{


template<>
void Serialize(ISerializer* _s, sr::Obj::Model& _model)
{
	Serialize(_s, _model.m_meshes);
	Serialize(_s, _model.m_materials);
}

template<>
void Serialize(ISerializer* _s, sr::Obj::Mesh& _mesh)
{
	Serialize(_s, _mesh.m_indexType);
	Serialize(_s, _mesh.m_indexData);
	Serialize(_s, _mesh.m_numIndicies);
	Serialize(_s, _mesh.m_vertexData);
	Serialize(_s, _mesh.m_matIdx);
}

template<>
void Serialize(ISerializer* _s, sr::Obj::Material& _mat)
{
	Serialize(_s, _mat.m_diffuse);
	Serialize(_s, _mat.m_name);
}


}

namespace sr
{


namespace Obj
{

struct TempFace
{
	uint32_t posIdx;
	uint32_t uvIdx;
	uint32_t normIdx;
};

bool operator==(TempFace const& _lhs, TempFace const& _rhs)
{
	return _lhs.normIdx == _rhs.normIdx && _lhs.posIdx == _rhs.posIdx && _lhs.uvIdx == _rhs.uvIdx;
}

template <typename T>
static void FlipMeshWinding(Mesh& _m, void* _idxBuff)
{
	T* buff = (T*)_idxBuff;
	for (uint32_t i = 0; i < _m.m_numIndicies; i += 3)
	{
		kt::Swap(buff[i + 1], buff[i + 2]);
	}
}

struct MeshParserState
{
	MeshParserState(kt::IAllocator* _tempAllocator)
		: m_faceMap(_tempAllocator)
		, m_tempVerticies(_tempAllocator)
		, m_tempIndexBuffer(_tempAllocator)
		, m_tempPos(_tempAllocator)
		, m_tempUv(_tempAllocator)
		, m_tempNormal(_tempAllocator)
	{
		m_faceMap.Reserve(4096);
		m_tempPos.Reserve(2048);
		m_tempUv.Reserve(2048);
		m_tempNormal.Reserve(2048);
		m_tempIndexBuffer.Reserve(2048);
		m_tempVerticies.Reserve(2048);
	}

	void FinalizeMesh(Mesh* m, uint32_t _matIdx)
	{
		m->Clear();
		m->m_matIdx = _matIdx;
		if (!m_tempVerticies.Size())
		{
			return;
		}

		m->m_indexType = m_tempVerticies.Size() > UINT16_MAX ? IndexType::u32 : IndexType::u16;

		if (m->m_indexType == IndexType::u32)
		{
			uint32_t const allocSz = sizeof(uint32_t) * m_tempIndexBuffer.Size();
			memcpy(m->m_indexData.PushBack_Raw(allocSz), m_tempIndexBuffer.Data(), allocSz);
		}
		else
		{
			uint32_t const allocSz = sizeof(uint16_t) * m_tempIndexBuffer.Size();
			uint16_t *pDst = (uint16_t*)m->m_indexData.PushBack_Raw(allocSz);
			for (uint32_t idx32 : m_tempIndexBuffer)
			{
				*pDst++ = (uint16_t)idx32;
			}
		}

		m->m_numIndicies = m_tempIndexBuffer.Size();
		Vertex* vertexWrite = m->m_vertexData.PushBack_Raw(m_tempVerticies.Size());
		memcpy(vertexWrite, m_tempVerticies.Data(), sizeof(Vertex) * m_tempVerticies.Size());

		m_faceMap.Clear();
		m_tempVerticies.Clear();
		m_tempIndexBuffer.Clear();
	}

	// data is index into m_tempVerticies
	kt::HashMap<TempFace, uint32_t> m_faceMap;

	kt::Array<Vertex> m_tempVerticies;
	kt::Array<uint32_t> m_tempIndexBuffer;

	kt::Array<kt::Vec3> m_tempPos;
	kt::Array<kt::Vec2> m_tempUv;
	kt::Array<kt::Vec3> m_tempNormal;
};

static char* StripWhiteSpaceAndNewLine(char* _buff)
{
	char* ret = _buff;
	// strip preceding whitespace
	while (*ret && (*ret == ' ' || *ret == '\t')) { ++ret; }

	char* temp = ret;

	if (size_t const len = strlen(temp))
	{
		// strip trailing whitespace/newline
		temp += (len - 1);
		while (temp != ret && (*temp == ' ' || *temp == '\t' || *temp == '\r' || *temp == '\n'))
		{
			*temp-- = '\0';
		}
	}

	return ret;
}

static int32_t FixupFaceIdx(int32_t _idx, int32_t _totalVerticies)
{
	if (!_idx) return 0;

	return _idx < 0 ? _totalVerticies + _idx : _idx - 1;
}

static bool ParseFace(MeshParserState& _parserState, char const* _line)
{
	char const* p = _line + 2;

	// 6 so we can expand quad into 2 tris
	int32_t posIndicies[6] = {};
	int32_t texIndicies[6] = {};
	int32_t normIndicies[6] = {};

	int32_t numFaceVerts = 0;

	while (*p && numFaceVerts < 4)
	{
		auto parseNum = [&p](int32_t& _output)
		{
			while (*p == ' ' || *p == '\t')
			{
				++p;
			}

			int32_t sign = 1;
			if (*p == '-')
			{
				sign = -1;
				++p;
			}

			int32_t i = 0;
			while (*p >= '0' && *p <= '9')
			{
				i = i * 10 + *p - '0';
				++p;
			}

			_output = i * sign;
		};

		int32_t vertIdx = numFaceVerts;

		parseNum(posIndicies[vertIdx]);

		if (posIndicies[vertIdx] == 0)
		{
			break;
		}

		++numFaceVerts;

		if (*p != '/')
		{
			continue;
		}
		++p;

		if (*p == '/')
		{
			++p;
			// Missing tex idx
			parseNum(normIndicies[vertIdx]);
			continue;
		}

		parseNum(texIndicies[vertIdx]);

		if (*p != '/')
		{
			continue;
		}
		++p;
		parseNum(normIndicies[vertIdx]);
	}

	// triangulate the quad
	if (numFaceVerts == 4)
	{
		// dumb triangulation
		posIndicies[5] = posIndicies[3];
		normIndicies[5] = normIndicies[3];
		texIndicies[5] = texIndicies[3];

		posIndicies[3] = posIndicies[0];
		normIndicies[3] = normIndicies[0];
		texIndicies[3] = texIndicies[0];
	
		posIndicies[4] = posIndicies[2];
		normIndicies[4] = normIndicies[2];
		texIndicies[4] = texIndicies[2];
		numFaceVerts = 6;
	}


	TempFace tempFace;
	for (int32_t vId = 0; vId < numFaceVerts; ++vId)
	{
		tempFace.posIdx = FixupFaceIdx(posIndicies[vId], (int32_t)_parserState.m_tempPos.Size());
		tempFace.uvIdx = FixupFaceIdx(texIndicies[vId], (int32_t)_parserState.m_tempUv.Size());
		tempFace.normIdx = FixupFaceIdx(normIndicies[vId], (int32_t)_parserState.m_tempNormal.Size());

		kt::HashMap<TempFace, uint32_t>::Iterator it = _parserState.m_faceMap.Find(tempFace);
		if (it != _parserState.m_faceMap.End())
		{
			_parserState.m_tempIndexBuffer.PushBack(it->m_val);
		}
		else
		{
			uint32_t const idx = _parserState.m_tempVerticies.Size();
			Vertex* v = _parserState.m_tempVerticies.PushBack_Raw();
			_parserState.m_tempIndexBuffer.PushBack(idx);

			_parserState.m_faceMap[tempFace] = idx;

			if (tempFace.posIdx >= _parserState.m_tempPos.Size())
			{
				KT_ASSERT(false);
				return false;
			}

			memcpy(&v->pos, &_parserState.m_tempPos[tempFace.posIdx], sizeof(float) * 3);

			if (_parserState.m_tempNormal.Size() != 0)
			{
				if (tempFace.normIdx >= _parserState.m_tempNormal.Size())
				{
					KT_ASSERT(false);
					return false;
				}

				memcpy(&v->norm, &_parserState.m_tempNormal[tempFace.normIdx], sizeof(float) * 3);
			}
			else
			{
				memset(&v->norm, 0, sizeof(float) * 3);
			}

			if (_parserState.m_tempUv.Size() != 0)
			{
				if (tempFace.uvIdx >= _parserState.m_tempUv.Size())
				{
					KT_ASSERT(false);
					return false;
				}

				memcpy(&v->uv, &_parserState.m_tempUv[tempFace.uvIdx], sizeof(float) * 2);
			}
			else
			{
				memset(&v->uv, 0, sizeof(float) * 2);
			}
		}
	}

	return true;
}

static void ParseMaterial(FILE* _file, Model& _model, kt::FilePath const& _rootPath)
{
	char lineBuff[2048];

	Material* curMat = nullptr;

	while (fgets(lineBuff, sizeof(lineBuff), _file))
	{
		char* line = StripWhiteSpaceAndNewLine(lineBuff);
		switch (*line)
		{
			case 'm':
			{
				static uint32_t const map_Kd_len = 6;
				if (strncmp(line, "map_Kd", map_Kd_len) == 0)
				{
					if (!curMat)
					{
						KT_LOG_INFO("No newmtl directive, can't parse mtl!");
						break;
					}

					char* fileName = StripWhiteSpaceAndNewLine(line + map_Kd_len);
					kt::FilePath diffusePath = _rootPath;

					diffusePath.Append(fileName);
					curMat->m_diffuse.CreateFromFile(diffusePath.Data());
				}
			} break;
		
			case 'n':
			{
				static uint32_t const newmtl_len = 6;
				if (strncmp(line, "newmtl", newmtl_len) == 0)
				{
					curMat = &_model.m_materials.PushBack();
					curMat->m_name = StripWhiteSpaceAndNewLine(line + newmtl_len);
				}
			} break;
		}
	}
}

Model::~Model()
{
	for (Mesh& m : m_meshes)
	{
		m.Clear();
	}

	for (Material& mat : m_materials)
	{
		mat.m_diffuse.Clear();
	}
}

bool Model::Load(char const* _path, kt::IAllocator* _tempAllocator, uint32_t const _flags)
{
	kt::String1024 binpath;
	binpath.AppendFmt("%s.bin", _path);
	
	if (kt::FileExists(binpath.Data()))
	{
		FILE* cacheFile = fopen(binpath.Data(), "rb");
		if (cacheFile)
		{
			KT_LOG_INFO("Found OBJ cache %s", binpath.Data());
			kt::FileReader reader(cacheFile);
			kt::ISerializer serializer(&reader, 0);
			kt::Serialize(&serializer, *this);
			return true; // todo: error checking
		}
		else
		{
			KT_LOG_ERROR("Failed to open obj cache file %s", binpath.Data());
		}
	}


	FILE* f = fopen(_path, "r");
	if (!f)
	{
		KT_LOG_ERROR("Failed to open obj file: %s", _path);
		return false;
	}
	KT_SCOPE_EXIT(fclose(f));

	KT_LOG_INFO("No cache found, parsing obj %s.", _path);

	Clear();

	char lineBuff[2048];

	MeshParserState parserState(_tempAllocator);

	kt::FilePath const rootPath(kt::FilePath(_path).GetPath());

	uint32_t curMatIdx = 0;

	while (fgets(lineBuff, sizeof(lineBuff), f))
	{
		char* line = StripWhiteSpaceAndNewLine(lineBuff);

		switch (line[0])
		{
			case 'v':
			{
				if (line[1] == ' ')
				{
					// pos				
					float* readVtx = (float*)parserState.m_tempPos.PushBack_Raw();
					int numRead = sscanf(line + 2, "%f %f %f", readVtx, readVtx + 1, readVtx + 2);
					if (numRead != 3)
					{
						KT_LOG_ERROR("Failed to parse obj, Bad vertex pos!");
						Clear();
						return false;
					}
				}
				else if (line[1] == 't')
				{
					// uv coord			
					float* readUv = (float*)parserState.m_tempUv.PushBack_Raw();
					int numRead = sscanf(line + 3, "%f %f", readUv, readUv + 1);

					if (numRead != 2)
					{
						KT_LOG_ERROR("Failed to parse obj, Bad uv coord!");
						Clear();
						return false;
					}

					if (_flags & LoadFlags::FlipUVs)
					{
						readUv[1] = 1.0f - readUv[1];
					}
				}
				else if (line[1] == 'n')
				{
					// normal	
					float* readNormal = (float*)parserState.m_tempNormal.PushBack_Raw();
					int numRead = sscanf(line + 3, "%f %f %f", readNormal, readNormal + 1, readNormal + 2);
					if (numRead != 3)
					{
						KT_LOG_ERROR("Failed to parse obj, Bad vertex normal!");
						Clear();
						return false;
					}
				}
			} break;

			case 'f':
			{
				if (!ParseFace(parserState, line))
				{
					KT_LOG_ERROR("Failed to parse obj, Bad vertex face!");
					Clear();
					return false;
				}
			} break;

			case 'g':
			{
				if (parserState.m_tempVerticies.Size())
				{
					parserState.FinalizeMesh(&m_meshes.PushBack(), curMatIdx);
				}
			} break;

			case 'u':
			{
				static uint32_t const usemtl_len = 6;
				if (strncmp(line, "usemtl", usemtl_len) == 0)
				{
					char* mtl = StripWhiteSpaceAndNewLine(line + usemtl_len);
					for (uint32_t mtlIdx = 0; mtlIdx < m_materials.Size(); ++mtlIdx)
					{
						if (m_materials[mtlIdx].m_name == mtl)
						{
							curMatIdx = mtlIdx;
							break;
						}
					}
				}
			} break;

			case 'm':
			{
				static uint32_t const mtlliblen = 6;
				if (strncmp("mtllib", line, mtlliblen) == 0)
				{
					char* mtlName = line + mtlliblen;
					while (*mtlName && (*mtlName == ' ' || *mtlName == '\t')) { ++mtlName; }

					if (!(*mtlName))
					{
						KT_LOG_ERROR("Invalid material name in obj file: %s", _path);
						break;
					}

					kt::FilePath mtlPath = rootPath;
					mtlPath.Append(mtlName);

					FILE* mtlFile = fopen(mtlPath.Data(), "r");
					KT_SCOPE_EXIT(fclose(mtlFile));
					if (!mtlFile)
					{
						KT_LOG_ERROR("Failed to open material file: %s", mtlPath.Data());
					}
					ParseMaterial(mtlFile, *this, rootPath);
				}
			} break;
		}
	}

	if (parserState.m_tempVerticies.Size())
	{
		parserState.FinalizeMesh(&m_meshes.PushBack(), curMatIdx);
	}

	if (_flags & LoadFlags::FlipWinding)
	{
		for (Mesh& mesh : m_meshes)
		{
			if (mesh.m_indexType == IndexType::u16)
			{
				FlipMeshWinding<uint16_t>(mesh, mesh.m_indexData.Data());
			}
			else
			{
				FlipMeshWinding<uint32_t>(mesh, mesh.m_indexData.Data());
			}
		}
	}

	// Now serialize to binary.
	FILE* cacheFile = fopen(binpath.Data(), "wb");
	if (cacheFile)
	{
		kt::FileWriter writer(cacheFile);
		kt::ISerializer serializer(&writer, 0);
		kt::Serialize(&serializer, *this);
	}
	else
	{
		KT_LOG_ERROR("Failed to write obj cache file %s.", binpath.Data());
	}

	return true;
}

void Model::Clear()
{
	for (Mesh& m : m_meshes)
	{
		m.Clear();
	}
	m_meshes.Clear();
}

void Mesh::Clear()
{
	m_vertexData.ClearAndFree();
	m_indexData.ClearAndFree();
}

}

}
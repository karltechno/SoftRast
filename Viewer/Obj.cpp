#include "Obj.h"

#include <stdio.h>
#include <kt/Memory.h>
#include <kt/HashMap.h>
#include <kt/Logging.h>
#include <kt/FilePath.h>

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
static void FlipMeshWinding(Mesh& _m, T* _idxBuff)
{
	for (uint32_t i = 0; i < _m.m_numIndicies; i += 3)
	{
		kt::Swap(_idxBuff[i + 1], _idxBuff[i + 2]);
	}
}

struct MeshParserState
{
	MeshParserState(kt::IAllocator* _tempAllocator)
		: m_faceMap(_tempAllocator)
		, m_tempPos(_tempAllocator)
		, m_tempUv(_tempAllocator)
		, m_tempNormal(_tempAllocator)
		, m_tempVerticies(_tempAllocator)
		, m_tempIndexBuffer(_tempAllocator)
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
		uint32_t const idxSize = m->m_indexType == IndexType::u32 ? 4 : 2;

		if (m->m_indexType == IndexType::u32)
		{
			uint32_t const allocSz = sizeof(uint32_t) * m_tempIndexBuffer.Size();
			m->m_indexData.index32 = (uint32_t*)kt::Malloc(allocSz);
			memcpy(m->m_indexData.index32, m_tempIndexBuffer.Data(), allocSz);
		}
		else
		{
			uint32_t const allocSz = sizeof(uint16_t) * m_tempIndexBuffer.Size();
			m->m_indexData.index16 = (uint16_t*)kt::Malloc(allocSz);
			uint16_t *pDst = m->m_indexData.index16;
			for (uint32_t idx32 : m_tempIndexBuffer)
			{
				*pDst++ = (uint16_t)idx32;
			}
		}

		m->m_numIndicies = m_tempIndexBuffer.Size();
		m->m_numVertices = m_tempVerticies.Size();
		m->m_vertexData = (Vertex*)kt::Malloc(sizeof(Vertex) * m->m_numVertices);
		memcpy(m->m_vertexData, m_tempVerticies.Data(), sizeof(Vertex) * m->m_numVertices);

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
	while (*ret && *ret == ' ' || *ret == '\t') { ++ret; }

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

	int32_t posIndicies[3] = {};
	int32_t texIndicies[3] = {};
	int32_t normIndicies[3] = {};

	int32_t numFaceVerts = 0;

	while (*p && numFaceVerts < 3)
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
				memset(&v->norm, 0, sizeof(float) * 2);
			}
		}
	}

	return true;
}

static void ParseMaterial(FILE* _file, Material* _mat, kt::FilePath const& _rootPath)
{
	char lineBuff[2048];

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
					char* fileName = StripWhiteSpaceAndNewLine(line + map_Kd_len);
					kt::FilePath diffusePath = _rootPath;

					diffusePath.Append(fileName);
					_mat->m_diffuse.CreateFromFile(diffusePath.Data());
				}
			} break;
		
			case 'n':
			{
				static uint32_t const newmtl_len = 6;
				if (strncmp(line, "newmtl", newmtl_len) == 0)
				{
					_mat->m_name = StripWhiteSpaceAndNewLine(line + newmtl_len);
				}
			} break;
		}
	}
}

bool Model::Load(char const* _path, kt::IAllocator* _tempAllocator, uint32_t  const _flags)
{
	FILE* f = fopen(_path, "r");
	if (!f)
	{
		KT_LOG_ERROR("Failed to open obj file: %s", _path);
		return false;
	}
	KT_SCOPE_EXIT(fclose(f));

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
					ParseMaterial(mtlFile, &m_materials.PushBack(), rootPath);
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
				FlipMeshWinding(mesh, mesh.m_indexData.index16);
			}
			else
			{
				FlipMeshWinding(mesh, mesh.m_indexData.index32);
			}
		}
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
	kt::Free(m_vertexData);
	kt::Free(m_indexData.index16);
}

}

}
#include <kt/Memory.h>

#define STBI_MALLOC(sz)           kt::Malloc(sz)
#define STBI_REALLOC(p,newsz)     kt::Realloc(p,newsz)
#define STBI_FREE(p)              kt::Free(p)

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
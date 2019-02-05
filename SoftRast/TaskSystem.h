#pragma once
#include <stdint.h>
#include "kt/LinearAllocator.h"

namespace sr
{

using ThreadScratchAllocator = kt::LinearAllocator<kt::LinearAllocatorThreadSafety::ThreadSafeAlloc>;

using TaskFn = void(*)(void const* _data, uint32_t _start, uint32_t _end);

struct Packet
{
	TaskFn m_fn;
	uint32_t m_granularity;
	
	// Number of partitions done
	uint32_t m_num;

	// Total partitions.
	uint32_t m_total;
};


class System
{
public:
	static const uint32_t c_maxTasks = 1 << 16;

	void Init(int32_t const _numThreads = -1);

private:
};

}
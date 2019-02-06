#pragma once
#include <stdint.h>
#include "kt/LinearAllocator.h"
#include "kt/Array.h"

namespace sr
{

using ThreadScratchAllocator = kt::LinearAllocator<kt::LinearAllocatorThreadSafety::ThreadSafeAlloc>;

struct Task;

using TaskFn = void(*)(Task const* _task, uint32_t _threadIdx, uint32_t _start, uint32_t _end);

struct Task
{
	Task() = default;

	Task(TaskFn _fn, uint32_t _numPartitions, uint32_t _granularity, void* _user)
		: m_fn(_fn), m_totalPartitions(_numPartitions), m_granularity(_granularity), m_userData(_user)
	{}

	// Task function
	TaskFn m_fn = nullptr;

	// Iterations per task
	uint32_t m_granularity = 0;

	// Number of partitions done
	uint32_t m_numCompletedPartitions = 0;

	// Total partitions.
	uint32_t m_totalPartitions = 0;

	int32_t* m_taskCounter = nullptr;

	// User defined data.
	void* m_userData = nullptr;
};

struct TaskPacket
{
	// The task
	Task* m_task = nullptr;

	// Begin index
	uint32_t m_begin = 0;

	// End index
	uint32_t m_end = 0;
};


class TaskSystem
{
public:
	static uint32_t const MAX_TASK_PACKETS = 1 << 16;
	static uint32_t const QUEUE_MASK = MAX_TASK_PACKETS - 1;

	static uint32_t TlsThreadIdx();

	void InitFromMainThread(uint32_t const _numWorkers);
	void WaitAndShutdown();

	void PushTask(Task* _task);

	void SyncAndWaitForAll();

	void WaitForCounter(int32_t* _counter);

	uint32_t TotalThreadsIncludingMainThread();

private:
	void WorkerLoop(uint32_t _threadId);

	kt::Thread* m_threads = nullptr;
	uint32_t m_numWorkers = 0;

	TaskPacket* m_packets = nullptr;
	
	uint32_t m_queueHead = 0;
	uint32_t m_queueTail = 0;
	uint32_t m_numEntriesInQueue = 0;

	// Todo: lock free
	kt::Mutex m_queueMutex;
	kt::Event m_queueSignal;

	int32_t m_keepRunning = 1;

	int32_t m_numActiveWorkers = 0;
};

}
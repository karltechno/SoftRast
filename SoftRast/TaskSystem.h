#pragma once
#include <stdint.h>

#include "kt/LinearAllocator.h"
#include "kt/Array.h"
#include "kt/Concurrency.h"

namespace sr
{

using ThreadScratchAllocator = kt::LinearAllocator;

struct Task;

using TaskFn = void(*)(Task const* _task, uint32_t _threadIdx, uint32_t _start, uint32_t _end);

struct Task
{
	Task() = default;

	Task(TaskFn _fn, uint32_t _numPartitions, uint32_t _granularity, void* _user, std::atomic<uint32_t>* _counter = nullptr)
		: m_fn(_fn), m_granularity(_granularity),  m_numCompletedPartitions(0), m_totalPartitions(_numPartitions),  m_taskCounter(_counter), m_userData(_user)
	{}

	void Set(TaskFn _fn, uint32_t _numPartitions, uint32_t _granularity, void* _user, std::atomic<uint32_t>* _counter = nullptr)
	{
		KT_ASSERT(!m_taskCounter || m_taskCounter->load() == 0);
		m_fn = _fn;
		m_numCompletedPartitions = _numPartitions;
		m_granularity = _granularity;
		m_userData = _user;
		m_taskCounter = _counter;
		m_totalPartitions = _numPartitions;
	}

	// Task function
	TaskFn m_fn = nullptr;

	// Iterations per task
	uint32_t m_granularity = 0;

	// Number of partitions done
	std::atomic<uint32_t> m_numCompletedPartitions;

	// Total partitions.
	uint32_t m_totalPartitions = 0;

	std::atomic<uint32_t>* m_taskCounter = nullptr;

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

	TaskSystem()
		: m_numEntriesInQueue(0)
		, m_keepRunning(1)
		, m_numActiveWorkers(0)
	{}

	void InitFromMainThread(uint32_t const _numWorkers);
	void WaitAndShutdown();

	void PushTask(Task* _task);

	void SyncAndWaitForAll();

	void WaitForCounter(std::atomic<uint32_t>* _counter);

	uint32_t TotalThreadsIncludingMainThread() const;

	ThreadScratchAllocator& ThreadAllocator() const;
	void ResetAllocators();

private:
	void WorkerLoop(uint32_t _threadId);

	bool TryRunOnePacket_NoLock();

	ThreadScratchAllocator* m_allocators = nullptr;

	kt::Thread* m_threads = nullptr;
	uint32_t m_numWorkers = 0;

	TaskPacket* m_packets = nullptr;
	
	uint32_t m_queueHead = 0;
	uint32_t m_queueTail = 0;
	std::atomic<uint32_t> m_numEntriesInQueue;

	// Todo: lock free
	kt::Mutex m_queueMutex;
	kt::Event m_queueSignal;

	std::atomic<uint32_t> m_keepRunning;

	std::atomic<uint32_t> m_numActiveWorkers;
};

}
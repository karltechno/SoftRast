#define WIN32_LEAN_AND_MEAN
#include <Windows.h>

#include "TaskSystem.h"
#include "kt/Strings.h"

#include "microprofile.h"

thread_local uint32_t tls_threadIndex;

namespace sr
{

uint32_t TaskSystem::TlsThreadIdx()
{
	return tls_threadIndex;
}

void TaskSystem::InitFromMainThread(uint32_t const _numWorkers)
{
	// main thread is 0

	m_packets = (TaskPacket*)kt::Malloc(sizeof(TaskPacket) * MAX_TASK_PACKETS);

	tls_threadIndex = 0;

	m_numWorkers = _numWorkers;

	// including main thread
	m_allocators = new ThreadScratchAllocator[TotalThreadsIncludingMainThread()];

	for (uint32_t i = 0; i < TotalThreadsIncludingMainThread(); ++i)
	{
		static const size_t MEM_SIZE = 128 * 1024 * 1024;

		ThreadScratchAllocator& alloc = m_allocators[i];
		// Todo: virtual mem kt wrapper
		void* ptr = ::VirtualAlloc(nullptr, MEM_SIZE, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
		KT_ASSERT(ptr);
		alloc.Init(ptr, MEM_SIZE);
	}

	if (!_numWorkers)
	{
		return;
	}

	m_threads = new kt::Thread[_numWorkers];
	
	std::atomic<uint32_t> initCounter{ _numWorkers };
	m_numWorkers = _numWorkers;

	struct ThreadInitData
	{
		std::atomic<uint32_t>* initCounter;
		uint32_t threadId;
		TaskSystem* sys;
		char const* name;
	};

	ThreadInitData* initData = (ThreadInitData*)KT_ALLOCA(sizeof(ThreadInitData) * _numWorkers);

	kt::String128* threadNames = (kt::String128*)KT_ALLOCA(sizeof(kt::String128) * _numWorkers);

	for (uint32_t i = 0; i < _numWorkers; ++i)
	{
		kt::Thread& t = m_threads[i];
		ThreadInitData& data = initData[i];
		threadNames[i].Clear();
		threadNames[i].AppendFmt("SoftRast Worker %u", i + 1);

		// main thread has 0
		data.threadId = i + 1;
		data.sys = this;
		data.initCounter = &initCounter;
		data.name = threadNames[i].Data();

		t.Run([](kt::Thread* _self) 
		{ 
			ThreadInitData* data = (ThreadInitData*)_self->GetUserData();
			TaskSystem* sys = data->sys;
			tls_threadIndex = data->threadId;
			MicroProfileOnThreadCreate(data->name);
			std::atomic_fetch_sub_explicit(data->initCounter, 1, std::memory_order_acquire);
			sys->WorkerLoop(data->threadId);
		}, 
		&data, threadNames[i].Data());
	}
	
	while(std::atomic_load(&initCounter) != 0) {}
}

void TaskSystem::WaitAndShutdown()
{
	std::atomic_store_explicit(&m_keepRunning, 0, std::memory_order_relaxed);

	for (uint32_t i = 0; i < m_numWorkers; ++i)
	{
		m_queueSignal.Signal();
	}

	for (uint32_t i = 0; i < m_numWorkers; ++i)
	{
		m_threads[i].Join();
	}

	ResetAllocators();

	delete[] m_threads;
	delete[] m_allocators;

	kt::Free(m_packets);

	m_threads = nullptr;
	m_packets = nullptr;
	m_allocators = nullptr;
}

void TaskSystem::PushTask(Task* _task)
{
	uint32_t totalTasks = 0;
	{
		kt::ScopedLock<kt::Mutex> lk(m_queueMutex);

		// split up task
		uint32_t lastEnd = 0;

		uint32_t tasksPushed = 0;

		totalTasks = (_task->m_totalPartitions + _task->m_granularity - 1) / _task->m_granularity;
		KT_ASSERT(totalTasks);

		if (_task->m_taskCounter)
		{
			std::atomic_fetch_add_explicit(_task->m_taskCounter, totalTasks, std::memory_order_acquire);
		}

		while (lastEnd < _task->m_totalPartitions)
		{
			++tasksPushed;
			if (m_numEntriesInQueue == MAX_TASK_PACKETS)
			{
				KT_ASSERT(false);
				return;
			}

			TaskPacket& p = m_packets[m_queueHead];
			m_queueHead = (m_queueHead + 1) & QUEUE_MASK;
			std::atomic_fetch_add_explicit(&m_numEntriesInQueue, 1, std::memory_order_acquire);
			p.m_task = _task;
			p.m_begin = lastEnd;
			lastEnd = kt::Min(lastEnd + _task->m_granularity, _task->m_totalPartitions);
			p.m_end = lastEnd;
		}
		KT_ASSERT(totalTasks == tasksPushed);
	}

	for (uint32_t i = 0; i < kt::Min(m_numWorkers, totalTasks); ++i)
	{
		m_queueSignal.Signal();
	}
}

void TaskSystem::SyncAndWaitForAll()
{
	KT_ASSERT(false); // ?
}

void TaskSystem::WaitForCounter(std::atomic<uint32_t>* _counter)
{
	for(;;)
	{
		if (std::atomic_load_explicit(_counter, std::memory_order_acquire) == 0)
		{
			return;
		}
		if (!TryRunOnePacket_NoLock())
		{
			break;
		}
	}

	MICROPROFILE_SCOPEI("TaskSystem", "IDLE WAIT FOR COUNTER", MP_RED);
	while (std::atomic_load_explicit(_counter, std::memory_order_acquire) > 0)
	{
		// dumb spin
		_mm_pause();
	}
}

uint32_t TaskSystem::TotalThreadsIncludingMainThread() const
{
	return m_numWorkers + 1;
}

ThreadScratchAllocator& TaskSystem::ThreadAllocator() const
{
	uint32_t const idx = tls_threadIndex;
	KT_ASSERT(idx < TotalThreadsIncludingMainThread());
	return m_allocators[idx];
}

void TaskSystem::ResetAllocators()
{
	for (uint32_t i = 0; i < TotalThreadsIncludingMainThread(); ++i)
	{
		m_allocators[i].Reset();
	}
}

void TaskSystem::WorkerLoop(uint32_t _threadId)
{
	while (std::atomic_load_explicit(&m_keepRunning, std::memory_order_acquire))
	{
		{
			MICROPROFILE_SCOPEI("TaskSystem", "IDLE", MP_RED);
			m_queueSignal.Wait();
		}

		std::atomic_fetch_add_explicit(&m_numActiveWorkers, 1, std::memory_order_acquire);

		while(TryRunOnePacket_NoLock()) {}

		std::atomic_fetch_sub_explicit(&m_numActiveWorkers, 1, std::memory_order_release);
	}
}


bool TaskSystem::TryRunOnePacket_NoLock()
{
	TaskPacket packet;
	bool popped = false;
	{
		kt::ScopedLock<kt::Mutex> mt(m_queueMutex);

		if (std::atomic_load_explicit(&m_numEntriesInQueue, std::memory_order_relaxed))
		{
			packet = m_packets[m_queueTail];
			m_queueTail = (m_queueTail + 1) & QUEUE_MASK;
			std::atomic_fetch_sub_explicit(&m_numEntriesInQueue, 1, std::memory_order_relaxed);
			popped = true;
		}
	}

	if (popped)
	{
		packet.m_task->m_fn(packet.m_task, tls_threadIndex, packet.m_begin, packet.m_end);
		if (packet.m_task->m_taskCounter)
		{
			std::atomic_fetch_sub_explicit(packet.m_task->m_taskCounter, 1, std::memory_order_release);
		}
	}

	return popped;
}

}
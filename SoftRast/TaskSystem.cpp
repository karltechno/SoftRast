#include "TaskSystem.h"
#include "kt/Strings.h"

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

	KT_ASSERT(_numWorkers);

	m_threads = (kt::Thread*)kt::Malloc(sizeof(kt::Thread) * _numWorkers);
	
	// Todo: we should make a kt ... new[]

	for (uint32_t i = 0; i < _numWorkers; ++i)
	{
		kt::PlacementNew(m_threads + i);
	}

	std::atomic<uint32_t> initCounter{ _numWorkers };
	m_numWorkers = _numWorkers;

	struct ThreadInitData
	{
		std::atomic<uint32_t>* initCounter;
		uint32_t threadId;
		TaskSystem* sys;
	};

	ThreadInitData* initData = (ThreadInitData*)KT_ALLOCA(sizeof(ThreadInitData) * _numWorkers);

	kt::String128* threadNames = (kt::String128*)KT_ALLOCA(sizeof(kt::String128) * _numWorkers);

	for (uint32_t i = 0; i < _numWorkers; ++i)
	{
		kt::Thread& t = m_threads[i];
		ThreadInitData& data = initData[i];
		threadNames[i].Clear();
		threadNames[i].AppendFmt("SoftRast Worker %u", i);

		// main thread has 0
		data.threadId = i + 1;
		data.sys = this;
		data.initCounter = &initCounter;

		t.Run([](kt::Thread* _self) 
		{ 
			ThreadInitData* data = (ThreadInitData*)_self->GetUserData();
			TaskSystem* sys = data->sys;
			tls_threadIndex = data->threadId;
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

	kt::Free(m_threads);
	kt::Free(m_packets);

	m_threads = nullptr;
	m_packets = nullptr;
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
	while (std::atomic_load_explicit(_counter, std::memory_order_acquire) > 0)
	{
		TaskPacket packet;

		// Todo: duplicate code
		bool popped = false;
		{
			kt::ScopedLock<kt::Mutex> lk(m_queueMutex);

			if (m_numEntriesInQueue)
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

			if (std::atomic_fetch_sub_explicit(&packet.m_task->m_numCompletedPartitions, 1, std::memory_order_release) == 0)
			{
				// do anything?
			}
		}
	}
}

uint32_t TaskSystem::TotalThreadsIncludingMainThread()
{
	return m_numWorkers + 1;
}

void TaskSystem::WorkerLoop(uint32_t _threadId)
{
	while (std::atomic_load_explicit(&m_keepRunning, std::memory_order_acquire))
	{
		m_queueSignal.Wait();

		std::atomic_fetch_add_explicit(&m_numActiveWorkers, 1, std::memory_order_acquire);

		for (;;)
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
				packet.m_task->m_fn(packet.m_task, _threadId, packet.m_begin, packet.m_end);
				if (packet.m_task->m_taskCounter)
				{
					std::atomic_fetch_sub_explicit(packet.m_task->m_taskCounter, 1, std::memory_order_release);
				}

				if (std::atomic_fetch_sub_explicit(&packet.m_task->m_numCompletedPartitions, 1, std::memory_order_release) == 0)
				{
					// do anything?
				}
			}
			else
			{
				break;
			}
		}
		std::atomic_fetch_sub_explicit(&m_numActiveWorkers, 1, std::memory_order_release);
	}
}




}
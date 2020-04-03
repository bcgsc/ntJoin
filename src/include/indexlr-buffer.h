#ifndef INDEXLR_BUFFER_H
#define INDEXLR_BUFFER_H

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

const size_t BLOCK_SIZE = 64;
const size_t BUFFER_SIZE = 16384;

// A block of data, so that threads don't work with one read at a time
template<typename T>
struct Block
{
	T data[BLOCK_SIZE];
	size_t dataCounter = 0;
	size_t num = 0;
};

// Surrounds pieces of data in the buffer with a busy mutex
// for exclusive access
template<typename T>
struct Slot
{
	T data;
	std::mutex busy;
	bool occupied = false;
	std::condition_variable occupancyChanged;
	int64_t lastTenant = -1; // Required to ensure read order
};

template<typename T>
class Buffer
{

  public:
	size_t elements() const { return elementCount; }

	void close()
	{
		closed = true;
		for (Slot<T>& slot : this->slots) {
			slot.occupancyChanged.notify_all();
		}
	}

	bool isClosed() const { return closed; }

  protected:
	std::vector<Slot<T>> slots{ BUFFER_SIZE };
	size_t readCounter = 0;
	std::atomic<size_t> elementCount{ 0 };
	std::atomic<bool> closed{ false };
};

template<typename T>
class InputBuffer : public Buffer<T>
{

  public:
	// For a direct write to a slot, instead of copying into it
	T& getWriteAccess(size_t num)
	{
		Slot<T>& target = this->slots[num / BLOCK_SIZE % BUFFER_SIZE];
		std::unique_lock<std::mutex> busyLock(target.busy);
		target.occupancyChanged.wait(busyLock, [&] { return !target.occupied; });
		busyLock.release();
		return target.data;
	}

	void releaseWriteAccess(size_t num)
	{
		Slot<T>& target = this->slots[num / BLOCK_SIZE % BUFFER_SIZE];
		target.occupied = true;
		target.occupancyChanged.notify_one();
		++(this->elementCount);
		target.busy.unlock();
	}

	void write(T& data)
	{
		size_t num = data.num / BLOCK_SIZE;
		Slot<T>& target = this->slots[num % BUFFER_SIZE];
		std::unique_lock<std::mutex> busyLock(target.busy);
		target.occupancyChanged.wait(busyLock, [&] { return !target.occupied; });
		target.data = std::move(data);
		target.occupied = true;
		target.occupancyChanged.notify_one();
		++(this->elementCount);
	}

	void read(T& data)
	{
		static std::mutex readMutex;
		std::unique_lock<std::mutex> readLock(readMutex);

		Slot<T>& target = this->slots[this->readCounter % BUFFER_SIZE];
		std::unique_lock<std::mutex> busyLock(target.busy);
		target.occupancyChanged.wait(busyLock, [&] { return target.occupied || this->closed; });
		if (this->closed) {
			return;
		}
		++(this->readCounter);

		readLock.unlock();

		data = std::move(target.data);
		target.occupied = false;
		target.occupancyChanged.notify_one();
		--(this->elementCount);
	}
};

template<typename T>
class OutputBuffer : public Buffer<T>
{

  public:
	void write(T& data)
	{
		size_t num = data.num / BLOCK_SIZE;
		Slot<T>& target = this->slots[num % BUFFER_SIZE];
		std::unique_lock<std::mutex> busyLock(target.busy);
		target.occupancyChanged.wait(
		    busyLock, [&] { return !target.occupied && (num - target.lastTenant <= BUFFER_SIZE); });
		target.data = std::move(data);
		target.occupied = true;
		target.lastTenant = num;
		target.occupancyChanged.notify_all();
		++(this->elementCount);
	}

	void read(T& data)
	{
		Slot<T>& target = this->slots[this->readCounter % BUFFER_SIZE];
		std::unique_lock<std::mutex> busyLock(target.busy);
		target.occupancyChanged.wait(busyLock, [&] { return target.occupied || this->closed; });
		if (this->closed) {
			return;
		}
		++(this->readCounter);
		data = std::move(target.data);
		target.occupied = false;
		target.occupancyChanged.notify_all();
		--(this->elementCount);
	}
};

#endif

// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#include "cas/allocator.h"
#include "cas/detail/config.h"
#include "cas/detail/memoryPool.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>

namespace Cas {
namespace {

/** @brief 内置后备池，与参考工程 defaultStorage 同为 64 MiB 静态区。 */
alignas(std::max_align_t) unsigned char defaultStorage[Detail::kDefaultPoolSize];

Detail::MemoryPool pool;
bool poolReady = false;

} // namespace

bool Allocator::initialize(void *storage, std::size_t capacity) {
    // 已就绪时不替换缓冲区：只接受"不指定"或"同一块缓冲区"的重复调用，与参考实现一致。
    if (poolReady) return storage == nullptr;
    if (!storage) {
        storage = defaultStorage;
        capacity = sizeof(defaultStorage);
    }
    poolReady = pool.initialize(storage, capacity);
    return poolReady;
}

void *Allocator::allocate(std::size_t bytes, std::size_t alignment) {
    if (!poolReady && !initialize()) return nullptr;
    return pool.allocate(bytes ? bytes : 1, alignment);
}

bool Allocator::release(void *pointer) {
    return !pointer || (poolReady && pool.release(pointer));
}

void *Allocator::reallocate(void *pointer, std::size_t bytes) {
    if (!poolReady && !initialize()) return nullptr;
    return pool.reallocate(pointer, bytes, alignof(std::max_align_t));
}

MemoryStats Allocator::statistics() {
    MemoryStats result;
    if (!poolReady) return result;
    const auto &s = pool.statistics();
    result.liveBytes = s.liveBytes;
    result.liveBlocks = s.liveBlocks;
    result.peakBytes = s.peakBytes;
    result.totalBlocks = s.totalBlocks;
    result.capacity = pool.capacity();
    return result;
}

bool Allocator::ready() { return poolReady; }

bool Allocator::quiescent() { return !poolReady || pool.liveBlocks() == 0; }

namespace Detail {

/** @brief 分配裸内存，空间不足转换为 std::bad_alloc。 */
void *poolMalloc(std::size_t bytes) {
    void *result = Allocator::allocate(bytes);
    if (!result) throw std::bad_alloc();
    return result;
}

/** @brief 释放裸内存。 */
void poolFree(void *pointer) { Allocator::release(pointer); }

/** @brief 池内 strdup。 */
char *poolStrdup(const char *text) {
    if (!text) return nullptr;
    const std::size_t size = std::strlen(text) + 1;
    auto *result = static_cast<char *>(Allocator::allocate(size, alignof(char)));
    if (result) std::memcpy(result, text, size);
    return result;
}

} // namespace Detail
} // namespace Cas

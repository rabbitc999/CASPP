// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#include "cas/detail/memoryPool.h"
#include "cas/detail/config.h"

#include <cstring>
#include <limits>

namespace Cas::Detail {
namespace {

/** @brief 向上对齐。 */
inline std::size_t alignUp(std::size_t value, std::size_t alignment) {
    return (value + alignment - 1) & ~(alignment - 1);
}

/** @brief 一个可拆分空块的最小尺寸：块头 + 偏移记录 + 一个对齐单位。 */
constexpr std::size_t kMinimumBlock = sizeof(PoolBlock) + sizeof(std::size_t) + MemoryPool::kAlignment;

} // namespace

bool MemoryPool::initialize(void *storage, std::size_t capacity) {
    m_stats = PoolStats{};
    m_freeHead = kNone;
    m_storage = nullptr;
    m_capacity = 0;
    m_rawStorage = storage;
    m_rawCapacity = capacity;
    if (!storage) return false;

    auto start = reinterpret_cast<std::uintptr_t>(storage);
    const auto alignment = static_cast<std::uintptr_t>(kAlignment);
    if (start > std::numeric_limits<std::uintptr_t>::max() - (alignment - 1)) return false;
    const auto aligned = (start + alignment - 1) & ~(alignment - 1);
    const auto skipped = static_cast<std::size_t>(aligned - start);
    if (capacity < skipped) return false;
    capacity = (capacity - skipped) & ~(kAlignment - 1);
    if (capacity < kMinimumBlock) return false;

    m_storage = static_cast<unsigned char *>(storage) + skipped;
    m_capacity = capacity;

    PoolBlock whole{};
    whole.size = capacity;
    addFree(0, whole);
    return true;
}

PoolBlock MemoryPool::loadBlock(std::size_t at) const {
    PoolBlock block{};
    if (at + sizeof(PoolBlock) <= m_capacity) std::memcpy(&block, m_storage + at, sizeof(block));
    return block;
}

void MemoryPool::storeBlock(std::size_t at, const PoolBlock &block) {
    if (at + sizeof(PoolBlock) <= m_capacity) std::memcpy(m_storage + at, &block, sizeof(block));
}

void MemoryPool::removeFree(std::size_t at, const PoolBlock &block) {
    if (block.previousFree == kNone) {
        m_freeHead = block.nextFree;
    } else {
        PoolBlock prev = loadBlock(block.previousFree);
        prev.nextFree = block.nextFree;
        storeBlock(block.previousFree, prev);
    }
    if (block.nextFree != kNone) {
        PoolBlock next = loadBlock(block.nextFree);
        next.previousFree = block.previousFree;
        storeBlock(block.nextFree, next);
    }
    (void)at;
}

void MemoryPool::addFree(std::size_t at, const PoolBlock &block) {
    PoolBlock head = block;
    head.previousFree = kNone;
    head.nextFree = m_freeHead;
    if (m_freeHead != kNone) {
        PoolBlock next = loadBlock(m_freeHead);
        next.previousFree = at;
        storeBlock(m_freeHead, next);
    }
    m_freeHead = at;
    storeBlock(at, head);
}

bool MemoryPool::owns(const void *pointer) const {
    if (!m_storage || !pointer) return false;
    const auto address = reinterpret_cast<std::uintptr_t>(pointer);
    const auto start = reinterpret_cast<std::uintptr_t>(m_storage);
    return address >= start && address - start < m_capacity;
}

void *MemoryPool::allocate(std::size_t bytes, std::size_t alignment) {
    if (!m_storage || !memoryPoolValidAlignment(alignment)) return nullptr;
    if (bytes == 0) bytes = 1;

    for (std::size_t at = m_freeHead; at != kNone;) {
        const PoolBlock block = loadBlock(at);
        const auto base = reinterpret_cast<std::uintptr_t>(m_storage + at);
        if (base > std::numeric_limits<std::uintptr_t>::max() - sizeof(PoolBlock) - sizeof(std::size_t)) return nullptr;
        if (base + sizeof(PoolBlock) + sizeof(std::size_t) > std::numeric_limits<std::uintptr_t>::max() - (alignment - 1)) return nullptr;
        const auto address = (base + sizeof(PoolBlock) + sizeof(std::size_t) + alignment - 1) & ~static_cast<std::uintptr_t>(alignment - 1);
        const auto payload = static_cast<std::size_t>(address - base);
        if (payload > std::numeric_limits<std::size_t>::max() - bytes) return nullptr;
        if (payload + bytes > std::numeric_limits<std::size_t>::max() - (kAlignment - 1)) return nullptr;
        const std::size_t wanted = alignUp(payload + bytes, kAlignment);
        if (wanted <= block.size) {
            removeFree(at, block);
            PoolBlock used = block;
            // 剩余部分不足以构成独立块时，把整块拨给调用者；否则块链会出现长度缺口。
            const bool splittable = block.size - wanted >= kMinimumBlock;
            used.size = splittable ? wanted : block.size;
            if (splittable) {
                PoolBlock rest{};
                rest.size = block.size - wanted;
                rest.previousSize = used.size;
                addFree(at + used.size, rest);
                const std::size_t next = at + used.size + rest.size;
                if (next < m_capacity) {
                    PoolBlock tail = loadBlock(next);
                    tail.previousSize = rest.size;
                    storeBlock(next, tail);
                }
            }
            used.payload = payload;
            used.requested = bytes;
            storeBlock(at, used);
            // 链接区紧贴负载起点之前：写入的是块起始的**绝对**偏移，lookup 才能反查回块头。
            const std::size_t linkAt = at;
            std::memcpy(m_storage + at + payload - sizeof(std::size_t), &linkAt, sizeof(linkAt));
            m_stats.liveBytes += used.size;
            ++m_stats.liveBlocks;
            ++m_stats.totalBlocks;
            if (m_stats.liveBytes > m_stats.peakBytes) m_stats.peakBytes = m_stats.liveBytes;
            return reinterpret_cast<void *>(address);
        }
        at = block.nextFree;
    }
    return nullptr;
}

bool MemoryPool::lookup(const void *pointer, std::size_t *at, PoolBlock *block) const {
    if (!owns(pointer)) return false;
    const auto offset = reinterpret_cast<std::uintptr_t>(pointer) - reinterpret_cast<std::uintptr_t>(m_storage);
    if (offset < sizeof(PoolBlock) + sizeof(std::size_t)) return false;
    std::size_t start = 0;
    std::memcpy(&start, m_storage + offset - sizeof(std::size_t), sizeof(start));
    if (start % kAlignment || start > m_capacity - sizeof(PoolBlock)) return false;
    const PoolBlock found = loadBlock(start);
    if (!found.payload || found.size < kMinimumBlock || found.size > m_capacity - start) return false;
    if (found.payload >= found.size) return false;
    if (m_storage + start + found.payload != pointer) return false;
    *at = start;
    *block = found;
    return true;
}

bool MemoryPool::release(void *pointer) {
    std::size_t at = 0;
    PoolBlock block{};
    if (!lookup(pointer, &at, &block)) return false;
    m_stats.liveBytes -= block.size;
    --m_stats.liveBlocks;
    block.payload = 0;
    block.requested = 0;
    storeBlock(at, block);

    std::size_t next = at + block.size;
    if (next < m_capacity) {
        const PoolBlock following = loadBlock(next);
        if (!following.payload) {
            removeFree(next, following);
            block.size += following.size;
        }
    }
    if (block.previousSize && block.previousSize <= at) {
        const std::size_t previousAt = at - block.previousSize;
        const PoolBlock previous = loadBlock(previousAt);
        if (!previous.payload) {
            removeFree(previousAt, previous);
            PoolBlock merged = previous;
            merged.size += block.size;
            block = merged;
            at = previousAt;
        }
    }
    next = at + block.size;
    if (next < m_capacity) {
        PoolBlock tail = loadBlock(next);
        tail.previousSize = block.size;
        storeBlock(next, tail);
    }
    addFree(at, block);
    return true;
}

void *MemoryPool::reallocate(void *pointer, std::size_t bytes, std::size_t alignment) {
    if (!pointer) return allocate(bytes, alignment);
    std::size_t at = 0;
    PoolBlock block{};
    if (!lookup(pointer, &at, &block)) return nullptr;
    if (block.requested > block.size - block.payload) return nullptr;
    if (bytes == 0) {
        release(pointer);
        return nullptr;
    }
    void *result = allocate(bytes, alignment);
    if (!result) return nullptr;
    std::memcpy(result, pointer, bytes < block.requested ? bytes : block.requested);
    release(pointer);
    return result;
}

} // namespace Cas::Detail

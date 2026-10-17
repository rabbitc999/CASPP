// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#ifndef CAS_ALLOCATOR_H
#define CAS_ALLOCATOR_H

/**
 * @file allocator.h
 * @brief 进程级自定义分配池，接口与参考工程 Cas::Allocator 对齐。
 */

#include <cstddef>

namespace Cas {

/** @brief 自定义内存池的占用统计。 */
struct MemoryStats {
    std::size_t liveBytes = 0;   ///< 活动块占用字节（含块头、对齐填充与块内碎片）。
    std::size_t liveBlocks = 0;  ///< 活动块数。
    std::size_t peakBytes = 0;   ///< liveBytes 历史峰值。
    std::size_t totalBlocks = 0; ///< 累计分配次数。
    std::size_t capacity = 0;    ///< 池容量。
};

/** @brief 自定义分配池适配层，支持外部缓冲区与对齐分配。 */
class Allocator {
public:
    /** @brief 首次初始化内存池；已初始化时拒绝替换缓冲区。
     * @param storage 外部缓冲区，空指针使用内置后备池。
     * @param capacity 外部缓冲区容量；使用内置池时忽略该参数。
     * @return 初始化成功（或已用同一缓冲区初始化过）返回 true。 */
    static bool initialize(void *storage = nullptr, std::size_t capacity = 0);
    /** @brief 从池中分配内存。
     * @param bytes 字节数。
     * @param alignment 对齐要求（2 的幂）。
     * @return 池内地址，空间不足返回 nullptr。 */
    static void *allocate(std::size_t bytes, std::size_t alignment = alignof(std::max_align_t));
    /** @brief 释放池中内存。
     * @param pointer 池内地址或 nullptr。
     * @return 地址属于本池且释放成功返回 true。 */
    static bool release(void *pointer);
    /** @brief 调整池内分配大小，失败保留旧块。
     * @param pointer 旧分配。
     * @param bytes 新字节数。
     * @return 新地址或 nullptr。 */
    static void *reallocate(void *pointer, std::size_t bytes);
    /** @brief 获取内存占用统计。
     * @return 当前池统计。 */
    static MemoryStats statistics();
    /** @brief 池是否已就绪。 */
    static bool ready();
    /** @brief 池内是否有活动块；会话结束前应清零。 */
    static bool quiescent();
};

} // namespace Cas

#endif // CAS_ALLOCATOR_H

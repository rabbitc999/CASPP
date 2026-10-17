#ifndef CAS_DETAIL_MEMORYPOOL_H
#define CAS_DETAIL_MEMORYPOOL_H

/**
 * @file memoryPool.h
 * @brief 首次适配链表内存池。
 *
 * 语义移植自参考工程 ParseEngine/library/allocator/pool.c：
 *  - 存储完全由调用者提供（默认 64 MiB 内部后备缓冲区），元数据内嵌在池内，不依赖宿主分配器；
 *  - 块头记录 size / payload / requested / 前后空闲链，分配时按首次适配拆分，释放时双向合并；
 *  - 返回地址前 8 字节保存块偏移，释放与扩容均可反查块头，并校验对齐、归属与重复释放；
 *  - liveBytes 统计包含块头、对齐填充与块内碎片，与参考库口径一致。
 *
 * 与参考库的差异：改为 C++17 实现，块头直接以结构体形式驻留池内（无 memcpy 搬运）；
 * 不支持线性池模式（本工程只用链表池）。
 */

#include <cstddef>
#include <cstdint>

namespace Cas::Detail {

/** @brief 池内块头。紧随其后是 size_t 的记录偏移与对齐填充，然后是负载。 */
struct PoolBlock {
    std::size_t size;         ///< 整块占用字节数（含块头与对齐填充）。
    std::size_t previousSize; ///< 前一块占用字节数，0 表示池首。
    std::size_t payload;      ///< 块首到负载起点的字节数，0 表示空闲块。
    std::size_t requested;    ///< 调用者请求的字节数。
    std::size_t nextFree;     ///< 下一个空闲块偏移，kNone 表示链尾。
    std::size_t previousFree; ///< 上一个空闲块偏移，kNone 表示链首。
};

/** @brief 池统计。 */
struct PoolStats {
    std::size_t liveBytes = 0;   ///< 活动块占用（含全部管理开销）。
    std::size_t liveBlocks = 0;  ///< 活动块数。
    std::size_t peakBytes = 0;   ///< liveBytes 历史峰值。
    std::size_t totalBlocks = 0; ///< 累计分配次数。
};

/** @brief 固定缓冲区上的首次适配内存池。 */
class MemoryPool {
public:
    /** @brief 空池；需再调用 initialize。 */
    MemoryPool() = default;
    /** @brief 连接外部缓冲区并初始化。 @param storage 缓冲区。 @param capacity 字节数。 */
    MemoryPool(void *storage, std::size_t capacity) { initialize(storage, capacity); }
    MemoryPool(const MemoryPool &) = delete;
    MemoryPool &operator=(const MemoryPool &) = delete;

    /** @brief 连接外部缓冲区。地址会先按指针宽度对齐，容量按对齐向下取整。
     * @param storage 缓冲区首址，空指针返回 false。
     * @param capacity 字节数。
     * @return 容量足够返回 true。 */
    bool initialize(void *storage, std::size_t capacity);
    /** @brief 分配内存。 @param bytes 字节数（0 按 1 处理）。 @param alignment 必须是 2 的幂。 @return 地址或空指针。 */
    void *allocate(std::size_t bytes, std::size_t alignment = kAlignment);
    /** @brief 归还内存。 @param pointer 分配地址或空指针。 @return 有效分配返回 true。 */
    bool release(void *pointer);
    /** @brief 扩容或缩容，内容按较小者复制。失败时旧块保持不变，0 字节等价于释放。
     * @return 新地址或空指针。 */
    void *reallocate(void *pointer, std::size_t bytes, std::size_t alignment = kAlignment);
    /** @brief 判断地址是否落在池范围内（不保证是分配起点）。 */
    bool owns(const void *pointer) const;
    /** @brief 获取统计。 */
    const PoolStats &statistics() const { return m_stats; }
    /** @brief 池容量（对齐后）。 */
    std::size_t capacity() const { return m_capacity; }
    /** @brief 是否已连接存储。 */
    bool ready() const { return m_storage != nullptr; }
    /** @brief 释放活动块数，用于会话复用前的一致性检查。 */
    std::size_t liveBlocks() const { return m_stats.liveBlocks; }

    /** @brief 指针宽度对齐量，与参考库 MEMORY_POOL_ALIGNMENT 一致。 */
    static constexpr std::size_t kAlignment = sizeof(void *);
    /** @brief 空链哨兵。 */
    static constexpr std::size_t kNone = static_cast<std::size_t>(-1);

private:
    PoolBlock loadBlock(std::size_t at) const;
    void storeBlock(std::size_t at, const PoolBlock &block);
    void removeFree(std::size_t at, const PoolBlock &block);
    void addFree(std::size_t at, const PoolBlock &block);
    bool lookup(const void *pointer, std::size_t *at, PoolBlock *block) const;

    unsigned char *m_storage = nullptr;   ///< 对齐后的池起点。
    std::size_t m_capacity = 0;           ///< 对齐后的池容量。
    void *m_rawStorage = nullptr;         ///< 调用者缓冲区原始地址。
    std::size_t m_rawCapacity = 0;        ///< 调用者声明的容量。
    std::size_t m_freeHead = kNone;       ///< 空闲链首偏移。
    PoolStats m_stats{};
};

/** @brief 查询对齐量是否合法（非零且为 2 的幂）。 */
constexpr bool memoryPoolValidAlignment(std::size_t alignment) {
    return alignment != 0 && (alignment & (alignment - 1)) == 0;
}

} // namespace Cas::Detail

#endif // CAS_DETAIL_MEMORYPOOL_H

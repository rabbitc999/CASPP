#ifndef CAS_DETAIL_CONFIG_H
#define CAS_DETAIL_CONFIG_H

/**
 * @file config.h
 * @brief 引擎与内存池的可配置参数。
 *
 * 参考工程的 config.h 把引擎、内存池、大数三类参数集中定义；本工程沿用这一做法，
 * 所有需要按目标平台调整的量都在这里，改完需整体重编译。
 */

#include <cstddef>
#include <cstdint>

namespace Cas::Detail {

/* ============================ 内存池参数 ============================ */

/** @brief 默认池容量（字节）。参考工程 Cas::Allocator 的后备池同为 64 MiB。 */
#ifndef CAS_POOL_SIZE
#define CAS_POOL_SIZE (64u * 1024u * 1024u)
#endif

/** @brief 内部后备缓冲区容量。 */
inline constexpr std::size_t kDefaultPoolSize = CAS_POOL_SIZE;

/** @brief 池粒度提示，仅用于文档与统计显示，块头实际按 8 字节对齐。 */
inline constexpr std::size_t kPoolBlockGranularity = 8;

/* ============================ 表达式树参数 ============================ */

/** @brief 规范化后 Add/Mul 最多保留的操作数个数，超出则打包嵌套（避免 n 元平铺的哈希开销）。 */
inline constexpr int kMaxNAryOperands = 64;

/** @brief 展开运算的操作数上限，防止 (x+y)^n 的规模失控。 */
inline constexpr long long kExpandTermLimit = 20000;

/** @brief 递归深度上限，超出抛出 RecursionError。 */
inline constexpr int kMaxRecursionDepth = 1024;

/** @brief 符号哈希表初始桶数。 */
inline constexpr std::size_t kSymbolTableBuckets = 4096;

/* ============================== 大整数 ============================== */

/** @brief 大整数转浮点时的有效位数（十进制）。 */
inline constexpr int kDefaultDigits = 15;

/** @brief NTT 模数（支持 2^24 长度的变换）。 */
inline constexpr std::uint32_t kNttPrime = 2013265921u;
/** @brief 上述模数的原根。 */
inline constexpr std::uint32_t kNttRoot = 31u;

/** @brief 低于该值（绝对值）的整数用 64 位内联存储，不分配大数缓冲区。 */
inline constexpr long long kSmallIntLimit = (1LL << 61);

} // namespace Cas::Detail

#endif // CAS_DETAIL_CONFIG_H

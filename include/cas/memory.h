#ifndef CAS_MEMORY_H
#define CAS_MEMORY_H

/**
 * @file memory.h
 * @brief 把标准容器接到 CAS 池上的适配层，对应参考工程的 containerAllocator.hpp。
 *
 * 表达式树的所有容器（Add/Mul/Pow 的操作数数组、符号名、输出缓冲区）都使用这里的别名，
 * 因此 memoryStats() 统计到的就是整个引擎的真实占用。
 */

#include "allocator.h"
#include "detail/config.h"
#include "detail/memoryPool.h"

#include <cstddef>
#include <deque>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <new>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace Cas::Memory {

/** @brief 记录分配归属的容器分配器；对象存活期间即使切换池也能正确归还。 */
template <class T>
class PoolAllocator {
public:
    using value_type = T;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using propagate_on_container_move_assignment = std::true_type;
    using propagate_on_container_copy_assignment = std::false_type;
    using propagate_on_container_swap = std::true_type;
    using is_always_equal = std::true_type;

    PoolAllocator() noexcept = default;
    template <class U>
    PoolAllocator(const PoolAllocator<U> &) noexcept {}

    /** @brief 从进程级 CAS 池分配；池未初始化时退回标准分配。 */
    T *allocate(std::size_t count) {
        if (count > static_cast<std::size_t>(-1) / sizeof(T)) throw std::bad_alloc();
        void *raw = Allocator::ready() ? Allocator::allocate(count * sizeof(T), alignof(T)) : nullptr;
        if (!raw) {
            // 池未就绪或容量不足：交给标准分配器。释放时按同一判据归还，配对不会错。
            raw = ::operator new(count * sizeof(T));
        }
        return static_cast<T *>(raw);
    }
    /** @brief 归还：池内地址回池，其余交给标准分配器。 */
    void deallocate(T *pointer, std::size_t) noexcept {
        if (!pointer) return;
        if (!Allocator::release(pointer)) ::operator delete(static_cast<void *>(pointer));
    }

    template <class U>
    bool operator==(const PoolAllocator<U> &) const noexcept { return true; }
    template <class U>
    bool operator!=(const PoolAllocator<U> &) const noexcept { return false; }
};

template <class T>
using Vector = std::vector<T, PoolAllocator<T>>;
template <class T>
using List = std::list<T, PoolAllocator<T>>;
template <class T>
using Deque = std::deque<T, PoolAllocator<T>>;
template <class K, class V, class C = std::less<K>>
using Map = std::map<K, V, C, PoolAllocator<std::pair<const K, V>>>;
template <class K, class C = std::less<K>>
using Set = std::set<K, C, PoolAllocator<K>>;
template <class K, class V, class H = std::hash<K>, class E = std::equal_to<K>>
using UnorderedMap = std::unordered_map<K, V, H, E, PoolAllocator<std::pair<const K, V>>>;
template <class K, class H = std::hash<K>, class E = std::equal_to<K>>
using UnorderedSet = std::unordered_set<K, H, E, PoolAllocator<K>>;

/** @brief 池内字符串。 */
using String = std::basic_string<char, std::char_traits<char>, PoolAllocator<char>>;

/** @brief 池内字符串流。 */
using StringStream = std::basic_stringstream<char, std::char_traits<char>, PoolAllocator<char>>;

/** @brief 由 std::string 转池内字符串。 */
inline String toPoolString(const std::string &text) { return String(text.begin(), text.end()); }

/** @brief 由池内字符串转 std::string。 */
inline std::string toStdString(const String &text) { return std::string(text.begin(), text.end()); }

/** @brief 以池分配构造共享对象。 */
template <class T, class... Args>
std::shared_ptr<T> makeShared(Args &&...args) {
    return std::allocate_shared<T>(PoolAllocator<T>(), std::forward<Args>(args)...);
}

} // namespace Cas::Memory

#endif // CAS_MEMORY_H

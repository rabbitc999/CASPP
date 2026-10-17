#ifndef CAS_DETAIL_NODE_H
#define CAS_DETAIL_NODE_H

/**
 * @file node.h
 * @brief 表达式树的节点定义。
 *
 * 参考工程用 mark-sweep GC + 全局求值栈管理 Atom 生命周期，其 doc/refactor-plan.md
 * 记录了栈弹出后裸指针悬垂的 4 类真实缺陷。本工程改为**引用计数 + 不可变节点**：
 * 每个 Node 一经创建就不再修改，任何变换都产生新节点，因此不存在"被别处改坏"的
 * 中间状态；引用计数让释放时机完全确定，不需要 stop-the-world 扫描。
 */

#include "cas/detail/bigInt.h"
#include "cas/detail/config.h"
#include "cas/memory.h"

#include <cstddef>
#include <cstdint>

namespace Cas::Detail {

/** @brief 节点种类。 */
enum class NodeKind : std::uint8_t {
    Integer,   ///< 任意精度整数
    Rational,  ///< 有理数 p/q，q > 0，已约分
    Symbol,    ///< 符号（变量名驻留在全局符号表）
    Add,       ///< n 元加法
    Mul,       ///< n 元乘法
    Pow,       ///< 二元幂 base^exponent
    Function,  ///< 命名函数调用（sin/cos/...）
    Constant,  ///< 命名常数（Pi/E/Infinity）
};

/** @brief 常数符号编号。 */
enum class ConstantId : std::uint8_t {
    Pi = 0,
    E = 1,
    EulerGamma = 2,
    Infinity = 3,
    NegativeInfinity = 4,
};

/** @brief 节点内部的子节点引用计数句柄。 */
class Node;
using NodePtr = Node *;

/**
 * @brief 表达式节点。
 *
 * 采用侵入式引用计数：`refs` 记录持有者数量，归零时整棵子树被释放回池。
 * 操作数数组按需分配（Add/Mul 用 operands，Pow 用 operands[0..1]，
 * Function 用 operands 表示实参，Symbol/Integer/Rational/Constant 不用）。
 */
class Node {
public:
    /** @brief 创建整数节点。 */
    static NodePtr makeInteger(const BigInt &value);
    /** @brief 创建有理数节点，自动约分并规范符号。 */
    static NodePtr makeRational(const BigInt &numerator, const BigInt &denominator);
    /** @brief 创建符号节点。 */
    static NodePtr makeSymbol(const std::string &name);
    /** @brief 创建符号节点（已有驻留编号）。 */
    static NodePtr makeSymbolById(std::uint32_t id);
    /** @brief 创建 n 元节点。 */
    static NodePtr makeNary(NodeKind kind, const Memory::Vector<NodePtr> &operands);
    /** @brief 创建二元幂节点。 */
    static NodePtr makePower(NodePtr base, NodePtr exponent);
    /** @brief 创建函数调用节点。 */
    static NodePtr makeFunction(const std::string &name, const Memory::Vector<NodePtr> &arguments);
    /** @brief 创建常数节点。 */
    static NodePtr makeConstant(ConstantId id);

    /** @brief 增加引用计数。 */
    static void retain(NodePtr node);
    /** @brief 减少引用计数，归零时销毁。 */
    static void release(NodePtr node);

    NodeKind kind() const { return m_kind; }
    /** @brief 操作数个数。 */
    std::size_t size() const { return m_operands.size(); }
    /** @brief 第 i 个操作数。 */
    NodePtr operand(std::size_t index) const { return m_operands[index]; }
    const Memory::Vector<NodePtr> &operands() const { return m_operands; }

    /** @brief 整数/有理数的分子。 */
    const BigInt &numerator() const { return m_numerator; }
    /** @brief 有理数的分母。 */
    const BigInt &denominator() const { return m_denominator; }
    /** @brief 符号名或函数名（驻留字符串指针）。 */
    const char *name() const { return m_name; }
    /** @brief 符号驻留编号。 */
    std::uint32_t symbolId() const { return m_symbolId; }
    /** @brief 常数编号。 */
    ConstantId constant() const { return m_constant; }

    /** @brief 小整数快速取值，仅在 isSmallInteger 为真时有效。 */
    bool isSmallInteger(long long *out = nullptr) const;
    /** @brief 是否为整数值（Integer 或分母为 1 的 Rational）。 */
    bool isInteger() const;
    /** @brief 是否为数值（Integer 或 Rational）。 */
    bool isNumber() const { return m_kind == NodeKind::Integer || m_kind == NodeKind::Rational; }
    /** @brief 是否为精确的 0 / 1 / -1。 */
    bool isZero() const;
    bool isOne() const;
    bool isNegativeOne() const;
    /** @brief 数值是否为负。 */
    bool isNegative() const;

    /** @brief 结构相等比较（不含自动化简）。 */
    bool equals(const NodePtr &other) const;
    /** @brief 建立哈希值（缓存）。 */
    std::uint64_t hash() const;

    /** @brief 排序键：用于把 Add/Mul 的操作数规范化排序。 */
    int compareTo(const NodePtr &other) const;

private:
    Node() = default;
    ~Node();
    Node(const Node &) = delete;
    Node &operator=(const Node &) = delete;

    NodeKind m_kind = NodeKind::Integer;
    /**
     * @brief 引用计数。
     *
     * 节点创建时即为 1，代表"创建者持有的那个引用"；随后必须立即交给 Ptr 接管，
     * 由 Ptr 负责在析构时 release。这样 makeXxx 与 Ptr 之间没有中间态。
     */
    std::uint32_t m_refs = 1;
    std::uint32_t m_symbolId = 0;
    ConstantId m_constant = ConstantId::Pi;
    const char *m_name = nullptr;
    BigInt m_numerator;
    BigInt m_denominator;
    Memory::Vector<NodePtr> m_operands;
    mutable std::uint64_t m_hash = 0;

    friend class NodeArena;
};

/* ------------------------------------------------------------------ 符号表 */

/** @brief 把名字驻留到全局符号表，返回编号；相同名字总是同一编号。 */
std::uint32_t internSymbol(const std::string &name);
/** @brief 由编号取回名字。 */
const char *symbolName(std::uint32_t id);
/** @brief 取全部已驻留符号（按编号顺序）。 */
Memory::Vector<std::uint32_t> allSymbols();
/** @brief 在符号表里按名字查找编号，未找到返回 UINT32_MAX。 */
std::uint32_t lookupSymbol(const std::string &name);

/* --------------------------------------------------------------- 引用工具 */

/**
 * @brief 引用计数守卫：让"创建后必须 release"的约定由 RAII 保证。
 *
 * 变换函数一律返回新节点（`Ptr`），调用方用完自动释放，避免了手工配对出错。
 */
class Ptr {
public:
    Ptr() = default;
    /**
     * @brief 接管一个**新产生的**引用（makeXxx 的返回值）。
     *
     * 契约：node 的引用计数里已经包含"这一次"所有权，Ptr 析构时消费它。
     * 不要传一个别处还持有着的裸指针进来，否则会双重释放。
     */
    explicit Ptr(NodePtr node) : m_node(node) {}
    /** @brief 借用一颗由别处持有的子树：额外 retain 一份，形成独立持有者。 */
    static Ptr borrow(NodePtr node) {
        Node::retain(node);
        return Ptr(node);
    }
    Ptr(const Ptr &other) : m_node(other.m_node) { Node::retain(m_node); }
    Ptr(Ptr &&other) noexcept : m_node(other.m_node) { other.m_node = nullptr; }
    Ptr &operator=(const Ptr &other) {
        if (this != &other) {
            Node::retain(other.m_node);
            Node::release(m_node);
            m_node = other.m_node;
        }
        return *this;
    }
    Ptr &operator=(Ptr &&other) noexcept {
        if (this != &other) {
            Node::release(m_node);
            m_node = other.m_node;
            other.m_node = nullptr;
        }
        return *this;
    }
    ~Ptr() { Node::release(m_node); }

    /** @brief 复制出一个新的持有者。 */
    Ptr share() const { return Ptr(*this); }
    /** @brief 转移所有权给裸指针调用方（不再释放）。 */
    NodePtr detach() {
        NodePtr result = m_node;
        m_node = nullptr;
        return result;
    }
    NodePtr get() const { return m_node; }
    NodePtr operator->() const { return m_node; }
    explicit operator bool() const { return m_node != nullptr; }
    NodePtr &ref() { return m_node; }

private:
    NodePtr m_node = nullptr;
};

} // namespace Cas::Detail

#endif // CAS_DETAIL_NODE_H

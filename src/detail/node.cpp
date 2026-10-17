// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#include "cas/detail/node.h"
#include "cas/allocator.h"

#include <algorithm>
#include <cstring>
#include <memory>
#include <mutex>
#include <new>
#include <string>
#include <unordered_map>

namespace Cas::Detail {
namespace {

/* ------------------------------------------------------------- 节点内存池 */

/**
 * @brief 节点分配走 CAS 内存池；池未初始化时退化为普通堆分配。
 *
 * 退化路径刻意用 malloc/free 而不是 operator new/delete：
 * 分配与释放必须是同一对函数，混用 nothrow new 与普通 delete 会破坏堆。
 */
void *nodeAllocate(std::size_t bytes) {
    if (!Allocator::ready()) return std::malloc(bytes);
    return Allocator::allocate(bytes, alignof(std::max_align_t));
}

void nodeFree(void *pointer) {
    if (!pointer) return;
    if (!Allocator::release(pointer)) std::free(pointer);
}

/* --------------------------------------------------------------- 符号表 */

struct SymbolTable {
    std::mutex guard;
    /**
     * @brief 符号名存储。
     *
     * 用 unique_ptr 逐个持有字符串，保证扩容时已有名字的地址不变——
     * Node 里缓存了 name 指针，若用 vector<string> 扩容会让所有缓存指针悬垂。
     */
    Memory::Vector<std::unique_ptr<std::string>> names;
    std::unordered_map<std::string, std::uint32_t> index;
};

SymbolTable &symbolTable() {
    // 故意泄漏：符号表要活到进程结束，避免静态析构顺序问题。
    static SymbolTable *table = new SymbolTable();
    return *table;
}

/* ------------------------------------------------------- 规范化排序键 */

/**
 * @brief 种类的排序权重。
 *
 * 这是**结构比较**用的基准顺序，不承担"降幂"或"系数在前"这类排版意图：
 * 那两种需求分别由 simplify.cpp 里的 AddOrder / MulOrder 处理。
 * 这里只保证同构表达式在任何构建路径下都得到同一顺序，从而结构相等可判定。
 */
int kindRank(NodeKind kind) {
    switch (kind) {
    case NodeKind::Integer:
    case NodeKind::Rational: return 0;
    case NodeKind::Symbol: return 1;
    case NodeKind::Constant: return 2;
    case NodeKind::Pow: return 3;
    case NodeKind::Function: return 4;
    case NodeKind::Mul: return 5;
    case NodeKind::Add: return 6;
    }
    return 7;
}

/** @brief 数值比较：先比大小关系，保证排序稳定。 */
int compareNumbers(const Node *left, const Node *right) {
    const BigInt l = left->numerator() * right->denominator();
    const BigInt r = right->numerator() * left->denominator();
    return l.compare(r);
}

/** @brief 字符串比较，用于符号名与函数名。 */
int compareText(const char *left, const char *right) {
    if (left == right) return 0;
    if (!left) return -1;
    if (!right) return 1;
    const int result = std::strcmp(left, right);
    return result < 0 ? -1 : (result > 0 ? 1 : 0);
}

} // namespace

/* --------------------------------------------------------------- 构造 */

Node::~Node() {
    for (NodePtr child : m_operands) Node::release(child);
}

NodePtr Node::makeInteger(const BigInt &value) {
    void *raw = nodeAllocate(sizeof(Node));
    if (!raw) return nullptr;
    NodePtr node = new (raw) Node();
    node->m_kind = NodeKind::Integer;
    node->m_numerator = value;
    node->m_denominator = BigInt(1);
    return node;
}

NodePtr Node::makeRational(const BigInt &numerator, const BigInt &denominator) {
    if (denominator.isZero()) return nullptr;
    BigInt p = numerator, q = denominator;
    if (q.isNegative()) {
        p = p.negated();
        q = q.negated();
    }
    const BigInt divisor = BigInt::gcd(p.abs(), q);
    if (!divisor.isOne() && !divisor.isZero()) {
        p = p / divisor;
        q = q / divisor;
    }
    if (q.isOne()) return makeInteger(p);

    void *raw = nodeAllocate(sizeof(Node));
    if (!raw) return nullptr;
    NodePtr node = new (raw) Node();
    node->m_kind = NodeKind::Rational;
    node->m_numerator = p;
    node->m_denominator = q;
    return node;
}

NodePtr Node::makeSymbol(const std::string &name) { return makeSymbolById(internSymbol(name)); }

NodePtr Node::makeSymbolById(std::uint32_t id) {
    void *raw = nodeAllocate(sizeof(Node));
    if (!raw) return nullptr;
    NodePtr node = new (raw) Node();
    node->m_kind = NodeKind::Symbol;
    node->m_symbolId = id;
    node->m_name = symbolName(id);
    node->m_denominator = BigInt(1);
    return node;
}

NodePtr Node::makeNary(NodeKind kind, const Memory::Vector<NodePtr> &operands) {
    void *raw = nodeAllocate(sizeof(Node));
    if (!raw) return nullptr;
    NodePtr node = new (raw) Node();
    node->m_kind = kind;
    node->m_denominator = BigInt(1);
    node->m_operands = operands;
    for (NodePtr child : node->m_operands) Node::retain(child);
    return node;
}

NodePtr Node::makePower(NodePtr base, NodePtr exponent) {
    Memory::Vector<NodePtr> operands;
    operands.push_back(base);
    operands.push_back(exponent);
    return makeNary(NodeKind::Pow, operands);
}

NodePtr Node::makeFunction(const std::string &name, const Memory::Vector<NodePtr> &arguments) {
    void *raw = nodeAllocate(sizeof(Node));
    if (!raw) return nullptr;
    NodePtr node = new (raw) Node();
    node->m_kind = NodeKind::Function;
    node->m_symbolId = internSymbol(name);
    node->m_name = symbolName(node->m_symbolId);
    node->m_denominator = BigInt(1);
    node->m_operands = arguments;
    for (NodePtr child : node->m_operands) Node::retain(child);
    return node;
}

NodePtr Node::makeConstant(ConstantId id) {
    void *raw = nodeAllocate(sizeof(Node));
    if (!raw) return nullptr;
    NodePtr node = new (raw) Node();
    node->m_kind = NodeKind::Constant;
    node->m_constant = id;
    node->m_denominator = BigInt(1);
    return node;
}

/* ---------------------------------------------------------- 引用计数 */

void Node::retain(NodePtr node) {
    if (node) ++node->m_refs;
}

void Node::release(NodePtr node) {
    if (!node) return;
    if (--node->m_refs) return;
    node->~Node();
    nodeFree(node);
}

/* ------------------------------------------------------------- 查询 */

bool Node::isSmallInteger(long long *out) const {
    if (m_kind != NodeKind::Integer) return false;
    if (!m_numerator.fitsLongLong()) return false;
    if (out) *out = m_numerator.toLongLong();
    return true;
}

bool Node::isInteger() const {
    if (m_kind == NodeKind::Integer) return true;
    if (m_kind == NodeKind::Rational) return m_denominator.isOne();
    return false;
}

bool Node::isZero() const {
    if (m_kind != NodeKind::Integer) return false;
    return m_numerator.isZero();
}

bool Node::isOne() const {
    if (m_kind != NodeKind::Integer) return false;
    return m_numerator.isOne();
}

bool Node::isNegativeOne() const {
    if (m_kind != NodeKind::Integer) return false;
    return m_numerator == BigInt(-1);
}

bool Node::isNegative() const {
    if (!isNumber()) return false;
    return m_numerator.isNegative();
}

/* ------------------------------------------------------------ 比较 */

bool Node::equals(const NodePtr &other) const {
    if (this == other) return true;
    if (!other || m_kind != other->m_kind) return false;
    switch (m_kind) {
    case NodeKind::Integer:
        return m_numerator == other->m_numerator;
    case NodeKind::Rational:
        return m_numerator == other->m_numerator && m_denominator == other->m_denominator;
    case NodeKind::Symbol:
        return m_symbolId == other->m_symbolId;
    case NodeKind::Constant:
        return m_constant == other->m_constant;
    default: break;
    }
    if (m_name != other->m_name) {
        const int order = compareText(m_name, other->m_name);
        if (order != 0) return false;
    }
    if (m_operands.size() != other->m_operands.size()) return false;
    for (std::size_t i = 0; i < m_operands.size(); ++i) {
        if (!m_operands[i]->equals(other->m_operands[i])) return false;
    }
    return true;
}

std::uint64_t Node::hash() const {
    if (m_hash) return m_hash;
    std::uint64_t seed = 1469598103934665603ull;
    auto mix = [&seed](std::uint64_t value) {
        seed ^= value;
        seed *= 1099511628211ull;
    };
    mix(static_cast<std::uint64_t>(m_kind) + 0x9e3779b97f4a7c15ull);
    switch (m_kind) {
    case NodeKind::Integer:
        for (std::size_t i = 0; i < m_numerator.limbCount(); ++i) mix(m_numerator.limbAt(i));
        mix(m_numerator.sign() + 1);
        break;
    case NodeKind::Rational:
        for (std::size_t i = 0; i < m_numerator.limbCount(); ++i) mix(m_numerator.limbAt(i));
        mix(m_numerator.sign() + 1);
        for (std::size_t i = 0; i < m_denominator.limbCount(); ++i) mix(m_denominator.limbAt(i));
        break;
    case NodeKind::Symbol:
        mix(m_symbolId);
        break;
    case NodeKind::Constant:
        mix(static_cast<std::uint64_t>(m_constant));
        break;
    default:
        if (m_name) {
            for (const char *p = m_name; *p; ++p) mix(static_cast<unsigned char>(*p));
        }
        for (NodePtr child : m_operands) mix(child->hash());
        break;
    }
    m_hash = seed ? seed : 1;
    return m_hash;
}

int Node::compareTo(const NodePtr &other) const {
    if (this == other) return 0;
    if (!other) return 1;
    const int leftRank = kindRank(m_kind), rightRank = kindRank(other->m_kind);
    if (leftRank != rightRank) return leftRank < rightRank ? -1 : 1;

    switch (m_kind) {
    case NodeKind::Integer:
    case NodeKind::Rational:
        return compareNumbers(this, other);
    case NodeKind::Symbol: {
        // 按名字字典序而不是驻留编号：驻留顺序依赖程序执行历史，会让输出顺序不可预期。
        if (m_symbolId == other->m_symbolId) return 0;
        return compareText(m_name, other->m_name);
    }
    case NodeKind::Constant: {
        if (m_constant == other->m_constant) return 0;
        return m_constant < other->m_constant ? -1 : 1;
    }
    default: break;
    }

    const int nameOrder = compareText(m_name, other->m_name);
    if (nameOrder != 0) return nameOrder;
    const std::size_t count = std::min(m_operands.size(), other->m_operands.size());
    for (std::size_t i = 0; i < count; ++i) {
        const int order = m_operands[i]->compareTo(other->m_operands[i]);
        if (order != 0) return order;
    }
    if (m_operands.size() == other->m_operands.size()) return 0;
    return m_operands.size() < other->m_operands.size() ? -1 : 1;
}

/* ------------------------------------------------------------ 符号表 */

std::uint32_t internSymbol(const std::string &name) {
    SymbolTable &table = symbolTable();
    std::lock_guard<std::mutex> lock(table.guard);
    auto it = table.index.find(name);
    if (it != table.index.end()) return it->second;
    const auto id = static_cast<std::uint32_t>(table.names.size());
    table.names.push_back(std::unique_ptr<std::string>(new std::string(name)));
    table.index.emplace(name, id);
    return id;
}

const char *symbolName(std::uint32_t id) {
    SymbolTable &table = symbolTable();
    std::lock_guard<std::mutex> lock(table.guard);
    if (id >= table.names.size()) return "";
    return table.names[id]->c_str();
}

Memory::Vector<std::uint32_t> allSymbols() {
    SymbolTable &table = symbolTable();
    std::lock_guard<std::mutex> lock(table.guard);
    Memory::Vector<std::uint32_t> result;
    result.reserve(table.names.size());
    for (std::size_t i = 0; i < table.names.size(); ++i) result.push_back(static_cast<std::uint32_t>(i));
    return result;
}

std::uint32_t lookupSymbol(const std::string &name) {
    SymbolTable &table = symbolTable();
    std::lock_guard<std::mutex> lock(table.guard);
    auto it = table.index.find(name);
    return it == table.index.end() ? static_cast<std::uint32_t>(-1) : it->second;
}

} // namespace Cas::Detail

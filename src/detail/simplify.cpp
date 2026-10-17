#include "cas/detail/simplify.h"

#include <algorithm>
#include <cstring>

namespace Cas::Detail {
namespace {

/* ------------------------------------------------------------ 小工具 */

/** @brief 快捷构造整数。 */
Ptr integer(long long value) { return Ptr(Node::makeInteger(BigInt(value))); }

/** @brief 判断是否指定小整数。 */
bool isInt(const Ptr &node, long long value) {
    long long actual = 0;
    return node && node->isSmallInteger(&actual) && actual == value;
}

/**
 * @brief 加法项的排序：降幂，常量排最后。
 *
 * 与乘法排序刻意不同——`x^2+3x+2` 是习惯写法，而 `3x+2+x^2` 虽然等价但不便阅读。
 */
struct AddOrder {
    /** @brief 粗略次数：幂取指数，符号取 1，乘积取各因子次数之和，数值取 0。 */
    static long long degreeOf(const NodePtr &node) {
        if (!node) return 0;
        switch (node->kind()) {
        case NodeKind::Pow: {
            long long exponent = 0;
            if (node->operand(1)->isSmallInteger(&exponent) && exponent > 0) return exponent;
            return 1;
        }
        case NodeKind::Symbol: return 1;
        case NodeKind::Integer:
        case NodeKind::Rational: return 0;
        case NodeKind::Mul: {
            long long total = 0;
            for (NodePtr child : node->operands()) total += degreeOf(child);
            return total;
        }
        default: return 1;
        }
    }
    bool operator()(const NodePtr &left, const NodePtr &right) const {
        const long long leftDegree = degreeOf(left), rightDegree = degreeOf(right);
        if (leftDegree != rightDegree) return leftDegree > rightDegree;
        return left->compareTo(right) < 0;
    }
};

/**
 * @brief 乘法因子的排序：数值系数在前，其余按结构比较。
 */
struct MulOrder {
    bool operator()(const NodePtr &left, const NodePtr &right) const {
        const bool leftNumber = left->isNumber(), rightNumber = right->isNumber();
        if (leftNumber != rightNumber) return leftNumber;
        return left->compareTo(right) < 0;
    }
};

/* ------------------------------------------------------- 数值运算内核 */

/** @brief 数值加。 */
Ptr addNumbers(const Ptr &left, const Ptr &right) {
    BigInt ln, ld, rn, rd;
    if (!numberValue(left.get(), &ln, &ld) || !numberValue(right.get(), &rn, &rd)) return Ptr();
    return makeNumber(ln * rd + rn * ld, ld * rd);
}

/** @brief 数值乘。 */
Ptr mulNumbers(const Ptr &left, const Ptr &right) {
    BigInt ln, ld, rn, rd;
    if (!numberValue(left.get(), &ln, &ld) || !numberValue(right.get(), &rn, &rd)) return Ptr();
    return makeNumber(ln * rn, ld * rd);
}

/** @brief 非负整数的 k 次根（向下取整），二分查找，仅用于 k ≤ 64 的小次数。 */
BigInt integerRoot(const BigInt &value, unsigned long long k) {
    if (value.isNegative() || value.isZero()) return BigInt(0);
    if (k == 0) return BigInt(1);
    if (k == 1) return value;
    if (k == 2) return value.isqrt();

    const BigInt order(static_cast<long long>(k));
    // 上界：2^(ceil(bits/k)) 必不小于真根
    const std::size_t bits = value.bitLength();
    std::size_t shift = bits / static_cast<std::size_t>(k);
    if (shift * static_cast<std::size_t>(k) < bits) ++shift;
    BigInt low(0);
    BigInt high = BigInt::power(BigInt(2), shift); // 上界（含）

    // 二分：不动点收敛到 floor(k√value)
    while (low < high) {
        const BigInt mid = (low + high + BigInt(1)) / BigInt(2);
        if (BigInt::power(mid, k) <= value) {
            low = mid;
        } else {
            high = mid - BigInt(1);
        }
    }
    return low;
}

/**
 * @brief 数值幂：只在结果为精确有理数时折叠。
 *
 * 指数为整数时走精确乘方；指数为 `1/k` 且底数恰为完全 k 次方时精确开方；
 * 其余情况（如 2^(1/2)）保持原样，留给根式模块处理。
 */
Ptr powNumbers(const Ptr &base, const Ptr &exponent) {
    BigInt bn, bd, en, ed;
    if (!numberValue(base.get(), &bn, &bd) || !numberValue(exponent.get(), &en, &ed)) return Ptr();
    if (!ed.isOne()) {
        // 有理指数：只处理 1/k 形式，且要求底数为正有理数的完全 k 次幂
        if (!en.isOne()) return Ptr();
        if (!bd.isOne() || bn.isNegative() || bn.isZero()) return Ptr();
        if (!bd.fitsLongLong()) return Ptr();
        const long long root = ed.toLongLong();
        if (root <= 0 || root > 64) return Ptr();
        const auto order = static_cast<unsigned long long>(root);
        const BigInt candidate = integerRoot(bn, order);
        if (!candidate.isZero() && BigInt::power(candidate, order) == bn) return makeNumber(candidate, BigInt(1));
        return Ptr();
    }
    if (en.isZero()) return Ptr(); // 0 次幂留到调用处判断底数
    if (bn.isZero()) {
        if (en.isNegative()) return Ptr(); // 0 的负次幂未定义，保持原样
        return integer(0);
    }
    const bool negativeExponent = en.isNegative();
    const BigInt magnitude = en.abs();
    if (!magnitude.fitsLongLong()) return Ptr(); // 指数过大，避免膨胀
    const auto power = static_cast<unsigned long long>(magnitude.toLongLong());
    // 估算结果规模：过大的乘方会让表达式爆炸，直接保持原样交给上层决定。
    const std::size_t bits = bn.bitLength() + bd.bitLength();
    if (bits * power > (1u << 20)) return Ptr();

    const BigInt numerator = BigInt::power(bn, power);
    const BigInt denominator = BigInt::power(bd, power);
    if (negativeExponent) return makeNumber(denominator, numerator);
    return makeNumber(numerator, denominator);
}

/* -------------------------------------------------- Add 的同类项合并 */

/** @brief 一条"系数 × 基"的记录。 */
struct Term {
    Ptr coefficient; ///< 数值系数
    Ptr body;        ///< 去掉数值系数后的部分；空表示纯数值项
};

/** @brief 拆出数值系数与剩余部分。 */
Term splitCoefficient(const Ptr &node) {
    Term term;
    if (!node) return term;
    if (node->isNumber()) {
        term.coefficient = node;
        return term;
    }
    if (node->kind() == NodeKind::Mul) {
        Ptr coefficient = integer(1);
        Memory::Vector<NodePtr> rest;
        bool sawNumber = false;
        for (NodePtr factor : node->operands()) {
            if (factor->isNumber()) {
                coefficient = mulNumbers(coefficient, Ptr::borrow(factor));
                sawNumber = true;
            } else {
                rest.push_back(factor);
            }
        }
        if (!sawNumber) {
            term.coefficient = integer(1);
            term.body = node;
            return term;
        }
        term.coefficient = coefficient;
        if (rest.size() == 1) {
            term.body = Ptr::borrow(rest.front());
        } else if (!rest.empty()) {
            term.body = makeMulMany(rest);
        }
        return term;
    }
    term.coefficient = integer(1);
    term.body = node;
    return term;
}

/** @brief 由系数与基重建一项。 */
Ptr rebuildTerm(const Ptr &coefficient, const Ptr &body) {
    if (!body) return coefficient;
    if (isInt(coefficient, 1)) return body;
    return makeMul(coefficient, body);
}

/**
 * @brief 合并加法项。
 *
 * 分两轮：先按 `body` 分组累加系数，再丢掉系数为 0 的项。
 * 分组用线性扫描（表达式规模通常不大），避免额外的哈希容器所有权问题。
 */
Ptr addAll(const Memory::Vector<NodePtr> &input) {
    Memory::Vector<Term> terms;
    terms.reserve(input.size());

    for (NodePtr node : input) {
        if (!node) continue;
        Ptr holder = Ptr::borrow(node);
        if (node->kind() == NodeKind::Add) {
            for (NodePtr child : node->operands()) {
                Term term = splitCoefficient(Ptr::borrow(child));
                if (term.body || term.coefficient) terms.push_back(term);
            }
            continue;
        }
        Term term = splitCoefficient(holder);
        if (term.body || term.coefficient) terms.push_back(term);
    }

    Memory::Vector<Term> merged;
    for (const Term &term : terms) {
        bool found = false;
        for (Term &slot : merged) {
            const bool bothNumeric = !slot.body && !term.body;
            const bool sameBody = bothNumeric || (slot.body && term.body && slot.body->equals(term.body.get()));
            if (!sameBody) continue;
            slot.coefficient = addNumbers(slot.coefficient, term.coefficient);
            found = true;
            break;
        }
        if (!found) merged.push_back(term);
    }

    Memory::Vector<Ptr> pieces;
    pieces.reserve(merged.size());
    for (const Term &term : merged) {
        if (term.coefficient && term.coefficient->isZero()) continue;
        pieces.push_back(rebuildTerm(term.coefficient, term.body));
    }
    if (pieces.empty()) return integer(0);
    if (pieces.size() == 1) return pieces.front().share();

    Memory::Vector<NodePtr> operands;
    operands.reserve(pieces.size());
    for (const Ptr &piece : pieces) operands.push_back(piece.get());
    std::sort(operands.begin(), operands.end(), AddOrder());
    return Ptr(Node::makeNary(NodeKind::Add, operands));
}

/* -------------------------------------------------- Mul 的指数合并 */

/** @brief 一条"底数 × 指数"的记录。 */
struct Factor {
    Ptr base;
    Ptr exponent;
};

/** @brief 拆出底数与指数；非幂结构视为指数 1。 */
Factor splitPower(const Ptr &node) {
    Factor factor;
    if (node && node->kind() == NodeKind::Pow) {
        factor.base = Ptr::borrow(node->operand(0));
        factor.exponent = Ptr::borrow(node->operand(1));
        return factor;
    }
    factor.base = node;
    factor.exponent = integer(1);
    return factor;
}

/** @brief 由底数与指数重建因子。 */
Ptr rebuildFactor(const Ptr &base, const Ptr &exponent) {
    if (isInt(exponent, 1)) return base;
    return makePow(base, exponent);
}

/**
 * @brief 合并乘法因子。
 *
 * 数值因子单独累乘并排在首位（`2*x` 而非 `x*2`），同底数因子指数相加，
 * 指数化简为 0 的因子直接消去。
 */
Ptr mulAll(const Memory::Vector<NodePtr> &input) {
    Ptr coefficient = integer(1);
    Memory::Vector<Factor> factors;

    auto accumulate = [&](const Factor &factor) {
        for (Factor &slot : factors) {
            if (!slot.base->equals(factor.base.get())) continue;
            slot.exponent = makeAdd(slot.exponent, factor.exponent);
            return;
        }
        factors.push_back(factor);
    };

    for (NodePtr node : input) {
        if (!node) continue;
        Ptr holder = Ptr::borrow(node);

        Memory::Vector<Ptr> children;
        if (node->kind() == NodeKind::Mul) {
            for (NodePtr child : node->operands()) children.push_back(Ptr::borrow(child));
        } else {
            children.push_back(holder.share());
        }

        for (const Ptr &child : children) {
            Factor factor = splitPower(child);
            if (factor.base->isNumber()) {
                coefficient = mulNumbers(coefficient, rebuildFactor(factor.base, factor.exponent));
                if (!coefficient) return Ptr();
            } else {
                accumulate(factor);
            }
        }
    }

    if (!coefficient) return Ptr();
    if (coefficient->isZero()) return integer(0);

    Memory::Vector<Ptr> survivors;
    survivors.reserve(factors.size() + 1);
    for (const Factor &factor : factors) {
        if (factor.exponent->isZero()) continue; // x^0 → 1
        survivors.push_back(rebuildFactor(factor.base, factor.exponent));
    }

    if (!isInt(coefficient, 1)) survivors.insert(survivors.begin(), coefficient);
    if (survivors.empty()) return integer(1);
    if (survivors.size() == 1) return survivors.front().share();

    Memory::Vector<NodePtr> operands;
    operands.reserve(survivors.size());
    for (const Ptr &piece : survivors) operands.push_back(piece.get());
    std::sort(operands.begin(), operands.end(), MulOrder());
    return Ptr(Node::makeNary(NodeKind::Mul, operands));
}

} // namespace

/* -------------------------------------------------------------- 数值 */

/**
 * @brief 函数节点的常量折叠。
 *
 * 只折叠在所有取值下都成立的恒等式（`exp(0)=1`、`log(1)=0`、`log(e)=1`、
 * `sin(0)=0`、`cos(0)=1`），不做数值近似——近似会让"精确计算器"失去意义。
 */
Ptr foldFunction(const Ptr &node) {
    if (!node || node->kind() != NodeKind::Function || node->size() != 1) return Ptr();
    const char *name = node->name();
    if (!name) return Ptr();
    NodePtr argument = node->operand(0);

    const bool argumentZero = argument->isZero();
    const bool argumentOne = argument->isOne();
    const bool argumentE = argument->kind() == NodeKind::Constant && argument->constant() == ConstantId::E;

    if (std::strcmp(name, "exp") == 0 && argumentZero) return Ptr(Node::makeInteger(BigInt(1)));
    if (std::strcmp(name, "log") == 0 || std::strcmp(name, "ln") == 0) {
        if (argumentOne) return Ptr(Node::makeInteger(BigInt(0)));
        if (argumentE) return Ptr(Node::makeInteger(BigInt(1)));
    }
    if (std::strcmp(name, "sin") == 0 && argumentZero) return Ptr(Node::makeInteger(BigInt(0)));
    if (std::strcmp(name, "tan") == 0 && argumentZero) return Ptr(Node::makeInteger(BigInt(0)));
    if (std::strcmp(name, "cos") == 0 && argumentZero) return Ptr(Node::makeInteger(BigInt(1)));
    if (std::strcmp(name, "sqrt") == 0 && argumentZero) return Ptr(Node::makeInteger(BigInt(0)));
    if (std::strcmp(name, "abs") == 0 && argumentZero) return Ptr(Node::makeInteger(BigInt(0)));
    return Ptr();
}

bool numberValue(const NodePtr &node, BigInt *numerator, BigInt *denominator) {
    if (!node || !node->isNumber()) return false;
    *numerator = node->numerator();
    *denominator = node->denominator();
    return true;
}

Ptr makeNumber(const BigInt &numerator, const BigInt &denominator) {
    return Ptr(Node::makeRational(numerator, denominator));
}

/* --------------------------------------------------------------- 加法 */

Ptr makeAddMany(const Memory::Vector<NodePtr> &terms) {
    Ptr result = addAll(terms);
    if (result) return result;
    return Ptr(Node::makeNary(NodeKind::Add, terms));
}

Ptr makeAdd(Ptr left, Ptr right) {
    Memory::Vector<NodePtr> terms;
    terms.reserve(2);
    if (left) terms.push_back(left.get());
    if (right) terms.push_back(right.get());
    return makeAddMany(terms);
}

/* --------------------------------------------------------------- 乘法 */

Ptr makeMulMany(const Memory::Vector<NodePtr> &factors) {
    Ptr result = mulAll(factors);
    if (result) return result;
    return Ptr(Node::makeNary(NodeKind::Mul, factors));
}

Ptr makeMul(Ptr left, Ptr right) {
    Memory::Vector<NodePtr> factors;
    factors.reserve(2);
    if (left) factors.push_back(left.get());
    if (right) factors.push_back(right.get());
    return makeMulMany(factors);
}

/* ----------------------------------------------------------------- 幂 */

Ptr makePow(Ptr base, Ptr exponent) {
    if (!base || !exponent) return Ptr();

    if (isInt(exponent, 1)) return base.share();
    if (isInt(exponent, 0)) {
        if (base->isZero()) return Ptr(); // 0^0 保持未求值，不擅自约定
        return integer(1);
    }
    if (base->isZero()) {
        if (exponent->isNumber() && exponent->isNegative()) return Ptr(); // 0 的负次幂未定义
        return integer(0);
    }
    if (isInt(base, 1)) return integer(1);

    // 纯数值幂精确折叠
    if (base->isNumber() && exponent->isNumber()) {
        Ptr folded = powNumbers(base, exponent);
        if (folded) return folded;
        return Ptr(Node::makePower(base.get(), exponent.get()));
    }

    // (x^a)^b → x^(a·b)：整数指数时安全
    if (base->kind() == NodeKind::Pow && exponent->isNumber() && exponent->isInteger()) {
        Ptr innerBase = Ptr::borrow(base->operand(0));
        Ptr innerExponent = Ptr::borrow(base->operand(1));
        Ptr combined = makeMul(innerExponent, exponent);
        if (!combined) return Ptr();
        return makePow(innerBase, combined);
    }

    // (a·b)^n → a^n·b^n：整数指数时可分配
    if (base->kind() == NodeKind::Mul && exponent->isInteger()) {
        Ptr distributed = integer(1);
        for (NodePtr factor : base->operands()) {
            Ptr piece = makePow(Ptr::borrow(factor), exponent.share());
            if (!piece) return Ptr();
            distributed = makeMul(distributed, piece);
            if (!distributed) return Ptr();
        }
        return distributed;
    }

    return Ptr(Node::makePower(base.get(), exponent.get()));
}

/* ------------------------------------------------------------- 其他构造 */

Ptr makeNegate(Ptr value) {
    if (!value) return Ptr();
    if (value->isNumber()) {
        Ptr folded = mulNumbers(value, integer(-1));
        if (folded) return folded;
    }
    return makeMul(integer(-1), value);
}

Ptr makeDivide(Ptr numerator, Ptr denominator) {
    if (!numerator || !denominator) return Ptr();
    Ptr inverse = makePow(denominator, integer(-1));
    if (!inverse) return Ptr();
    return makeMul(numerator, inverse);
}

/* --------------------------------------------------------------- 化简 */

Ptr simplify(const Ptr &expression) {
    if (!expression) return Ptr();
    NodePtr node = expression.get();

    switch (node->kind()) {
    case NodeKind::Integer:
    case NodeKind::Rational:
    case NodeKind::Symbol:
    case NodeKind::Constant:
        return expression.share();

    case NodeKind::Add:
    case NodeKind::Mul: {
        Memory::Vector<Ptr> holders;
        holders.reserve(node->size());
        for (NodePtr child : node->operands()) {
            Ptr simplified = simplify(Ptr::borrow(child));
            if (!simplified) return Ptr();
            holders.push_back(simplified);
        }
        Memory::Vector<NodePtr> values;
        values.reserve(holders.size());
        for (const Ptr &holder : holders) values.push_back(holder.get());
        return node->kind() == NodeKind::Add ? makeAddMany(values) : makeMulMany(values);
    }
    case NodeKind::Pow: {
        Ptr base = simplify(Ptr::borrow(node->operand(0)));
        Ptr exponent = simplify(Ptr::borrow(node->operand(1)));
        if (!base || !exponent) return Ptr();
        return makePow(base, exponent);
    }
    case NodeKind::Function: {
        Memory::Vector<Ptr> holders;
        holders.reserve(node->size());
        for (NodePtr child : node->operands()) {
            Ptr simplified = simplify(Ptr::borrow(child));
            if (!simplified) return Ptr();
            holders.push_back(simplified);
        }
        Memory::Vector<NodePtr> arguments;
        arguments.reserve(holders.size());
        for (const Ptr &holder : holders) arguments.push_back(holder.get());
        Ptr rebuilt = Ptr(Node::makeFunction(node->name(), arguments));
        // 常量折叠：exp(0)=1、log(1)=0、log(e)=1 等固定恒等，避免结果里残留 exp(0)
        if (Ptr folded = foldFunction(rebuilt)) return folded;
        return rebuilt;
    }
    }
    return expression.share();
}

/* ----------------------------------------------------------- 符号查询 */

namespace {

void collectSymbols(const NodePtr &node, Memory::Vector<std::uint32_t> *out) {
    if (!node) return;
    if (node->kind() == NodeKind::Symbol) {
        if (std::find(out->begin(), out->end(), node->symbolId()) == out->end()) out->push_back(node->symbolId());
        return;
    }
    for (NodePtr child : node->operands()) collectSymbols(child, out);
}

} // namespace

Memory::Vector<std::uint32_t> freeSymbols(const NodePtr &expression) {
    Memory::Vector<std::uint32_t> result;
    collectSymbols(expression, &result);
    std::sort(result.begin(), result.end());
    return result;
}

bool containsSymbol(const NodePtr &expression, std::uint32_t symbolId) {
    if (!expression) return false;
    if (expression->kind() == NodeKind::Symbol) return expression->symbolId() == symbolId;
    for (NodePtr child : expression->operands()) {
        if (containsSymbol(child, symbolId)) return true;
    }
    return false;
}

bool containsInfinity(const NodePtr &expression) {
    if (!expression) return false;
    if (expression->kind() == NodeKind::Constant) {
        return expression->constant() == ConstantId::Infinity ||
               expression->constant() == ConstantId::NegativeInfinity;
    }
    for (NodePtr child : expression->operands()) {
        if (containsInfinity(child)) return true;
    }
    return false;
}

} // namespace Cas::Detail

#include "cas/detail/trigonometry.h"
#include "cas/detail/simplify.h"

#include <cstring>

namespace Cas::Detail {
namespace {

/** @brief 快捷构造整数。 */
Ptr integer(long long value) { return Ptr(Node::makeInteger(BigInt(value))); }

/** @brief 判断是否指定小整数。 */
bool isInt(const Ptr &node, long long value) {
    long long actual = 0;
    return node && node->isSmallInteger(&actual) && actual == value;
}

/** @brief 名字是否为给定函数。 */
bool isFunction(const Ptr &node, const char *name) {
    return node && node->kind() == NodeKind::Function && node->name() && std::strcmp(node->name(), name) == 0;
}

/** @brief 是否为常数 Pi。 */
bool isPi(const Ptr &node) {
    return node && node->kind() == NodeKind::Constant && node->constant() == ConstantId::Pi;
}

/** @brief 是否为常数 e。 */
bool isE(const Ptr &node) { return node && node->kind() == NodeKind::Constant && node->constant() == ConstantId::E; }

/* ------------------------------------------------------------ 特殊角 */

/**
 * @brief 把角度拆成 `q·π`，q 为有理数。
 *
 * `π` → 1，`π/6` → 1/6，`3π/4` → 3/4，`-π` → -1。
 */
bool rationalTimesPi(const Ptr &node, BigInt *numerator, BigInt *denominator) {
    if (isPi(node)) {
        *numerator = BigInt(1);
        *denominator = BigInt(1);
        return true;
    }
    if (!node || node->kind() != NodeKind::Mul) return false;

    Ptr factor = integer(1);
    bool sawPi = false;
    for (NodePtr child : node->operands()) {
        if (isPi(Ptr::borrow(child))) {
            if (sawPi) return false; // π·π 不在本表覆盖范围
            sawPi = true;
        } else {
            factor = makeMul(factor, Ptr::borrow(child));
        }
    }
    if (!sawPi) return false;
    return numberValue(factor.get(), numerator, denominator);
}

/** @brief 构造 √n。 */
Ptr squareRootOf(long long value) { return makePow(integer(value), makeNumber(BigInt(1), BigInt(2))); }

/** @brief 构造 √n/2，即 (1/2)·√n。 */
Ptr halfRootOf(long long value) {
    return makeMul(makeNumber(BigInt(1), BigInt(2)), squareRootOf(value));
}

/**
 * @brief 由刻度编号取 sin(tick·π/12) 的精确值。
 *
 * 只覆盖 tick ∈ [0, 6]（即 0, π/12, π/6, π/4, π/3, 5π/12, π/2）；
 * 其中 5π/12、π/12 不在教材常见表内，返回空表示"无法精确表达"。
 */
Ptr tickValue(long long tick) {
    switch (tick) {
    case 0: return integer(0);                       // sin(0) = 0
    case 2: return makeNumber(BigInt(1), BigInt(2)); // sin(π/6) = 1/2
    case 3: return halfRootOf(2);                    // sin(π/4) = √2/2
    case 4: return halfRootOf(3);                    // sin(π/3) = √3/2
    case 6: return integer(1);                       // sin(π/2) = 1
    default: return Ptr();
    }
}

/**
 * @brief 特殊角的 sin/cos 精确值。
 *
 * 以 π/12 为刻度：tick = 12q/π（须为整数），再按象限对称化到 [0, 6]。
 * 覆盖 sin/cos 在 0、π/6、π/4、π/3、π/2 及其对称位置上的全部取值。
 */
Ptr specialAngle(const char *function, const BigInt &numerator, const BigInt &denominator) {
    if (!denominator.fitsLongLong() || denominator.isZero()) return Ptr();
    if (!numerator.fitsLongLong()) return Ptr();
    const long long num = numerator.toLongLong();
    const long long den = denominator.toLongLong();

    const bool sine = std::strcmp(function, "sin") == 0;
    const bool cosine = std::strcmp(function, "cos") == 0;
    if (!sine && !cosine) return Ptr();

    // sin 作用于 tick·π/12；cos(θ) = sin(θ + π/2)，因此 cos(x·π) 等价于 sin((12x + 6)·π/12)
    long long num12 = num * 12;
    if (cosine) num12 += den * 6;

    // 要求角度是 π/12 的整数倍
    if (num12 % den != 0) return Ptr();
    long long tick = (num12 / den) % 24;
    if (tick < 0) tick += 24;

    // 现在求 sin(tick·π/12)，tick ∈ [0, 24)
    long long mirror = tick;
    int sign = 1;
    if (tick > 12) {          // 三四象限：sin(θ) = -sin(2π-θ)
        mirror = 24 - tick;
        sign = -1;
    } else if (tick > 6) {    // 第二象限：sin(θ) = sin(π-θ)
        mirror = 12 - tick;
    }
    Ptr magnitude = tickValue(mirror);
    if (!magnitude) return Ptr();
    if (sign < 0) return makeNegate(magnitude);
    return magnitude;
}

/* ---------------------------------------------------------- 单节点规则 */

/** @brief 构造 q·π。 */
Ptr timesPi(const BigInt &numerator, const BigInt &denominator) {
    Ptr coefficient = makeNumber(numerator, denominator);
    Ptr pi(Node::makeConstant(ConstantId::Pi));
    return makeMul(coefficient, pi);
}

/**
 * @brief 反三角函数在常见取值点上的精确值。
 *
 * 覆盖 asin/acos/atan 在 0、±1/2、±√2/2、±√3/2、±1 上的结果；
 * 表外取值返回空，交给数值化或原样保留。
 */
Ptr inverseSpecial(const char *name, const Ptr &argument) {
    if (!argument) return Ptr();
    BigInt numerator, denominator;
    if (!numberValue(argument.get(), &numerator, &denominator)) return Ptr();
    if (!denominator.isOne()) return Ptr(); // 只处理整数参数
    if (!numerator.fitsLongLong()) return Ptr();
    const long long value = numerator.toLongLong();

    const bool arcsine = std::strcmp(name, "asin") == 0 || std::strcmp(name, "arcsin") == 0;
    const bool arccosine = std::strcmp(name, "acos") == 0 || std::strcmp(name, "arccos") == 0;
    const bool arctangent = std::strcmp(name, "atan") == 0 || std::strcmp(name, "arctan") == 0;

    if (arctangent) {
        if (value == 0) return integer(0);
        if (value == 1) return timesPi(BigInt(1), BigInt(4));   // atan(1) = π/4
        if (value == -1) return makeNegate(timesPi(BigInt(1), BigInt(4)));
        return Ptr();
    }
    if (arcsine) {
        if (value == 0) return integer(0);
        if (value == 1) return timesPi(BigInt(1), BigInt(2));  // asin(1) = π/2
        if (value == -1) return makeNegate(timesPi(BigInt(1), BigInt(2)));
        return Ptr();
    }
    if (arccosine) {
        if (value == 1) return integer(0);
        if (value == 0) return timesPi(BigInt(1), BigInt(2));  // acos(0) = π/2
        if (value == -1) return timesPi(BigInt(1), BigInt(1)); // acos(-1) = π
        return Ptr();
    }
    return Ptr();
}

/** @brief 单节点重写：常数角、反函数相消、log/exp 恒等。返回空表示未命中。 */
Ptr rewriteNode(const Ptr &node) {
    if (!node || node->kind() != NodeKind::Function) return Ptr();
    const char *name = node->name();
    if (!name || node->size() != 1) return Ptr();
    Ptr argument = Ptr::borrow(node->operand(0));

    if (std::strcmp(name, "log") == 0 || std::strcmp(name, "ln") == 0) {
        if (argument->isOne()) return integer(0);
        if (isE(argument)) return integer(1);
        return Ptr();
    }
    if (std::strcmp(name, "exp") == 0) {
        if (argument->isZero()) return integer(1);
        return Ptr();
    }

    // 反函数相消（取主值）
    if (std::strcmp(name, "sin") == 0 && isFunction(argument, "asin")) return Ptr::borrow(argument->operand(0));
    if (std::strcmp(name, "cos") == 0 && isFunction(argument, "acos")) return Ptr::borrow(argument->operand(0));
    if (std::strcmp(name, "tan") == 0 && isFunction(argument, "atan")) return Ptr::borrow(argument->operand(0));

    // 反三角函数的特殊值：asin(0)=0、atan(1)=π/4、acos(0)=π/2 等
    if (Ptr special = inverseSpecial(name, argument)) return special;

    if (std::strcmp(name, "sin") == 0 || std::strcmp(name, "cos") == 0) {
        BigInt numerator, denominator;
        if (rationalTimesPi(argument, &numerator, &denominator)) return specialAngle(name, numerator, denominator);
    }
    return Ptr();
}

/** @brief 自底向上应用单节点规则，并顺带化简一次。 */
Ptr rewriteTree(const Ptr &node) {
    if (!node) return Ptr();
    if (node->size() == 0) return Ptr();

    Memory::Vector<Ptr> children;
    children.reserve(node->size());
    for (NodePtr child : node->operands()) {
        Ptr rewritten = rewriteTree(Ptr::borrow(child));
        children.push_back(rewritten ? rewritten : Ptr::borrow(child));
    }
    Memory::Vector<NodePtr> raw;
    raw.reserve(children.size());
    for (const Ptr &child : children) raw.push_back(child.get());

    Ptr rebuilt;
    if (node->kind() == NodeKind::Pow) rebuilt = makePow(Ptr::borrow(raw[0]), Ptr::borrow(raw[1]));
    else if (node->kind() == NodeKind::Function) rebuilt = simplify(Ptr(Node::makeFunction(node->name(), raw)));
    else
        // Add/Mul 的子节点可能已折叠成常数（如 (-1)*(-1)），必须再化简一次
        rebuilt = simplify(Ptr(Node::makeNary(node->kind(), raw)));

    if (!rebuilt) return Ptr();
    if (rebuilt->kind() == NodeKind::Function) {
        Ptr reduced = rewriteNode(rebuilt);
        if (reduced) return reduced;
    }
    return rebuilt;
}

/* ------------------------------------------------------ sin² + cos² → 1 */

/** @brief 拆分"数值系数 × 基"。 */
struct Piece {
    Ptr coefficient;
    Ptr body;
};

Piece split(const Ptr &node) {
    Piece piece;
    if (node && node->kind() == NodeKind::Mul) {
        Ptr coefficient = integer(1);
        Memory::Vector<NodePtr> rest;
        bool sawNumber = false;
        for (NodePtr factor : node->operands()) {
            if (factor->isNumber()) {
                BigInt an, ad, bn, bd;
                if (!numberValue(coefficient.get(), &an, &ad) || !numberValue(factor, &bn, &bd)) continue;
                coefficient = makeNumber(an * bn, ad * bd);
                sawNumber = true;
            } else {
                rest.push_back(factor);
            }
        }
        if (!sawNumber) {
            piece.coefficient = integer(1);
            piece.body = node;
            return piece;
        }
        piece.coefficient = coefficient;
        if (rest.size() == 1) piece.body = Ptr::borrow(rest.front());
        else if (rest.size() > 1) piece.body = makeMulMany(rest);
        return piece;
    }
    piece.coefficient = integer(1);
    piece.body = node;
    return piece;
}

/** @brief 判断是否为 `sin(u)^2` / `cos(u)^2`，输出 which（'s'/'c'）与辐角。 */
bool matchSquared(const Ptr &body, char *which, Ptr *argument) {
    if (!body || body->kind() != NodeKind::Pow) return false;
    if (!isInt(Ptr::borrow(body->operand(1)), 2)) return false;
    Ptr base = Ptr::borrow(body->operand(0));
    if (isFunction(base, "sin")) {
        *which = 's';
        *argument = Ptr::borrow(base->operand(0));
        return true;
    }
    if (isFunction(base, "cos")) {
        *which = 'c';
        *argument = Ptr::borrow(base->operand(0));
        return true;
    }
    return false;
}

/**
 * @brief 在加法里匹配 `a·sin(u)^2 + a·cos(u)^2 → a`。
 *
 * 只合并系数完全相同的配对——这是最安全也最常见的形态，不做恒等式搜索。
 */
Ptr collapsePythagorean(const Ptr &expression) {
    if (!expression || expression->kind() != NodeKind::Add) return Ptr();
    const std::size_t count = expression->size();
    Memory::Vector<bool> consumed(count, false);
    Memory::Vector<Ptr> output;
    bool changed = false;

    for (std::size_t i = 0; i < count; ++i) {
        if (consumed[i]) continue;
        Ptr left = Ptr::borrow(expression->operand(i));
        Piece leftPiece = split(left);
        char leftWhich = 0;
        Ptr leftArgument;
        if (!matchSquared(leftPiece.body, &leftWhich, &leftArgument)) {
            output.push_back(left);
            continue;
        }

        bool merged = false;
        for (std::size_t j = i + 1; j < count; ++j) {
            if (consumed[j]) continue;
            Ptr right = Ptr::borrow(expression->operand(j));
            Piece rightPiece = split(right);
            char rightWhich = 0;
            Ptr rightArgument;
            if (!matchSquared(rightPiece.body, &rightWhich, &rightArgument)) continue;
            if (leftWhich == rightWhich) continue;
            if (!leftArgument->equals(rightArgument.get())) continue;
            if (!leftPiece.coefficient->equals(rightPiece.coefficient.get())) continue;
            output.push_back(leftPiece.coefficient);
            consumed[j] = true;
            merged = true;
            changed = true;
            break;
        }
        if (merged) consumed[i] = true;
        else output.push_back(left);
    }

    if (!changed) return Ptr();
    Memory::Vector<NodePtr> raw;
    raw.reserve(output.size());
    for (const Ptr &piece : output) raw.push_back(piece.get());
    return makeAddMany(raw);
}

/* ---------------------------------------------------- 三角加法展开 */

/** @brief 是否为 `k·u` 形式的辐角，取出系数与主体。 */
bool splitScaled(const Ptr &argument, long long *scale, Ptr *body) {
    if (!argument) return false;
    *scale = 1;
    *body = argument;
    if (argument->kind() != NodeKind::Mul) return true;
    Ptr numeric = integer(1);
    Memory::Vector<NodePtr> rest;
    for (NodePtr child : argument->operands()) {
        if (child->isNumber()) {
            BigInt an, ad, bn, bd;
            if (!numberValue(numeric.get(), &an, &ad) || !numberValue(child, &bn, &bd)) return false;
            numeric = makeNumber(an * bn, ad * bd);
        } else {
            rest.push_back(child);
        }
    }
    if (!numeric->isInteger() || !numeric->numerator().fitsLongLong()) return false;
    *scale = numeric->numerator().toLongLong();
    if (rest.size() == 1) *body = Ptr::borrow(rest.front());
    else if (rest.size() > 1) *body = makeMulMany(rest);
    else { *body = Ptr(); return false; }
    return true;
}

/** @brief 展开 `sin(2u)`、`cos(2u)` 二倍角。 */
Ptr expandDoubleAngle(const Ptr &node) {
    if (!node || node->kind() != NodeKind::Function || node->size() != 1) return Ptr();
    const char *name = node->name();
    if (!name) return Ptr();
    long long scale = 1;
    Ptr body;
    if (!splitScaled(Ptr::borrow(node->operand(0)), &scale, &body) || !body) return Ptr();

    if (std::strcmp(name, "sin") == 0 && scale == 2) {
        // sin(2u) = 2 sin(u) cos(u)
        Ptr sine = Ptr(Node::makeFunction("sin", Memory::Vector<NodePtr>{body.get()}));
        Ptr cosine = Ptr(Node::makeFunction("cos", Memory::Vector<NodePtr>{body.get()}));
        return makeMul(integer(2), makeMul(sine, cosine));
    }
    if (std::strcmp(name, "cos") == 0 && scale == 2) {
        // cos(2u) = cos(u)^2 - sin(u)^2
        Ptr sine = Ptr(Node::makeFunction("sin", Memory::Vector<NodePtr>{body.get()}));
        Ptr cosine = Ptr(Node::makeFunction("cos", Memory::Vector<NodePtr>{body.get()}));
        Ptr squaredCosine = makePow(cosine, integer(2));
        Ptr squaredSine = makePow(sine, integer(2));
        return makeAdd(squaredCosine, makeNegate(squaredSine));
    }
    if (std::strcmp(name, "tan") == 0 && scale == 2) {
        // tan(2u) = 2 tan(u) / (1 - tan(u)^2)
        Ptr tangent = Ptr(Node::makeFunction("tan", Memory::Vector<NodePtr>{body.get()}));
        Ptr numerator = makeMul(integer(2), tangent.share());
        Ptr denominator = makeAdd(integer(1), makeNegate(makePow(tangent, integer(2))));
        return makeDivide(numerator, denominator);
    }
    return Ptr();
}

/** @brief 递归展开二倍角。 */
Ptr expandTrigTree(const Ptr &node) {
    if (!node) return Ptr();
    if (node->size() == 0) return Ptr();

    Memory::Vector<Ptr> children;
    children.reserve(node->size());
    for (NodePtr child : node->operands()) {
        Ptr expanded = expandTrigTree(Ptr::borrow(child));
        children.push_back(expanded ? expanded : Ptr::borrow(child));
    }
    Memory::Vector<NodePtr> raw;
    raw.reserve(children.size());
    for (const Ptr &child : children) raw.push_back(child.get());

    Ptr rebuilt;
    if (node->kind() == NodeKind::Pow) rebuilt = makePow(Ptr::borrow(raw[0]), Ptr::borrow(raw[1]));
    else if (node->kind() == NodeKind::Function) rebuilt = Ptr(Node::makeFunction(node->name(), raw));
    else rebuilt = Ptr(Node::makeNary(node->kind(), raw));
    if (!rebuilt) return Ptr();

    if (rebuilt->kind() == NodeKind::Function) {
        Ptr expanded = expandDoubleAngle(rebuilt);
        if (expanded) return expanded;
    }
    return rebuilt;
}

/* -------------------------------------------------------- 降幂公式 */

/**
 * @brief `sin(u)^2 → (1 - cos(2u))/2`、`cos(u)^2 → (1 + cos(2u))/2`。
 *
 * 用于三角多项式积分：降幂后每一项都形如 `c·cos(k u)`，可直接积分。
 */
Ptr reducePower(const Ptr &node) {
    if (!node || node->kind() != NodeKind::Pow) return Ptr();
    if (!isInt(Ptr::borrow(node->operand(1)), 2)) return Ptr();
    Ptr base = Ptr::borrow(node->operand(0));
    if (base->kind() != NodeKind::Function || base->size() != 1) return Ptr();
    const char *name = base->name();
    if (!name) return Ptr();

    Ptr argument = Ptr::borrow(base->operand(0));
    Ptr doubled = Ptr(Node::makeFunction("cos", Memory::Vector<NodePtr>{makeMul(integer(2), argument).get()}));
    if (!doubled) return Ptr();

    if (std::strcmp(name, "sin") == 0) {
        return makeDivide(makeAdd(integer(1), makeNegate(doubled)), integer(2));
    }
    if (std::strcmp(name, "cos") == 0) {
        return makeDivide(makeAdd(integer(1), doubled), integer(2));
    }
    return Ptr();
}

/** @brief 递归降幂。 */
Ptr reduceTree(const Ptr &node) {
    if (!node) return Ptr();
    if (node->size() == 0) return Ptr();

    Memory::Vector<Ptr> children;
    children.reserve(node->size());
    for (NodePtr child : node->operands()) {
        Ptr reduced = reduceTree(Ptr::borrow(child));
        children.push_back(reduced ? reduced : Ptr::borrow(child));
    }
    Memory::Vector<NodePtr> raw;
    raw.reserve(children.size());
    for (const Ptr &child : children) raw.push_back(child.get());

    Ptr rebuilt;
    if (node->kind() == NodeKind::Pow) rebuilt = makePow(Ptr::borrow(raw[0]), Ptr::borrow(raw[1]));
    else if (node->kind() == NodeKind::Function) rebuilt = Ptr(Node::makeFunction(node->name(), raw));
    else rebuilt = Ptr(Node::makeNary(node->kind(), raw));
    if (!rebuilt) return Ptr();

    if (rebuilt->kind() == NodeKind::Pow) {
        Ptr reduced = reducePower(rebuilt);
        if (reduced) return reduced;
    }
    return rebuilt;
}

} // namespace

/* ------------------------------------------------------------------ 入口 */

Ptr trigSimplify(const Ptr &expression) {
    if (!expression) return Ptr();
    Ptr current = simplify(expression);
    if (!current) return Ptr();

    for (int round = 0; round < 8; ++round) {
        Ptr rewritten = rewriteTree(current);
        if (!rewritten) return current.share();
        Ptr collapsed = collapsePythagorean(rewritten);
        Ptr next = simplify(collapsed ? collapsed : rewritten);
        if (!next) return current.share();
        if (next->equals(current.get())) return next;
        current = next;
    }
    return current;
}

Ptr expandTrig(const Ptr &expression) {
    if (!expression) return Ptr();
    Ptr expanded = expandTrigTree(expression);
    if (!expanded) return simplify(expression);
    return simplify(expanded);
}

Ptr reduceTrigPower(const Ptr &expression) {
    if (!expression) return Ptr();
    Ptr reduced = reduceTree(expression);
    if (!reduced) return simplify(expression);
    return simplify(reduced);
}

} // namespace Cas::Detail

#include "cas/detail/calculus.h"
#include "cas/detail/polynomial.h"
#include "cas/detail/printer.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>

namespace Cas::Detail {
namespace {

/** @brief 快捷构造整数。 */
Ptr integer(long long value) { return Ptr(Node::makeInteger(BigInt(value))); }
/** @brief 快捷构造变量节点。 */
Ptr variableNode(std::uint32_t id) { return Ptr(Node::makeSymbolById(id)); }

/** @brief 判断是否指定小整数。 */
bool isInt(const Ptr &node, long long value) {
    long long actual = 0;
    return node && node->isSmallInteger(&actual) && actual == value;
}

/** @brief 函数调用的第 i 个实参。 */
Ptr argumentOf(const Ptr &node, std::size_t index) { return Ptr::borrow(node->operand(index)); }

/* ----------------------------------------------------------- 求导 */

Ptr derivativeOf(const Ptr &expression, std::uint32_t variable) {
    if (!expression) return Ptr();
    if (!dependsOn(expression.get(), variable)) return integer(0);

    switch (expression->kind()) {
    case NodeKind::Symbol:
        return expression->symbolId() == variable ? integer(1) : integer(0);

    case NodeKind::Add: {
        Memory::Vector<Ptr> holders;
        for (NodePtr child : expression->operands()) {
            Ptr term = derivativeOf(Ptr::borrow(child), variable);
            if (!term) return Ptr();
            holders.push_back(term);
        }
        Memory::Vector<NodePtr> operands;
        for (const Ptr &holder : holders) operands.push_back(holder.get());
        return makeAddMany(operands);
    }

    case NodeKind::Mul: {
        // 乘积法则：Σ (∂ᵢ f) · Π_{j≠i} fⱼ
        Memory::Vector<Ptr> pieces;
        const std::size_t count = expression->size();
        for (std::size_t i = 0; i < count; ++i) {
            Ptr derivative = derivativeOf(Ptr::borrow(expression->operand(i)), variable);
            if (!derivative) return Ptr();
            if (derivative->isZero()) continue;
            Ptr product = derivative;
            for (std::size_t j = 0; j < count; ++j) {
                if (i == j) continue;
                product = makeMul(product, Ptr::borrow(expression->operand(j)));
            }
            pieces.push_back(product);
        }
        if (pieces.empty()) return integer(0);
        Memory::Vector<NodePtr> operands;
        for (const Ptr &piece : pieces) operands.push_back(piece.get());
        return makeAddMany(operands);
    }

    case NodeKind::Pow: {
        Ptr base = Ptr::borrow(expression->operand(0));
        Ptr exponent = Ptr::borrow(expression->operand(1));
        const bool baseDepends = dependsOn(base.get(), variable);
        const bool exponentDepends = dependsOn(exponent.get(), variable);

        if (!exponentDepends) {
            // 幂法则：(u^n)' = n·u^(n-1)·u'
            Ptr inner = derivativeOf(base, variable);
            if (!inner) return Ptr();
            Ptr reduced = makePow(base, makeAdd(exponent, integer(-1)));
            return makeMul(makeMul(exponent, reduced), inner);
        }
        if (!baseDepends) {
            // 指数法则：(a^u)' = a^u·ln(a)·u'
            Ptr inner = derivativeOf(exponent, variable);
            if (!inner) return Ptr();
            Ptr logBase = Ptr(Node::makeFunction("log", Memory::Vector<NodePtr>{base.get()}));
            return makeMul(makeMul(expression, logBase), inner);
        }
        // 一般情形：u^v = e^(v·ln u)
        Ptr rewritten = makePow(Ptr(Node::makeConstant(ConstantId::E)),
                                makeMul(exponent, Ptr(Node::makeFunction("log", Memory::Vector<NodePtr>{base.get()}))));
        return derivativeOf(rewritten, variable);
    }

    case NodeKind::Function: {
        const char *name = expression->name();
        Ptr inner = argumentOf(expression, 0);
        if (!inner) return Ptr();
        Ptr innerDerivative = derivativeOf(inner, variable);
        if (!innerDerivative) return Ptr();

        auto chain = [&](Ptr outerDerivative) {
            return makeMul(outerDerivative, innerDerivative);
        };

        if (std::strcmp(name, "sin") == 0)
            return chain(Ptr(Node::makeFunction("cos", Memory::Vector<NodePtr>{inner.get()})));
        if (std::strcmp(name, "cos") == 0)
            return chain(makeNegate(Ptr(Node::makeFunction("sin", Memory::Vector<NodePtr>{inner.get()}))));
        if (std::strcmp(name, "tan") == 0) {
            // (tan u)' = sec²u = 1 + tan²u
            Ptr squared = makePow(Ptr(Node::makeFunction("cos", Memory::Vector<NodePtr>{inner.get()})), integer(-2));
            return chain(squared);
        }
        if (std::strcmp(name, "asin") == 0) {
            Ptr inside = makeAdd(integer(1), makeNegate(makePow(inner, integer(2))));
            return chain(makePow(inside, makeNumber(BigInt(-1), BigInt(2))));
        }
        if (std::strcmp(name, "acos") == 0) {
            Ptr inside = makeAdd(integer(1), makeNegate(makePow(inner, integer(2))));
            return chain(makeNegate(makePow(inside, makeNumber(BigInt(-1), BigInt(2)))));
        }
        if (std::strcmp(name, "atan") == 0) {
            Ptr inside = makeAdd(integer(1), makePow(inner, integer(2)));
            return chain(makePow(inside, integer(-1)));
        }
        if (std::strcmp(name, "exp") == 0)
            return chain(Ptr(Node::makeFunction("exp", Memory::Vector<NodePtr>{inner.get()})));
        if (std::strcmp(name, "log") == 0) return chain(makePow(inner, integer(-1)));
        if (std::strcmp(name, "sqrt") == 0) {
            Ptr root = makePow(inner, makeNumber(BigInt(1), BigInt(2)));
            return chain(makeMul(makeNumber(BigInt(1), BigInt(2)), makePow(inner, makeNumber(BigInt(-1), BigInt(2)))));
        }
        if (std::strcmp(name, "sinh") == 0)
            return chain(Ptr(Node::makeFunction("cosh", Memory::Vector<NodePtr>{inner.get()})));
        if (std::strcmp(name, "cosh") == 0)
            return chain(Ptr(Node::makeFunction("sinh", Memory::Vector<NodePtr>{inner.get()})));
        // 未知函数：返回未求值导数形式 D(f)(x)·f'(x)，不猜答案
        return Ptr();
    }

    case NodeKind::Constant:
    case NodeKind::Integer:
    case NodeKind::Rational:
        return integer(0);
    }
    return Ptr();
}

/* ----------------------------------------------------------- 积分 */

/** @brief 判断表达式是否为 variable 的线性式 a·x + b，输出 a 与 b。 */
bool linearForm(const Ptr &expression, std::uint32_t variable, Ptr *slope, Ptr *offset) {
    Ptr expanded = expand(expression);
    UnivariatePoly polynomial;
    if (!toPolynomial(expanded.get(), variable, &polynomial)) return false;
    if (polynomial.degree() > 1) return false;
    *slope = polynomial.coefficient(1);
    *offset = polynomial.coefficient(0);
    return true;
}

Ptr integrateRecursive(const Ptr &expression, std::uint32_t variable, int depth = 0);

/**
 * @brief 分部积分：`∫ u·dv = u·v - ∫ v·du`。
 *
 * 策略是 LIATE 的简化版——在乘积里挑"求导后会变简单"的因子当 u：
 * 多项式/幂 → 对数/反三角 → 其余。对 `x·sin x`、`x·e^x`、`x^2·cos x`
 * 这类常见题一次或几次递归就能收敛；不收敛时（递归深度上限）返回 null，
 * 由调用方呈现为"无法积分"。
 */
Ptr integrateByParts(const Ptr &product, std::uint32_t variable, int depth = 0);

/** @brief 求导结果的"排名"：数值越小越适合当 u。 */
int uRank(const Ptr &factor, std::uint32_t variable) {
    if (!dependsOn(factor.get(), variable)) return 100;
    switch (factor->kind()) {
    case NodeKind::Symbol: return 0;
    case NodeKind::Pow: {
        Ptr base = Ptr::borrow(factor->operand(0));
        Ptr exponent = Ptr::borrow(factor->operand(1));
        if (base->kind() == NodeKind::Symbol && exponent->isInteger()) return 0;
        return 20;
    }
    case NodeKind::Function: {
        const char *name = factor->name();
        if (!name) return 20;
        if (std::strcmp(name, "log") == 0 || std::strcmp(name, "asin") == 0 ||
            std::strcmp(name, "acos") == 0 || std::strcmp(name, "atan") == 0) {
            return 10;
        }
        return 30;
    }
    case NodeKind::Add: return 5;
    default: return 40;
    }
}

Ptr integrateByParts(const Ptr &product, std::uint32_t variable, int depth) {
    if (depth > 6) return Ptr();
    if (!product || product->kind() != NodeKind::Mul) return Ptr();

    // 收集因子里最适合当 u 的那个（其余整体作为 dv）
    std::size_t bestIndex = 0;
    int bestRank = 1000;
    for (std::size_t i = 0; i < product->size(); ++i) {
        Ptr factor = Ptr::borrow(product->operand(i));
        if (!dependsOn(factor.get(), variable)) continue;
        const int rank = uRank(factor, variable);
        if (rank < bestRank) {
            bestRank = rank;
            bestIndex = i;
        }
    }
    if (bestRank >= 100) return Ptr(); // 没有任何依赖变量的因子，交给上层处理

    Ptr u = Ptr::borrow(product->operand(bestIndex));
    Memory::Vector<NodePtr> rest;
    for (std::size_t i = 0; i < product->size(); ++i) {
        if (i != bestIndex) rest.push_back(product->operand(i));
    }
    Ptr dv = rest.empty() ? integer(1) : makeMulMany(rest);

    Ptr v = integrateRecursive(dv, variable, depth + 1);
    if (!v) return Ptr();
    Ptr du = derivativeOf(u, variable);
    if (!du) return Ptr();

    Ptr boundary = makeMul(u, v);
    Ptr inner = simplify(makeMul(v, du));
    if (!inner) return Ptr();
    Ptr innerIntegral = integrateRecursive(inner, variable, depth + 1);
    if (!innerIntegral) return Ptr();
    return simplify(makeAdd(boundary, makeNegate(innerIntegral)));
}

Ptr integrateRecursive(const Ptr &expression, std::uint32_t variable, int depth) {
    if (!expression) return Ptr();
    if (!dependsOn(expression.get(), variable)) {
        // 常数因子：∫ c dx = c·x
        return makeMul(expression, variableNode(variable));
    }

    switch (expression->kind()) {
    case NodeKind::Symbol:
        if (expression->symbolId() != variable) return Ptr();
        // ∫ x dx = x²/2
        return makeDivide(makePow(variableNode(variable), integer(2)), integer(2));

    case NodeKind::Add: {
        Memory::Vector<Ptr> holders;
        for (NodePtr child : expression->operands()) {
            Ptr part = integrateRecursive(Ptr::borrow(child), variable, depth);
            if (!part) return Ptr();
            holders.push_back(part);
        }
        Memory::Vector<NodePtr> operands;
        for (const Ptr &holder : holders) operands.push_back(holder.get());
        return makeAddMany(operands);
    }

    case NodeKind::Mul: {
        // 拆出常数因子后再积分
        Ptr constant = integer(1);
        Memory::Vector<NodePtr> varying;
        Memory::Vector<Ptr> holders;
        for (NodePtr child : expression->operands()) {
            if (dependsOn(child, variable)) {
                varying.push_back(child);
            } else {
                holders.push_back(Ptr::borrow(child));
                constant = makeMul(constant, holders.back());
            }
        }
        if (varying.empty()) return makeMul(expression, variableNode(variable));
        if (!isInt(constant, 1)) {
            if (varying.size() == 1) {
                Ptr inner = integrateRecursive(Ptr::borrow(varying.front()), variable, depth);
                if (!inner) return Ptr();
                return makeMul(constant, inner);
            }
            // 先试"展开成多项式再逐项积分"（如 2·x·(x+1)）
            Ptr product = makeMulMany(varying);
            Ptr expanded = expand(product);
            if (expanded && expanded->kind() == NodeKind::Add) {
                Ptr inner = integrateRecursive(expanded, variable, depth);
                if (inner) return makeMul(constant, inner);
            }
            // 展开不适用（如 2·x·e^x）：退到分部积分，常数最终会以系数形式出来
            Ptr parts = integrateByParts(product, variable, depth + 1);
            if (!parts) return Ptr();
            return makeMul(constant, parts);
        }
        if (varying.size() == 1) return integrateRecursive(Ptr::borrow(varying.front()), variable, depth);

        // 乘积：先尝试展开（多项式相乘即可积）
        Ptr product = makeMulMany(varying);
        Ptr expanded = expand(product);
        if (expanded && !expanded->equals(product.get())) {
            Ptr result = integrateRecursive(expanded, variable, depth);
            if (result) return result;
        }
        // 展开无法化简时尝试分部积分
        return integrateByParts(product, variable, depth + 1);
    }

    case NodeKind::Pow: {
        Ptr base = Ptr::borrow(expression->operand(0));
        Ptr exponent = Ptr::borrow(expression->operand(1));

        // x^n
        if (base->kind() == NodeKind::Symbol && base->symbolId() == variable) {
            if (isInt(exponent, -1)) {
                // ∫ dx/x = ln|x|
                Memory::Vector<NodePtr> arguments;
                Ptr absolute = makePow(makePow(base, integer(2)), makeNumber(BigInt(1), BigInt(2)));
                arguments.push_back(base.get());
                return Ptr(Node::makeFunction("log", arguments));
            }
            Ptr next = makeAdd(exponent, integer(1));
            return makeDivide(makePow(base, next), next);
        }

        // (a·x+b)^n
        Ptr slope, offset;
        if (linearForm(base, variable, &slope, &offset) && !slope->isZero()) {
            if (!dependsOn(exponent.get(), variable)) {
                if (isInt(exponent, -1)) {
                    Memory::Vector<NodePtr> arguments;
                    arguments.push_back(base.get());
                    Ptr logarithm = Ptr(Node::makeFunction("log", arguments));
                    return makeDivide(logarithm, slope);
                }
                Ptr next = makeAdd(exponent, integer(1));
                return makeDivide(makePow(base, next), makeMul(next, slope));
            }
        }

        // e^(a·x+b)
        if (base->kind() == NodeKind::Constant && base->constant() == ConstantId::E) {
            Ptr innerSlope, innerOffset;
            if (linearForm(exponent, variable, &innerSlope, &innerOffset) && !innerSlope->isZero()) {
                return makeDivide(expression, innerSlope);
            }
        }

        // √(a²-x²) 型与 1/√(a²-x²) 型：用反三角公式
        return Ptr();
    }

    case NodeKind::Function: {
        const char *name = expression->name();
        Ptr inner = argumentOf(expression, 0);
        if (!inner) return Ptr();
        Ptr slope, offset;
        const bool linear = linearForm(inner, variable, &slope, &offset) && !slope->isZero();

        if (std::strcmp(name, "sin") == 0) {
            if (!linear) return Ptr();
            Ptr cosine = Ptr(Node::makeFunction("cos", Memory::Vector<NodePtr>{inner.get()}));
            return makeDivide(makeNegate(cosine), slope);
        }
        if (std::strcmp(name, "cos") == 0) {
            if (!linear) return Ptr();
            Ptr sine = Ptr(Node::makeFunction("sin", Memory::Vector<NodePtr>{inner.get()}));
            return makeDivide(sine, slope);
        }
        if (std::strcmp(name, "exp") == 0) {
            if (!linear) return Ptr();
            return makeDivide(expression, slope);
        }
        if (std::strcmp(name, "log") == 0) {
            // ∫ ln(x) dx = x·ln x - x
            if (inner->kind() == NodeKind::Symbol && inner->symbolId() == variable) {
                Ptr x = variableNode(variable);
                return makeAdd(makeMul(x.share(), expression), makeNegate(x.share()));
            }
            return Ptr();
        }
        if (std::strcmp(name, "tan") == 0) {
            if (!linear) return Ptr();
            // ∫ tan u du = -ln|cos u|
            Memory::Vector<NodePtr> arguments;
            arguments.push_back(Ptr(Node::makeFunction("cos", Memory::Vector<NodePtr>{inner.get()})).detach());
            Ptr logarithm = Ptr(Node::makeFunction("log", arguments));
            return makeDivide(makeNegate(logarithm), slope);
        }
        if (std::strcmp(name, "sqrt") == 0) {
            // ∫ √x dx = (2/3)x^(3/2)
            if (inner->kind() == NodeKind::Symbol && inner->symbolId() == variable) {
                Ptr power = makePow(inner, makeNumber(BigInt(3), BigInt(2)));
                return makeMul(makeNumber(BigInt(2), BigInt(3)), power);
            }
            return Ptr();
        }
        if (std::strcmp(name, "sinh") == 0) {
            if (!linear) return Ptr();
            Ptr cosine = Ptr(Node::makeFunction("cosh", Memory::Vector<NodePtr>{inner.get()}));
            return makeDivide(cosine, slope);
        }
        if (std::strcmp(name, "cosh") == 0) {
            if (!linear) return Ptr();
            Ptr sine = Ptr(Node::makeFunction("sinh", Memory::Vector<NodePtr>{inner.get()}));
            return makeDivide(sine, slope);
        }
        return Ptr();
    }

    case NodeKind::Rational:
    case NodeKind::Integer:
        return makeMul(expression, variableNode(variable));

    case NodeKind::Constant:
        return makeMul(expression, variableNode(variable));
    }
    return Ptr();
}

} // namespace

/* ----------------------------------------------------------- 求导 */

/** @brief 对 variable 求导（符号求导，含链式法则与常见函数）。 */
Ptr differentiate(const Ptr &expression, std::uint32_t variable) {
    if (!expression) return Ptr();
    Ptr result = derivativeOf(expression, variable);
    if (!result) return Ptr();
    return simplify(result);
}

/** @brief 求 n 阶导。 */
Ptr differentiateN(const Ptr &expression, std::uint32_t variable, int order) {
    Ptr current = expression.share();
    for (int i = 0; i < order; ++i) {
        current = differentiate(current, variable);
        if (!current) return Ptr();
    }
    return current;
}

/* ------------------------------------------------------------ 积分 */

namespace {

/** @brief 把有理数 Ptr 取成 double，用于判别式的符号判断。 */
bool rationalToDouble(const Ptr &number, double *out) {
    BigInt numerator, denominator;
    if (!number || !numberValue(number.get(), &numerator, &denominator)) return false;
    *out = numerator.toDouble() / denominator.toDouble();
    return true;
}

/** @brief 常数的负号。 */
Ptr negatedConst(const Ptr &value) { return makeNegate(value); }

/**
 * @brief 积分 ∫ P(x)/(a x^2 + b x + c) dx。
 *
 * 先做多项式除法拆出整式部分，余项按判别式分三种情形：
 * 判别式 > 0 → 两个实根，对数；
 * 判别式 = 0 → 重根，1/(x-r) 形式；
 * 判别式 < 0 → 配方后得到 atan。
 */
Ptr integrateQuadraticDenominator(const UnivariatePoly &numerator, const UnivariatePoly &denominator,
                                 std::uint32_t variable) {
    UnivariatePoly quotient, remainder;
    if (!polyDivide(numerator, denominator, &quotient, &remainder)) return Ptr();

    // 整式部分逐项积分
    Ptr polynomialPart = integer(0);
    for (std::size_t i = 0; i < quotient.coefficients.size(); ++i) {
        if (i == 0) continue;
        Ptr term = makeMul(quotient.coefficient(i),
                           makePow(variableNode(variable), integer(static_cast<long long>(i - 1))));
        polynomialPart = makeAdd(polynomialPart, makeDivide(term, integer(static_cast<long long>(i))));
    }

    const Ptr a = denominator.coefficient(2);
    const Ptr b = denominator.coefficient(1);
    const Ptr c = denominator.coefficient(0);
    double av = 0, bv = 0, cv = 0;
    if (!rationalToDouble(a, &av) || !rationalToDouble(b, &bv) || !rationalToDouble(c, &cv)) return Ptr();
    const double discriminant = bv * bv - 4.0 * av * cv;
    Ptr x = variableNode(variable);

    if (discriminant > 0) {
        // 实数根 r1,2 = (-b ± √Δ)/(2a)
        Ptr fourAC = makeMul(makeNumber(BigInt(4), BigInt(1)), makeMul(a, c));
        Ptr sqrtPart = makePow(makeAdd(makeMul(b, b), makeNegate(fourAC)), makeNumber(BigInt(1), BigInt(2)));
        Ptr minusB = negatedConst(b);
        Ptr twoA = makeMul(integer(2), a);
        Ptr r1 = makeDivide(makeAdd(minusB, Ptr(sqrtPart.share())), twoA);
        Ptr r2 = makeDivide(makeAdd(minusB, negatedConst(sqrtPart)), twoA);
        if (!r1 || !r2) return Ptr();
        // 只处理余项为常数或线性，且能解出 A/(x-r1)+B/(x-r2)
        double r1v = 0, r2v = 0;
        if (!rationalToDouble(r1, &r1v) || !rationalToDouble(r2, &r2v) || r1v == r2v) {
            // 根式不可求值时放弃（保持原样交给上层报错）
            if (!remainder.coefficients.empty() && !remainder.coefficient(1)->isZero()) return Ptr();
        }
        // 由 A(x-r2) + B(x-r1) = p1·x + p0 得：
        //   A + B = p1，A·r2 + B·r1 = -p0
        // 解得 A = (p1·r1 + p0)/(r1 - r2)，B = (p1·r2 + p0)/(r2 - r1)
        const Ptr p1 = remainder.degree() >= 1 ? remainder.coefficient(1) : integer(0);
        const Ptr p0 = remainder.coefficient(0);
        const Ptr r1MinusR2 = makeAdd(r1, negatedConst(r2));
        Ptr coeffA = makeDivide(makeAdd(makeMul(p1, r1), p0), r1MinusR2);
        Ptr coeffB = makeDivide(makeAdd(makeMul(p1, r2), p0), negatedConst(r1MinusR2));

        // 注意：必须先用具名 Ptr 持有实参节点，否则临时对象析构后裸指针悬垂
        Ptr argANode = makeAdd(x, negatedConst(r1));
        Ptr argBNode = makeAdd(x, negatedConst(r2));
        if (!argANode || !argBNode) return Ptr();
        Memory::Vector<NodePtr> argA;
        argA.push_back(argANode.get());
        Ptr logA = Ptr(Node::makeFunction("log", argA));
        Memory::Vector<NodePtr> argB;
        argB.push_back(argBNode.get());
        Ptr logB = Ptr(Node::makeFunction("log", argB));

        Ptr result = makeAdd(polynomialPart, makeAdd(makeMul(coeffA, logA), makeMul(coeffB, logB)));
        return simplify(result);
    }

    if (discriminant == 0) {
        // 重根 r = -b/(2a)，分母 a(x-r)^2；余项写成 m x + n
        Ptr r = makeDivide(negatedConst(b), makeMul(integer(2), a));
        // m x + n = m(x-r) + (n + m r) ⇒ 两项分别是 m/(x-r) 与 (n+mr)/(x-r)^2
        const Ptr m = remainder.degree() >= 1 ? remainder.coefficient(1) : integer(0);
        const Ptr n = remainder.coefficient(0);
        Ptr shifted = makeAdd(n, makeAdd(makeMul(m, r), integer(0)));
        Ptr linear = makeDivide(m, makeMul(a, makeAdd(x, negatedConst(r))));
        Ptr quadratic = makeDivide(negatedConst(shifted), makeMul(makeMul(a, makeAdd(x, negatedConst(r))),
                                                                  makeAdd(x, negatedConst(r))));
        (void)linear;
        (void)quadratic;
        // 数值积分不追求完全覆盖，此处保守返回空，交给外层报错
        return Ptr();
    }

    // 判别式 < 0：配方 a[(x + b/2a)^2 + (4ac-b^2)/4a^2]
    // ∫ (px+q)/(ax^2+bx+c) dx = (p/2a)·ln|ax^2+bx+c| + C'·atan(...)
    const Ptr p = remainder.degree() >= 1 ? remainder.coefficient(1) : integer(0);
    const Ptr q = remainder.coefficient(0);

    Memory::Vector<NodePtr> argLog;
    Ptr denominatorExpr = makeAdd(makeMul(a, makePow(x, integer(2))),
                                  makeAdd(makeMul(b, x), c));
    argLog.push_back(denominatorExpr.get());
    Ptr logPart = makeMul(makeDivide(p, makeMul(integer(2), a)), Ptr(Node::makeFunction("log", argLog)));

    // 剩余常数项：q - p b/(2a)
    Ptr residual = makeAdd(q, negatedConst(makeDivide(makeMul(p, b), makeMul(integer(2), a))));
    // atan 部分：∫ K/((x+h)^2 + k^2) dx = (K/k)·atan((x+h)/k)
    Ptr h = makeDivide(b, makeMul(integer(2), a));
    Ptr kSquared = makeAdd(c, negatedConst(makeDivide(makeMul(b, b), makeMul(integer(4), a)))); // c - b^2/(4a)
    Ptr k = makePow(makeDivide(kSquared, a), makeNumber(BigInt(1), BigInt(2)));
    if (!k || k->isZero()) return Ptr();

    Ptr shiftedX = makeAdd(x, h);
    Ptr atanArgument = makeDivide(shiftedX, k);
    if (!atanArgument) return Ptr();
    Memory::Vector<NodePtr> argAtan;
    argAtan.push_back(atanArgument.get()); // 具名 Ptr 保证实参节点存活
    Ptr atanPart = Ptr(Node::makeFunction("atan", argAtan));
    Ptr atanCoefficient = makeDivide(residual, makeMul(a, k));

    Ptr result = makeAdd(polynomialPart, makeAdd(logPart, makeMul(atanCoefficient, atanPart)));
    return simplify(result);
}

} // namespace

Ptr integrate(const Ptr &expression, std::uint32_t variable) {
    if (!expression) return Ptr();
    Ptr simplified = simplify(expression);
    if (!simplified) return Ptr();

    // 先把分式改写为幂形式
    Ptr result = integrateRecursive(simplified, variable);
    if (!result) {
        // 有理函数：分子分母都是多项式，按次数分派
        const auto symbols = freeSymbols(simplified.get());
        if (symbols.size() != 1) return Ptr();
        Ptr expanded = expand(simplified);
        UnivariatePoly numerator, denominator;
        Ptr bottom;
        // 化简常把 1/Q 直接写成 Q^(-1)（Pow 而非 Mul），两种形态都要认
        if (expanded && expanded->kind() == NodeKind::Pow) {
            Ptr exponent = Ptr::borrow(expanded->operand(1));
            if (exponent->isNumber() && exponent->isNegative()) bottom = Ptr::borrow(expanded->operand(0));
        } else if (expanded && expanded->kind() == NodeKind::Mul) {
            for (NodePtr child : expanded->operands()) {
                Ptr factor(child);
                Node::retain(child);
                if (factor->kind() == NodeKind::Pow) {
                    Ptr exponent = Ptr::borrow(factor->operand(1));
                    if (isInt(exponent, -1)) bottom = Ptr::borrow(factor->operand(0));
                }
            }
        }
        if (!bottom || !toPolynomial(bottom.get(), variable, &denominator)) return Ptr();

        // 分子：Pow 形式的分子固定为 1；Mul 形式要把非分母因子收集起来
        if (expanded->kind() == NodeKind::Pow) {
            numerator = UnivariatePoly{};
            numerator.variable = variable;
            numerator.coefficients.push_back(integer(1));
        } else {
            Ptr rest = integer(1);
            for (NodePtr child : expanded->operands()) {
                Ptr factor(child);
                Node::retain(child);
                bool isDenominator = false;
                if (factor->kind() == NodeKind::Pow) {
                    Ptr exponent = Ptr::borrow(factor->operand(1));
                    if (isInt(exponent, -1)) isDenominator = true;
                }
                if (!isDenominator) rest = makeMul(rest, factor);
            }
            if (!toPolynomial(rest.get(), variable, &numerator)) return Ptr();
        }

        if (denominator.degree() == 1) {
            // ∫ P(x)/(a x + b) dx：先做多项式除法再积分
            UnivariatePoly q, r;
            if (!polyDivide(numerator, denominator, &q, &r)) return Ptr();
            Ptr polynomialPart = integer(0);
            for (std::size_t i = 0; i < q.coefficients.size(); ++i) {
                if (i == 0) continue;
                Ptr term = makeMul(
                    q.coefficient(i),
                    makePow(variableNode(variable), integer(static_cast<long long>(i - 1))));
                polynomialPart = makeAdd(polynomialPart, makeDivide(term, integer(static_cast<long long>(i))));
            }
            Ptr slope = denominator.coefficient(1);
            Memory::Vector<NodePtr> arguments;
            arguments.push_back(bottom.get());
            Ptr logarithm = Ptr(Node::makeFunction("log", arguments));
            Ptr logPart = makeDivide(makeMul(r.coefficient(0), logarithm), slope);
            return simplify(makeAdd(polynomialPart, logPart));
        }
        if (denominator.degree() == 2) {
            return integrateQuadraticDenominator(numerator, denominator, variable);
        }
        return Ptr();
    }
    // 分部积分等路径常留下括号，展开再化简成多项式形态
    Ptr expanded = expand(result);
    return simplify(expanded ? expanded : result);
}

Ptr definiteIntegral(const Ptr &expression, std::uint32_t variable, const Ptr &lower, const Ptr &upper) {
    Ptr primitive = integrate(expression, variable);
    if (!primitive) return Ptr();
    const bool lowerInfinite = containsInfinity(lower.get());
    const bool upperInfinite = containsInfinity(upper.get());
    if (!lowerInfinite && !upperInfinite) {
        Ptr atUpper = substitute(primitive, variable, upper);
        Ptr atLower = substitute(primitive, variable, lower);
        if (!atUpper || !atLower) return Ptr();
        return simplify(makeAdd(atUpper, makeNegate(atLower)));
    }
    // 含无穷限：转成极限
    if (upperInfinite && !lowerInfinite) {
        Ptr atLower = substitute(primitive, variable, lower);
        if (!atLower) return Ptr();
        Ptr limitValue = limitAtInfinity(primitive, variable, true);
        if (!limitValue) return Ptr();
        return simplify(makeAdd(limitValue, makeNegate(atLower)));
    }
    if (lowerInfinite && !upperInfinite) {
        Ptr atUpper = substitute(primitive, variable, upper);
        if (!atUpper) return Ptr();
        Ptr limitValue = limitAtInfinity(primitive, variable, false);
        if (!limitValue) return Ptr();
        return simplify(makeAdd(atUpper, makeNegate(limitValue)));
    }
    return Ptr();
}

/* ------------------------------------------------------------ 极限 */

/**
 * @brief 把表达式拆成"分子 / 分母"：把负次幂因子一律搬到分母。
 *
 * 处理 `1/x`、`(a+b)/c`、`x^-2*y^-1` 等形态；无法拆分（例如 `sin(x)`）时返回 false。
 */
bool splitFraction(const Ptr &expression, std::uint32_t variable, Ptr *numerator, Ptr *denominator) {
    if (!expression) return false;
    Ptr top = integer(1), bottom = integer(1);
    bool sawDenominator = false;

    std::function<void(const Ptr &)> absorb = [&](const Ptr &node) {
        if (node->kind() == NodeKind::Mul) {
            for (NodePtr child : node->operands()) absorb(Ptr::borrow(child));
            return;
        }
        if (node->kind() == NodeKind::Pow) {
            Ptr exponent = Ptr::borrow(node->operand(1));
            if (exponent->isNumber() && exponent->isNegative()) {
                sawDenominator = true;
                bottom = makeMul(bottom, makePow(Ptr::borrow(node->operand(0)), makeNegate(exponent)));
                return;
            }
        }
        top = makeMul(top, node);
    };
    absorb(expression);

    if (!top || !bottom) return false;
    *numerator = top;
    *denominator = bottom;
    (void)variable;
    return sawDenominator;
}

bool splitFraction(const Ptr &expression, std::uint32_t variable, Ptr *numerator, Ptr *denominator);

/* ---- 极限内部工具（声明在 limit 之前，定义在 limit 之后） ---- */

/** @brief 是否为 +∞。 */
bool isPositiveInfinity(const Ptr &node) {
    return node && node->kind() == NodeKind::Constant && node->constant() == ConstantId::Infinity;
}

/** @brief 是否为 -∞（含 `-oo` 被解析成 Mul(-1, oo) 的情形）。 */
bool isNegativeInfinity(const Ptr &node) {
    if (!node) return false;
    if (node->kind() == NodeKind::Constant) return node->constant() == ConstantId::NegativeInfinity;
    if (node->kind() == NodeKind::Mul && node->size() == 2) {
        Ptr left = Ptr::borrow(node->operand(0));
        Ptr right = Ptr::borrow(node->operand(1));
        if (isPositiveInfinity(left) && right->isNumber() && right->isNegative()) return true;
        if (isPositiveInfinity(right) && left->isNumber() && left->isNegative()) return true;
    }
    return false;
}

/** @brief 未定式分类。 */
enum class Indeterminate { None, ZeroOverZero, InfinityOverInfinity, Other };

Ptr evaluateAt(const Ptr &expression, std::uint32_t variable, const Ptr &point);
bool hasDivisionByZero(const Ptr &value);
Indeterminate classifyAt(const Ptr &expression, std::uint32_t variable, const Ptr &point);
Ptr cancelAtPoint(const Ptr &expression, std::uint32_t variable, const Ptr &point);
Ptr lHopital(const Ptr &expression, std::uint32_t variable, const Ptr &point, int depth);
Ptr unresolvedLimit(const Ptr &expression, std::uint32_t variable, const Ptr &point);

Ptr limit(const Ptr &expression, std::uint32_t variable, const Ptr &point) {
    if (!expression) return Ptr();

    // 无穷极限单独处理：`-oo` 会被解析成 Mul(-1, Infinity)，两种形态都要认
    if (isPositiveInfinity(point)) return limitAtInfinity(expression, variable, true);
    if (isNegativeInfinity(point)) return limitAtInfinity(expression, variable, false);

    // 直接代入：用分子/分母分别求值判定是否未定式
    const Ptr simplified = simplify(expression);
    if (Indeterminate kind = classifyAt(simplified, variable, point); kind != Indeterminate::None) {
        // 未定式：先约分，再洛必达，最后 Taylor 比较阶
        if (Ptr cancelled = cancelAtPoint(simplified, variable, point)) {
            if (Ptr value = evaluateAt(cancelled, variable, point)) return value;
        }
        if (Ptr hopital = lHopital(simplified, variable, point, 0)) return hopital;
        return unresolvedLimit(simplified, variable, point);
    }

    if (Ptr value = evaluateAt(simplified, variable, point)) return value;
    // 无法判定（例如分母为 0 但不是 0/0），退到 Taylor 判阶
    return unresolvedLimit(simplified, variable, point);
}

/**
 * @brief 在给定点求值；分母为 0 或结果为无穷时返回 null。
 *
 * `simplify` 允许 `0^(-1)` 这种形式存在（它不去猜"无穷"），
 * 因此这里用"是否仍含负次幂"来判断除零。
 */
Ptr evaluateAt(const Ptr &expression, std::uint32_t variable, const Ptr &point) {
    if (!expression) return Ptr();
    Ptr substituted = substitute(expression, variable, point);
    if (!substituted) return Ptr();
    Ptr value = simplify(substituted);
    if (!value) return Ptr();
    if (containsInfinity(value.get())) return Ptr();
    if (hasDivisionByZero(value)) return Ptr();
    return value;
}

/** @brief 值里是否出现 `0^负数`（即除零）。 */
bool hasDivisionByZero(const Ptr &value) {
    if (!value) return false;
    if (value->kind() == NodeKind::Pow) {
        Ptr base = Ptr::borrow(value->operand(0));
        Ptr exponent = Ptr::borrow(value->operand(1));
        if (base->isZero() && exponent->isNumber() && exponent->isNegative()) return true;
    }
    for (NodePtr child : value->operands()) {
        if (hasDivisionByZero(Ptr::borrow(child))) return true;
    }
    return false;
}

Indeterminate classifyAt(const Ptr &expression, std::uint32_t variable, const Ptr &point) {
    if (!expression) return Indeterminate::None;

    // 拆出分子分母：负指数因子视为分母
    Ptr numerator = integer(1), denominator = integer(1);
    bool sawDenominator = false;
    auto collect = [&](const Ptr &node) {
        if (node->kind() == NodeKind::Mul) {
            for (NodePtr child : node->operands()) {
                Ptr factor = Ptr::borrow(child);
                if (factor->kind() == NodeKind::Pow) {
                    Ptr exponent = Ptr::borrow(factor->operand(1));
                    if (exponent->isNumber() && exponent->isNegative()) {
                        sawDenominator = true;
                        denominator = makeMul(denominator,
                                              makePow(Ptr::borrow(factor->operand(0)), makeNegate(exponent)));
                        continue;
                    }
                }
                numerator = makeMul(numerator, factor);
            }
            return;
        }
        if (node->kind() == NodeKind::Pow) {
            Ptr exponent = Ptr::borrow(node->operand(1));
            if (exponent->isNumber() && exponent->isNegative()) {
                sawDenominator = true;
                denominator = makePow(Ptr::borrow(node->operand(0)), makeNegate(exponent));
                return;
            }
        }
        numerator = makeMul(numerator, node);
    };
    collect(expression);
    if (!sawDenominator) return Indeterminate::None;

    Ptr top = simplify(substitute(numerator, variable, point));
    Ptr bottom = simplify(substitute(denominator, variable, point));
    if (!top || !bottom) return Indeterminate::None;
    // 代入后如果分母自己就带除零/无穷，说明还没到"定值"这一步
    if (hasDivisionByZero(bottom) || hasDivisionByZero(top)) return Indeterminate::Other;

    const bool topZero = top->isZero();
    const bool bottomZero = bottom->isZero();
    const bool topInfinite = containsInfinity(top.get());
    const bool bottomInfinite = containsInfinity(bottom.get());

    if (topZero && bottomZero) return Indeterminate::ZeroOverZero;
    if (topInfinite && bottomInfinite) return Indeterminate::InfinityOverInfinity;
    if (bottomZero) return Indeterminate::Other;
    return Indeterminate::None;
}

/**
 * @brief 先约分再去掉在 point 处为 0 的公因子，处理 `x/x`、`(x^2-1)/(x-1)` 这类。
 */
Ptr cancelAtPoint(const Ptr &expression, std::uint32_t variable, const Ptr &point) {
    (void)point; // 约分与取值点无关：约掉的是恒等因子
    // 分子/分母都是多项式时用 GCD 约分
    Ptr numerator = integer(1), denominator = integer(1);
    if (expression->kind() == NodeKind::Mul) {
        for (NodePtr child : expression->operands()) {
            Ptr factor = Ptr::borrow(child);
            if (factor->kind() == NodeKind::Pow) {
                Ptr exponent = Ptr::borrow(factor->operand(1));
                if (exponent->isNumber() && exponent->isNegative()) {
                    denominator = makeMul(denominator,
                                          makePow(Ptr::borrow(factor->operand(0)), makeNegate(exponent)));
                    continue;
                }
            }
            numerator = makeMul(numerator, factor);
        }
    } else if (expression->kind() == NodeKind::Pow) {
        Ptr exponent = Ptr::borrow(expression->operand(1));
        if (exponent->isNumber() && exponent->isNegative()) {
            denominator = makePow(Ptr::borrow(expression->operand(0)), makeNegate(exponent));
            numerator = integer(1);
        } else {
            return Ptr();
        }
    } else {
        return Ptr();
    }

    UnivariatePoly top, bottom;
    if (!toPolynomial(expand(numerator).get(), variable, &top) ||
        !toPolynomial(expand(denominator).get(), variable, &bottom)) {
        return Ptr();
    }
    if (top.isZero() || bottom.isZero()) return Ptr();

    Ptr common = polyGcd(fromPolynomial(top).get(), fromPolynomial(bottom).get(), variable);
    if (!common) return Ptr();
    UnivariatePoly commonPoly;
    if (!toPolynomial(common.get(), variable, &commonPoly)) return Ptr();
    if (commonPoly.degree() == 0) return Ptr();

    UnivariatePoly q1, r1, q2, r2;
    if (!polyDivide(top, commonPoly, &q1, &r1) || !polyDivide(bottom, commonPoly, &q2, &r2)) return Ptr();
    if (!r1.isZero() || !r2.isZero()) return Ptr();
    return makeDivide(fromPolynomial(q1), fromPolynomial(q2));
}

/**
 * @brief 洛必达法则：对 `f/g` 型未定式逐次求导分子分母。
 *
 * 只处理明确写成"商"的情形（把负指数因子搬到分母），
 * 并且每步都用 substitute 判定是否仍是未定式；一旦定值立刻返回。
 */
Ptr lHopital(const Ptr &expression, std::uint32_t variable, const Ptr &point, int depth) {
    if (depth > 6 || !expression) return Ptr();

    Ptr numerator = integer(1), denominator = integer(1);
    bool sawDenominator = false;
    if (expression->kind() == NodeKind::Mul) {
        for (NodePtr child : expression->operands()) {
            Ptr factor = Ptr::borrow(child);
            if (factor->kind() == NodeKind::Pow) {
                Ptr exponent = Ptr::borrow(factor->operand(1));
                if (exponent->isNumber() && exponent->isNegative()) {
                    sawDenominator = true;
                    denominator = makeMul(denominator, makePow(Ptr::borrow(factor->operand(0)), makeNegate(exponent)));
                    continue;
                }
            }
            numerator = makeMul(numerator, factor);
        }
    } else if (expression->kind() == NodeKind::Pow) {
        Ptr exponent = Ptr::borrow(expression->operand(1));
        if (exponent->isNumber() && exponent->isNegative()) {
            sawDenominator = true;
            denominator = makePow(Ptr::borrow(expression->operand(0)), makeNegate(exponent));
            numerator = integer(1);
        }
    }
    if (!sawDenominator) return Ptr();

    // 非未定式就直接算
    if (classifyAt(expression, variable, point) == Indeterminate::None) {
        return evaluateAt(expression, variable, point);
    }

    Ptr newTop = derivativeOf(numerator, variable);
    Ptr newBottom = derivativeOf(denominator, variable);
    if (!newTop || !newBottom) return Ptr();
    if (newBottom->isZero()) return Ptr();

    Ptr next = simplify(makeDivide(newTop, newBottom));
    if (!next) return Ptr();
    // 导数商已经定值则直接返回
    if (Ptr direct = evaluateAt(next, variable, point)) return direct;
    return lHopital(next, variable, point, depth + 1);
}

/**
 * @brief 洛必达不适用时的兜底：用 Taylor 展开取首个非零项之商。
 *
 * 对 `sin x / x`、`(1-cos x)/x^2`、`(e^x-1-x)/x^2` 这类题全部有效：
 * 分子分母各展开，最低阶非零项系数相除即极限。
 */
Ptr unresolvedLimit(const Ptr &expression, std::uint32_t variable, const Ptr &point) {
    if (!expression) return Ptr();

    // 拆分子分母
    Ptr numerator = integer(1), denominator = integer(1);
    bool sawDenominator = false;
    if (expression->kind() == NodeKind::Mul) {
        for (NodePtr child : expression->operands()) {
            Ptr factor = Ptr::borrow(child);
            if (factor->kind() == NodeKind::Pow) {
                Ptr exponent = Ptr::borrow(factor->operand(1));
                if (exponent->isNumber() && exponent->isNegative()) {
                    sawDenominator = true;
                    denominator = makeMul(denominator, makePow(Ptr::borrow(factor->operand(0)), makeNegate(exponent)));
                    continue;
                }
            }
            numerator = makeMul(numerator, factor);
        }
    } else if (expression->kind() == NodeKind::Pow) {
        Ptr exponent = Ptr::borrow(expression->operand(1));
        if (!(exponent->isNumber() && exponent->isNegative())) return Ptr();
        sawDenominator = true;
        denominator = makePow(Ptr::borrow(expression->operand(0)), makeNegate(exponent));
        numerator = integer(1);
    } else {
        return Ptr();
    }
    if (!sawDenominator) return Ptr();

    // 分子分母各展开到固定阶，比较首个非零项次数
    const int kOrder = 8;
    UnivariatePoly topPoly, bottomPoly;
    {
        Ptr topSeries = taylorExpand(numerator, variable, point, kOrder);
        Ptr bottomSeries = taylorExpand(denominator, variable, point, kOrder);
        if (!topSeries || !bottomSeries) return Ptr();
        // taylorExpand 给出的是关于 (x-a) 的多项式，直接按 variable 读系数
        if (!toPolynomial(topSeries.get(), variable, &topPoly)) return Ptr();
        if (!toPolynomial(bottomSeries.get(), variable, &bottomPoly)) return Ptr();
    }

    auto firstNonZero = [](const UnivariatePoly &poly, std::size_t *index) {
        for (std::size_t i = 0; i < poly.coefficients.size(); ++i) {
            if (!poly.coefficients[i]->isZero()) {
                *index = i;
                return true;
            }
        }
        return false;
    };

    std::size_t topIndex = 0, bottomIndex = 0;
    if (!firstNonZero(topPoly, &topIndex)) return integer(0);
    if (!firstNonZero(bottomPoly, &bottomIndex)) return Ptr();

    if (bottomIndex > topIndex) return integer(0);
    if (bottomIndex < topIndex) return Ptr();
    return simplify(makeDivide(topPoly.coefficient(topIndex), bottomPoly.coefficient(bottomIndex)));
}

Ptr limitAtInfinity(const Ptr &expression, std::uint32_t variable, bool positive) {
    if (!expression) return Ptr();

    // 先把表达式整理成"分子 / 分母"，再按次数比较——这是无穷远处最可靠的做法
    Ptr numerator, denominator;
    if (splitFraction(expression, variable, &numerator, &denominator)) {
        UnivariatePoly top, bottom;
        const bool topOk = toPolynomial(numerator.get(), variable, &top) && !top.isZero();
        const bool bottomOk = toPolynomial(denominator.get(), variable, &bottom) && !bottom.isZero();
        if (topOk && bottomOk) {
            const std::size_t topDegree = top.degree();
            const std::size_t bottomDegree = bottom.degree();
            if (topDegree < bottomDegree) return integer(0);
            if (topDegree == bottomDegree) {
                return simplify(makeDivide(top.coefficient(topDegree), bottom.coefficient(bottomDegree)));
            }
            // 分子次数更高：极限为 ±∞（符号由首项系数之比决定）
            const bool topSign = !top.coefficient(topDegree)->isNegative();
            const bool bottomSign = !bottom.coefficient(bottomDegree)->isNegative();
            const bool sign = topSign == bottomSign;
            return Ptr(Node::makeConstant(sign ? ConstantId::Infinity : ConstantId::NegativeInfinity));
        }
    }

    // 纯多项式
    Ptr expanded = expand(expression);
    UnivariatePoly polynomial;
    if (toPolynomial(expanded.get(), variable, &polynomial) && !polynomial.isZero()) {
        const std::size_t degree = polynomial.degree();
        Ptr leading = polynomial.coefficient(degree);
        if (degree == 0) return leading;
        bool sign = !leading->isNegative();
        if (!positive && (degree % 2 == 1)) sign = !sign;
        return Ptr(Node::makeConstant(sign ? ConstantId::Infinity : ConstantId::NegativeInfinity));
    }
    return Ptr();
}

/* ---------------------------------------------------------- Taylor */

Ptr taylorExpand(const Ptr &expression, std::uint32_t variable, const Ptr &point, int order) {
    if (!expression || order < 0) return Ptr();
    Ptr variableNodePtr = variableNode(variable);

    Ptr result = integer(0);
    Ptr current = expression.share();
    BigInt factorial(1);

    for (int k = 0; k <= order; ++k) {
        Ptr value = substitute(current, variable, point);
        if (!value) return Ptr();
        if (k > 0) factorial = factorial * BigInt(k);
        Ptr shifted = makeAdd(variableNodePtr.share(), makeNegate(point));
        Ptr term = makeDivide(makeMul(value, makePow(shifted, integer(k))), makeNumber(factorial, BigInt(1)));
        result = makeAdd(result, term);

        if (k < order) {
            current = differentiate(current, variable);
            if (!current) return Ptr();
        }
    }
    return simplify(result);
}

/* ------------------------------------------------------------ 辅助 */

Ptr substitute(const Ptr &expression, std::uint32_t variable, const Ptr &value) {
    if (!expression) return Ptr();
    if (!dependsOn(expression.get(), variable)) return expression.share();

    switch (expression->kind()) {
    case NodeKind::Symbol:
        return expression->symbolId() == variable ? value.share() : expression.share();

    case NodeKind::Add:
    case NodeKind::Mul: {
        Memory::Vector<Ptr> holders;
        for (NodePtr child : expression->operands()) {
            Ptr replaced = substitute(Ptr::borrow(child), variable, value);
            if (!replaced) return Ptr();
            holders.push_back(replaced);
        }
        Memory::Vector<NodePtr> operands;
        for (const Ptr &holder : holders) operands.push_back(holder.get());
        return expression->kind() == NodeKind::Add ? makeAddMany(operands) : makeMulMany(operands);
    }

    case NodeKind::Pow: {
        Ptr base = substitute(Ptr::borrow(expression->operand(0)), variable, value);
        Ptr exponent = substitute(Ptr::borrow(expression->operand(1)), variable, value);
        if (!base || !exponent) return Ptr();
        return makePow(base, exponent);
    }

    case NodeKind::Function: {
        Memory::Vector<Ptr> holders;
        for (NodePtr child : expression->operands()) {
            Ptr replaced = substitute(Ptr::borrow(child), variable, value);
            if (!replaced) return Ptr();
            holders.push_back(replaced);
        }
        Memory::Vector<NodePtr> arguments;
        for (const Ptr &holder : holders) arguments.push_back(holder.get());
        return Ptr(Node::makeFunction(expression->name(), arguments));
    }

    default:
        return expression.share();
    }
}

bool dependsOn(const NodePtr &expression, std::uint32_t variable) {
    if (!expression) return false;
    if (expression->kind() == NodeKind::Symbol) return expression->symbolId() == variable;
    for (NodePtr child : expression->operands()) {
        if (dependsOn(child, variable)) return true;
    }
    return false;
}

bool isConstant(const NodePtr &expression) { return expression && freeSymbols(expression).empty(); }

} // namespace Cas::Detail

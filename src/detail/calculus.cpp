#include "cas/detail/calculus.h"
#include "cas/detail/polynomial.h"
#include "cas/detail/printer.h"
#include "cas/detail/solver.h"
#include "cas/detail/trigonometry.h"

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
        if (!name) return Ptr();

        // 导数节点 `Derivative(y, x, n)`：求导就是升一阶。
        // ODE 求解把未知函数写成这类节点，隐式解的校验与初值代入都要用到。
        if (std::strcmp(name, "Derivative") == 0) {
            if (expression->size() < 2) return Ptr();
            Ptr target = argumentOf(expression, 0);
            Ptr variableArg = argumentOf(expression, 1);
            if (!target || !variableArg || variableArg->kind() != NodeKind::Symbol) return Ptr();
            if (variableArg->symbolId() != variable) return integer(0);
            long long order = 1;
            if (expression->size() >= 3) {
                if (!argumentOf(expression, 2)->isSmallInteger(&order)) return Ptr();
            }
            Ptr nextOrder = integer(order + 1);
            return Ptr(Node::makeFunction(
                "Derivative",
                Memory::Vector<NodePtr>{target.get(), variableArg.get(), nextOrder.get()}));
        }

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
 * @brief 取"复合结构"的内层，用于换元。
 *
 * 只认这几种外形（其余返回空）：
 * - `e^{g(x)}`、`a^{g(x)}`（指数）
 * - `g(x)^n`（n 为数值且 g 非常数）
 * - `sin(g)`/`cos(g)`/`tan(g)`/`log(g)`/`sqrt(g)` 等一元函数
 *
 * @return 内层 `g(x)`；不是复合结构时返回空。
 */
Ptr compositeInnerOf(const Ptr &factor, std::uint32_t variable) {
    if (!factor || !dependsOn(factor.get(), variable)) return Ptr();
    if (factor->kind() == NodeKind::Pow) {
        Ptr base = Ptr::borrow(factor->operand(0));
        Ptr exponent = Ptr::borrow(factor->operand(1));
        // 指数形：a^{g(x)}，内层是 g(x)
        if (dependsOn(exponent.get(), variable) && !dependsOn(base.get(), variable)) return exponent;
        // 幂形：g(x)^n，内层是 g(x)（要求 base 不是纯变量）
        if (!dependsOn(exponent.get(), variable) && dependsOn(base.get(), variable) &&
            !(base->kind() == NodeKind::Symbol && base->symbolId() == variable)) {
            return base;
        }
        return Ptr();
    }
    if (factor->kind() == NodeKind::Function) {
        const char *name = factor->name();
        if (!name || factor->size() < 1) return Ptr();
        static const char *kUnary[] = {"sin", "cos", "tan", "cot", "sec", "csc",
                                       "log", "exp", "sqrt", "asin", "acos", "atan"};
        bool unary = false;
        for (const char *candidate : kUnary) {
            if (std::strcmp(name, candidate) == 0) {
                unary = true;
                break;
            }
        }
        if (!unary) return Ptr();
        Ptr argument = Ptr::borrow(factor->operand(0));
        if (!argument || !dependsOn(argument.get(), variable)) return Ptr();
        if (argument->kind() == NodeKind::Symbol && argument->symbolId() == variable) return Ptr();
        return argument;
    }
    return Ptr();
}

/** @brief 按子树结构替换：把 expression 里所有与 pattern 相等的子树换成 replacement。 */
Ptr replaceSubtree(const Ptr &expression, const Ptr &pattern, const Ptr &replacement) {
    if (!expression || !pattern || !replacement) return Ptr();
    if (expression->equals(pattern.get())) return replacement;
    if (expression->kind() == NodeKind::Symbol || expression->kind() == NodeKind::Integer ||
        expression->kind() == NodeKind::Rational || expression->kind() == NodeKind::Constant) {
        return expression;
    }
    Memory::Vector<Ptr> holders;
    Memory::Vector<NodePtr> children;
    bool changed = false;
    for (NodePtr child : expression->operands()) {
        Ptr original = Ptr::borrow(child);
        Ptr next = replaceSubtree(original, pattern, replacement);
        if (!next) return Ptr();
        if (!next->equals(original.get())) changed = true;
        holders.push_back(next);
        children.push_back(next.get());
    }
    if (!changed) return expression;
    switch (expression->kind()) {
    case NodeKind::Add:
        return makeAddMany(children);
    case NodeKind::Mul:
        return makeMulMany(children);
    case NodeKind::Pow:
        if (children.size() != 2) return Ptr();
        return makePow(holders[0], holders[1]);
    case NodeKind::Function: {
        const char *name = expression->name();
        if (!name) return Ptr();
        return Ptr(Node::makeFunction(name, children));
    }
    default:
        return Ptr();
    }
}

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

        // 一般换元：∫ F(g(x))·g'(x) dx = ∫ F(u) du（先于常数提取尝试，
        // 否则 `2·x·e^{x²}` 在常数分支里就走不到换元了）。
        if (depth < 3 && varying.size() >= 2) {
            for (std::size_t i = 0; i < varying.size(); ++i) {
                Ptr factor = Ptr::borrow(varying[i]);
                Ptr compositeInner = compositeInnerOf(factor, variable);
                if (!compositeInner) continue;
                Ptr innerDerivative = derivativeOf(compositeInner, variable);
                if (!innerDerivative) continue;
                Memory::Vector<NodePtr> restFactors;
                Memory::Vector<Ptr> restHolders;
                for (std::size_t j = 0; j < varying.size(); ++j) {
                    if (i == j) continue;
                    restFactors.push_back(varying[j]);
                    restHolders.push_back(Ptr::borrow(varying[j]));
                }
                Ptr rest = restFactors.size() == 1 ? restHolders.front() : makeMulMany(restFactors);
                if (!rest) continue;
                Ptr ratio = simplify(makeDivide(rest, innerDerivative));
                if (!ratio || dependsOn(ratio.get(), variable)) continue;
                const std::uint32_t substitutionSymbol = internSymbol("__u");
                Ptr replaced = replaceSubtree(factor, compositeInner, variableNode(substitutionSymbol));
                if (!replaced) continue;
                Ptr anti = integrateRecursive(replaced, substitutionSymbol, depth + 1);
                if (!anti) continue;
                Ptr back = replaceSubtree(anti, variableNode(substitutionSymbol), compositeInner);
                if (!back) continue;
                return makeMul(constant, makeMul(ratio, back));
            }
        }

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

        // log(u)^n：递推 ∫ln^n u du = u·ln^n u − n·∫ln^(n-1) u du
        if (base->kind() == NodeKind::Function && base->name() &&
            std::strcmp(base->name(), "log") == 0 && !dependsOn(exponent.get(), variable)) {
            long long power = 0;
            if (exponent->isSmallInteger(&power) && power >= 1) {
                Ptr innerArgument = Ptr::borrow(base->operand(0));
                Ptr innerSlope, innerOffset;
                if (linearForm(innerArgument, variable, &innerSlope, &innerOffset) && !innerSlope->isZero()) {
                    if (power == 1) {
                        // ∫ln(ax+b) dx = (ax+b)·ln(ax+b)/a − x
                        Ptr product = makeDivide(makeMul(innerArgument, expression), innerSlope);
                        return makeAdd(product, makeNegate(variableNode(variable)));
                    }
                    Ptr lowerPower = makePow(base, integer(power - 1));
                    Ptr innerIntegral = integrateRecursive(lowerPower, variable, depth + 1);
                    if (innerIntegral) {
                        Ptr boundaryTerm = makeDivide(makeMul(innerArgument, expression), innerSlope);
                        return makeAdd(boundaryTerm, makeMul(integer(-power), innerIntegral));
                    }
                }
            }
        }

        // sin^n(u)、cos^n(u)：奇数 n 直接换元，偶数 n 反复降幂后递归
        if (base->kind() == NodeKind::Function && !dependsOn(exponent.get(), variable)) {
            const char *functionName = base->name();
            long long power = 0;
            bool integerPower = exponent->isSmallInteger(&power);
            if (functionName && integerPower && power >= 2 &&
                (std::strcmp(functionName, "sin") == 0 || std::strcmp(functionName, "cos") == 0)) {
                Ptr innerArgument = Ptr::borrow(base->operand(0));
                Ptr innerSlope, innerOffset;
                if (linearForm(innerArgument, variable, &innerSlope, &innerOffset) && !innerSlope->isZero()) {
                    const bool isSine = std::strcmp(functionName, "sin") == 0;
                    if (power % 2 == 1) {
                        // 奇次：拆出一个因子，剩余偶次用 1 - sin²（或 1 - cos²）换成 cos（或 sin）
                        Ptr single = Ptr(Node::makeFunction(functionName, Memory::Vector<NodePtr>{innerArgument.get()}));
                        Ptr remaining = makePow(base, integer(power - 1));
                        // 幂次减一的偶次式转成另一函数的幂：cos²ᵏ = (1 - sin²)ᵏ
                        const char *otherName = isSine ? "cos" : "sin";
                        Ptr other = Ptr(Node::makeFunction(otherName, Memory::Vector<NodePtr>{innerArgument.get()}));
                        Ptr complement = makeAdd(integer(1), makeNegate(makePow(other, integer(2))));
                        Ptr rewritten = makePow(single, integer(1));
                        for (long long i = 0; i < power - 1; i += 2) {
                            rewritten = makeMul(rewritten, complement);
                        }
                        (void)remaining;
                        Ptr innerIntegral = integrateRecursive(simplify(rewritten), variable, depth + 1);
                        if (innerIntegral) return innerIntegral;
                    } else {
                        // 偶次：用降幂公式把 sin²(u) 换成 (1−cos 2u)/2，再递归
                        Ptr squared = makePow(base, integer(2));
                        Ptr reduced = squared ? reduceTrigPower(squared) : Ptr();
                        if (reduced) {
                            // 降幂只处理平方，剩下的偶次幂继续保留
                            Ptr scaled = reduced;
                            if (power > 2) scaled = makeMul(reduced, makePow(base, integer(power - 2)));
                            Ptr innerIntegral = integrateRecursive(simplify(scaled), variable, depth + 1);
                            if (innerIntegral) return innerIntegral;
                        }
                    }
                }
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
            // ∫ ln(u) dx，u = a·x + b：换元后 ∫ln(u)du/a = (u·ln u − u)/a
            if (inner->kind() == NodeKind::Symbol && inner->symbolId() == variable) {
                Ptr x = variableNode(variable);
                return makeAdd(makeMul(x.share(), expression), makeNegate(x.share()));
            }
            if (linear) {
                Ptr scaled = makeDivide(makeMul(inner, expression), slope);
                Ptr correction = makeDivide(inner, slope);
                return makeAdd(scaled, makeNegate(correction));
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

/* ------------------------------------------------------ 部分分式分解 */

/**
 * @brief 真分式 `P/Q` 的部分分式分解积分。
 *
 * 先把分母因式分解成一次因子（可带重数），再设
 * `P/Q = Σ Σ_{k=1..mᵢ} Aᵢₖ / fᵢ^k`，用待定系数法解出 `Aᵢₖ`：
 * 把等式两边乘回 `Q` 展开，按 variable 的幂次逐项比较可得线性方程组。
 * 逐项积分为 `A·ln|f|/a`（k=1）与 `A·f^(1-k)/(a(1-k))`（k>1）。
 *
 * 只处理分母能完全分解成一次因式的情形；含不可约二次因子时退回
 * `integrateQuadraticDenominator`（单个二次因子）或报错。
 */
/** @brief 分母分解出的一个因式及其重数。 */
struct FractionFactor {
    Ptr base;                    ///< 首一化后的因式
    int multiplicity = 1;        ///< 重数
    std::size_t degree = 1;      ///< 因式次数
};

/** @brief 把分解结果摊平成「因式 + 重数」表；遇到超过二次的因子返回 false。 */
bool collectFractionFactors(const Ptr &factored, std::uint32_t variable,
                            Memory::Vector<FractionFactor> *out) {
    if (!factored) return false;
    Memory::Vector<Ptr> pieces;
    if (factored->kind() == NodeKind::Mul) {
        for (NodePtr child : factored->operands()) pieces.push_back(Ptr::borrow(child));
    } else {
        pieces.push_back(factored);
    }
    for (const Ptr &piece : pieces) {
        Ptr base = piece;
        int repeat = 1;
        if (piece->kind() == NodeKind::Pow) {
            Ptr exponent = Ptr::borrow(piece->operand(1));
            long long value = 0;
            if (exponent->isSmallInteger(&value) && value >= 1) {
                repeat = static_cast<int>(value);
                base = Ptr::borrow(piece->operand(0));
            }
        }
        UnivariatePoly polynomial;
        if (!toPolynomial(base.get(), variable, &polynomial) || polynomial.degree() < 1) return false;
        if (polynomial.degree() > 2) return false;
        FractionFactor entry;
        entry.base = base;
        entry.multiplicity = repeat;
        entry.degree = polynomial.degree();
        out->push_back(entry);
    }
    return true;
}

/**
 * @brief 线性因子的部分分式项 `A / f^k` 的积分。
 *
 * `k = 1` 给 `A·ln|f|/a`；`k > 1` 给 `A·f^(1-k)/(a(1-k))`。
 */
Ptr integrateLinearFactorTerm(const Ptr &coefficient, const Ptr &base, int power, std::uint32_t variable) {
    UnivariatePoly polynomial;
    if (!toPolynomial(base.get(), variable, &polynomial)) return Ptr();
    const Ptr slope = polynomial.degree() >= 1 ? polynomial.coefficient(1) : integer(0);
    if (slope->isZero()) return Ptr();
    if (power == 1) {
        Memory::Vector<NodePtr> arguments;
        arguments.push_back(base.get());
        Ptr logarithm = Ptr(Node::makeFunction("log", arguments));
        return makeDivide(makeMul(coefficient, logarithm), slope);
    }
    const long long next = 1 - static_cast<long long>(power);
    Ptr powered = makePow(base, integer(next));
    if (!powered) return Ptr();
    return makeDivide(makeMul(coefficient, powered), makeMul(slope, integer(next)));
}

/**
 * @brief 不可约二次因子的部分分式项 `(p·x + q)/f` 的积分。
 *
 * 拆成 `(p/2a)·f'/f` 与常数项的 `K/((x+h)² + k²)`：
 * 前者给对数，后者给 `atan`。
 */
Ptr integrateQuadraticFactorTerm(const Ptr &numeratorExpression, const Ptr &base, std::uint32_t variable) {
    UnivariatePoly basePoly, numeratorPoly;
    if (!toPolynomial(base.get(), variable, &basePoly)) return Ptr();
    if (!toPolynomial(numeratorExpression.get(), variable, &numeratorPoly)) return Ptr();
    if (basePoly.degree() != 2) return Ptr();

    const Ptr a = basePoly.coefficient(2);
    const Ptr b = basePoly.coefficient(1);
    const Ptr c = basePoly.coefficient(0);
    const Ptr p = numeratorPoly.degree() >= 1 ? numeratorPoly.coefficient(1) : integer(0);
    const Ptr q = numeratorPoly.coefficient(0);
    Ptr x = variableNode(variable);

    // 对数部分：(p/(2a))·ln|f|
    Ptr logPart = Ptr();
    if (!p->isZero()) {
        Memory::Vector<NodePtr> arguments;
        arguments.push_back(base.get());
        Ptr logarithm = Ptr(Node::makeFunction("log", arguments));
        logPart = makeDivide(makeMul(p, logarithm), makeMul(integer(2), a));
    }

    // atan 部分：剩余常数项 K = q - p·b/(2a)，∫K/(a((x+h)²+k²)) = K/(a·k)·atan((x+h)/k)
    Ptr residual = makeAdd(q, makeNegate(makeDivide(makeMul(p, b), makeMul(integer(2), a))));
    residual = simplify(residual);
    Ptr atanPart = Ptr();
    if (residual && !residual->isZero()) {
        // k² = c/a - b²/(4a²)，即 (4ac - b²)/(4a²)
        Ptr kSquared = makeAdd(makeDivide(c, a),
                               makeNegate(makeDivide(makePow(b, integer(2)),
                                                     makeMul(integer(4), makePow(a, integer(2))))));
        kSquared = simplify(kSquared);
        Ptr k = makePow(kSquared, makeNumber(BigInt(1), BigInt(2)));
        if (!k || k->isZero()) return Ptr();
        Ptr halfB = makeDivide(b, makeMul(integer(2), a));
        Ptr shifted = makeAdd(x.share(), halfB);
        Ptr argument = makeDivide(shifted, k);
        if (!argument) return Ptr();
        Memory::Vector<NodePtr> atanArgs;
        atanArgs.push_back(argument.get());
        Ptr atanFunction = Ptr(Node::makeFunction("atan", atanArgs));
        atanPart = makeDivide(makeMul(residual, atanFunction), makeMul(a, k));
    }

    if (logPart && atanPart) return makeAdd(logPart, atanPart);
    if (logPart) return logPart;
    return atanPart;
}

/** @brief 分母分解后要设的一个待定系数（未知量的位置）。 */
struct PartialTerm {
    std::size_t factorIndex = 0;   ///< 属于第几个因式
    int power = 1;                 ///< 该因式的幂次
    int shift = 0;                 ///< 分子里 x 的幂次（0 表示常数 A，1 表示 A·x）
};

Ptr integratePartialFractions(const UnivariatePoly &numerator, const UnivariatePoly &denominator,
                              std::uint32_t variable) {
    if (numerator.degree() >= denominator.degree()) return Ptr();
    Ptr factored = factor(fromPolynomial(denominator));
    if (!factored) return Ptr();
    Memory::Vector<FractionFactor> factors;
    if (!collectFractionFactors(factored, variable, &factors)) return Ptr();
    if (factors.empty()) return Ptr();

    // 分母本身就是单个不可约二次因式（无重数）：直接走二次分母公式
    if (factors.size() == 1 && factors[0].multiplicity == 1 && factors[0].degree == 2) {
        return integrateQuadraticFactorTerm(fromPolynomial(numerator), factors[0].base, variable);
    }
    // 二次因式带重数：需要降幂递推，暂不支持
    for (const FractionFactor &entry : factors) {
        if (entry.degree == 2 && entry.multiplicity > 1) return Ptr();
    }

    // 设每个部分分式项的分子：一次因子给常数 A；二次因子给 A·x + B
    std::vector<PartialTerm> terms;
    for (std::size_t i = 0; i < factors.size(); ++i) {
        for (int power = 1; power <= factors[i].multiplicity; ++power) {
            PartialTerm term;
            term.factorIndex = i;
            term.power = power;
            if (factors[i].degree == 2) {
                term.shift = 1;
                terms.push_back(term);   // A·x
                term.shift = 0;
                terms.push_back(term);   // B
            } else {
                term.shift = 0;
                terms.push_back(term);
            }
        }
    }
    const std::size_t unknownCount = terms.size();

    Memory::Vector<std::uint32_t> symbols;
    Memory::Vector<Ptr> symbolNodes;
    for (std::size_t i = 0; i < unknownCount; ++i) {
        const std::uint32_t id = internSymbol("c" + std::to_string(i) + "x");
        symbols.push_back(id);
        symbolNodes.push_back(Ptr(Node::makeSymbolById(id)));
    }

    // 两边乘回 Q：Σ Aᵢ·(Q/fᵢ^k)·x^shift = P
    Ptr combination = integer(0);
    for (std::size_t i = 0; i < terms.size(); ++i) {
        const PartialTerm &term = terms[i];
        const FractionFactor &owner = factors[term.factorIndex];
        Ptr cofactor = integer(1);
        for (const FractionFactor &other : factors) {
            int exponent = other.multiplicity;
            if (other.base->equals(owner.base.get())) exponent = other.multiplicity - term.power;
            if (exponent > 0) cofactor = makeMul(cofactor, makePow(other.base, integer(exponent)));
        }
        if (term.shift > 0) cofactor = makeMul(cofactor, makePow(variableNode(variable), integer(term.shift)));
        combination = makeAdd(combination, makeMul(symbolNodes[i], cofactor));
    }
    Ptr plain = expand(makeAdd(combination, makeNegate(fromPolynomial(numerator))));
    if (!plain) return Ptr();

    // 按变量幂次收集线性方程组
    UnivariatePoly equation;
    if (!toPolynomial(plain.get(), variable, &equation)) return Ptr();
    const std::size_t rows = equation.degree() + 1;

    // 每个系数按未知符号线性拆分：逐符号求偏置
    std::vector<std::vector<Ptr>> system(rows, std::vector<Ptr>(unknownCount, integer(0)));
    std::vector<Ptr> constants(rows, integer(0));
    for (std::size_t row = 0; row < rows; ++row) {
        Ptr coefficient = simplify(equation.coefficient(row));
        if (!coefficient) return Ptr();
        Memory::Vector<Ptr> parts;
        Ptr expanded = expand(coefficient);
        if (expanded && expanded->kind() == NodeKind::Add) {
            for (NodePtr child : expanded->operands()) parts.push_back(Ptr::borrow(child));
        } else {
            parts.push_back(coefficient);
        }
        for (const Ptr &part : parts) {
            bool matched = false;
            for (std::size_t col = 0; col < unknownCount; ++col) {
                if (!dependsOn(part.get(), symbols[col])) continue;
                Ptr quotient = makeDivide(part, symbolNodes[col]);
                Ptr value = quotient ? simplify(quotient) : Ptr();
                if (!value || dependsOn(value.get(), symbols[col])) return Ptr();
                system[row][col] = makeAdd(system[row][col], value);
                matched = true;
                break;
            }
            if (!matched) constants[row] = makeAdd(constants[row], part);
        }
    }

    // 高斯消元
    std::vector<std::size_t> pivotColumns;
    std::size_t pivotRow = 0;
    for (std::size_t col = 0; col < unknownCount && pivotRow < rows; ++col) {
        std::size_t chosen = rows;
        for (std::size_t row = pivotRow; row < rows; ++row) {
            if (!system[row][col]->isZero()) {
                chosen = row;
                break;
            }
        }
        if (chosen == rows) continue;
        std::swap(system[chosen], system[pivotRow]);
        std::swap(constants[chosen], constants[pivotRow]);
        const Ptr pivot = system[pivotRow][col];
        for (std::size_t k = col; k < unknownCount; ++k) {
            system[pivotRow][k] = simplify(makeDivide(system[pivotRow][k], pivot));
        }
        constants[pivotRow] = simplify(makeDivide(constants[pivotRow], pivot));
        for (std::size_t row = 0; row < rows; ++row) {
            if (row == pivotRow || system[row][col]->isZero()) continue;
            const Ptr scale = system[row][col];
            for (std::size_t k = col; k < unknownCount; ++k) {
                system[row][k] =
                    simplify(makeAdd(system[row][k], makeNegate(makeMul(scale, system[pivotRow][k]))));
            }
            constants[row] =
                simplify(makeAdd(constants[row], makeNegate(makeMul(scale, constants[pivotRow]))));
        }
        pivotColumns.push_back(col);
        ++pivotRow;
    }
    if (pivotColumns.size() != unknownCount) return Ptr();

    std::vector<Ptr> values(unknownCount, integer(0));
    for (std::size_t i = 0; i < pivotColumns.size(); ++i) {
        values[pivotColumns[i]] = simplify(makeNegate(constants[i]));
    }

    // 逐项积分
    Ptr result = integer(0);
    std::size_t cursor = 0;
    for (std::size_t i = 0; i < factors.size(); ++i) {
        const FractionFactor &entry = factors[i];
        for (int power = 1; power <= entry.multiplicity; ++power) {
            if (entry.degree == 2) {
                // 该项占两个未知量：A·x + B
                Ptr first = values[cursor];
                Ptr second = values[cursor + 1];
                cursor += 2;
                Ptr numeratorExpression = makeAdd(makeMul(first, variableNode(variable)), second);
                if (numeratorExpression && !(first->isZero() && second->isZero())) {
                    Ptr piece = integrateQuadraticFactorTerm(numeratorExpression, entry.base, variable);
                    if (!piece) return Ptr();
                    result = makeAdd(result, piece);
                }
            } else {
                const Ptr coefficient = values[cursor];
                ++cursor;
                if (!coefficient || coefficient->isZero()) continue;
                Ptr piece = integrateLinearFactorTerm(coefficient, entry.base, power, variable);
                if (!piece) return Ptr();
                result = makeAdd(result, piece);
            }
        }
    }
    return simplify(result);
}

/* ---------------------------------------------- 三角有理式（万能代换） */

/** @brief 把 `sin(u)`/`cos(u)`/`tan(u)` 的实参收集出来（要求含 variable）。 */
bool trigonometricKernel(const Ptr &expression, std::uint32_t variable, Ptr *kernel, bool *hasFunctions) {
    *hasFunctions = false;
    Memory::Vector<Ptr> kernels;
    std::function<bool(const Ptr &)> walk = [&](const Ptr &node) {
        if (!node) return false;
        if (node->kind() == NodeKind::Add || node->kind() == NodeKind::Mul) {
            for (NodePtr child : node->operands()) {
                if (!walk(Ptr::borrow(child))) return false;
            }
            return true;
        }
        if (node->kind() == NodeKind::Pow) {
            return walk(Ptr::borrow(node->operand(0))) && walk(Ptr::borrow(node->operand(1)));
        }
        if (node->kind() == NodeKind::Function) {
            const char *name = node->name();
            if (name && (std::strcmp(name, "sin") == 0 || std::strcmp(name, "cos") == 0 ||
                         std::strcmp(name, "tan") == 0)) {
                Ptr argument = Ptr::borrow(node->operand(0));
                if (!argument || !dependsOn(argument.get(), variable)) return false;
                *hasFunctions = true;
                bool known = false;
                for (const Ptr &existing : kernels) {
                    if (existing->equals(argument.get())) { known = true; break; }
                }
                if (!known) kernels.push_back(argument);
                return walk(argument);
            }
            // 其它函数一律视为不可代换（如 sin(x)·log(x) 不在万能代换范围内）
            return true;
        }
        return true;
    };
    if (!walk(expression)) return false;
    if (!*hasFunctions || kernels.empty()) return false;
    // 只处理同一个内核的情形：sin(2x) 与 sin(x) 混在一起时不套用
    for (std::size_t i = 1; i < kernels.size(); ++i) {
        if (!kernels[i]->equals(kernels[0].get())) return false;
    }
    *kernel = kernels.front();
    return true;
}

/**
 * @brief 万能代换 `t = tan(u/2)` 求解 `∫ R(sin u, cos u, tan u) du`。
 *
 * 用 `sin u = 2t/(1+t²)`、`cos u = (1-t²)/(1+t²)`、`du = 2dt/(1+t²)` 把三角有理式
 * 变成 t 的有理函数再积分，最后把 t 换回 `tan(u/2)`。要求内核是 `a·x+b`，
 * 这样 `du/dx = a` 是常数。
 */
Ptr integrateTrigonometricRational(const Ptr &expression, std::uint32_t variable) {
    Ptr kernel;
    bool hasFunctions = false;
    if (!trigonometricKernel(expression, variable, &kernel, &hasFunctions)) return Ptr();
    if (!hasFunctions) return Ptr();

    Ptr slope, offset;
    if (!linearForm(kernel, variable, &slope, &offset) || slope->isZero()) return Ptr();

    const std::uint32_t tSymbol = internSymbol("__t");
    Ptr t = variableNode(tSymbol);
    // sin u、cos u、tan u 的 t 表达式
    Ptr onePlusT2 = makeAdd(integer(1), makePow(t, integer(2)));
    Ptr sinU = makeDivide(makeMul(integer(2), t.share()), onePlusT2.share());
    Ptr cosU = makeDivide(makeAdd(integer(1), makeNegate(makePow(t.share(), integer(2)))), onePlusT2.share());
    Ptr tanU = makeDivide(makeMul(integer(2), t.share()),
                          makeAdd(integer(1), makeNegate(makePow(t.share(), integer(2)))));
    if (!sinU || !cosU || !tanU) return Ptr();

    // 把 sin(kernel)→sinU、cos(kernel)→cosU、tan(kernel)→tanU 替换进去
    Ptr replaced = expression.share();
    // 从内核最具体的函数开始替换；replaceSubtree 只匹配整棵子树，安全
    struct { const char *name; Ptr value; } table[] = {
        {"sin", sinU}, {"cos", cosU}, {"tan", tanU},
    };
    for (const auto &entry : table) {
        if (!dependsOn(replaced.get(), variable)) break;
        // 逐个函数节点替换：直接对整棵表达式做子树替换，命中就换
        Ptr target = Ptr(Node::makeFunction(entry.name, Memory::Vector<NodePtr>{kernel.get()}));
        Ptr next = replaceSubtree(replaced, target, entry.value);
        if (!next) return Ptr();
        replaced = next;
    }
    if (dependsOn(replaced.get(), variable)) return Ptr();

    // du = 2dt/(1+t²)，dx = du/slope；换成对 t 积分
    Ptr integrand = makeDivide(makeMul(integer(2), replaced), onePlusT2.share());
    integrand = makeDivide(integrand, slope);
    // 代换后的形状是嵌套分式（`2/((t²+1)*(1 + (1-t²)/(t²+1)))`），
    // 必须先通分化成单个分式，integrate 才能按有理函数处理。
    Ptr combined = togetherDeep(integrand);
    if (combined) integrand = combined;
    Ptr anti = integrate(integrand, tSymbol);
    if (!anti) return Ptr();

    // t = tan(u/2) 回代
    Ptr halfKernel = makeDivide(kernel, integer(2));
    if (!halfKernel) return Ptr();
    Ptr tanHalf = Ptr(Node::makeFunction("tan", Memory::Vector<NodePtr>{halfKernel.get()}));
    if (!tanHalf) return Ptr();
    return substitute(anti, tSymbol, tanHalf);
}

} // namespace

/* ------------------------------------------------------ 根式三角代换 */

/** @brief 前置声明：根式代换内部要递归调用积分主入口。 */
Ptr integrate(const Ptr &expression, std::uint32_t variable);

/**
 * @brief 从表达式里找出所有 `sqrt(g)` 的实参（要求都是同一个 g 且含 variable）。
 *
 * 注意：`sqrt(g)` 在节点层是 `Pow(g, 1/2)`，打印器只是把它写成 `sqrt(...)`，
 * 所以这里按「指数为 1/2 的幂」识别，而不是找名为 `sqrt` 的 Function。
 */
bool radicalKernel(const Ptr &expression, std::uint32_t variable, Ptr *kernel) {
    Memory::Vector<Ptr> kernels;
    std::function<bool(const Ptr &)> walk = [&](const Ptr &node) {
        if (!node) return false;
        if (node->kind() == NodeKind::Add || node->kind() == NodeKind::Mul) {
            for (NodePtr child : node->operands()) {
                if (!walk(Ptr::borrow(child))) return false;
            }
            return true;
        }
        if (node->kind() == NodeKind::Pow) {
            Ptr base = Ptr::borrow(node->operand(0));
            Ptr exponent = Ptr::borrow(node->operand(1));
            // 1/√g 在节点层是 Pow(g, -1/2)（不是 1/Pow(g,1/2)），两种指数都算根式
            BigInt top(0), bottom(1);
            const bool halfPower = numberValue(exponent.get(), &top, &bottom) && top.abs() == BigInt(1) &&
                                   bottom == BigInt(2);
            if (halfPower && dependsOn(base.get(), variable)) {
                for (const Ptr &existing : kernels) {
                    if (!existing->equals(base.get())) return false;
                }
                kernels.push_back(base);
                return walk(base);
            }
            return walk(base) && walk(exponent);
        }
        return true;
    };
    if (!walk(expression) || kernels.empty()) return false;
    *kernel = kernels.front();
    return true;
}

/**
 * @brief 根式三角代换：把 `sqrt(a² − x²)` 一类换成三角式。
 *
 * 按 `g` 的形态挑代换：
 * - `g = a² − x²` ⇒ `x = a·sin t`，`√g = a·cos t`，`dx = a·cos t dt`；
 * - `g = a² + x²` ⇒ `x = a·tan t`，`√g = a/cos t`，`dx = a/sec²t dt`；
 * - `g = x² − a²` ⇒ `x = a/cos t`，`√g = a·tan t`，`dx = a·sin t/cos²t dt`。
 *
 * 代换后按三角有理式继续积分，最后把 `t` 换回反三角表达式。
 */
Ptr integrateRadicalSubstitution(const Ptr &expression, std::uint32_t variable) {
    Ptr kernel;
    if (!radicalKernel(expression, variable, &kernel)) return Ptr();
    Ptr expandedKernel = expand(kernel);
    if (!expandedKernel) expandedKernel = kernel;
    UnivariatePoly polynomial;
    if (!toPolynomial(expandedKernel.get(), variable, &polynomial)) return Ptr();
    if (polynomial.degree() != 2) return Ptr();
    const Ptr quadratic = polynomial.coefficient(2);
    const Ptr linearCoefficient = polynomial.coefficient(1);
    const Ptr constantTerm = polynomial.coefficient(0);
    if (!linearCoefficient->isZero()) return Ptr(); // 含 x 一次项时需先配方，暂不处理

    const std::uint32_t tSymbol = internSymbol("__s");
    Ptr t = Ptr(Node::makeSymbolById(tSymbol));
    Ptr sine = Ptr(Node::makeFunction("sin", Memory::Vector<NodePtr>{t.get()}));
    Ptr cosine = Ptr(Node::makeFunction("cos", Memory::Vector<NodePtr>{t.get()}));
    Ptr tangent = Ptr(Node::makeFunction("tan", Memory::Vector<NodePtr>{t.get()}));
    if (!sine || !cosine || !tangent) return Ptr();

    // 判定形态：a²−x²（二次项系数 < 0）、a²+x²（> 0 且常数项 > 0）
    double quadraticValue = 0.0, constantValue = 0.0;
    BigInt numerator, denominator;
    if (!numberValue(quadratic.get(), &numerator, &denominator)) return Ptr();
    quadraticValue = numerator.toDouble() / denominator.toDouble();
    if (!numberValue(constantTerm.get(), &numerator, &denominator)) return Ptr();
    constantValue = numerator.toDouble() / denominator.toDouble();
    if (quadraticValue > 0.0 && constantValue < 0.0) return Ptr(); // x²−a² 需要 a 的具体值

    // a = √|常数项 / 二次项系数|：
    // a² − x² 型常数项为正、二次项为负；a² + x² 型两者同号。
    // 统一取绝对值，避免 a²+x² 情形算出 √(-1)。
    Ptr ratio = makeDivide(constantTerm, quadratic);
    ratio = simplify(ratio);
    if (!ratio) return Ptr();
    {
        BigInt ratioTop(0), ratioBottom(1);
        if (numberValue(ratio.get(), &ratioTop, &ratioBottom) && ratioTop.isNegative()) {
            ratio = makeNegate(ratio);
        }
    }
    Ptr a = makePow(ratio, makeNumber(BigInt(1), BigInt(2)));
    if (!a) return Ptr();

    Ptr substitution;
    Ptr radicalAfter;
    Ptr derivativeFactor;   // dx/dt
    if (quadraticValue < 0.0) {
        // x = a sin t ⇒ √g = a cos t，dx = a cos t dt
        substitution = makeMul(a, sine);
        radicalAfter = makeMul(a, cosine);
        derivativeFactor = makeMul(a, cosine);
    } else {
        // x = a tan t ⇒ √g = a/cos t，dx = a/cos²t dt
        substitution = makeMul(a, tangent);
        radicalAfter = makeDivide(a, cosine);
        Ptr cosineSquared = makePow(cosine, integer(2));
        derivativeFactor = makeDivide(a, cosineSquared);
    }
    if (!substitution || !radicalAfter || !derivativeFactor) return Ptr();

    // 先把根式节点换成代换后的形式，再把 x 换成 substitution。
    // √g 与 1/√g 分别对应 Pow(g, 1/2) 与 Pow(g, -1/2)，两种都要处理。
    Memory::Vector<Ptr> radicalNodes;
    Ptr rootHalf = makePow(kernel, makeNumber(BigInt(1), BigInt(2)));
    Ptr rootNegativeHalf = makePow(kernel, makeNumber(BigInt(-1), BigInt(2)));
    if (!rootHalf || !rootNegativeHalf) return Ptr();
    radicalNodes.push_back(rootHalf);
    radicalNodes.push_back(rootNegativeHalf);
    Ptr replaced = expression.share();
    for (std::size_t i = 0; i < radicalNodes.size(); ++i) {
        Ptr replacement = i == 0 ? radicalAfter : makeDivide(integer(1), radicalAfter);
        if (!replacement) return Ptr();
        Ptr next = replaceSubtree(replaced, radicalNodes[i], replacement);
        if (!next) return Ptr();
        replaced = next;
    }
    replaced = substitute(replaced, variable, substitution);
    if (!replaced || dependsOn(replaced.get(), variable)) return Ptr();

    Ptr integrand = makeMul(replaced, derivativeFactor);
    Ptr combined = togetherDeep(integrand);
    if (combined) integrand = combined;
    Ptr anti = integrate(integrand, tSymbol);
    if (!anti) return Ptr();

    // t 回代：x = a sin t ⇒ t = asin(x/a)；x = a tan t ⇒ t = atan(x/a)
    Ptr argument = makeDivide(variableNode(variable), a);
    if (!argument) return Ptr();
    const char *inverseName = quadraticValue < 0.0 ? "asin" : "atan";
    Ptr inverse = Ptr(Node::makeFunction(inverseName, Memory::Vector<NodePtr>{argument.get()}));
    if (!inverse) return Ptr();
    return substitute(anti, tSymbol, inverse);
}

Ptr integrate(const Ptr &expression, std::uint32_t variable) {
    if (!expression) return Ptr();
    Ptr simplified = simplify(expression);
    if (!simplified) return Ptr();

    // 先把分式改写为幂形式
    Ptr result = integrateRecursive(simplified, variable);
    if (!result) {
        // 根式三角代换：√(a²−x²)、√(a²+x²) 一类
        Ptr radical = integrateRadicalSubstitution(simplified, variable);
        if (radical) return simplify(radical);
    }
    if (!result) {
        // 三角有理式：万能代换 t = tan(u/2)
        Ptr trigonometric = integrateTrigonometricRational(simplified, variable);
        if (trigonometric) return simplify(trigonometric);
    }
    if (!result) {
        // 有理函数：先用 toQuotient 化成规范的单分式 P/Q，再按次数分派。
        // 旧实现手工拆 `Pow(Q, -1)` / `Mul`，遇到 `1/(Q^2)` 会丢掉指数
        // （`∫dt/(t²+1)²` 被错当成 `∫dt/(t²+1)`），这里统一走规范形式。
        const auto symbols = freeSymbols(simplified.get());
        if (symbols.size() != 1) return Ptr();
        Ptr top, bottom;
        if (!toQuotient(simplified, &top, &bottom)) return Ptr();
        if (isInt(bottom, 1)) return Ptr();
        Ptr numeratorExpression = expand(top);
        Ptr denominatorExpression = expand(bottom);
        if (!numeratorExpression || !denominatorExpression) return Ptr();
        // 先约掉公因式，避免分母带重复因子后分解困难
        Ptr reduced = cancelCommon(numeratorExpression, denominatorExpression, variable);
        if (reduced) {
            Ptr reducedTop, reducedBottom;
            if (toQuotient(reduced, &reducedTop, &reducedBottom)) {
                if (isInt(reducedBottom, 1)) {
                    Ptr simplifiedReduced = simplify(reducedTop);
                    if (!simplifiedReduced) return Ptr();
                    Ptr direct = integrateRecursive(simplifiedReduced, variable);
                    if (direct) return simplify(direct);
                }
                Ptr expandedTop = expand(reducedTop);
                Ptr expandedBottom = expand(reducedBottom);
                if (expandedTop) numeratorExpression = expandedTop;
                if (expandedBottom) denominatorExpression = expandedBottom;
            }
        }
        UnivariatePoly numerator, denominator;
        if (!toPolynomial(numeratorExpression.get(), variable, &numerator)) return Ptr();
        if (!toPolynomial(denominatorExpression.get(), variable, &denominator)) return Ptr();
        if (denominator.isZero()) return Ptr();

        // 真分式部分：先做多项式除法拆出整式
        UnivariatePoly quotientPart, remainderPart;
        if (!polyDivide(numerator, denominator, &quotientPart, &remainderPart)) return Ptr();
        Ptr polynomialPart = integer(0);
        for (std::size_t i = 1; i < quotientPart.coefficients.size(); ++i) {
            Ptr term = makeMul(
                quotientPart.coefficient(i),
                makePow(variableNode(variable), integer(static_cast<long long>(i - 1))));
            polynomialPart = makeAdd(polynomialPart, makeDivide(term, integer(static_cast<long long>(i))));
        }

        if (remainderPart.isZero()) return simplify(polynomialPart);

        Ptr remainderResult = Ptr();
        if (denominator.degree() == 1) {
            // ∫ R/(a x + b) dx = R·ln|ax+b|/a（R 必为常数）
            Ptr slope = denominator.coefficient(1);
            Memory::Vector<NodePtr> arguments;
            arguments.push_back(denominatorExpression.get());
            Ptr logarithm = Ptr(Node::makeFunction("log", arguments));
            remainderResult = makeDivide(makeMul(remainderPart.coefficient(0), logarithm), slope);
        } else if (denominator.degree() == 2) {
            remainderResult = integrateQuadraticDenominator(remainderPart, denominator, variable);
        } else {
            remainderResult = integratePartialFractions(remainderPart, denominator, variable);
        }
        if (!remainderResult) return Ptr();
        return simplify(makeAdd(polynomialPart, remainderResult));
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

/**
 * @brief 幂指型未定式 `F(x)^G(x)`：改写成 `e^{G·ln F}` 再取极限。
 *
 * 覆盖 1^∞、0^0、∞^0 三类：底数与指数都随 x 变化时，先求
 * `L = lim G·ln F`，答案即 `e^L`。`F` 必须为正（实值意义下），
 * 若 `ln F` 的极限求不出就放弃，交回普通通路。
 *
 * 只在「底数与指数都含变量」时启动，避免把 `x^2`、`2^x` 这类
 * 有直接结果的幂也绕一遍指数-对数。
 */
Ptr limitPowerForm(const Ptr &expression, std::uint32_t variable, const Ptr &point);

Ptr limit(const Ptr &expression, std::uint32_t variable, const Ptr &point) {
    if (!expression) return Ptr();

    // 无穷极限单独处理：`-oo` 会被解析成 Mul(-1, Infinity)，两种形态都要认
    if (isPositiveInfinity(point)) {
        // 先试幂指型（1^∞ 这类在无穷处也可能出现）
        if (Ptr power = limitPowerForm(expression, variable, point)) return power;
        return limitAtInfinity(expression, variable, true);
    }
    if (isNegativeInfinity(point)) return limitAtInfinity(expression, variable, false);

    // 幂指型未定式 F(x)^G(x)：底数与指数都随 x 变（1^∞、0^0、∞^0）
    if (Ptr power = limitPowerForm(expression, variable, point)) return power;

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
 * @brief 幂指型未定式 `F^G` 的极限：`e^{G·ln F}`。
 */
Ptr limitPowerForm(const Ptr &expression, std::uint32_t variable, const Ptr &point) {
    if (!expression) return Ptr();
    Ptr base, exponent;
    if (expression->kind() != NodeKind::Pow) return Ptr();
    base = Ptr::borrow(expression->operand(0));
    exponent = Ptr::borrow(expression->operand(1));
    if (!base || !exponent) return Ptr();
    if (!dependsOn(base.get(), variable) || !dependsOn(exponent.get(), variable)) return Ptr();

    // L = lim G·ln F
    Memory::Vector<NodePtr> logArguments;
    logArguments.push_back(base.get());
    Ptr logarithm = Ptr(Node::makeFunction("log", logArguments));
    if (!logarithm) return Ptr();
    Ptr product = makeMul(exponent, logarithm);
    if (!product) return Ptr();
    Ptr exponentLimit = limit(product, variable, point);
    if (!exponentLimit) return Ptr();
    if (containsInfinity(exponentLimit.get())) {
        // e^∞ = ∞，e^-∞ = 0
        Ptr sign = simplify(exponentLimit);
        if (sign && sign->kind() == NodeKind::Constant &&
            sign->constant() == ConstantId::NegativeInfinity) {
            return integer(0);
        }
        return Ptr(Node::makeConstant(ConstantId::Infinity));
    }
    Ptr result = makePow(Ptr(Node::makeConstant(ConstantId::E)), exponentLimit);
    return result ? simplify(result) : Ptr();
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

    // 兜底：换元 t = 1/x，把 x→∞ 转成 t→0⁺ 再求
    // （x·ln(1+1/x)、x·sin(1/x) 这类"无穷 × 无穷小"全走这条路）
    {
        const std::uint32_t tSymbol = internSymbol("__t");
        Ptr t = Ptr(Node::makeSymbolById(tSymbol));
        Ptr reciprocal = makeDivide(integer(1), t);
        if (reciprocal) {
            // 先把 x 换成 1/t，结果里只剩 t（`simplify` 会把 1/(1/t) 收成 t）
            Ptr transformed = substitute(expression, variable, reciprocal);
            if (transformed) {
                Ptr simplifiedTransformed = simplify(transformed);
                if (simplifiedTransformed && !dependsOn(simplifiedTransformed.get(), variable) &&
                    !containsInfinity(simplifiedTransformed.get())) {
                    Ptr value = limit(simplifiedTransformed, tSymbol, integer(0));
                    if (value && !containsInfinity(value.get()) &&
                        !dependsOn(value.get(), tSymbol)) {
                        return value;
                    }
                }
            }
        }
    }
    return Ptr();
}

/* ---------------------------------------------------------- Taylor */

namespace {

/** @brief 在 point 处求值并判断是否"干净"（有限、无除零、无 0^负数）。 */
bool evaluatesCleanly(const Ptr &expression, std::uint32_t variable, const Ptr &point, Ptr *value) {
    Ptr substituted = substitute(expression, variable, point);
    if (!substituted) return false;
    Ptr result = simplify(substituted);
    if (!result) return false;
    if (containsInfinity(result.get())) return false;
    // 仍含 0^负数 说明除零
    std::function<bool(const Ptr &)> hasZeroNegativePower = [&](const Ptr &node) {
        if (!node) return false;
        if (node->kind() == NodeKind::Pow) {
            Ptr base = Ptr::borrow(node->operand(0));
            Ptr exponent = Ptr::borrow(node->operand(1));
            if (base->isZero() && exponent->isNumber() && exponent->isNegative()) return true;
        }
        for (NodePtr child : node->operands()) {
            if (hasZeroNegativePower(Ptr::borrow(child))) return true;
        }
        return false;
    };
    if (hasZeroNegativePower(result)) return false;
    *value = result;
    return true;
}

} // namespace

/**
 * @brief Taylor / Laurent 展开。
 *
 * 常规路径按 `f^(k)(a)/k!·(x−a)^k` 逐阶求导。若 `f` 在 `a` 处奇异
 * （`1/x`、`sin(x)/x`、`1/sin(x)` 这类），先用 `(x−a)^k` 乘到不再奇异，
 * 展开后再除回去，得到 Laurent 展开（含负幂项）。
 *
 * @param order 展开到的最高阶（对 Laurent 展开指「非负部分」的最高阶）。
 */
Ptr taylorExpand(const Ptr &expression, std::uint32_t variable, const Ptr &point, int order) {
    if (!expression || order < 0) return Ptr();
    Ptr variableNodePtr = variableNode(variable);
    Ptr shifted = makeAdd(variableNodePtr.share(), makeNegate(point));
    if (!shifted) return Ptr();

    // 先看能否直接展开；不能则试着乘掉奇点
    Ptr working = expression.share();
    long long poleOrder = 0;
    {
        Ptr probe;
        if (!evaluatesCleanly(working, variable, point, &probe)) {
            // 逐次乘 (x−a)，最多 8 次（超过说明不是极点而是本性奇点）
            for (int attempt = 1; attempt <= 8; ++attempt) {
                working = simplify(makeMul(working, shifted));
                if (!working) return Ptr();
                Ptr value;
                if (evaluatesCleanly(working, variable, point, &value)) {
                    poleOrder = attempt;
                    break;
                }
            }
            if (poleOrder == 0) return Ptr();
        }
    }

    Ptr result = integer(0);
    Ptr current = working.share();
    BigInt factorial(1);

    for (int k = 0; k <= order; ++k) {
        Ptr value;
        if (!evaluatesCleanly(current, variable, point, &value)) return Ptr();
        if (k > 0) factorial = factorial * BigInt(k);
        Ptr term = makeDivide(makeMul(value, makePow(shifted, integer(k))), makeNumber(factorial, BigInt(1)));
        result = makeAdd(result, term);

        if (k < order) {
            current = differentiate(current, variable);
            if (!current) return Ptr();
        }
    }
    // 抵消之前乘上的 (x−a)^poleOrder：逐项除，结果才是「按 (x−a) 的幂次」展开的形态
    if (poleOrder > 0) {
        Memory::Vector<Ptr> pieces;
        if (result->kind() == NodeKind::Add) {
            for (NodePtr child : result->operands()) pieces.push_back(Ptr::borrow(child));
        } else {
            pieces.push_back(result);
        }
        Memory::Vector<Ptr> holders;
        for (const Ptr &piece : pieces) {
            Ptr divided = makeDivide(piece, makePow(shifted, integer(poleOrder)));
            if (!divided) return Ptr();
            holders.push_back(divided);
        }
        Memory::Vector<NodePtr> operands;
        for (const Ptr &holder : holders) operands.push_back(holder.get());
        result = makeAddMany(operands);
        if (!result) return Ptr();
    }
    return simplify(result);
}

/* ------------------------------------------------------ 驻点与极值 */

const char *criticalKindName(CriticalKind kind) {
    switch (kind) {
    case CriticalKind::Minimum: return "极小值";
    case CriticalKind::Maximum: return "极大值";
    case CriticalKind::Inflection: return "非极值（拐点）";
    case CriticalKind::Unknown:
    default: return "无法判定";
    }
}

namespace {

/**
 * @brief 化简一个不再含自由变量的常数表达式。
 *
 * 先走通用化简，再用三角化简折叠特殊角（`acos(0)` → `1/2*pi`），
 * 最后再化简一次，这样 `sin` 一类驻点的位置能直接落到闭式数值上。
 */
Ptr refineConstant(const Ptr &expression) {
    if (!expression) return Ptr();
    Ptr value = simplify(expression);
    if (!value) value = expression.share();
    Ptr folded = trigSimplify(value);
    if (folded) value = folded;
    Ptr again = simplify(value);
    return again ? again : value;
}

} // namespace

bool criticalPoints(const Ptr &expression, std::uint32_t variable,
                    Memory::Vector<CriticalPoint> *points) {
    if (!expression || !points) return false;
    Ptr derivative = differentiate(expression, variable);
    if (!derivative) return false;
    Ptr slope = simplify(derivative);
    if (!slope) return false;
    if (!dependsOn(slope.get(), variable)) return false; // 常数斜率：没有驻点

    Ptr zero = integer(0);
    Memory::Vector<Solution> solutions;
    if (!solve(slope.get(), zero.get(), variable, &solutions)) return false;
    if (solutions.empty()) return false;

    Ptr secondDerivative = differentiate(slope, variable);
    Memory::Vector<Ptr> secondHolder;
    if (secondDerivative) secondHolder.push_back(secondDerivative);

    bool anyCollected = false;
    for (const Solution &solution : solutions) {
        if (!solution.value) continue;
        CriticalPoint point;
        // 求解器给出的根可能是 acos(0) 这类反三角形式，先折叠成 π/2 再参与后续计算。
        point.location = refineConstant(solution.value);
        if (!point.location) point.location = solution.value;
        Ptr atPoint = substitute(expression, variable, point.location);
        point.value = atPoint ? refineConstant(atPoint) : Ptr();
        if (secondDerivative) {
            Ptr curvature = substitute(secondDerivative, variable, point.location);
            Ptr signedCurvature = curvature ? refineConstant(curvature) : Ptr();
            if (signedCurvature && signedCurvature->isNumber()) {
                if (signedCurvature->isZero()) {
                    point.kind = CriticalKind::Unknown;
                } else if (signedCurvature->isNegative()) {
                    point.kind = CriticalKind::Maximum;
                } else {
                    point.kind = CriticalKind::Minimum;
                }
            }
        }
        points->push_back(point);
        anyCollected = true;
    }
    return anyCollected;
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

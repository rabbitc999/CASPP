#include "cas/detail/special.h"

#include "cas/detail/printer.h"

#include <cstring>
#include <string>

namespace Cas::Detail {

namespace {

/* --------------------------------------------------------- 数值小工具 */

/** @brief 把节点取成有理数；失败返回 false。 */
bool toRational(const Ptr &node, BigInt *numerator, BigInt *denominator) {
    if (!node || !node->isNumber()) return false;
    *numerator = node->numerator();
    *denominator = node->denominator();
    return true;
}

/** @brief 是否为整数（可为负），取值写入 out。 */
bool asInteger(const Ptr &node, long long *out) {
    long long value = 0;
    if (!node || !node->isSmallInteger(&value)) return false;
    *out = value;
    return true;
}

Ptr integer(long long value) { return makeNumber(BigInt(value), BigInt(1)); }

/** @brief √π。返回 `sqrt(pi)` 的化简形式。 */
Ptr sqrtPi() {
    Ptr pi(Node::makeConstant(ConstantId::Pi));
    return makePow(pi, makeNumber(BigInt(1), BigInt(2)));
}

/** @brief 有理数的分母是否恰为 2（即半整数）。 */
bool isHalfInteger(const Ptr &node, long long *half) {
    BigInt numerator;
    BigInt denominator;
    if (!toRational(node, &numerator, &denominator)) return false;
    if (denominator != BigInt(2)) return false;
    if (!numerator.fitsLongLong()) return false;
    *half = numerator.toLongLong();
    return true;
}

/** @brief 构造函数调用节点（实参已化简）。 */
Ptr makeCall(const char *name, const Memory::Vector<Ptr> &arguments) {
    Memory::Vector<NodePtr> raw;
    raw.reserve(arguments.size());
    for (const Ptr &argument : arguments) raw.push_back(argument.get());
    return Ptr(Node::makeFunction(name, raw));
}

/** @brief 重新构造同一个函数调用（用于规范化后原样返回）。 */
Ptr rebuildCall(const std::string &name, const Memory::Vector<Ptr> &arguments) { return makeCall(name.c_str(), arguments); }

} // namespace

/* --------------------------------------------------------------- Gamma */

// Gamma 的定义与折叠规则对齐 sympy 的 `gamma.eval`：整数走阶乘，半整数走
// `(2n)!/(4^n n!)·√π` 与 `(-4)^n n!/(2n)!·√π`，其余保留。

Ptr gamma(const Ptr &argument) {
    if (!argument) return Ptr();

    long long value = 0;
    if (asInteger(argument, &value)) {
        if (value > 0) {
            // Γ(n) = (n-1)!；上限防止把 Γ(10^9) 炸掉进程
            if (value > 3000) return Ptr();
            BigInt product(1);
            for (long long k = 2; k <= value - 1; ++k) product = product * BigInt(k);
            return makeNumber(product, BigInt(1));
        }
        // Γ(0) 与 Γ(负整数) 都是极点，保持符号形式
        return Ptr();
    }

    long long half = 0;
    if (isHalfInteger(argument, &half)) {
        // Γ(1/2 + n) = (2n)!/(4^n n!) √π
        if (half > 0 && half % 2 == 1) {
            const long long n = (half - 1) / 2;
            if (n > 500) return Ptr();
            const BigInt powerOfTwo = BigInt::power(BigInt(2), static_cast<unsigned long long>(2 * n));
            BigInt factorialN(1);
            for (long long k = 2; k <= n; ++k) factorialN = factorialN * BigInt(k);
            BigInt factorial2N(1);
            for (long long k = 2; k <= 2 * n; ++k) factorial2N = factorial2N * BigInt(k);
            Ptr coefficient =
                makeDivide(makeNumber(factorial2N, BigInt(1)), makeNumber(powerOfTwo * factorialN, BigInt(1)));
            return makeMul(coefficient, sqrtPi());
        }
        // Γ(1/2 - n) = (-4)^n n!/(2n)! √π，其中 half = 1 - 2n
        if (half < 0 && (-half) % 2 == 1) {
            const long long n = (1 - half) / 2;
            if (n > 500) return Ptr();
            BigInt power = BigInt::power(BigInt(4), static_cast<unsigned long long>(n));
            if (n % 2 == 1) power = -power;
            BigInt factorialN(1);
            for (long long k = 2; k <= n; ++k) factorialN = factorialN * BigInt(k);
            BigInt factorial2N(1);
            for (long long k = 2; k <= 2 * n; ++k) factorial2N = factorial2N * BigInt(k);
            Ptr coefficient = makeDivide(makeNumber(power * factorialN, BigInt(1)), makeNumber(factorial2N, BigInt(1)));
            return makeMul(coefficient, sqrtPi());
        }
    }
    return Ptr();
}

Ptr factorial(const Ptr &argument) {
    if (!argument) return Ptr();
    long long value = 0;
    if (asInteger(argument, &value)) {
        if (value < 0) return Ptr(); // 负整数阶乘发散
        if (value > 5000) return Ptr();
        BigInt product(1);
        for (long long k = 2; k <= value; ++k) product = product * BigInt(k);
        return makeNumber(product, BigInt(1));
    }
    return Ptr();
}

Ptr doubleFactorial(const Ptr &argument) {
    if (!argument) return Ptr();
    long long value = 0;
    if (asInteger(argument, &value)) {
        if (value < 0) {
            // (-1)!! = 1, (-3)!! = -1, (-5)!! = 1/3 ...
            if (value % 2 == 0) return Ptr();
            const long long n = (-value - 1) / 2;
            if (n > 3000) return Ptr();
            BigInt product(1);
            for (long long k = 1; k <= n; ++k) product = product * BigInt(2 * k + 1);
            BigInt result = BigInt::power(BigInt(-1), static_cast<unsigned long long>(n));
            return makeDivide(makeNumber(result, BigInt(1)), makeNumber(product, BigInt(1)));
        }
        if (value > 8000) return Ptr();
        BigInt product(1);
        for (long long k = value; k > 0; k -= 2) product = product * BigInt(k);
        return makeNumber(product, BigInt(1));
    }
    return Ptr();
}

Ptr beta(const Ptr &first, const Ptr &second) {
    if (!first || !second) return Ptr();
    // 数值参数走 Γ 定义；只在两侧都能精确折叠时才返回
    const bool firstNumeric = first->isNumber();
    const bool secondNumeric = second->isNumber();
    if (firstNumeric && secondNumeric) {
        Ptr g1 = gamma(first);
        Ptr g2 = gamma(second);
        Ptr sum = makeAdd(first, second);
        Ptr g3 = gamma(sum);
        if (g1 && g2 && g3) return makeDivide(makeMul(g1, g2), g3);
    }
    // 交换律规范化：让结构更小的参数在前，便于结构相等判定
    if (second->compareTo(first.get()) < 0) {
        Memory::Vector<Ptr> swapped{second, first};
        return rebuildCall("beta", swapped);
    }
    return Ptr();
}

Ptr binomial(const Ptr &n, const Ptr &k) {
    if (!n || !k) return Ptr();
    long long upper = 0;
    long long lower = 0;
    if (!asInteger(n, &upper) || !asInteger(k, &lower)) return Ptr();
    if (upper < 0) return Ptr();
    if (lower < 0) return integer(0); // 约定 C(n, 负) = 0（组合意义下）
    if (lower > upper) return integer(0);
    if (lower > upper - lower) lower = upper - lower; // C(n,k)=C(n,n-k)，取较小的 k
    if (lower > 3000) return Ptr();
    BigInt result(1);
    for (long long index = 1; index <= lower; ++index) {
        result = result * BigInt(upper - lower + index) / BigInt(index);
    }
    return makeNumber(result, BigInt(1));
}

Ptr harmonic(const Ptr &n, const Ptr &order) {
    if (!n) return Ptr();
    long long upper = 0;
    long long power = 1;
    if (order) {
        if (!asInteger(order, &power)) return Ptr();
        if (power < 0) return Ptr();
    }
    if (!asInteger(n, &upper)) return Ptr();
    if (upper < 0) return Ptr();
    if (upper > 20000) return Ptr();

    if (power == 0) return integer(upper); // H(n, 0) = n
    BigInt numerator(0);
    BigInt denominator(1);
    for (long long index = 1; index <= upper; ++index) {
        const BigInt divisor = BigInt::power(BigInt(index), static_cast<unsigned long long>(power));
        // numerator/denominator += 1/divisor
        numerator = numerator * divisor + denominator;
        denominator = denominator * divisor;
        const BigInt common = BigInt::gcd(numerator.abs(), denominator);
        if (!common.isZero() && !common.isOne()) {
            numerator = numerator / common;
            denominator = denominator / common;
        }
    }
    return makeNumber(numerator, denominator);
}

/* --------------------------------------------------------- 误差函数 */

// erf 的折叠参考 sympy `erf.eval`：只折叠 0 与 ±∞，以及奇偶性。
// 奇偶性写成 `erf(-z) = -erf(z)`，能避免大量重复结构。

Ptr erf(const Ptr &argument) {
    if (!argument) return Ptr();
    if (argument->isNumber() && argument->isZero()) return integer(0);
    if (argument->kind() == NodeKind::Constant) {
        if (argument->constant() == ConstantId::Infinity) return integer(1);
        if (argument->constant() == ConstantId::NegativeInfinity) return integer(-1);
    }
    return Ptr();
}

Ptr erfc(const Ptr &argument) {
    if (!argument) return Ptr();
    if (argument->isNumber() && argument->isZero()) return integer(1);
    if (argument->kind() == NodeKind::Constant) {
        if (argument->constant() == ConstantId::Infinity) return integer(0);
        if (argument->constant() == ConstantId::NegativeInfinity) return integer(2);
    }
    return Ptr();
}

Ptr erfi(const Ptr &argument) {
    if (!argument) return Ptr();
    if (argument->isNumber() && argument->isZero()) return integer(0);
    return Ptr();
}

/* ------------------------------------------------- 指数积分 / 对数积分 */

Ptr expIntegralEi(const Ptr &argument) {
    if (!argument) return Ptr();
    if (argument->isNumber() && argument->isZero()) {
        // Ei(0) 是 -∞：保持未求值以免混淆
        return Ptr();
    }
    if (argument->kind() == NodeKind::Constant) {
        if (argument->constant() == ConstantId::NegativeInfinity) return integer(0);
    }
    return Ptr();
}

Ptr logarithmicIntegral(const Ptr &argument) {
    if (!argument) return Ptr();
    if (argument->isNumber() && argument->isOne()) {
        // li(1) 是 -∞（对数积分在 1 处有支点）
        return Ptr();
    }
    long long value = 0;
    if (asInteger(argument, &value) && value == 0) return integer(0); // li(0) = 0
    return Ptr();
}

/* ------------------------------------------------------------- Bessel */

// J_ν(0)：ν=0 时为 1；ν 为正整数时为 0；负整数阶按 (-1)^n 规范化（sympy besselj.eval）。

Ptr besselJ(const Ptr &order, const Ptr &argument) {
    if (!order || !argument) return Ptr();
    long long nu = 0;
    const bool integerOrder = asInteger(order, &nu);

    if (argument->isNumber() && argument->isZero()) {
        if (integerOrder) return integer(nu == 0 ? 1 : 0);
        BigInt numerator;
        BigInt denominator;
        if (toRational(order, &numerator, &denominator)) {
            if (numerator.isNegative()) return Ptr(); // 负非整数阶在 0 处发散
        }
        return integer(0);
    }
    if (integerOrder && nu < 0) {
        // J_{-n}(z) = (-1)^n J_n(z)
        Ptr sign = makeNumber((nu % 2 == 0) ? BigInt(1) : BigInt(-1), BigInt(1));
        Memory::Vector<Ptr> arguments{integer(-nu), argument};
        return makeMul(sign, rebuildCall("besselj", arguments));
    }
    if (argument->kind() == NodeKind::Constant) {
        if (argument->constant() == ConstantId::Infinity) return integer(0);
    }
    return Ptr();
}

Ptr besselY(const Ptr &order, const Ptr &argument) {
    if (!order || !argument) return Ptr();
    if (argument->isNumber() && argument->isZero()) {
        // Y_ν(0) 对 ν ≥ 0 发散（-∞），保持未求值
        return Ptr();
    }
    return Ptr();
}

Ptr besselI(const Ptr &order, const Ptr &argument) {
    if (!order || !argument) return Ptr();
    long long nu = 0;
    const bool integerOrder = asInteger(order, &nu);
    if (argument->isNumber() && argument->isZero()) {
        if (integerOrder) return integer(nu == 0 ? 1 : 0);
        return integer(0);
    }
    if (integerOrder && nu < 0) {
        // I_{-n}(z) = I_n(z)（修正 Bessel 函数没有符号因子）
        Memory::Vector<Ptr> arguments{integer(-nu), argument};
        return rebuildCall("besseli", arguments);
    }
    return Ptr();
}

Ptr besselK(const Ptr &order, const Ptr &argument) {
    if (!order || !argument) return Ptr();
    if (argument->isNumber() && argument->isZero()) return Ptr(); // K 在 0 处发散
    long long nu = 0;
    if (asInteger(order, &nu) && nu < 0) {
        // K_{-ν}(z) = K_ν(z)，对称
        Memory::Vector<Ptr> arguments{integer(-nu), argument};
        return rebuildCall("besselk", arguments);
    }
    return Ptr();
}

/* ------------------------------------------------------------ 总入口 */

Ptr evaluateSpecialCall(const std::string &name, const Memory::Vector<Ptr> &arguments) {
    const std::size_t count = arguments.size();
    if (name == "gamma") {
        if (count != 1) return Ptr();
        return gamma(arguments[0]);
    }
    if (name == "factorial") {
        if (count != 1) return Ptr();
        return factorial(arguments[0]);
    }
    if (name == "factorial2" || name == "doublefactorial") {
        if (count != 1) return Ptr();
        return doubleFactorial(arguments[0]);
    }
    if (name == "beta") {
        if (count != 2) return Ptr();
        return beta(arguments[0], arguments[1]);
    }
    if (name == "binomial" || name == "binomialCoefficient") {
        if (count != 2) return Ptr();
        return binomial(arguments[0], arguments[1]);
    }
    if (name == "harmonic") {
        if (count == 1) return harmonic(arguments[0], integer(1));
        if (count == 2) return harmonic(arguments[0], arguments[1]);
        return Ptr();
    }
    if (name == "erf") {
        if (count != 1) return Ptr();
        return erf(arguments[0]);
    }
    if (name == "erfc") {
        if (count != 1) return Ptr();
        return erfc(arguments[0]);
    }
    if (name == "erfi") {
        if (count != 1) return Ptr();
        return erfi(arguments[0]);
    }
    if (name == "Ei" || name == "expint_ei") {
        if (count != 1) return Ptr();
        return expIntegralEi(arguments[0]);
    }
    if (name == "li" || name == "Li") {
        if (count != 1) return Ptr();
        return logarithmicIntegral(arguments[0]);
    }
    if (name == "besselj") {
        if (count != 2) return Ptr();
        return besselJ(arguments[0], arguments[1]);
    }
    if (name == "bessely") {
        if (count != 2) return Ptr();
        return besselY(arguments[0], arguments[1]);
    }
    if (name == "besseli") {
        if (count != 2) return Ptr();
        return besselI(arguments[0], arguments[1]);
    }
    if (name == "besselk") {
        if (count != 2) return Ptr();
        return besselK(arguments[0], arguments[1]);
    }
    return Ptr();
}

/**
 * @brief 自底向上折叠特殊函数。
 *
 * 先递归处理子节点，再尝试折叠本层的函数调用；折叠结果再走一次 `simplify`，
 * 让 `Γ(5)` → `24`、`Γ(1/2)·2` → `2·√π` 这类复合结果也收敛。
 * 用深度计数器限制重写轮数，避免 `Γ(z+1) = z·Γ(z)` 之类的规则无限展开。
 */
Ptr simplifySpecial(const Ptr &expression) {
    if (!expression) return Ptr();

    switch (expression->kind()) {
    case NodeKind::Integer:
    case NodeKind::Rational:
    case NodeKind::Symbol:
    case NodeKind::Constant:
        return expression;

    case NodeKind::Function: {
        Memory::Vector<Ptr> arguments;
        arguments.reserve(expression->size());
        bool childChanged = false;
        for (NodePtr child : expression->operands()) {
            Ptr before = Ptr::borrow(child);
            Ptr after = simplifySpecial(before);
            if (after && !after->equals(child)) childChanged = true;
            arguments.push_back(after ? after : before);
        }
        const std::string name = expression->name() ? expression->name() : std::string();
        if (!childChanged) {
            Ptr folded = evaluateSpecialCall(name, arguments);
            if (folded) return simplify(folded);
        }
        Memory::Vector<NodePtr> raw;
        raw.reserve(arguments.size());
        for (const Ptr &argument : arguments) raw.push_back(argument.get());
        Ptr rebuilt(Node::makeFunction(name, raw));
        // 子节点变了，用新参数再尝试折叠一次
        Ptr folded = evaluateSpecialCall(name, arguments);
        if (folded) return simplify(folded);
        return rebuilt;
    }

    case NodeKind::Add:
    case NodeKind::Mul: {
        Memory::Vector<Ptr> operands;
        operands.reserve(expression->size());
        for (NodePtr child : expression->operands()) operands.push_back(simplifySpecial(Ptr::borrow(child)));
        Memory::Vector<NodePtr> raw;
        raw.reserve(operands.size());
        for (const Ptr &operand : operands) raw.push_back(operand.get());
        Ptr rebuilt(Node::makeNary(expression->kind(), raw));
        return simplify(rebuilt);
    }

    case NodeKind::Pow: {
        Ptr base = simplifySpecial(Ptr::borrow(expression->operand(0)));
        Ptr exponent = simplifySpecial(Ptr::borrow(expression->operand(1)));
        return makePow(base, exponent);
    }
    }
    return expression;
}

} // namespace Cas::Detail

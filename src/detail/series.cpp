#include "cas/detail/series.h"

#include "cas/detail/calculus.h"
#include "cas/detail/polynomial.h"
#include "cas/detail/printer.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace Cas::Detail {

namespace {

Ptr integer(long long value) { return makeNumber(BigInt(value), BigInt(1)); }

/** @brief 是否为无穷常数节点。 */
bool isInfinity(const Ptr &node, bool *positive = nullptr) {
    if (!node || node->kind() != NodeKind::Constant) return false;
    if (node->constant() == ConstantId::Infinity) {
        if (positive) *positive = true;
        return true;
    }
    if (node->constant() == ConstantId::NegativeInfinity) {
        if (positive) *positive = false;
        return true;
    }
    return false;
}

/** @brief 表达式是否含该变量。 */
bool hasVariable(const Ptr &expression, std::uint32_t variable) {
    return expression && dependsOn(expression.get(), variable);
}

/* ------------------------------------------------------ 伯努利数 */

/**
 * @brief 伯努利数 B_n，约定 B_1 = -1/2。
 *
 * 用递推 `B_m = -1/(m+1) · Σ_{k=0}^{m-1} C(m+1, k) B_k`：
 * B_0 = 1, B_1 = -1/2, B_2 = 1/6, B_3 = 0, B_4 = -1/30, ...
 * 奇数下标（除 1）恒为 0，这里按下标缓存避免重复计算。
 */
Memory::Vector<Ptr> &bernoulliCache() {
    static Memory::Vector<Ptr> cache;
    return cache;
}

Ptr bernoulliAt(std::size_t n) {
    Memory::Vector<Ptr> &cache = bernoulliCache();
    if (n < cache.size() && cache[n]) return cache[n];
    // 逐级补齐，保证前驱已就绪
    if (cache.empty()) cache.push_back(integer(1));

    while (cache.size() <= n) {
        const std::size_t m = cache.size(); // 要计算的 B_m
        BigInt sum(0);
        BigInt sumDenominator(1);
        for (std::size_t k = 0; k < m; ++k) {
            // C(m+1, k) * B_k
            BigInt binom(1);
            for (std::size_t j = 1; j <= k; ++j) {
                binom = binom * BigInt(static_cast<long long>(m + 1 - k + j)) / BigInt(static_cast<long long>(j));
            }
            const Ptr &bk = cache[k];
            const BigInt termNumerator = binom * bk->numerator();
            const BigInt termDenominator = bk->denominator();
            // sum += term
            sum = sum * termDenominator + termNumerator * sumDenominator;
            sumDenominator = sumDenominator * termDenominator;
            const BigInt common = BigInt::gcd(sum.abs(), sumDenominator);
            if (!common.isZero() && !common.isOne()) {
                sum = sum / common;
                sumDenominator = sumDenominator / common;
            }
        }
        // B_m = -sum / (m+1)
        sum = -sum;
        sumDenominator = sumDenominator * BigInt(static_cast<long long>(m + 1));
        cache.push_back(makeNumber(sum, sumDenominator));
    }
    return cache[n];
}

/* -------------------------------------------------- 多项式求和 */

} // namespace

/* ------------------------------------------------------------ 伯努利 */

Ptr bernoulliNumber(long long n) {
    if (n < 0) return Ptr();
    if (n > 400) return Ptr(); // 再大就没意义了，计算量也爆炸
    return bernoulliAt(static_cast<std::size_t>(n));
}

namespace {

/**
 * @brief 幂和 `Σ_{i=0}^{n-1} i^p`（Faulhaber）。
 *
 * 闭式：`Σ_{i=0}^{n-1} i^p = (B_{p+1}(n) - B_{p+1}(0))/(p+1)`，
 * 其中 `B_{p+1}(x) = Σ_k C(p+1,k) B_k x^{p+1-k}` 是伯努利多项式。
 * `B_{p+1}(0) = B_{p+1}`，而 p ≥ 1 时 B_{p+1} 为 0（奇数下标），故多数情形可以省略。
 */
Ptr powerSum(const BigInt &p, const Ptr &upperInclusive) {
    if (p.isNegative()) return Ptr();
    const std::size_t order = static_cast<std::size_t>(p.toLongLong());
    if (order > 200) return Ptr();

    // Σ_{i=1}^{n} 1 = n（i^0 的特殊情形，伯努利公式在这里会给出 B_1 相关的错误结果）
    if (order == 0) return simplify(upperInclusive);

    // Σ_{i=1}^{upper} i^p = Σ_{i=0}^{upper} i^p（0^p 对 p>0 为 0）
    // = (B_{p+1}(N+1) - B_{p+1}) / (p+1)，N = upper
    Ptr shifted = makeAdd(upperInclusive, integer(1));
    Memory::Vector<NodePtr> terms;
    // 持有者：terms 里存的是裸指针，生成 term 的 Ptr 一旦离开作用域就会释放节点，
    // 所以必须由 holders 把这些中间节点保活到 makeAddMany 之后。
    Memory::Vector<Ptr> holders;

    for (std::size_t k = 0; k <= order + 1; ++k) {
        const Ptr &bk = bernoulliAt(k);
        if (!bk) return Ptr();
        if (bk->isNumber() && bk->numerator().isZero()) continue; // 跳过 0 项

        // C(p+1, k) * B_k * N^(p+1-k)
        BigInt binom(1);
        for (std::size_t j = 1; j <= k; ++j) {
            binom = binom * BigInt(static_cast<long long>(order + 1 - k + j)) / BigInt(static_cast<long long>(j));
        }
        const std::size_t exponent = order + 1 - k;
        Ptr term = makeMul(bk, makeNumber(binom, BigInt(1)));
        if (exponent > 0) term = makeMul(term, makePow(shifted, integer(static_cast<long long>(exponent))));
        holders.push_back(term);
        terms.push_back(term.get());
    }
    Ptr numerator = makeAddMany(terms);
    // 减去常数项 B_{p+1}
    numerator = makeAdd(numerator, makeNegate(bernoulliAt(order + 1)));
    return makeDivide(numerator, integer(static_cast<long long>(order + 1)));
}

/* ------------------------------------------------ 几何级数辅助 */

/**
 * @brief 把二次多项式 `c2·i² + c1·i + c0` 分解成 `(i + a)(i + b)` 的整数偏移。
 *
 * 只处理首项系数为 1、判别式为完全平方的情形（教科书里的裂项题都在这个范围）。
 * @param shifts 输出 `a` 与 `b`，并保证 `shifts[0] < shifts[1]`。
 * @return 成功返回 true。
 */
bool integerLinearFactors(const UnivariatePoly &polynomial, long long *shifts) {
    if (polynomial.degree() != 2) return false;
    long long constant = 0, linear = 0, quadratic = 0;
    if (!polynomial.coefficients[0]->isSmallInteger(&constant)) return false;
    if (!polynomial.coefficients[1]->isSmallInteger(&linear)) return false;
    if (!polynomial.coefficients[2]->isSmallInteger(&quadratic)) return false;
    if (quadratic != 1) return false;

    // 展开 `(i + a)(i + b) = i² + (a+b)·i + a·b`：
    // 一次项系数就是 a+b（不是 -(a+b)），常数项是 a·b。
    const long long discriminant = linear * linear - 4 * constant;
    if (discriminant < 0) return false;
    const long long root = static_cast<long long>(std::llround(std::sqrt(static_cast<double>(discriminant))));
    if (root * root != discriminant) return false;
    if ((linear + root) % 2 != 0) return false;

    const long long first = (linear + root) / 2;
    const long long second = (linear - root) / 2;
    shifts[0] = std::min(first, second);
    shifts[1] = std::max(first, second);
    return shifts[0] != shifts[1];
}

/** @brief 有限上限时 `Σ_{i=a}^{b} 1/(i + shift)`，用谐波数表示。 */
Ptr sumShiftedReciprocal(long long shift, const SumLimits &limits, std::uint32_t variable) {
    (void)variable;
    // Σ 1/(i + s) = H(b + s) - H(a + s - 1)
    Ptr upperShifted = makeAdd(limits.upper, integer(shift));
    Ptr lowerShifted = makeAdd(limits.lower, integer(shift - 1));
    Memory::Vector<NodePtr> upperArgs{upperShifted.get()};
    Memory::Vector<NodePtr> lowerArgs{lowerShifted.get()};
    Ptr harmonicUpper(Node::makeFunction("harmonic", upperArgs));
    Ptr harmonicLower(Node::makeFunction("harmonic", lowerArgs));
    if (!harmonicUpper || !harmonicLower) return Ptr();
    return makeAdd(harmonicUpper, makeNegate(harmonicLower));
}

/**
 * @brief 无穷上限的裂项求和：`Σ_{i=a}^{∞} [1/(i+s0) - 1/(i+s1)]`。
 *
 * 两个调和尾项各自发散，但相减后恰好只剩前 `gap = s1 - s0` 项：
 * `Σ_{k=0}^{gap-1} 1/(a + s0 + k)`。这里逐项展开成有理数。
 * @param lower 求和下限（须为小整数，才能得到数值答案）。
 * @return 成功返回 true 并写入 `value`。
 */
bool sumShiftedReciprocalInfinite(long long shiftLow, long long shiftHigh, const SumLimits &limits, Ptr *value) {
    long long start = 0;
    if (!limits.lower->isSmallInteger(&start)) return false;
    const long long gap = shiftHigh - shiftLow;
    if (gap <= 0 || gap > 4096) return false;
    // Σ_{i=a}^{∞} [1/(i+s0) - 1/(i+s1)]：只剩 i = a .. a+gap-1 这些首项，写成
    // Σ_{k=0}^{gap-1} 1/(a + s0 + k)。
    // 注意下标必须是 `i + s0`（s0 为较小偏移），首项是 1/(a+s0)。
    Memory::Vector<NodePtr> pieces;
    Memory::Vector<Ptr> holders;
    for (long long k = 0; k < gap; ++k) {
        Ptr piece = makeDivide(integer(1), integer(start + shiftLow + k));
        holders.push_back(piece);
        pieces.push_back(piece.get());
    }
    Ptr sum = makeAddMany(pieces);
    *value = simplify(sum);
    return true;
}


} // namespace

namespace {

/**
 * @brief 把被加项整理成"单个分数 1/分母"的形式，输出分母多项式。
 *
 * 化简会把 `1/(i*(i+1))` 拆成 `1/i * 1/(i+1)`，只认 `Pow(...,-1)` 的裂项规则
 * 就看不见整个分母。这里统一两种外形：
 * - `Pow(p, -1)`      → 分母 p
 * - `Mul[Pow(f_k,-1)]` → 分母 ∏ f_k
 *
 * @return 成功返回 true，`*denominator` 为分母多项式。
 */
bool extractReciprocalDenominator(const Ptr &term, std::uint32_t variable, UnivariatePoly *denominator) {
    if (!term) return false;

    // 形如 Pow(p, -1)
    if (term->kind() == NodeKind::Pow && term->operand(1) && term->operand(1)->isNegativeOne()) {
        return toPolynomial(Ptr::borrow(term->operand(0)).get(), variable, denominator);
    }

    // 形如 Mul[f_1^(-1) · f_2^(-1) · ...]
    if (term->kind() == NodeKind::Mul && term->size() >= 2) {
        Memory::Vector<NodePtr> bases;
        Memory::Vector<Ptr> holders;
        for (NodePtr child : term->operands()) {
            Ptr factor = Ptr::borrow(child);
            if (factor->kind() != NodeKind::Pow || !factor->operand(1) ||
                !factor->operand(1)->isNegativeOne()) {
                return false;
            }
            Ptr base = Ptr::borrow(factor->operand(0));
            holders.push_back(base);
            bases.push_back(base.get());
        }
        if (bases.empty()) return false;
        return toPolynomial(makeMulMany(bases).get(), variable, denominator);
    }
    return false;
}

} // namespace

/* -------------------------------------------------------------- 求和 */

namespace {

/**
 * @brief 求和主实现。
 */
Ptr summationImpl(const Ptr &term, const SumLimits &limits) {
    if (!term || !limits.lower || !limits.upper) return Ptr();
    const std::uint32_t variable = limits.variable;


    // 1. 被加项与变量无关
    if (!hasVariable(term, variable)) {
        Ptr count = makeAdd(makeAdd(limits.upper, makeNegate(limits.lower)), integer(1));
        return makeMul(term, count);
    }

    // 上下限相等
    if (limits.lower->equals(limits.upper.get())) return substitute(term, variable, limits.lower);

    bool upperInfinite = false;
    bool lowerInfinite = false;
    const bool infiniteUpper = isInfinity(limits.upper, &upperInfinite) && upperInfinite;
    isInfinity(limits.lower, &lowerInfinite);

    long long lowerValue = 0;
    const bool lowerIsInteger = limits.lower->isSmallInteger(&lowerValue);
    long long upperValue = 0;
    const bool upperIsInteger = !infiniteUpper && limits.upper->isSmallInteger(&upperValue);

    // 2. 上限为具体整数、项数不多 → 直接展开（最可靠，也能兜住各种怪式子）
    //    注意：terms 存的是裸指针，必须由 holders 保活到 makeAddMany 返回之后，
    //    否则每轮循环结束 value 析构、节点被归还池，terms 里全是悬垂指针。
    if (upperIsInteger && upperValue - lowerValue < 64 && upperValue >= lowerValue) {
        Memory::Vector<NodePtr> terms;
        Memory::Vector<Ptr> holders;
        for (long long index = lowerValue; index <= upperValue; ++index) {
            Ptr value = substitute(term, variable, integer(index));
            holders.push_back(value);
            terms.push_back(value.get());
        }
        return makeAddMany(terms);
    }

    // 3. 多项式通项：能写成 i 的多项式就逐次幂求和（Faulhaber）。
    //    必须排在 Add/Mul 拆分之前——先拆会把 `3i^2-2i+1` 变成三项分别求和，
    //    结果再被求和函数打包成 `3*(...)`，失去了展开成标准多项式的机会。
    if (!infiniteUpper) {
        UnivariatePoly polynomial;
        if (toPolynomial(term.get(), variable, &polynomial) && !polynomial.isZero()) {
            Memory::Vector<Ptr> pieces;
            bool allOk = true;
            for (std::size_t k = 0; k < polynomial.coefficients.size(); ++k) {
                const Ptr &coefficient = polynomial.coefficients[k];
                if (!coefficient || (coefficient->isNumber() && coefficient->numerator().isZero())) continue;
                Ptr partial = powerSum(BigInt(static_cast<long long>(k)), limits.upper);
                if (!partial) {
                    allOk = false;
                    break;
                }
                pieces.push_back(makeMul(coefficient, partial));
            }
            if (allOk && !pieces.empty()) {
                Memory::Vector<NodePtr> raw;
                raw.reserve(pieces.size());
                for (const Ptr &piece : pieces) raw.push_back(piece.get());
                // 只展开、不再 simplify：simplify 会把展开后的和重新组合成
                // `系数*(...)` 形式，等于把展开白做了。
                return expand(makeAddMany(raw), 20000);
            }
        }
    }

    // 4. 线性拆分：Add 逐项求和
    if (term->kind() == NodeKind::Add) {
        Memory::Vector<NodePtr> parts;
        bool allOk = true;
        Memory::Vector<Ptr> holders;
        for (NodePtr child : term->operands()) {
            Ptr part = summationImpl(Ptr::borrow(child), limits);
            if (!part) {
                allOk = false;
                break;
            }
            holders.push_back(part);
            parts.push_back(part.get());
        }
        if (allOk) return makeAddMany(parts);
    }

    // 4.5 算术-几何级数：Σ (i 的多项式)·q^i。
    //     必须排在 Mul 拆分之前——先拆会把多项式因子单独拿走，几何部分算完之后
    //     变成 "Σ i · (几何级数的和)"，值就不对了。
    //     Σ_{i=a}^∞ q^i = q^a/(1-q)，Σ_{i=a}^∞ i·q^i = q^a·[a/(1-q) + q/(1-q)²]
    //     这里只处理一次多项式（覆盖 i、i+1、2i-3 这类常见题）。
    if (infiniteUpper && term->kind() == NodeKind::Mul) {
        Ptr polygon = Ptr();
        Ptr geoPower = Ptr();
        Ptr ratio = Ptr();          // q
        bool matched = false;
        Memory::Vector<NodePtr> others;
        Memory::Vector<Ptr> holders;
        for (NodePtr child : term->operands()) {
            Ptr factor = Ptr::borrow(child);
            // 形如 q^i 或 q^(-i)：指数必须是 i 的线性函数且系数为 ±1。
            // 指数可能写成 `i`、`-i`、`(-1)*i` 等多种外形，这里统一用多项式判定。
            if (!geoPower && factor->kind() == NodeKind::Pow && factor->operand(1)) {
                Ptr base = Ptr::borrow(factor->operand(0));
                Ptr exponent = Ptr::borrow(factor->operand(1));
                UnivariatePoly exponentPoly;
                if (toPolynomial(exponent.get(), variable, &exponentPoly) && exponentPoly.degree() == 1) {
                    long long slope = 0;
                    if (exponentPoly.coefficients[1]->isSmallInteger(&slope) &&
                        (slope == 1 || slope == -1)) {
                        ratio = slope == 1 ? base : makeDivide(integer(1), base);
                        holders.push_back(ratio);
                        geoPower = factor;
                        continue;
                    }
                }
            }
            if (!polygon && hasVariable(factor, variable)) {
                UnivariatePoly candidate;
                if (toPolynomial(factor.get(), variable, &candidate) && candidate.degree() <= 1) {
                    polygon = factor;
                    continue;
                }
            }
            others.push_back(factor.get());
        }
        matched = polygon && geoPower && ratio && ratio->isNumber();
        if (matched) {
            const double magnitude = (ratio->numerator().toDouble()) / (ratio->denominator().toDouble());
            matched = magnitude > -1.0 && magnitude < 1.0;
        }
        if (matched) {
            UnivariatePoly polynomial;
            if (toPolynomial(polygon.get(), variable, &polynomial)) {
                const Ptr a = limits.lower;
                Ptr oneMinusQ = makeAdd(integer(1), makeNegate(ratio));
                Ptr oneMinusQSquared = makePow(oneMinusQ, integer(2));
                Ptr qPowerA = makePow(ratio, a);
                // Σ_{i=a}^∞ q^i = q^a/(1-q)
                Ptr geo = makeDivide(qPowerA, oneMinusQ);
                Ptr second = makeDivide(makeMul(qPowerA, ratio), oneMinusQSquared);
                Ptr weighted = makeAdd(makeDivide(makeMul(a, qPowerA), oneMinusQ), second);
                // 多项式在 i = a 处的值，配合线性项 i 分解：c1·i + c0
                Ptr linearCoefficient = polynomial.degree() >= 1 ? polynomial.coefficients[1] : integer(0);
                Ptr constantTerm = polynomial.coefficients[0];
                Ptr total = makeAdd(makeMul(linearCoefficient, weighted), makeMul(constantTerm, geo));
                // 剩余与变量无关的因子
                if (!others.empty()) {
                    Ptr rest = makeMulMany(others);
                    total = makeMul(rest, total);
                }
                return simplify(total);
            }
        }
    }

    // Mul：把与变量无关的因子提到求和号外面
    if (term->kind() == NodeKind::Mul) {
        Memory::Vector<Ptr> constantFactors;
        Memory::Vector<NodePtr> varyingFactors;
        Memory::Vector<Ptr> holders;
        for (NodePtr child : term->operands()) {
            Ptr factor = Ptr::borrow(child);
            if (hasVariable(factor, variable)) {
                varyingFactors.push_back(child);
            } else {
                constantFactors.push_back(factor);
            }
        }
        if (!constantFactors.empty() && !varyingFactors.empty()) {
            Memory::Vector<NodePtr> holder;
            holder.reserve(constantFactors.size());
            for (const Ptr &factor : constantFactors) holder.push_back(factor.get());
            Ptr inner = summationImpl(makeMulMany(varyingFactors), limits);
            if (inner) return makeMul(makeMulMany(holder), inner);
        }
    }


    // 5. 幂项：i^p 与 i^(-p)（多项式分支未命中的情形，如负指数）
    if (term->kind() == NodeKind::Pow) {
        Ptr base = Ptr::borrow(term->operand(0));
        Ptr exponent = Ptr::borrow(term->operand(1));
        if (base->kind() == NodeKind::Symbol && base->symbolId() == variable && exponent->isNumber()) {
            const BigInt &numerator = exponent->numerator();
            if (exponent->denominator().isOne()) {
                if (!numerator.isNegative()) {
                    // Σ i^p，p ≥ 0
                    if (infiniteUpper && !lowerIsInteger) return Ptr();
                    if (infiniteUpper) return Ptr(); // Σ i^p 在 p ≥ 0 时发散
                    return powerSum(numerator, limits.upper);
                }
                // Σ 1/i^r
                const BigInt order = -numerator;
                if (infiniteUpper) {
                    if (lowerValue < 1) return Ptr();
                    // ζ(r)：只有 r 为偶数或 r=1（发散）时有闭式
                    if (order.isOne()) return Ptr(); // 调和级数发散
                    if (!order.fitsLongLong()) return Ptr();
                    const long long r = order.toLongLong();
                    if (r > 200 || r % 2 != 0) return Ptr();
                    // ζ(2n) = (-1)^(n+1) B_{2n} (2π)^(2n) / (2 (2n)!)
                    const long long n = r / 2;
                    Ptr bk = bernoulliAt(static_cast<std::size_t>(r));
                    Ptr pi(Node::makeConstant(ConstantId::Pi));
                    Ptr piPower = makePow(pi, integer(r));
                    BigInt factorial(1);
                    for (long long k = 2; k <= r; ++k) factorial = factorial * BigInt(k);
                    Ptr sign = integer((n % 2 == 1) ? 1 : -1);
                    Ptr result = makeMul(sign, bk);
                    result = makeMul(result, piPower);
                    result = makeMul(result, makeNumber(BigInt::power(BigInt(2), static_cast<unsigned long long>(r)),
                                                        BigInt(1)));
                    result = makeDivide(result, makeNumber(BigInt(2) * factorial, BigInt(1)));
                    if (lowerValue > 1) {
                        // 从 2 开始的 p 级数：扣掉前 lower-1 项
                        Memory::Vector<Ptr> head;
                        for (long long index = 1; index < lowerValue; ++index) {
                            head.push_back(makeDivide(integer(1), makePow(integer(index), integer(r))));
                        }
                        Memory::Vector<NodePtr> headRaw;
                        headRaw.reserve(head.size());
                        for (const Ptr &item : head) headRaw.push_back(item.get());
                        result = makeAdd(result, makeNegate(makeAddMany(headRaw)));
                    }
                    return simplify(result);
                }
                // 有限和：Σ 1/i^r = H(上限, r) - H(下限-1, r)
                Ptr harmonicUpper = Ptr();
                Memory::Vector<NodePtr> args{limits.upper.get(), integer(static_cast<long long>(order.toLongLong())).get()};
                harmonicUpper = Ptr(Node::makeFunction("harmonic", args));
                Ptr lowerShifted = makeAdd(limits.lower, integer(-1));
                Memory::Vector<NodePtr> args2{lowerShifted.get(), integer(static_cast<long long>(order.toLongLong())).get()};
                Ptr harmonicLower(Node::makeFunction("harmonic", args2));
                return makeAdd(harmonicUpper, makeNegate(harmonicLower));
            }
        }

        // 几何级数：q^i 或 c^(a*i+b)
        if (hasVariable(exponent, variable)) {
            // 尝试把指数写成 i 的线性函数：目前只处理恰好等于 i 或 -i 的情形
            const bool exponentIsVariable =
                exponent->kind() == NodeKind::Symbol && exponent->symbolId() == variable;
            const bool exponentIsNegated = [&]() {
                if (exponent->kind() != NodeKind::Mul || exponent->size() != 2) return false;
                Ptr first = Ptr::borrow(exponent->operand(0));
                Ptr second = Ptr::borrow(exponent->operand(1));
                Ptr sym = first->isNumber() ? second : first;
                Ptr num = first->isNumber() ? first : second;
                return sym->kind() == NodeKind::Symbol && sym->symbolId() == variable && num->isNumber() &&
                       num->isNegativeOne();
            }();
            if (exponentIsVariable || exponentIsNegated) {
                // q = base 或 1/base
                Ptr ratio = exponentIsVariable ? base : makeDivide(integer(1), base);
                if (infiniteUpper) {
                    // Σ_{i=a}^∞ q^i = q^a/(1-q)，要求 |q|<1；数值 q 直接判断
                    if (!ratio->isNumber()) return Ptr();
                    const double magnitude = (ratio->numerator().toDouble()) / (ratio->denominator().toDouble());
                    if (!(magnitude > -1.0 && magnitude < 1.0)) return Ptr();
                    Ptr head = makePow(ratio, limits.lower);
                    return makeDivide(head, makeAdd(integer(1), makeNegate(ratio)));
                }
                // 有限：q^a (1 - q^(b-a+1))/(1-q)；q=1 时退化为项数
                Ptr count = makeAdd(makeAdd(limits.upper, makeNegate(limits.lower)), integer(1));
                if (ratio->isNumber() && ratio->isOne()) return count;
                Ptr head = makePow(ratio, limits.lower);
                Ptr remaining = makeAdd(integer(1), makeNegate(makePow(ratio, count)));
                return simplify(makeDivide(makeMul(head, remaining), makeAdd(integer(1), makeNegate(ratio))));
            }
        }
    }

    // 6. 裂项（telescoping）：1/((i+a)(i+b)) 这类分母可拆成两个相差常数的线性因式。
    //    1/((i+a)(i+b)) = 1/(b-a) · (1/(i+a) - 1/(i+b))，求和只剩首尾少数几项。
    //    只处理分母为二次、系数全为小整数的情形——够覆盖教科书级的裂项题。
    //    注意分母可能被化简拆成 `1/i * 1/(i+1)`，所以用 extractReciprocalDenominator
    //    统一还原，而不是只认单个 Pow 节点。
    {
        UnivariatePoly denominator;
        if (extractReciprocalDenominator(term, variable, &denominator) && denominator.degree() == 2) {
            long long shifts[2] = {0, 0};
            if (integerLinearFactors(denominator, shifts)) {
                // gap = |a - b|，系数 1/gap
                const long long gap = shifts[1] - shifts[0];
                if (gap != 0) {
                    Ptr scale = makeNumber(BigInt(1), BigInt(gap));
                    bool infiniteUpper = false;
                    if (isInfinity(limits.upper, &infiniteUpper) && infiniteUpper) {
                        Ptr tailSum;
                        if (sumShiftedReciprocalInfinite(shifts[0], shifts[1], limits, &tailSum)) {
                            return simplify(makeMul(scale, tailSum));
                        }
                    } else {
                        Ptr head = sumShiftedReciprocal(shifts[0], limits, variable);
                        Ptr tail = sumShiftedReciprocal(shifts[1], limits, variable);
                        if (head && tail) return simplify(makeMul(scale, makeAdd(head, makeNegate(tail))));
                    }
                }
            }
        }
    }

    // 7. 直接展开：上限为整数时逐项算（项数受控）
    if (upperIsInteger && upperValue >= lowerValue) {
        const long long span = upperValue - lowerValue;
        if (span <= 4096) {
            Memory::Vector<NodePtr> terms;
            Memory::Vector<Ptr> holders;
            for (long long index = lowerValue; index <= upperValue; ++index) {
                Ptr value = simplify(substitute(term, variable, integer(index)));
                holders.push_back(value);
                terms.push_back(value.get());
            }
            return makeAddMany(terms);
        }
    }

    return Ptr();
}

} // namespace

Ptr summation(const Ptr &term, const SumLimits &limits) { return summationImpl(term, limits); }

/* ---------------------------------------------------------- 敛散性 */

bool seriesConvergence(const Ptr &term, const SumLimits &limits, Convergence *verdict, Ptr *value) {
    if (!term || !verdict) return false;
    *verdict = Convergence::Unknown;
    bool positive = true;
    if (!isInfinity(limits.upper, &positive) || !positive) return false;

    // 先尝试闭式求和：能求出有限值就是收敛，能得到 Infinity 就是发散
    Ptr closed = summationImpl(term, limits);
    if (closed) {
        if (closed->kind() == NodeKind::Constant &&
            (closed->constant() == ConstantId::Infinity || closed->constant() == ConstantId::NegativeInfinity)) {
            *verdict = Convergence::Divergent;
            if (value) *value = Ptr();
            return true;
        }
        if (!isConstant(closed.get())) return false;
        *verdict = Convergence::Convergent;
        if (value) *value = closed;
        return true;
    }

    // p 级数：1/i^p
    if (term->kind() == NodeKind::Pow) {
        Ptr base = Ptr::borrow(term->operand(0));
        Ptr exponent = Ptr::borrow(term->operand(1));
        if (base->kind() == NodeKind::Symbol && base->symbolId() == limits.variable && exponent->isNumber() &&
            exponent->isNegative()) {
            const BigInt order = -exponent->numerator();
            if (!order.fitsLongLong()) return false;
            *verdict = (order.toLongLong() > 1) ? Convergence::Convergent : Convergence::Divergent;
            return true;
        }
    }

    // 几何级数：q^i，|q|<1 收敛
    if (term->kind() == NodeKind::Pow) {
        Ptr exponent = Ptr::borrow(term->operand(1));
        if (exponent->kind() == NodeKind::Symbol && exponent->symbolId() == limits.variable) {
            Ptr base = Ptr::borrow(term->operand(0));
            if (base->isNumber()) {
                const double magnitude =
                    base->numerator().toDouble() / base->denominator().toDouble();
                if (magnitude > -1.0 && magnitude < 1.0) {
                    *verdict = Convergence::Convergent;
                    if (value) *value = summationImpl(term, limits);
                } else {
                    *verdict = Convergence::Divergent;
                }
                return true;
            }
        }
    }
    return false;
}

/* -------------------------------------------------------- 级数展开 */

Ptr seriesExpand(const Ptr &expression, std::uint32_t variable, const Ptr &point, int order) {
    if (!expression || !point || order < 0) return Ptr();
    return taylorExpand(expression, variable, point, order);
}

/* -------------------------------------------------------- 序列识别 */

Ptr guessSequence(const Memory::Vector<Ptr> &terms, std::uint32_t index) {
    if (terms.size() < 2) return Ptr();

    // 全部项必须是数值，否则无法比较相邻项
    for (const Ptr &term : terms) {
        if (!term || !term->isNumber()) return Ptr();
    }
    Ptr indexSymbol(Node::makeSymbolById(index));

    // 等比：相邻项比值恒定（先判，因为 0 项会让差分法误判）
    bool geometric = true;
    for (std::size_t i = 1; i < terms.size(); ++i) {
        if (terms[i - 1]->isZero()) {
            geometric = false;
            break;
        }
        Ptr ratio = makeDivide(terms[i], terms[i - 1]);
        if (!ratio || !ratio->isNumber()) {
            geometric = false;
            break;
        }
        if (i == 1) continue; // 记下第一个比值，后面逐个比较
        Ptr previous = makeDivide(terms[i - 1], terms[i - 2]);
        if (!previous || !previous->isNumber() || !previous->equals(ratio.get())) {
            geometric = false;
            break;
        }
    }
    if (geometric) {
        Ptr ratio = terms[1]->isZero() ? integer(0) : makeDivide(terms[1], terms[0]);
        if (ratio && ratio->isNumber()) {
            if (ratio->isOne()) return terms.front();      // 常数序列
            if (ratio->isZero()) return integer(0);        // 0,0,0,…
            // a_i = a_0 · q^i
            return simplify(makeMul(terms.front(), makePow(ratio, indexSymbol)));
        }
    }

    // 等差：一阶差分恒定
    Memory::Vector<Ptr> firstDifference;
    for (std::size_t i = 1; i < terms.size(); ++i) {
        Ptr difference = simplify(makeAdd(terms[i].share(), makeNegate(terms[i - 1].share())));
        if (!difference || !difference->isNumber()) return Ptr();
        firstDifference.push_back(difference);
    }
    bool constantDifference = true;
    for (std::size_t i = 1; i < firstDifference.size(); ++i) {
        if (!firstDifference[i]->equals(firstDifference[i - 1].get())) constantDifference = false;
    }
    if (constantDifference) {
        // a_i = a_0 + d·i（i 从 0 计）
        Ptr difference = firstDifference.empty() ? integer(0) : firstDifference.front();
        if (difference->isZero()) return terms.front();
        return simplify(makeAdd(terms.front().share(), makeMul(difference, indexSymbol)));
    }

    // 有限差分表：若 n 阶差分恒定，通项为 n 次多项式（Newton 前向插值）
    // 这里只做低阶（≤ 4 阶），避免误判与数值噪声。
    Memory::Vector<Ptr> table = firstDifference;
    Ptr binomialWeight = integer(1);
    Ptr result = terms.front().share();
    for (int order = 1; order <= 4; ++order) {
        if (table.empty()) break;
        // 牛顿插值项：Δ^order · C(i, order)
        Ptr coefficient = table.front();
        // C(i, order) = i(i-1)…(i-order+1)/order!
        Ptr falling = integer(1);
        BigInt factorial(1);
        for (int k = 0; k < order; ++k) {
            falling = makeMul(falling, makeAdd(indexSymbol, integer(-k)));
            factorial = factorial * BigInt(k + 1);
        }
        Ptr weight = makeDivide(falling, makeNumber(factorial, BigInt(1)));
        result = makeAdd(result, makeMul(coefficient, weight));

        // 计算下一阶差分
        Memory::Vector<Ptr> next;
        for (std::size_t i = 1; i < table.size(); ++i) {
            next.push_back(simplify(makeAdd(table[i].share(), makeNegate(table[i - 1].share()))));
        }
        bool constant = true;
        for (std::size_t i = 1; i < next.size(); ++i) {
            if (!next[i]->equals(next[i - 1].get())) constant = false;
        }
        if (constant && !next.empty() && next.front()->isZero()) {
            // 牛顿插值天然产出 `n(n-1)/2` 这类降阶乘形式，展开成标准多项式更易读。
            return simplify(expand(result, 20000));
        }
        table = next;
    }
    return Ptr();
}

} // namespace Cas::Detail

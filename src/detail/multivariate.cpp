// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#include "cas/detail/multivariate.h"
#include "cas/detail/polynomial.h"

#include <algorithm>

namespace Cas::Detail {
namespace {

/** @brief 快捷构造整数。 */
Ptr integer(long long value) { return Ptr(Node::makeInteger(BigInt(value))); }

/* --------------------------------------------------- 单项式序与运算 */

/** @brief 单项式的总次数。 */
long long totalDegree(const Monomial &monomial) {
    long long total = 0;
    for (long long exponent : monomial.exponents) total += exponent;
    return total;
}

/**
 * @brief 单项式比较，返回 >0 表示 left 在序中更大。
 *
 * Lex：从左到右比较指数（最左边的变量最"重要"）；
 * Grevlex：先比总次数，再比"最右非零差"，差的符号取反（1/x 顺序）。
 */
int compareMonomial(const Monomial &left, const Monomial &right, MonomialOrder order) {
    const std::size_t size = std::min(left.exponents.size(), right.exponents.size());
    if (order == MonomialOrder::Lex) {
        for (std::size_t i = 0; i < size; ++i) {
            if (left.exponents[i] != right.exponents[i]) return left.exponents[i] < right.exponents[i] ? -1 : 1;
        }
        return 0;
    }

    const long long leftDegree = totalDegree(left), rightDegree = totalDegree(right);
    if (leftDegree != rightDegree) return leftDegree < rightDegree ? -1 : 1;
    // Grevlex 的平局规则：从最右（最次要的变量）往左找第一个不同的指数，
    // **指数小的那个反而更大**。也就是比较"逆字典序"。
    // 之前写成 `left > right ? -1 : 1`（指数大的更大），方向正好相反，
    // 于是 x^2*y 会排在 x*y^2 后面，`terms.front()` 拿到的不是真正的首项，
    // 约化就会约错项，Groebner 基永远算不对。
    for (std::size_t i = size; i-- > 0;) {
        if (left.exponents[i] != right.exponents[i]) return left.exponents[i] < right.exponents[i] ? 1 : -1;
    }
    return 0;
}

/** @brief 单项式乘法：指数相加、系数相乘。 */
Monomial multiplyMonomial(const Monomial &left, const Monomial &right) {
    Monomial result;
    result.exponents.reserve(left.exponents.size());
    for (std::size_t i = 0; i < left.exponents.size(); ++i) result.exponents.push_back(left.exponents[i] + right.exponents[i]);
    result.numerator = left.numerator * right.numerator;
    result.denominator = left.denominator * right.denominator;
    return result;
}

/** @brief 单项式除法：返回 false 表示不能整除（有负指数）。 */
bool divideMonomial(const Monomial &left, const Monomial &right, Monomial *out) {
    out->exponents.clear();
    out->exponents.reserve(left.exponents.size());
    for (std::size_t i = 0; i < left.exponents.size(); ++i) {
        const long long exponent = left.exponents[i] - right.exponents[i];
        if (exponent < 0) return false;
        out->exponents.push_back(exponent);
    }
    out->numerator = left.numerator * right.denominator;
    out->denominator = left.denominator * right.numerator;
    return true;
}

/** @brief 单项式是否整除（不含系数的整数性）。 */
bool dividesMonomial(const Monomial &divisor, const Monomial &dividend) {
    for (std::size_t i = 0; i < divisor.exponents.size(); ++i) {
        if (divisor.exponents[i] > dividend.exponents[i]) return false;
    }
    return true;
}

/* ----------------------------------------------------- 多项式基本运算 */

/** @brief 插入单项式并合并同类项，再把有序保持住。 */
void addTerm(MultivariatePoly *polynomial, const Monomial &term, MonomialOrder order) {
    for (Monomial &existing : polynomial->terms) {
        bool same = existing.exponents.size() == term.exponents.size();
        for (std::size_t i = 0; same && i < existing.exponents.size(); ++i) {
            if (existing.exponents[i] != term.exponents[i]) same = false;
        }
        if (!same) continue;
        // 系数相加
        const BigInt numerator = existing.numerator * term.denominator + term.numerator * existing.denominator;
        const BigInt denominator = existing.denominator * term.denominator;
        const BigInt gcd = BigInt::gcd(numerator.abs(), denominator.abs());
        if (!gcd.isZero() && !gcd.isOne()) {
            existing.numerator = numerator / gcd;
            existing.denominator = denominator / gcd;
        } else {
            existing.numerator = numerator;
            existing.denominator = denominator;
        }
        return;
    }
    polynomial->terms.push_back(term);
    std::sort(polynomial->terms.begin(), polynomial->terms.end(),
              [order](const Monomial &a, const Monomial &b) { return compareMonomial(a, b, order) > 0; });
}

/** @brief 去掉系数为零的项。 */
void normalizePoly(MultivariatePoly *polynomial, MonomialOrder order) {
    Memory::Vector<Monomial> kept;
    kept.reserve(polynomial->terms.size());
    for (const Monomial &term : polynomial->terms) {
        if (!term.numerator.isZero()) kept.push_back(term);
    }
    polynomial->terms = kept;
    std::sort(polynomial->terms.begin(), polynomial->terms.end(),
              [order](const Monomial &a, const Monomial &b) { return compareMonomial(a, b, order) > 0; });
}

/** @brief 多项式加法。 */
MultivariatePoly polySum(const MultivariatePoly &left, const MultivariatePoly &right, MonomialOrder order) {
    MultivariatePoly result;
    result.variables = left.variables;
    for (const Monomial &term : left.terms) addTerm(&result, term, order);
    for (const Monomial &term : right.terms) addTerm(&result, term, order);
    normalizePoly(&result, order);
    return result;
}

/** @brief 多项式取负。 */
MultivariatePoly polyNegate(const MultivariatePoly &polynomial) {
    MultivariatePoly result = polynomial;
    for (Monomial &term : result.terms) term.numerator = term.numerator.negated();
    return result;
}

/** @brief 多项式乘法。当前供外部/后续算法使用，暂未在 Groebner 流程中调用。 */
[[maybe_unused]] MultivariatePoly polyProduct(const MultivariatePoly &left, const MultivariatePoly &right,
                                              MonomialOrder order) {
    MultivariatePoly result;
    result.variables = left.variables;
    for (const Monomial &a : left.terms) {
        for (const Monomial &b : right.terms) {
            addTerm(&result, multiplyMonomial(a, b), order);
        }
    }
    normalizePoly(&result, order);
    return result;
}

/** @brief 单项式乘多项式。 */
MultivariatePoly polyScale(const MultivariatePoly &polynomial, const Monomial &monomial, MonomialOrder order) {
    MultivariatePoly result;
    result.variables = polynomial.variables;
    for (const Monomial &term : polynomial.terms) addTerm(&result, multiplyMonomial(term, monomial), order);
    normalizePoly(&result, order);
    return result;
}

/** @brief 首一化：把首项系数变成 1（保持调用方的单项式序）。 */
void makeMonic(MultivariatePoly *polynomial, MonomialOrder order = MonomialOrder::Grevlex) {
    if (polynomial->terms.empty()) return;
    Monomial inverse = polynomial->terms.front();
    inverse.exponents.assign(inverse.exponents.size(), 0);
    inverse.numerator = polynomial->terms.front().denominator;
    inverse.denominator = polynomial->terms.front().numerator;
    if (inverse.denominator.isNegative()) {
        inverse.numerator = inverse.numerator.negated();
        inverse.denominator = inverse.denominator.negated();
    }
    // 必须沿用调用方的序：写死 Grevlex 会让 Lex 下的首项被重新排到后面，
    // 于是基元素的首项不再是真正的首项，polyReduce 与 sPolynomial 全部失效。
    *polynomial = polyScale(*polynomial, inverse, order);
}

/**
 * @brief 带余除法：f = Σ q_i · g_i + r（关于给定单项式序）。
 */
bool polyReduce(const MultivariatePoly &dividend, const Memory::Vector<MultivariatePoly> &divisors,
                const MonomialOrder order, Memory::Vector<MultivariatePoly> *quotients,
                MultivariatePoly *remainder) {
    quotients->clear();
    quotients->reserve(divisors.size());
    for (std::size_t i = 0; i < divisors.size(); ++i) {
        MultivariatePoly zero;
        zero.variables = dividend.variables;
        quotients->push_back(zero);
    }
    MultivariatePoly rest = dividend;
    remainder->variables = dividend.variables;
    remainder->terms.clear();

    std::size_t guard = 0;
    while (!rest.terms.empty()) {
        if (++guard > 200000) return false; // 防御：异常数据下不死循环
        const Monomial leading = rest.terms.front();
        bool reduced = false;
        for (std::size_t i = 0; i < divisors.size(); ++i) {
            if (divisors[i].terms.empty()) continue;
            const Monomial &divisorLeading = divisors[i].terms.front();
            if (!dividesMonomial(divisorLeading, leading)) continue;

            Monomial factor;
            if (!divideMonomial(leading, divisorLeading, &factor)) continue;
            // 系数也要除：factor = leading / divisorLeading
            addTerm(&quotients->at(i), factor, order);
            MultivariatePoly product = polyScale(divisors[i], factor, order);
            rest = polySum(rest, polyNegate(product), order);
            reduced = true;
            break;
        }
        if (reduced) continue;
        // 首项无法被约去，移入余式
        addTerm(remainder, leading, order);
        rest.terms.erase(rest.terms.begin());
    }
    return true;
}

/** @brief S-多项式。 */
MultivariatePoly sPolynomial(const MultivariatePoly &left, const MultivariatePoly &right, MonomialOrder order) {
    const Monomial &a = left.terms.front();
    const Monomial &b = right.terms.front();

    // 单项式的最小公倍：指数逐位取最大，系数取 1。
    // 系数不能留空（Monomial 默认系数是 0/0），否则下面 divideMonomial
    // 会除出 0/0，整个 S-多项式恒等于 0，Buchberger 直接退化成"原样返回"。
    Monomial lcm;
    lcm.exponents.reserve(a.exponents.size());
    for (std::size_t i = 0; i < a.exponents.size(); ++i) lcm.exponents.push_back(std::max(a.exponents[i], b.exponents[i]));
    lcm.numerator = BigInt(1);
    lcm.denominator = BigInt(1);

    Monomial scaleA, scaleB;
    if (!divideMonomial(lcm, a, &scaleA)) return MultivariatePoly();
    if (!divideMonomial(lcm, b, &scaleB)) return MultivariatePoly();
    // S = (lcm/lt(f))·f - (lcm/lt(g))·g。
    // divideMonomial 已经把系数一并除掉（scaleA = lcm / a = 1/lc(f)），
    // 两项的首项因此都化成 lcm 的系数，相加即相消。不能再额外取一次倒数。
    MultivariatePoly first = polyScale(left, scaleA, order);
    MultivariatePoly second = polyScale(right, scaleB, order);
    return polySum(first, polyNegate(second), order);
}

/**
 * @brief Buchberger 第一判据：lcm(lt(f), lt(g)) 恰好等于两者的乘积时，
 * 该对的 S-多项式必然约化为 0，可以跳过。
 *
 * 注意不要和"首项能否被基中某元素整除"混淆——后者是化简判据，
 * 在 S-多项式尚未约化时套用会把有信息量的对丢掉。
 */
bool coprimeLeads(const MultivariatePoly &left, const MultivariatePoly &right) {
    if (left.terms.empty() || right.terms.empty()) return true;
    const Monomial &a = left.terms.front();
    const Monomial &b = right.terms.front();
    const std::size_t size = std::min(a.exponents.size(), b.exponents.size());
    for (std::size_t i = 0; i < size; ++i) {
        if (a.exponents[i] != 0 && b.exponents[i] != 0) return false;
    }
    return true;
}

/* ------------------------------------------------------ 表达式 → 多元 */

/** @brief 由表达式构造单项式；失败返回 false。 */
bool buildMonomial(const NodePtr &node, const Memory::Vector<std::uint32_t> &variables, Monomial *out) {
    out->exponents.assign(variables.size(), 0);
    out->numerator = BigInt(1);
    out->denominator = BigInt(1);
    if (!node) return false;

    switch (node->kind()) {
    case NodeKind::Integer:
    case NodeKind::Rational:
        out->numerator = node->numerator();
        out->denominator = node->denominator();
        return true;
    case NodeKind::Symbol: {
        for (std::size_t i = 0; i < variables.size(); ++i) {
            if (variables[i] == node->symbolId()) {
                out->exponents[i] = 1;
                return true;
            }
        }
        return false; // 未声明的符号
    }
    case NodeKind::Pow: {
        NodePtr base = node->operand(0);
        NodePtr exponent = node->operand(1);
        long long power = 0;
        if (!exponent->isSmallInteger(&power) || power < 0) return false;
        if (base->kind() != NodeKind::Symbol) return false;
        for (std::size_t i = 0; i < variables.size(); ++i) {
            if (variables[i] == base->symbolId()) {
                out->exponents[i] = power;
                return true;
            }
        }
        return false;
    }
    case NodeKind::Mul: {
        for (NodePtr factor : node->operands()) {
            // 注意各因子的系数要**乘**上去，不能共用一个 out：
            // buildMonomial 会先把 out 清零，若对每个因子复用 out，
            // 前一个因子的系数就被后一个覆盖了（-12 → 1），
            // 于是 S-多项式算出来是常数 0，Groebner 基永远消不出元。
            Monomial piece;
            if (!buildMonomial(factor, variables, &piece)) return false;
            for (std::size_t i = 0; i < variables.size(); ++i) out->exponents[i] += piece.exponents[i];
            const BigInt numerator = out->numerator * piece.numerator;
            const BigInt denominator = out->denominator * piece.denominator;
            const BigInt gcd = BigInt::gcd(numerator.abs(), denominator.abs());
            if (!gcd.isZero() && !gcd.isOne()) {
                out->numerator = numerator / gcd;
                out->denominator = denominator / gcd;
            } else {
                out->numerator = numerator;
                out->denominator = denominator;
            }
        }
        return true;
    }
    default: return false;
    }
}

} // namespace

/* ---------------------------------------------------------------- 公开接口 */

bool toMultivariate(const NodePtr &expression, const Memory::Vector<std::uint32_t> &variables,
                    MultivariatePoly *out) {
    out->variables = variables;
    out->terms.clear();
    if (!expression) return false;

    if (expression->kind() == NodeKind::Add) {
        for (NodePtr term : expression->operands()) {
            Monomial monomial;
            if (!buildMonomial(term, variables, &monomial)) return false;
            addTerm(out, monomial, MonomialOrder::Grevlex);
        }
        normalizePoly(out, MonomialOrder::Grevlex);
        return true;
    }
    Monomial monomial;
    if (!buildMonomial(expression, variables, &monomial)) return false;
    addTerm(out, monomial, MonomialOrder::Grevlex);
    return true;
}

/** @brief 按指定序重排多项式（把首项放到最前面）。 */
MultivariatePoly resort(const MultivariatePoly &polynomial, MonomialOrder order) {
    MultivariatePoly result = polynomial;
    std::sort(result.terms.begin(), result.terms.end(),
              [order](const Monomial &a, const Monomial &b) { return compareMonomial(a, b, order) > 0; });
    return result;
}

Ptr fromMultivariate(const MultivariatePoly &polynomial) {
    Memory::Vector<NodePtr> terms;
    terms.reserve(polynomial.terms.size());
    Memory::Vector<Ptr> holders;
    holders.reserve(polynomial.terms.size());

    for (const Monomial &term : polynomial.terms) {
        Ptr coefficient = makeNumber(term.numerator, term.denominator);
        Memory::Vector<Ptr> factors;
        factors.push_back(coefficient);
        for (std::size_t i = 0; i < term.exponents.size(); ++i) {
            if (term.exponents[i] == 0) continue;
            Ptr variable = Ptr(Node::makeSymbolById(polynomial.variables[i]));
            if (term.exponents[i] == 1) {
                factors.push_back(variable);
            } else {
                factors.push_back(makePow(variable, integer(term.exponents[i])));
            }
        }
        Memory::Vector<NodePtr> raw;
        raw.reserve(factors.size());
        for (const Ptr &factor : factors) raw.push_back(factor.get());
        holders.push_back(factors.size() == 1 ? factors.front().share() : makeMulMany(raw));
        terms.push_back(holders.back().get());
    }
    return terms.empty() ? integer(0) : makeAddMany(terms);
}

Ptr expandMultivariate(const Ptr &expression, std::size_t termLimit) {
    if (!expression) return Ptr();
    Memory::Vector<std::uint32_t> variables = freeSymbols(expression.get());
    if (variables.size() <= 1) return expand(expression, termLimit);

    Ptr expanded = expand(expression, termLimit);
    MultivariatePoly polynomial;
    if (!toMultivariate(expanded.get(), variables, &polynomial)) return expanded;
    if (polynomial.terms.size() > termLimit) return expanded;

    // 按降序重新输出，让结果读起来是"字典序降幂"
    std::sort(polynomial.terms.begin(), polynomial.terms.end(), [](const Monomial &a, const Monomial &b) {
        return compareMonomial(a, b, MonomialOrder::Lex) > 0;
    });
    return fromMultivariate(polynomial);
}

Memory::Vector<Ptr> groebnerBasis(const Memory::Vector<Ptr> &generators,
                                  const Memory::Vector<std::uint32_t> &variables, const MonomialOrder order) {
    Memory::Vector<MultivariatePoly> basis;
    basis.reserve(generators.size());

    for (const Ptr &generator : generators) {
        if (!generator) continue;
        MultivariatePoly polynomial;
        if (!toMultivariate(generator.get(), variables, &polynomial)) continue;
        // toMultivariate 固定按 Grevlex 排；这里必须按本次请求的序重排，
        // 否则 `terms.front()` 取到的不是该序下的首项，
        // polyReduce / sPolynomial 就会约错项，基永远算不对。
        polynomial = resort(polynomial, order);
        normalizePoly(&polynomial, order);
        if (polynomial.isZero()) continue;
        makeMonic(&polynomial, order);
        // 去重必须比较**整个多项式**，不能只比首项指数：
        // `x^2+y^2-10` 与 `x^2-y^2-8` 的项数和首项完全相同，
        // 只比首项会把后者当重复丢掉，基里就少了独立的一条方程。
        bool duplicate = false;
        for (const MultivariatePoly &existing : basis) {
            if (existing.terms.size() != polynomial.terms.size()) continue;
            bool same = true;
            for (std::size_t i = 0; i < existing.terms.size() && same; ++i) {
                const Monomial &left = existing.terms[i];
                const Monomial &right = polynomial.terms[i];
                if (left.exponents != right.exponents || left.numerator != right.numerator ||
                    left.denominator != right.denominator) {
                    same = false;
                }
            }
            if (same) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) basis.push_back(polynomial);
    }

    if (basis.size() < 2) {
        Memory::Vector<Ptr> result;
        for (const MultivariatePoly &element : basis) result.push_back(fromMultivariate(element));
        return result;
    }

    // Buchberger：枚举对，用第一判据跳过，把非零余式加入基
    Memory::Vector<std::pair<std::size_t, std::size_t>> pairs;
    for (std::size_t i = 0; i < basis.size(); ++i) {
        for (std::size_t j = i + 1; j < basis.size(); ++j) pairs.push_back({i, j});
    }

    std::size_t processed = 0;
    while (processed < pairs.size()) {
        const auto pair = pairs[processed++];
        if (pair.first >= basis.size() || pair.second >= basis.size()) continue;
        // Buchberger 第一判据：首项互素的对，S-多项式必定约化为 0，跳过
        if (coprimeLeads(basis[pair.first], basis[pair.second])) continue;

        MultivariatePoly s = sPolynomial(basis[pair.first], basis[pair.second], order);
        normalizePoly(&s, order);
        if (s.isZero()) continue;
        // 注意：这里**不能**用"首项被基中某元素整除就跳过"。
        // 那是化简判据，会在 S-多项式尚未约化时把它丢掉——
        // 例如 [x^2+y^2-25, x*y-12] 的 S = y^3-25y+12x，
        // 首项 y^3 并不被 x^2 / x*y 整除，但真正的消元多项式
        // y^4-25y^2+144 恰恰要靠它继续约化出来。
        // Buchberger 第一判据是 lcm(lt(f),lt(g)) == lt(f)*lt(g)，
        // 与"能否整除"无关。

        Memory::Vector<MultivariatePoly> quotients;
        MultivariatePoly remainder;
        if (!polyReduce(s, basis, order, &quotients, &remainder)) continue;
        normalizePoly(&remainder, order);
        if (remainder.isZero()) continue;

        makeMonic(&remainder, order);
        const std::size_t newIndex = basis.size();
        basis.push_back(remainder);
        for (std::size_t i = 0; i < newIndex; ++i) pairs.push_back({i, newIndex});

        if (basis.size() > 64) break; // 防御：极端输入下限制规模
    }

    // 约化：先删掉首项能被其它元素首项整除的（这些是多余的），
    // 再把剩下每个元素的首项也约掉，得到唯一的约化 Groebner 基。
    //
    // 顺序很重要：不能直接对每个元素做 reduce 再判断余式是否为零。
    // 当两个生成元互为倍数（如 `x+y-1` 与 `2x+2y-2`）时，
    // 各自用对方约化都会得到零余式，照原逻辑两个都被丢掉，基就成了空集，
    // 整个方程组被误判成矛盾。
    Memory::Vector<bool> redundant(basis.size(), false);
    for (std::size_t i = 0; i < basis.size(); ++i) {
        if (basis[i].terms.empty()) {
            redundant[i] = true;
            continue;
        }
        for (std::size_t j = 0; j < basis.size(); ++j) {
            if (i == j || basis[j].terms.empty()) continue;
            if (!dividesMonomial(basis[j].terms.front(), basis[i].terms.front())) continue;
            // 首项互相整除（互为倍数）：只保留索引较小的那个，避免两个都被删掉
            const bool mutual = dividesMonomial(basis[i].terms.front(), basis[j].terms.front());
            if (mutual && j > i) continue;
            redundant[i] = true;
            break;
        }
    }

    Memory::Vector<MultivariatePoly> kept;
    for (std::size_t i = 0; i < basis.size(); ++i) {
        if (!redundant[i]) kept.push_back(basis[i]);
    }

    Memory::Vector<MultivariatePoly> reduced;
    for (std::size_t i = 0; i < kept.size(); ++i) {
        Memory::Vector<MultivariatePoly> others;
        for (std::size_t j = 0; j < kept.size(); ++j) {
            if (j != i) others.push_back(kept[j]);
        }
        if (others.empty()) {
            reduced.push_back(kept[i]);
            continue;
        }
        Memory::Vector<MultivariatePoly> quotients;
        MultivariatePoly remainder;
        if (!polyReduce(kept[i], others, order, &quotients, &remainder)) {
            reduced.push_back(kept[i]);
            continue;
        }
        normalizePoly(&remainder, order);
        if (remainder.isZero()) continue;
        makeMonic(&remainder, order);
        reduced.push_back(remainder);
    }

    Memory::Vector<Ptr> result;
    result.reserve(reduced.size());
    for (const MultivariatePoly &element : reduced) result.push_back(fromMultivariate(element));
    return result;
}

Ptr multivariateGcd(const Ptr &left, const Ptr &right, const Memory::Vector<std::uint32_t> &variables) {
    if (!left || !right) return integer(1);
    if (variables.size() == 1) return polyGcd(left.get(), right.get(), variables[0]);

    // 多元情形：用 Groebner 基判断公因子不现实，退化为"提取数值内容"
    MultivariatePoly a, b;
    if (!toMultivariate(left.get(), variables, &a) || !toMultivariate(right.get(), variables, &b)) {
        return integer(1);
    }
    BigInt gcd = BigInt(1);
    for (const Monomial &term : a.terms) gcd = BigInt::gcd(gcd, term.numerator.abs());
    for (const Monomial &term : b.terms) gcd = BigInt::gcd(gcd, term.numerator.abs());
    return integer(gcd.fitsLongLong() ? gcd.toLongLong() : 1);
}

} // namespace Cas::Detail

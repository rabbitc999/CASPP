#include "cas/detail/polynomial.h"

#include <algorithm>

namespace Cas::Detail {
namespace {

/** @brief 快捷构造整数。 */
Ptr integer(long long value) { return Ptr(Node::makeInteger(BigInt(value))); }

/** @brief 判断是否指定小整数。 */
bool isInt(const Ptr &node, long long value) {
    long long actual = 0;
    return node && node->isSmallInteger(&actual) && actual == value;
}

/** @brief 把有理数值节点拆成 (分子, 分母)。 */
bool rationalOf(const Ptr &node, BigInt *numerator, BigInt *denominator) {
    if (!node || !node->isNumber()) return false;
    if (numerator) *numerator = node->numerator();
    if (denominator) *denominator = node->denominator();
    return true;
}

/** @brief 多项式系数是否全为数值。 */
bool allNumeric(const UnivariatePoly &polynomial) {
    for (const Ptr &coefficient : polynomial.coefficients) {
        if (!coefficient || !coefficient->isNumber()) return false;
    }
    return true;
}

/* ---------------------------------------------------------- 展开内核 */

/** @brief 把表达式的加法项拆成列表（不是加法则单元素）。 */
Memory::Vector<Ptr> asTerms(const Ptr &expression);

/**
 * @brief 两个项列表的卷积（多项式乘法）。
 *
 * 每步用化简器合并同类项，因此系数始终是规范形式，后续比较可以直接用结构相等。
 */
Memory::Vector<Ptr> convolve(const Memory::Vector<Ptr> &left, const Memory::Vector<Ptr> &right,
                             std::size_t termLimit, bool *overflow) {
    if (left.empty() || right.empty()) return {};
    if (left.size() > termLimit / right.size()) {
        *overflow = true;
        return {};
    }
    // 逐对相乘后交给 makeAddMany 统一合并同类项；
    // 注意不能按 i+j 索引写回固定槽位，否则不同乘积会被错误合并或丢失
    Memory::Vector<Ptr> products;
    products.reserve(left.size() * right.size());
    for (std::size_t i = 0; i < left.size(); ++i) {
        for (std::size_t j = 0; j < right.size(); ++j) {
            Ptr term = makeMul(left[i], right[j]);
            if (!term) continue;
            products.push_back(term);
        }
    }
    if (products.empty()) return {};
    Memory::Vector<NodePtr> raw;
    raw.reserve(products.size());
    for (const Ptr &term : products) raw.push_back(term.get());
    Ptr combined = makeAddMany(raw);
    return asTerms(combined);
}

/** @brief 把表达式的加法项拆成列表（不是加法则单元素）。 */
Memory::Vector<Ptr> asTerms(const Ptr &expression) {
    Memory::Vector<Ptr> terms;
    if (!expression) return terms;
    if (expression->kind() == NodeKind::Add) {
        // 必须持有引用：调用方可能在本函数返回后就不再持有 expression
        for (NodePtr child : expression->operands()) terms.push_back(Ptr::borrow(child).share());
    } else {
        terms.push_back(expression.share());
    }
    return terms;
}

/** @brief 由项列表重建表达式。 */
Ptr fromTerms(const Memory::Vector<Ptr> &terms) {
    if (terms.empty()) return integer(0);
    if (terms.size() == 1) return terms.front().share();
    Memory::Vector<NodePtr> operands;
    operands.reserve(terms.size());
    for (const Ptr &term : terms) operands.push_back(term.get());
    return makeAddMany(operands);
}

/* -------------------------------------------------- 单变量多项式解析 */

/**
 * @brief 递归解析表达式为 variable 的多项式。
 *
 * 只接受加、乘、非负整数幂，以及不含 variable 的任意常数因子。
 */
bool parseInto(const Ptr &expression, std::uint32_t variable, UnivariatePoly *out) {
    if (!expression) return false;

    if (!containsSymbol(expression.get(), variable)) {
        out->coefficients.assign(1, expression);
        return true;
    }

    if (expression->kind() == NodeKind::Symbol) {
        if (expression->symbolId() != variable) return false;
        out->coefficients.assign(2, integer(0));
        out->coefficients[1] = integer(1);
        return true;
    }

    if (expression->kind() == NodeKind::Pow) {
        Ptr base = Ptr::borrow(expression->operand(0));
        Ptr exponent = Ptr::borrow(expression->operand(1));
        if (base->kind() != NodeKind::Symbol || base->symbolId() != variable) return false;
        long long power = 0;
        if (!exponent->isSmallInteger(&power) || power < 0) return false;
        out->coefficients.assign(static_cast<std::size_t>(power) + 1, integer(0));
        out->coefficients[static_cast<std::size_t>(power)] = integer(1);
        return true;
    }

    if (expression->kind() == NodeKind::Mul) {
        out->coefficients.assign(1, integer(1));
        for (NodePtr child : expression->operands()) {
            UnivariatePoly part;
            part.variable = variable;
            if (!parseInto(Ptr::borrow(child), variable, &part)) return false;
            Memory::Vector<Ptr> product(out->coefficients.size() + part.coefficients.size() - 1);
            for (std::size_t i = 0; i < out->coefficients.size(); ++i) {
                for (std::size_t j = 0; j < part.coefficients.size(); ++j) {
                    Ptr term = makeMul(out->coefficients[i], part.coefficients[j]);
                    product[i + j] = product[i + j] ? makeAdd(product[i + j], term) : term;
                }
            }
            out->coefficients = product;
        }
        out->trim();
        return true;
    }

    if (expression->kind() == NodeKind::Add) {
        out->coefficients.clear();
        for (NodePtr child : expression->operands()) {
            UnivariatePoly part;
            part.variable = variable;
            if (!parseInto(Ptr::borrow(child), variable, &part)) return false;
            if (out->coefficients.size() < part.coefficients.size())
                out->coefficients.resize(part.coefficients.size(), integer(0));
            for (std::size_t i = 0; i < part.coefficients.size(); ++i) {
                if (!part.coefficients[i]) continue;
                out->coefficients[i] = out->coefficients[i]
                                           ? makeAdd(out->coefficients[i], part.coefficients[i])
                                           : part.coefficients[i];
            }
        }
        out->trim();
        return true;
    }

    return false;
}

/** @brief 展开递归体。 */
Ptr expandRecursive(const Ptr &expression, std::size_t termLimit, bool *overflow) {
    if (!expression) return Ptr();
    if (*overflow) return expression.share();

    switch (expression->kind()) {
    case NodeKind::Integer:
    case NodeKind::Rational:
    case NodeKind::Symbol:
    case NodeKind::Constant:
        return expression.share();

    case NodeKind::Add: {
        Memory::Vector<Ptr> holders;
        for (NodePtr child : expression->operands()) {
            Ptr expanded = expandRecursive(Ptr::borrow(child), termLimit, overflow);
            if (*overflow) return expression.share();
            holders.push_back(expanded);
        }
        Memory::Vector<NodePtr> operands;
        for (const Ptr &holder : holders) operands.push_back(holder.get());
        return makeAddMany(operands);
    }

    case NodeKind::Mul: {
        Memory::Vector<Memory::Vector<Ptr>> blocks;
        for (NodePtr child : expression->operands()) {
            Ptr expanded = expandRecursive(Ptr::borrow(child), termLimit, overflow);
            if (*overflow) return expression.share();
            blocks.push_back(asTerms(expanded));
        }
        Memory::Vector<Ptr> accumulator;
        if (blocks.empty()) {
            accumulator.push_back(integer(1));
        } else {
            accumulator = blocks.front();
            for (std::size_t i = 1; i < blocks.size(); ++i) {
                accumulator = convolve(accumulator, blocks[i], termLimit, overflow);
                if (*overflow) return expression.share();
            }
        }
        return fromTerms(accumulator);
    }

    case NodeKind::Pow: {
        Ptr base = Ptr::borrow(expression->operand(0));
        Ptr exponent = Ptr::borrow(expression->operand(1));
        long long power = 0;
        if (!exponent->isSmallInteger(&power) || power < 0 || power > 64) {
            Ptr expandedBase = expandRecursive(base, termLimit, overflow);
            if (*overflow) return expression.share();
            return makePow(expandedBase, exponent);
        }
        Ptr expandedBase = expandRecursive(base, termLimit, overflow);
        if (*overflow) return expression.share();
        Memory::Vector<Ptr> terms = asTerms(expandedBase);
        if (terms.size() == 1) return makePow(expandedBase, exponent);

        Memory::Vector<Ptr> result;
        result.push_back(integer(1));
        Memory::Vector<Ptr> squarer = terms;
        long long remaining = power;
        while (remaining > 0) {
            if (remaining & 1) {
                result = convolve(result, squarer, termLimit, overflow);
                if (*overflow) return expression.share();
            }
            remaining >>= 1;
            if (remaining) {
                squarer = convolve(squarer, squarer, termLimit, overflow);
                if (*overflow) return expression.share();
            }
        }
        return fromTerms(result);
    }

    case NodeKind::Function: {
        Memory::Vector<Ptr> holders;
        for (NodePtr child : expression->operands()) {
            Ptr expanded = expandRecursive(Ptr::borrow(child), termLimit, overflow);
            if (*overflow) return expression.share();
            holders.push_back(expanded);
        }
        Memory::Vector<NodePtr> arguments;
        for (const Ptr &holder : holders) arguments.push_back(holder.get());
        return Ptr(Node::makeFunction(expression->name(), arguments));
    }
    }
    return expression.share();
}

/* ------------------------------------------------------------ GCD 内核 */

/**
 * @brief 数值系数的多项式 GCD（有理数域上，结果首一化）。
 *
 * 用伪除法避免中途出现分式：每步把余式乘以 `lc(b)^(deg a - deg r + 1)`。
 * 由于全程是精确有理运算，伪余式的零点与原余式一致，不影响 GCD。
 */
bool numericGcd(UnivariatePoly a, UnivariatePoly b, UnivariatePoly *out) {
    if (!allNumeric(a) || !allNumeric(b)) return false;
    a.trim();
    b.trim();
    while (!b.isZero()) {
        UnivariatePoly q, r;
        if (!polyDivide(a, b, &q, &r)) return false;
        a = b;
        b = r;
    }
    if (a.isZero()) {
        *out = a;
        return true;
    }
    // 首一化：除以首项系数
    Ptr leading = a.coefficient(a.degree());
    if (leading && leading->isNumber() && !leading->isZero()) {
        for (auto &coefficient : a.coefficients) coefficient = makeDivide(coefficient, leading);
        a.trim();
    }
    *out = a;
    return true;
}

/** @brief 从整系数多项式里提取所有系数的公共有理因子。 */
Ptr commonNumericFactor(const UnivariatePoly &polynomial) {
    BigInt numeratorGcd(0), denominatorLcm(1);
    bool any = false;
    for (const Ptr &coefficient : polynomial.coefficients) {
        BigInt n, d;
        if (!rationalOf(coefficient, &n, &d)) return integer(1);
        numeratorGcd = any ? BigInt::gcd(numeratorGcd, n) : n;
        denominatorLcm = denominatorLcm / BigInt::gcd(denominatorLcm, d) * d;
        any = true;
    }
    if (!any || numeratorGcd.isZero()) return integer(1);
    if (numeratorGcd.isOne() && denominatorLcm.isOne()) return integer(1);
    return makeNumber(numeratorGcd, denominatorLcm);
}

/** @brief 把多项式系数除以一个数值因子。 */
void divideCoefficients(UnivariatePoly *polynomial, const Ptr &divisor) {
    if (isInt(divisor, 1)) return;
    for (auto &coefficient : polynomial->coefficients) coefficient = makeDivide(coefficient, divisor);
    polynomial->trim();
}

/**
 * @brief 有理根定理：剥离所有一次因子。
 *
 * 对整系数多项式 a_n x^n + … + a_0，任何有理根 p/q 满足 p | a_0、q | a_n。
 * 逐个候选用 Horner 精确求值验证，命中后用综合除法降次。
 * 二次式额外用求根公式兜底，避免候选搜索范围设限导致漏根。
 */
Memory::Vector<Ptr> extractLinearFactors(UnivariatePoly *polynomial) {
    Memory::Vector<Ptr> factors;
    Ptr variable(Node::makeSymbolById(polynomial->variable));

    for (int guard = 0; guard < 256; ++guard) {
        if (polynomial->degree() < 1 || !allNumeric(*polynomial)) break;

        // 常数项为 0：x 本身是因子，直接降次
        if (polynomial->coefficient(0)->isZero()) {
            factors.push_back(variable.share());
            Memory::Vector<Ptr> shifted(polynomial->coefficients.begin() + 1, polynomial->coefficients.end());
            polynomial->coefficients = shifted;
            polynomial->trim();
            continue;
        }

        // 通分成整系数：乘以所有分母的最小公倍数
        BigInt lcm(1);
        for (const Ptr &coefficient : polynomial->coefficients) {
            BigInt d;
            if (!rationalOf(coefficient, nullptr, &d)) return factors;
            lcm = lcm / BigInt::gcd(lcm, d) * d;
        }
        Memory::Vector<BigInt> values;
        values.reserve(polynomial->coefficients.size());
        for (const Ptr &coefficient : polynomial->coefficients) {
            BigInt n, d;
            if (!rationalOf(coefficient, &n, &d)) return factors;
            values.push_back(n * (lcm / d));
        }

        const BigInt constant = values.front();
        const BigInt leading = values.back();
        BigInt rootNumerator, rootDenominator;
        bool found = false;

        // 候选：p 整除常数项、q 整除首项系数；用整数运算精确判零
        for (long long pValue = 1; pValue <= 2048 && !found; ++pValue) {
            const BigInt numerator(pValue);
            if (numerator > constant.abs()) break;
            if (!(constant % numerator).isZero()) continue;
            for (long long qValue = 1; qValue <= 2048 && !found; ++qValue) {
                const BigInt denominator(qValue);
                if (denominator > leading.abs()) break;
                if (!(leading % denominator).isZero()) continue;
                for (int sign : {1, -1}) {
                    const BigInt candidateNumerator = sign > 0 ? numerator : numerator.negated();
                    // Σ a_i·p^i·q^(n-i)，逐项累加避免分母
                    BigInt total(0);
                    BigInt pPower(1);
                    const std::size_t n = values.size();
                    // 先算 q^(n-1-i)
                    Memory::Vector<BigInt> qPowers(n, BigInt(1));
                    for (std::size_t i = 1; i < n; ++i) qPowers[i] = qPowers[i - 1] * denominator;
                    for (std::size_t i = 0; i < n; ++i) {
                        total = total + values[i] * pPower * qPowers[n - 1 - i];
                        pPower = pPower * candidateNumerator;
                    }
                    if (!total.isZero()) continue;
                    rootNumerator = candidateNumerator;
                    rootDenominator = denominator;
                    found = true;
                    break;
                }
            }
        }
        if (!found) break;

        // 因子取整数化写法 (q·x - p)，阅读与后续展开都方便
        Ptr factor = makeAdd(makeMul(makeNumber(rootDenominator, BigInt(1)), variable.share()),
                             makeNumber(rootNumerator.negated(), BigInt(1)));
        factors.push_back(factor);

        // 用带余除法精确剥离因子，不依赖手写综合除法的缩放
        UnivariatePoly divisorPart;
        divisorPart.variable = polynomial->variable;
        divisorPart.coefficients.clear();
        divisorPart.coefficients.push_back(makeNumber(rootNumerator.negated(), BigInt(1)));
        divisorPart.coefficients.push_back(makeNumber(rootDenominator, BigInt(1)));
        divisorPart.trim();

        UnivariatePoly quotientPart, remainderPart;
        if (!polyDivide(*polynomial, divisorPart, &quotientPart, &remainderPart)) return factors;
        if (!remainderPart.isZero()) return factors;
        *polynomial = quotientPart;
    }
    return factors;
}

} // namespace

/* ---------------------------------------------------- UnivariatePoly */

Ptr UnivariatePoly::coefficient(std::size_t index) const {
    if (index >= coefficients.size() || !coefficients[index]) return integer(0);
    return coefficients[index];
}

void UnivariatePoly::trim() {
    while (!coefficients.empty() && coefficients.back() && coefficients.back()->isZero()) coefficients.pop_back();
    if (coefficients.empty()) return;
    for (auto &coefficient : coefficients) {
        if (!coefficient) coefficient = integer(0);
    }
}

bool toPolynomial(const NodePtr &expression, std::uint32_t variable, UnivariatePoly *out) {
    out->variable = variable;
    out->coefficients.clear();
    Ptr simplified = simplify(Ptr::borrow(expression));
    if (!simplified) return false;
    if (!parseInto(simplified, variable, out)) return false;
    out->trim();
    return true;
}

Ptr fromPolynomial(const UnivariatePoly &polynomial) {
    Ptr result = integer(0);
    Ptr variable(Node::makeSymbolById(polynomial.variable));
    for (std::size_t i = polynomial.coefficients.size(); i-- > 0;) {
        Ptr coefficient = polynomial.coefficients[i];
        if (!coefficient || coefficient->isZero()) continue;
        Ptr term = coefficient;
        if (i > 0) {
            Ptr power = i == 1 ? variable.share() : makePow(variable, integer(static_cast<long long>(i)));
            term = makeMul(coefficient, power);
        }
        result = makeAdd(result, term);
    }
    return result;
}

/* ------------------------------------------------------------ 展开 */

Ptr expand(const Ptr &expression, std::size_t termLimit) {
    if (!expression) return Ptr();
    bool overflow = false;
    Ptr result = expandRecursive(expression, termLimit, &overflow);
    if (overflow || !result) return simplify(expression);
    return result;
}

/* -------------------------------------------------------- 多项式运算 */

Ptr polyAdd(const UnivariatePoly &left, const UnivariatePoly &right) {
    UnivariatePoly sum;
    sum.variable = left.variable;
    const std::size_t size = std::max(left.coefficients.size(), right.coefficients.size());
    sum.coefficients.resize(size);
    for (std::size_t i = 0; i < size; ++i) sum.coefficients[i] = makeAdd(left.coefficient(i), right.coefficient(i));
    sum.trim();
    return fromPolynomial(sum);
}

Ptr polyMultiply(const UnivariatePoly &left, const UnivariatePoly &right) {
    if (left.isZero() || right.isZero()) return integer(0);
    UnivariatePoly product;
    product.variable = left.variable;
    product.coefficients.assign(left.coefficients.size() + right.coefficients.size() - 1, integer(0));
    for (std::size_t i = 0; i < left.coefficients.size(); ++i) {
        for (std::size_t j = 0; j < right.coefficients.size(); ++j) {
            Ptr term = makeMul(left.coefficients[i], right.coefficients[j]);
            product.coefficients[i + j] = makeAdd(product.coefficients[i + j], term);
        }
    }
    product.trim();
    return fromPolynomial(product);
}

bool polyDivide(const UnivariatePoly &dividend, const UnivariatePoly &divisor, UnivariatePoly *quotient,
                UnivariatePoly *remainder) {
    *remainder = dividend;
    *quotient = UnivariatePoly();
    quotient->variable = dividend.variable;
    remainder->variable = dividend.variable;
    if (divisor.isZero()) return false;
    if (dividend.isZero() || dividend.degree() < divisor.degree()) {
        quotient->coefficients.assign(1, integer(0));
        return true;
    }
    if (!allNumeric(dividend) || !allNumeric(divisor)) return false;

    const std::size_t count = dividend.degree() - divisor.degree() + 1;
    quotient->coefficients.assign(count, integer(0));
    const Ptr divisorLead = divisor.coefficient(divisor.degree());

    for (std::size_t i = count; i-- > 0;) {
        if (remainder->isZero() || remainder->degree() < i + divisor.degree()) {
            quotient->coefficients[i] = integer(0);
            continue;
        }
        const Ptr leading = remainder->coefficient(i + divisor.degree());
        if (leading->isZero()) {
            quotient->coefficients[i] = integer(0);
            continue;
        }
        const Ptr factor = makeDivide(leading, divisorLead);
        quotient->coefficients[i] = factor;
        for (std::size_t j = 0; j <= divisor.degree(); ++j) {
            const std::size_t at = i + j;
            if (at >= remainder->coefficients.size()) continue;
            remainder->coefficients[at] = makeAdd(remainder->coefficients[at], makeNegate(makeMul(factor, divisor.coefficient(j))));
        }
        remainder->trim();
    }
    quotient->trim();
    return true;
}

Ptr polyGcd(const NodePtr &left, const NodePtr &right, std::uint32_t variable) {
    UnivariatePoly a, b;
    if (!toPolynomial(left, variable, &a) || !toPolynomial(right, variable, &b)) return integer(1);
    if (a.isZero()) return fromPolynomial(b);
    if (b.isZero()) return fromPolynomial(a);

    UnivariatePoly result;
    if (!numericGcd(a, b, &result)) return integer(1);
    return fromPolynomial(result);
}

/**
 * @brief Resultant，用 Sylvester 矩阵行列式定义（对低次情形直接算，高次用余式递归）。
 *
 * 递推式：`Res(A,B) = (-1)^(deg A·deg B) · lc(B)^(deg A - deg R) · Res(B, R)`，
 * 其中 R = A mod B。递推在 deg R == 0 时收敛：`Res(A, c) = c^deg(A)`。
 */
Ptr polyResultant(const NodePtr &left, const NodePtr &right, std::uint32_t variable) {
    UnivariatePoly a, b;
    if (!toPolynomial(left, variable, &a) || !toPolynomial(right, variable, &b)) return Ptr();
    if (a.isZero() || b.isZero()) return integer(0);

    Ptr accumulator = integer(1);
    bool negative = false;

    UnivariatePoly currentA = a, currentB = b;
    for (int guard = 0; guard < 512; ++guard) {
        currentA.trim();
        currentB.trim();
        if (currentB.isZero()) return integer(0);

        if (currentB.degree() == 0) {
            // Res(A, c) = c^deg(A)
            const long long exponent = static_cast<long long>(currentA.degree());
            accumulator = makeMul(accumulator, makePow(currentB.coefficient(0), integer(exponent)));
            return negative ? makeNegate(accumulator) : accumulator;
        }

        const std::size_t degA = currentA.degree(), degB = currentB.degree();
        if (degA < degB) {
            // 交换一对：Res(A,B) = (-1)^(degA·degB) Res(B,A)
            if ((degA * degB) % 2 == 1) negative = !negative;
            std::swap(currentA, currentB);
            continue;
        }

        UnivariatePoly q, r;
        if (!polyDivide(currentA, currentB, &q, &r)) return Ptr();
        // 符号修正：交换带来的 (-1)^(degA·degB)
        if ((degA * degB) % 2 == 1) negative = !negative;
        // 首项系数因子：lc(B)^(degA - degR)
        const std::size_t degR = r.isZero() ? 0 : r.degree();
        const long long power = static_cast<long long>(degA) - static_cast<long long>(degR);
        if (power > 0) accumulator = makeMul(accumulator, makePow(currentB.coefficient(degB), integer(power)));

        currentA = currentB;
        currentB = r;
    }
    return Ptr();
}

/* ------------------------------------------------------------ 因式分解 */

Ptr factor(const Ptr &expression) {
    Ptr simplified = simplify(expression);
    if (!simplified) return Ptr();

    const auto symbols = freeSymbols(simplified.get());
    if (symbols.size() != 1) return simplified;
    const std::uint32_t variable = symbols.front();

    // 已有乘积结构：逐因子分解
    if (simplified->kind() == NodeKind::Mul) {
        Memory::Vector<Ptr> pieces;
        for (NodePtr child : simplified->operands()) {
            Ptr part = factor(Ptr::borrow(child));
            if (!part) continue;
            // 若子因子本身是可分解的加法，递归结果可能又是乘积，直接摊平
            if (part->kind() == NodeKind::Mul) {
                for (NodePtr inner : part->operands()) pieces.push_back(Ptr::borrow(inner));
            } else {
                pieces.push_back(part);
            }
        }
        Memory::Vector<NodePtr> operands;
        for (const Ptr &piece : pieces) operands.push_back(piece.get());
        return makeMulMany(operands);
    }

    Ptr expanded = expand(simplified);
    UnivariatePoly polynomial;
    if (!toPolynomial(expanded.get(), variable, &polynomial)) return simplified;
    if (polynomial.degree() < 1) return simplified;

    Memory::Vector<Ptr> factors;

    // 1. 数值公因子
    Ptr numeric = commonNumericFactor(polynomial);
    divideCoefficients(&polynomial, numeric);
    if (!isInt(numeric, 1)) factors.push_back(numeric);

    // 2. 一次因子（有理根定理 + 综合除法）
    const auto linearFactors = extractLinearFactors(&polynomial);
    for (const Ptr &linear : linearFactors) factors.push_back(linear);

    // 4. 剩余部分
    if (polynomial.degree() >= 1) factors.push_back(fromPolynomial(polynomial));

    if (factors.empty()) return simplified;
    if (factors.size() == 1) return factors.front().share();
    Memory::Vector<NodePtr> operands;
    for (const Ptr &piece : factors) operands.push_back(piece.get());
    return makeMulMany(operands);
}

/* -------------------------------------------------------- 通分与整理 */

/** @brief `1/Q` 或 `top/Q` 形态拆成 (分子, 分母)；否则返回 false。 */
bool splitFraction(const Ptr &node, Ptr *numerator, Ptr *denominator) {
    if (!node) return false;
    if (node->kind() == NodeKind::Pow) {
        Ptr exponent = Ptr::borrow(node->operand(1));
        if (exponent->isNumber() && exponent->isNegative()) {
            *numerator = integer(1);
            *denominator = makePow(Ptr::borrow(node->operand(0)), makeNegate(exponent));
            return true;
        }
        return false;
    }
    if (node->kind() == NodeKind::Mul) {
        Ptr top = integer(1), bottom = integer(1);
        bool anyFraction = false;
        for (NodePtr child : node->operands()) {
            Ptr factor = Ptr::borrow(child);
            Ptr innerTop, innerBottom;
            if (splitFraction(factor, &innerTop, &innerBottom)) {
                top = makeMul(top, innerTop);
                bottom = makeMul(bottom, innerBottom);
                anyFraction = true;
            } else {
                top = makeMul(top, factor);
            }
        }
        if (!anyFraction) return false;
        *numerator = top;
        *denominator = bottom;
        return true;
    }
    return false;
}

/**
 * @brief 把任意表达式递归化成 `P/Q` 的规范形式（P、Q 里都不再含分式）。
 *
 * 递归规则：
 * - `Mul`：逐因子拆成 `Pᵢ/Qᵢ`，分子连乘、分母连乘；
 * - `Add`：各项统一到公分母，分子按公分母补齐；
 * - `Pow`：指数为整数 n 时取 `(P/Q)^n`，n 为负则分子分母互换。
 *
 * 这条通路专门服务万能代换：`2/((t²+1)·((1−t²)/(t²+1)+1))` 的最外层是一个
 * `Mul`，第二个因子是 `Pow(分母, -1)`，而分母自身还是 `Mul` + `Add` 的嵌套分式。
 * 普通 `together` 只合并最外层的加法项，看不到里层结构。
 */
bool toQuotient(const Ptr &expression, Ptr *numerator, Ptr *denominator);

bool toQuotient(const Ptr &expression, Ptr *numerator, Ptr *denominator) {
    if (!expression) return false;

    if (expression->kind() == NodeKind::Mul) {
        Ptr top = integer(1), bottom = integer(1);
        for (NodePtr child : expression->operands()) {
            Ptr subTop, subBottom;
            if (!toQuotient(Ptr::borrow(child), &subTop, &subBottom)) return false;
            top = makeMul(top, subTop);
            bottom = makeMul(bottom, subBottom);
        }
        *numerator = top;
        *denominator = bottom;
        return true;
    }

    if (expression->kind() == NodeKind::Add) {
        Memory::Vector<Ptr> tops, bottoms;
        for (NodePtr child : expression->operands()) {
            Ptr subTop, subBottom;
            if (!toQuotient(Ptr::borrow(child), &subTop, &subBottom)) return false;
            tops.push_back(subTop);
            bottoms.push_back(subBottom);
        }
        if (tops.empty()) {
            *numerator = integer(0);
            *denominator = integer(1);
            return true;
        }
        Ptr common = integer(1);
        for (const Ptr &piece : bottoms) common = makeMul(common, piece);
        Ptr sum = integer(0);
        for (std::size_t i = 0; i < tops.size(); ++i) {
            Ptr scaled = tops[i];
            for (std::size_t j = 0; j < bottoms.size(); ++j) {
                if (i != j) scaled = makeMul(scaled, bottoms[j]);
            }
            sum = makeAdd(sum, scaled);
        }
        *numerator = sum;
        *denominator = common;
        return true;
    }

    if (expression->kind() == NodeKind::Pow) {
        Ptr base = Ptr::borrow(expression->operand(0));
        Ptr exponent = Ptr::borrow(expression->operand(1));
        long long power = 0;
        if (exponent->isSmallInteger(&power)) {
            Ptr baseTop, baseBottom;
            if (!toQuotient(base, &baseTop, &baseBottom)) return false;
            const long long magnitude = power < 0 ? -power : power;
            Ptr topPower = makePow(baseTop, integer(magnitude));
            Ptr bottomPower = makePow(baseBottom, integer(magnitude));
            if (!topPower || !bottomPower) return false;
            if (power >= 0) {
                *numerator = topPower;
                *denominator = bottomPower;
            } else {
                *numerator = bottomPower;
                *denominator = topPower;
            }
            return true;
        }
        // 非整数指数：底数若含分式则放弃（无法保持单分式形态）
        Ptr baseTop, baseBottom;
        if (!toQuotient(base, &baseTop, &baseBottom)) return false;
        if (!isInt(baseBottom, 1)) return false;
        *numerator = expression;
        *denominator = integer(1);
        return true;
    }

    *numerator = expression;
    *denominator = integer(1);
    return true;
}

Ptr togetherDeep(const Ptr &expression) {
    Ptr simplified = simplify(expression);
    if (!simplified) return Ptr();
    Ptr numerator, denominator;
    if (!toQuotient(simplified, &numerator, &denominator)) return together(simplified);
    if (isInt(denominator, 1)) return simplified;
    Ptr merged = makeDivide(numerator, denominator);
    return merged ? merged : simplified;
}

Ptr together(const Ptr &expression) {
    Ptr simplified = simplify(expression);
    if (!simplified) return Ptr();
    const auto terms = asTerms(simplified);

    Memory::Vector<Ptr> numerators, denominators;
    for (const Ptr &term : terms) {
        Ptr numerator = integer(1), denominator = integer(1);
        Memory::Vector<Ptr> pieces;
        if (term->kind() == NodeKind::Mul) {
            for (NodePtr child : term->operands()) pieces.push_back(Ptr::borrow(child));
        } else {
            pieces.push_back(term.share());
        }
        for (const Ptr &piece : pieces) {
            if (piece->kind() == NodeKind::Pow) {
                Ptr exponent = Ptr::borrow(piece->operand(1));
                if (exponent->isNumber() && exponent->isNegative()) {
                    Ptr base = Ptr::borrow(piece->operand(0));
                    denominator = makeMul(denominator, makePow(base, makeNegate(exponent)));
                    continue;
                }
            }
            numerator = makeMul(numerator, piece);
        }
        numerators.push_back(numerator);
        denominators.push_back(denominator);
    }
    if (numerators.empty()) return integer(0);

    Ptr common = integer(1);
    for (const Ptr &denominator : denominators) {
        if (!isInt(denominator, 1)) common = makeMul(common, denominator);
    }
    if (isInt(common, 1)) return simplified;

    Ptr numerator = integer(0);
    for (std::size_t i = 0; i < numerators.size(); ++i) {
        Ptr scaled = numerators[i];
        for (std::size_t j = 0; j < denominators.size(); ++j) {
            if (i != j) scaled = makeMul(scaled, denominators[j]);
        }
        numerator = makeAdd(numerator, scaled);
    }
    return makeDivide(numerator, common);
}

Ptr collect(const Ptr &expression, const Memory::Vector<std::uint32_t> &variables) {
    Ptr expanded = expand(expression);
    if (variables.size() != 1) return expanded;
    UnivariatePoly polynomial;
    if (!toPolynomial(expanded.get(), variables.front(), &polynomial)) return expanded;
    return fromPolynomial(polynomial);
}

/**
 * @brief 有理函数化简：分子分母取多项式 GCD 后约分。
 */
Ptr cancelCommon(const Ptr &numerator, const Ptr &denominator, std::uint32_t variable) {
    UnivariatePoly top, bottom;
    if (!toPolynomial(numerator.get(), variable, &top)) return Ptr();
    if (!toPolynomial(denominator.get(), variable, &bottom)) return Ptr();
    if (top.isZero() || bottom.isZero()) return Ptr();

    UnivariatePoly divisor;
    if (!numericGcd(top, bottom, &divisor) || divisor.degree() < 1) return Ptr();
    Ptr gcdExpression = fromPolynomial(divisor);
    UnivariatePoly qTop, rTop, qBottom, rBottom;
    if (!polyDivide(top, divisor, &qTop, &rTop)) return Ptr();
    if (!polyDivide(bottom, divisor, &qBottom, &rBottom)) return Ptr();
    return makeDivide(fromPolynomial(qTop), fromPolynomial(qBottom));
}

} // namespace Cas::Detail

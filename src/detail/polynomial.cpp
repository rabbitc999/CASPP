// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#include "cas/detail/polynomial.h"
#include "cas/detail/calculus.h"
#include "cas/detail/printer.h"
#include <cstdio>

#include <algorithm>
#include <functional>
#include <utility>

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

/* ------------------------------------------------ 高次分解：结构识别 */

/**
 * @brief 偶次多项式在 `y = x²` 下可约时，把分解结果换回 x。
 *
 * `x^4+x^2+1`、`x^6-1`、`x^8-1` 的共同结构是"只含偶次项，且作为 y 的多项式能拆开"。
 *
 * 做法不做符号替换，而是直接在指数层面升次：y 上的因子 `Σ bᵢ·yⁱ`
 * 对应 x 上的 `Σ bᵢ·x^(2i)`；升次后若还能继续分解，就继续拆。
 *
 * @param polynomial 待分解的多项式（成功时会被清空）
 * @return 成功返回至少两个因子，否则返回空表示"这条规则不适用"。
 */
Memory::Vector<Ptr> factorByEvenSubstitution(UnivariatePoly *polynomial) {
    Memory::Vector<Ptr> factors;
    const std::size_t degree = polynomial->degree();
    if (degree < 4 || degree % 2 != 0) return factors;

    // 必须只含偶次项，才能写成 y = x² 的多项式
    for (std::size_t index = 0; index <= degree; ++index) {
        if (index % 2 != 0 && !polynomial->coefficient(index)->isZero()) return factors;
    }

    // 换元：把多项式按 `y = x²` 重写成 y 的多项式。
    // 这里必须真的用独立符号 y 来构造，不能借用原变量——
    // 否则后续无法区分"y 的幂"和"x 的幂"。
    const std::uint32_t yId = internSymbol("__cas_y");
    UnivariatePoly inY;
    inY.variable = yId;
    for (std::size_t index = 0; index <= degree / 2; ++index) {
        inY.coefficients.push_back(polynomial->coefficient(index * 2));
    }
    inY.trim();
    if (inY.degree() < 2) return factors;

    // 先在 y 上分解。
    Ptr originalInY = fromPolynomial(inY);
    Ptr factoredY = factor(originalInY);
    if (!factoredY) return factors;
    {
        PrintOptions po;
        po.useUnicode = false;
        if (toText(factoredY.get(), po) == toText(originalInY.get(), po)) return factors;
    }

    // 把 y 上的分解结果整体升次：y → x²。
    // 注意这里不能再调用 factor（升次后仍可能是偶次式，会无限递归），
    // 也不要无脑 expand——`(y+1)^2` 升次后是 `(x^2+1)^2`，
    // 那本身就是一份完整的分解结果，展开反而把它揉回成加法。
    Ptr replaced = substitute(factoredY, yId, makePow(Ptr(Node::makeSymbolById(polynomial->variable)),
                                                      integer(2)));
    if (!replaced) return Memory::Vector<Ptr>();

    // 递归收集「乘积因子」：Mul 摊平，Pow 保留为整体（形如 (x^2+1)^2）
    std::function<void(const Ptr &)> collect = [&](const Ptr &node) {
        if (!node) return;
        if (node->kind() == NodeKind::Mul) {
            for (NodePtr inner : node->operands()) collect(Ptr::borrow(inner));
            return;
        }
        // 单个加法因子（如 x^2-1）还能拆一次线性因子，这里只做局部提取，
        // 不回调 factor，避免重新进入本函数。
        UnivariatePoly piecePoly;
        if (node->kind() == NodeKind::Add &&
            toPolynomial(node.get(), polynomial->variable, &piecePoly) && piecePoly.degree() >= 2) {
            const auto linear = extractLinearFactors(&piecePoly);
            if (!linear.empty()) {
                for (const Ptr &one : linear) factors.push_back(one);
                if (piecePoly.degree() >= 1) factors.push_back(fromPolynomial(piecePoly));
                return;
            }
        }
        factors.push_back(node);
    };
    collect(replaced);

    // 判据不是"因子数量 ≥ 2"——`(x^2+1)^2` 只有一个因子但确实是更细的形式。
    // 改为比较升次结果与原式是否不同。
    {
        PrintOptions po;
        po.useUnicode = false;
        Ptr self = fromPolynomial(*polynomial);
        if (toText(replaced.get(), po) == toText(self.get(), po)) return Memory::Vector<Ptr>();
    }
    if (factors.empty()) return Memory::Vector<Ptr>();

    polynomial->coefficients.clear();
    return factors;
}

/**
 * @brief 识别符号表达式的完全平方，返回其平方根。
 *
 * 只处理几类常见形状，够多元二次分解用：
 * - 单个 `Pow(base, 2)` → base
 * - 乘积 `k · base²`（k 为完全平方有理数）→ √k · base
 * - 已经是平方的乘法链 `a²·b²·…` → a·b·…
 * 其余返回空，交给调用方当作"不可约"。
 */
Ptr sqrtOfSymbolicSquare(const Ptr &value) {
    if (!value) return Ptr();

    // 纯数值：交给有理数开方
    if (value->isNumber()) {
        BigInt numerator(0), denominator(1);
        if (!rationalOf(value, &numerator, &denominator)) return Ptr();
        if (numerator.isNegative()) return Ptr();
        if (!numerator.isPerfectSquare() || !denominator.isPerfectSquare()) return Ptr();
        return makeNumber(numerator.isqrt(), denominator.isqrt());
    }

    if (value->kind() == NodeKind::Pow) {
        Ptr base = Ptr::borrow(value->operand(0));
        Ptr exponent = Ptr::borrow(value->operand(1));
        long long power = 0;
        if (exponent->isSmallInteger(&power) && power >= 2 && power % 2 == 0) {
            return makePow(base, integer(power / 2));
        }
        return Ptr();
    }

    if (value->kind() == NodeKind::Mul) {
        // 逐个因子开方：每个因子都得是完全平方
        Ptr product = integer(1);
        for (NodePtr child : value->operands()) {
            Ptr piece = sqrtOfSymbolicSquare(Ptr::borrow(child));
            if (!piece) return Ptr();
            Ptr merged = makeMul(product, piece);
            if (!merged) return Ptr();
            product = merged;
        }
        return product;
    }
    return Ptr();
}

/**
 * @brief 把二次式拆成两个一次因子（需要判别式为完全平方数）。
 * @return 成功返回两个因子，否则为空。
 */
Memory::Vector<Ptr> splitQuadratic(UnivariatePoly *polynomial) {
    Memory::Vector<Ptr> factors;
    if (polynomial->degree() != 2) return factors;
    Ptr a = polynomial->coefficient(2);
    Ptr b = polynomial->coefficient(1);
    Ptr c = polynomial->coefficient(0);

    // 判别式 b²−4ac。数值判别式要求是非负完全平方；
    // 符号判别式（如 x^2−y^2 里的 4y²）则尝试识别成 `(表达式)²`。
    Ptr discriminant = simplify(makeAdd(makePow(b, integer(2)),
                                        makeNegate(makeMul(makeNumber(BigInt(4), BigInt(1)),
                                                           makeMul(a, c)))));
    if (!discriminant) return factors;

    Ptr root;          // 判别式的平方根（有理数或符号表达式）
    bool perfectSquare = false;

    BigInt numerator(0), denominator(1);
    if (rationalOf(discriminant, &numerator, &denominator)) {
        if (numerator.isNegative()) return factors;
        if (!numerator.isPerfectSquare() || !denominator.isPerfectSquare()) return factors;
        const BigInt rootNumerator = numerator.isqrt();
        const BigInt rootDenominator = denominator.isqrt();
        // 判别式为 0：二次式本身是完全平方 (x − r)²，
        // 两个重根要各出一个因子，否则会丢掉平方（x^2+2xy+y^2 拆成 x+y）。
        perfectSquare = rootNumerator.isZero();
        root = makeNumber(rootNumerator, rootDenominator);
    } else {
        // 符号判别式：找形如 `k²·(...)²` 或 `±(...)²` 的平方根
        Ptr squared = sqrtOfSymbolicSquare(discriminant);
        if (!squared) return factors;
        root = squared;
        perfectSquare = root->isZero();
    }
    Ptr twiceA = makeNumber(BigInt(2), BigInt(1));
    twiceA = makeMul(twiceA, a);

    Ptr rootOne = simplify(makeDivide(makeAdd(makeNegate(b), root), twiceA));
    Ptr rootTwo = simplify(makeDivide(makeAdd(makeNegate(b), makeNegate(root)), twiceA));
    if (!rootOne || !rootTwo) return factors;

    // 注意不能因为根含符号分母就放弃：`x²y²−1` 按 x 读时根是 ±1/y，
    // 而 `a·(x − 1/y) = y²x − y = y(xy − 1)` 通分后就变成整系数了。
    // 是否有残留分式交给组装阶段判定。

    Ptr variable(Node::makeSymbolById(polynomial->variable));
    // 根含分母 2a，直接写成 `(x − root)` 会带出分式。
    // 这里统一用 `a·(x − root₁)` 与 `(x − root₂)`，再分别通分成整系数，
    // 分母并入数值因子——这样 4x²−9y² 得到 4·(x−3/2·y)(x+3/2·y) → (2x−3y)(2x+3y)。
    auto makeLinear = [&](const Ptr &rootValue) {
        return makeAdd(variable.share(), makeNegate(rootValue));
    };
    Ptr first = simplify(makeMul(a, makeLinear(rootOne)));
    Ptr second = perfectSquare ? makeLinear(rootOne) : makeLinear(rootTwo);
    if (!first) return factors;
    if (!second) second = makeLinear(rootTwo);

    // 把因子整理成整系数形式：逐项取出有理系数，求出所有分母的最小公倍数，
    // 分子留下做新因子，分母并入总的数值系数。
    // （toQuotient 只在整体是单个分式时有效，逐项近似分式它看不出来。）
    auto integralize = [](const Ptr &factor, Ptr *numerator, Ptr *denominator) {
        Ptr value = expand(factor);
        if (!value) value = factor;

        // 拆成项：Add 逐项，否则单项
        Memory::Vector<Ptr> terms;
        if (value->kind() == NodeKind::Add) {
            for (NodePtr child : value->operands()) terms.push_back(Ptr::borrow(child));
        } else {
            terms.push_back(value);
        }

        // 每项拆成「系数 × 基」，求系数分母的最小公倍数
        BigInt lcm(1);
        Memory::Vector<Ptr> coefficients;
        Memory::Vector<Ptr> bodies;
        bool allRational = true;
        for (const Ptr &term : terms) {
            Ptr coefficient = integer(1);
            Ptr body = term;
            if (term->kind() == NodeKind::Mul) {
                Ptr rest = integer(1);
                for (NodePtr child : term->operands()) {
                    Ptr piece = Ptr::borrow(child);
                    if (piece->isNumber()) {
                        Ptr merged = makeMul(coefficient, piece);
                        if (merged) coefficient = merged;
                    } else {
                        Ptr merged = makeMul(rest, piece);
                        if (merged) rest = merged;
                    }
                }
                body = rest;
            } else if (term->isNumber()) {
                coefficient = term;
                body = integer(1);
            }
            BigInt numeratorValue(0), denominatorValue(1);
            if (!rationalOf(coefficient, &numeratorValue, &denominatorValue)) {
                allRational = false;
                break;
            }
            lcm = lcm / BigInt::gcd(lcm, denominatorValue) * denominatorValue;
            coefficients.push_back(coefficient);
            bodies.push_back(body);
        }
        if (!allRational) {
            *numerator = value;
            *denominator = integer(1);
            return;
        }

        Ptr scaled = integer(0);
        for (std::size_t index = 0; index < terms.size(); ++index) {
            BigInt numeratorValue(0), denominatorValue(1);
            if (!rationalOf(coefficients[index], &numeratorValue, &denominatorValue)) {
                *numerator = value;
                *denominator = integer(1);
                return;
            }
            Ptr newCoefficient = makeNumber(numeratorValue * (lcm / denominatorValue), BigInt(1));
            Ptr term = bodies[index]->isOne() ? newCoefficient : makeMul(newCoefficient, bodies[index]);
            Ptr sum = makeAdd(scaled, term);
            if (!sum) {
                *numerator = value;
                *denominator = integer(1);
                return;
            }
            scaled = sum;
        }
        *numerator = simplify(scaled);
        *denominator = integer(lcm.fitsLongLong() ? lcm.toLongLong() : 1);
    };

    Ptr topOne, bottomOne, topTwo, bottomTwo;
    integralize(first, &topOne, &bottomOne);
    integralize(second, &topTwo, &bottomTwo);

    Ptr denominators = makeMul(bottomOne, bottomTwo);
    Ptr coefficient = simplify(makeDivide(integer(1), denominators));

    if (coefficient && !coefficient->isOne()) factors.push_back(coefficient);
    factors.push_back(topOne);
    factors.push_back(topTwo);

    // 分解成功，多项式按定义已经全部拆完
    polynomial->coefficients.clear();
    return factors;
}

/**
 * @brief 把首一四次式拆成两个首一二次因子。
 *
 * 对 `x^4 + A·x^3 + B·x^2 + C·x + D`，若存在 `(x^2+p·x+q)(x^2+s·x+r)`，
 * 比较系数得：
 *   q + r + p·s = B，p·r + q·s = C，p + s = A，q·r = D。
 * 这里只处理最常见的对称情形 `s = -p`（此时 A = 0，即缺少三次项），
 * 于是 `q + r = B + p²`、`p·(r - q) = C`、`q·r = D`。
 * 对 p = ±1, ±2, … 逐个试（p 必须整除 C 或 C = 0 才有解），能覆盖
 * `x^4+x^2+1`、`x^4+4`（Sophie Germain）、`x^4+4x^2+4` 一类。
 *
 * @return 成功返回两个因子，否则为空。
 */
Memory::Vector<Ptr> splitQuartic(UnivariatePoly *polynomial) {
    Memory::Vector<Ptr> factors;
    if (polynomial->degree() != 4) return factors;

    // 只处理首一且无三次项的情况
    Ptr leading = polynomial->coefficient(4);
    if (!isInt(leading, 1)) return factors;
    if (!polynomial->coefficient(3)->isZero()) return factors;

    Ptr a2 = polynomial->coefficient(2);
    Ptr a1 = polynomial->coefficient(1);
    Ptr a0 = polynomial->coefficient(0);

    Ptr variable(Node::makeSymbolById(polynomial->variable));

    // p 取整数试探：结果里会出 p²、p，因此整数 p 覆盖了绝大多数可分解情形
    for (long long pValue = 1; pValue <= 64; ++pValue) {
        for (int sign : {1, -1}) {
            Ptr p = integer(sign * pValue);
            Ptr pSquared = simplify(makePow(p, integer(2)));

            // q + r = B + p²，q·r = D，配合 p·(r − q) = C
            // r − q 必须等于 C/p（p ≠ 0 时可直接算出）
            Ptr difference = simplify(makeDivide(a1, p));
            Ptr sum = simplify(makeAdd(a2, pSquared));

            // q = (sum − difference)/2，r = (sum + difference)/2
            Ptr half = makeNumber(BigInt(1), BigInt(2));
            Ptr q = simplify(makeMul(half, makeAdd(sum, makeNegate(difference))));
            Ptr r = simplify(makeMul(half, makeAdd(sum, difference)));
            if (!q || !r) continue;

            // 校验 q·r = D，避免浮点式误判
            Ptr product = simplify(makeMul(q, r));
            Ptr residual = product ? simplify(makeAdd(product, makeNegate(a0))) : Ptr();
            if (!residual || !residual->isZero()) continue;

            Ptr first = makeAdd(makeAdd(makePow(variable.share(), integer(2)),
                                        makeMul(p, variable.share())),
                                q);
            Ptr second = makeAdd(makeAdd(makePow(variable.share(), integer(2)),
                                         makeMul(makeNegate(p), variable.share())),
                                 r);
            factors.push_back(simplify(first));
            factors.push_back(simplify(second));
            polynomial->coefficients.clear();
            return factors;
        }
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

/**
 * @brief 把表达式整理成 `系数 × 因子连乘` 的形式，因子按文本去重合并成幂。
 *
 * 多元分解会在不同变量上分别尝试，得到的结果经常是同一份因子的不同排列，
 * 这里统一展开、按文本排序、相同因子并成幂次，保证输出稳定且紧凑。
 */
Ptr assembleFactors(const Memory::Vector<Ptr> &input) {
    if (input.empty()) return Ptr();

    // 摊平嵌套乘积，按文本归类计数
    Memory::Vector<Ptr> flat;
    std::function<void(const Ptr &)> push = [&](const Ptr &node) {
        if (!node) return;
        if (node->kind() == NodeKind::Mul) {
            for (NodePtr inner : node->operands()) push(Ptr::borrow(inner));
            return;
        }
        flat.push_back(node);
    };
    for (const Ptr &item : input) push(item);
    if (flat.empty()) return Ptr();

    PrintOptions options;
    options.useUnicode = false;

    Ptr coefficient = integer(1);
    Memory::Vector<std::pair<std::string, std::pair<Ptr, long long>>> groups;
    for (const Ptr &item : flat) {
        if (!item) continue;
        // 纯数值因子并入总系数
        if (item->isNumber()) {
            Ptr merged = makeMul(coefficient, item);
            if (merged) coefficient = merged;
            continue;
        }
        const std::string key = toText(item.get(), options);
        bool found = false;
        for (auto &group : groups) {
            if (group.first == key) {
                group.second.second += 1;
                found = true;
                break;
            }
        }
        if (!found) groups.push_back({key, {item.share(), 1}});
    }

    // 按文本排序，保证同一表达式每次输出一致
    std::sort(groups.begin(), groups.end(),
              [](const auto &left, const auto &right) { return left.first < right.first; });

    // 把各因子里的数值公因子提到总系数上，形式更干净：
    // `1/2·(4x−6y)(2x+3y)` → `2·(2x−3y)(2x+3y)`。
    // 这里不要求因子是单变量的——直接对展开式逐项的数值系数求 gcd，
    // 多变量因子（如 4x−6y）同样适用。
    for (auto &group : groups) {
        Ptr factorValue = group.second.first;
        Ptr expandedFactor = expand(factorValue);
        if (!expandedFactor) continue;

        Memory::Vector<Ptr> terms;
        if (expandedFactor->kind() == NodeKind::Add) {
            for (NodePtr child : expandedFactor->operands()) terms.push_back(Ptr::borrow(child));
        } else {
            terms.push_back(expandedFactor);
        }

        // 每项取出数值系数，求 gcd；遇到非有理系数就跳过整个因子
        BigInt common(0);
        Memory::Vector<std::pair<Ptr, Ptr>> pieces; // (系数, 基)
        bool rational = true;
        for (const Ptr &term : terms) {
            Ptr coefficient = integer(1);
            Ptr body = term;
            if (term->kind() == NodeKind::Mul) {
                Ptr rest = integer(1);
                for (NodePtr child : term->operands()) {
                    Ptr piece = Ptr::borrow(child);
                    if (piece->isNumber()) {
                        Ptr merged = makeMul(coefficient, piece);
                        if (merged) coefficient = merged;
                    } else {
                        Ptr merged = makeMul(rest, piece);
                        if (merged) rest = merged;
                    }
                }
                body = rest;
            } else if (term->isNumber()) {
                coefficient = term;
                body = integer(1);
            }
            BigInt numeratorValue(0), denominatorValue(1);
            if (!rationalOf(coefficient, &numeratorValue, &denominatorValue)) {
                rational = false;
                break;
            }
            common = common.isZero() ? numeratorValue.abs() : BigInt::gcd(common, numeratorValue.abs());
            pieces.push_back({coefficient, body});
        }
        if (!rational || common.isZero() || common.isOne()) continue;

        // 用 gcd 重新构造因子：系数除以 common
        Ptr scaled = integer(0);
        bool ok = true;
        for (const auto &piece : pieces) {
            Ptr reduced = simplify(makeDivide(piece.first, integer(common.fitsLongLong() ? common.toLongLong() : 1)));
            if (!reduced) {
                ok = false;
                break;
            }
            Ptr term = piece.second->isOne() ? reduced : makeMul(reduced, piece.second);
            Ptr sum = makeAdd(scaled, term);
            if (!sum) {
                ok = false;
                break;
            }
            scaled = sum;
        }
        if (!ok) continue;

        Ptr divided = simplify(scaled);
        if (!divided) continue;
        group.first = toText(divided.get(), options);
        group.second.first = divided;
        Ptr merged = makeMul(coefficient, integer(common.fitsLongLong() ? common.toLongLong() : 1));
        if (merged) coefficient = merged;
    }

    // 重新去重（提取公因子后可能整体重复）
    {
        Memory::Vector<std::pair<std::string, std::pair<Ptr, long long>>> dedup;
        for (const auto &group : groups) {
            bool found = false;
            for (auto &other : dedup) {
                if (other.first == group.first) {
                    other.second.second += group.second.second;
                    found = true;
                    break;
                }
            }
            if (!found) dedup.push_back(group);
        }
        groups = dedup;
        std::sort(groups.begin(), groups.end(),
                  [](const auto &left, const auto &right) { return left.first < right.first; });
    }

    // 必须用 Ptr 持有每个操作数直到 makeMulMany 调用结束——
    // 只存裸 NodePtr 的话，生成它们的临时 Ptr 会在语句结束时析构，
    // 引用计数归零后节点被回收，留下的就是悬垂指针。
    Memory::Vector<Ptr> holders;
    if (!isInt(coefficient, 1)) holders.push_back(coefficient);
    for (const auto &group : groups) {
        if (group.second.second <= 1) {
            holders.push_back(group.second.first);
        } else {
            Ptr power = makePow(group.second.first, integer(group.second.second));
            if (power) holders.push_back(power);
        }
    }
    if (holders.empty()) return integer(1);
    if (holders.size() == 1) return holders.front().share();

    Memory::Vector<NodePtr> operands;
    operands.reserve(holders.size());
    for (const Ptr &holder : holders) operands.push_back(holder.get());
    return makeMulMany(operands);
}

/**
 * @brief 多元因式分解：逐个变量当作一元多项式尝试。
 *
 * 做法是对每个自由符号 `v`，把表达式视作 `v` 的一元多项式（其余符号作为系数），
 * 走一遍一元分解；只要某个 `v` 上能拆开，就把结果当作新的候选，继续对因子递归。
 * 这条策略覆盖了教科书层面的多元分解：
 *
 *   x^2 - y^2          （按 x 看，常数项 -y^2，判别式 4y^2 是平方）
 *   x^2 + 2xy + y^2    （按 x 看，判别式 0）
 *   x^2*y^2 - 1        （按 x 看，是 x^2 减常数）
 *   x*y + x + y + 1    （按 x 看，一次式，系数含 y）
 *
 * 子因子数明显多于输入时放弃，避免组合爆炸。
 */
Ptr factorMultivariate(const Ptr &expression, const Memory::Vector<std::uint32_t> &symbols) {
    PrintOptions options;
    options.useUnicode = false;
    const std::string originalText = toText(expression.get(), options);

    Ptr best; // 当前找到的最优结果
    for (std::uint32_t variable : symbols) {
        UnivariatePoly polynomial;
        if (!toPolynomial(expression.get(), variable, &polynomial)) continue;
        if (polynomial.degree() < 1) continue;

        // 复用一元路径：把系数当作不透明的表达式，走同一套结构规则。
        // 这里只做"提取数值公因子 + 线性因子 + 二次拆分"这几条，
        // 不对符号系数做有理根搜索（那是无意义的）。
        Memory::Vector<Ptr> pieces;

        if (polynomial.degree() == 2) {
            UnivariatePoly copy = polynomial;
            const auto quadratics = splitQuadratic(&copy);
            for (const Ptr &piece : quadratics) pieces.push_back(piece);
        } else if (polynomial.degree() == 4) {
            UnivariatePoly copy = polynomial;
            const auto quartics = splitQuartic(&copy);
            for (const Ptr &piece : quartics) pieces.push_back(piece);
            if (!quartics.empty()) polynomial = copy;
        }

        if (pieces.empty()) continue;

        Ptr candidate = assembleFactors(pieces);
        if (!candidate) continue;
        const std::string candidateText = toText(candidate.get(), options);
        if (candidateText == originalText) continue;
        // 注意不要用"文本更短"来过滤：`(x - y)*(x + y)` 比 `x^2 - y^2` 长，
        // 但它才是分解结果。只要拆开了就接受。

        // 对每个因子继续递归，把还能拆的部分拆干净
        Memory::Vector<Ptr> refined;
        if (candidate->kind() == NodeKind::Mul) {
            for (NodePtr child : candidate->operands()) {
                Ptr one = Ptr::borrow(child);
                const auto childSymbols = freeSymbols(one.get());
                if (childSymbols.size() <= 1) {
                    Ptr deeper = factor(one);
                    refined.push_back(deeper ? deeper : one);
                } else {
                    refined.push_back(one);
                }
            }
        } else {
            refined.push_back(candidate);
        }

        Ptr result = assembleFactors(refined);
        if (!result) continue;
        // 含负指数的候选直接排除：`(x + 1/y)(x*y^2 − y)` 展开就是原式，
        // 属于"表面上拆开、实际只是把分母挪来挪去"，不是真正的分解。
        PrintOptions local;
        local.useUnicode = false;
        auto hasNegativePower = [](const Ptr &node) {
            std::function<bool(const Ptr &)> walk = [&](const Ptr &current) -> bool {
                if (!current) return false;
                if (current->kind() == NodeKind::Pow) {
                    Ptr exponent = Ptr::borrow(current->operand(1));
                    if (exponent->isNumber() && exponent->isNegative()) return true;
                }
                for (std::size_t index = 0; index < current->size(); ++index) {
                    if (walk(Ptr::borrow(current->operand(index)))) return true;
                }
                return false;
            };
            return walk(node);
        };
        if (hasNegativePower(result)) continue;

        // 剩下比因子数：多的更好（拆得更彻底）。
        const std::size_t candidateCount = result->kind() == NodeKind::Mul ? result->size() : 1;
        const std::size_t bestCount = best && best->kind() == NodeKind::Mul ? best->size() : 1;
        if (!best || candidateCount > bestCount) best = result;
    }
    return best;
}

Ptr factor(const Ptr &expression) {
    Ptr simplified = simplify(expression);
    if (!simplified) return Ptr();

    const auto symbols = freeSymbols(simplified.get());
    if (symbols.size() > 1) {
        // 多元路径要把表达式摊成标准和式，才能按变量读系数
        Ptr expanded = expand(simplified);
        Ptr factored = factorMultivariate(expanded ? expanded : simplified, symbols);
        // 拆不开（例如 x^4+y^4 在有理域上不可约）就原样返回，不要当失败
        return factored ? factored : (expanded ? expanded : simplified);
    }
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
    if (polynomial.degree() >= 1) {
        const auto linearFactors = extractLinearFactors(&polynomial);
        for (const Ptr &linear : linearFactors) factors.push_back(linear);
    }

    // 3. 剩余部分仍是二次：判别式为完全平方时可直接拆成两个一次因子
    if (polynomial.degree() == 2) {
        const auto quadratics = splitQuadratic(&polynomial);
        for (const Ptr &piece : quadratics) factors.push_back(piece);
    }

    // 4. 四次式拆两个首一二次因子：x^4+x^2+1 → (x^2+x+1)(x^2-x+1)
    if (polynomial.degree() == 4) {
        const auto quarticFactors = splitQuartic(&polynomial);
        for (const Ptr &piece : quarticFactors) factors.push_back(piece);
    }

    // 5. 偶次代换：x^8-1 → (x^2+1)(x^4+1)，x^6-1 之类也走这条
    if (polynomial.degree() >= 4) {
        const auto evenFactors = factorByEvenSubstitution(&polynomial);
        for (const Ptr &piece : evenFactors) factors.push_back(piece);
    }

    // 6. 剩余部分
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

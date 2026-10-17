#ifndef CAS_DETAIL_POLYNOMIAL_H
#define CAS_DETAIL_POLYNOMIAL_H

/**
 * @file polynomial.h
 * @brief 单变量多项式（系数为任意表达式）与多元多项式工具。
 *
 * 设计取舍：
 *  - 单变量多项式用稠密系数数组 `coeff[i] * x^i`，便于 GCD、Resultant、求根；
 *    系数是通用 Node，因此支持 `(a+1)x^2 + b` 这类符号系数。
 *  - 展开与多元处理直接在表达式树上做（`expand`），比先转成多元多项式环再转回来
 *    更简单，也不会在中间表示上丢精度。
 */

#include "cas/detail/node.h"
#include "cas/detail/simplify.h"

namespace Cas::Detail {

/** @brief 单变量稠密多项式：coeffs[i] 是 x^i 的系数，最高位非零（除零多项式）。 */
struct UnivariatePoly {
    std::uint32_t variable = 0;      ///< 变量符号编号
    Memory::Vector<Ptr> coefficients; ///< 从低次到高次

    /** @brief 是否为零多项式。 */
    bool isZero() const { return coefficients.empty(); }
    /** @brief 次数；零多项式返回 0。 */
    std::size_t degree() const { return coefficients.empty() ? 0 : coefficients.size() - 1; }
    /** @brief 取第 i 个系数（越界返回 0 节点）。 */
    Ptr coefficient(std::size_t index) const;
    /** @brief 去掉高阶零系数。 */
    void trim();
};

/** @brief 把表达式当作 variable 的多项式解析；失败（含非多项式项）返回 false。 */
bool toPolynomial(const NodePtr &expression, std::uint32_t variable, UnivariatePoly *out);

/** @brief 多项式转回表达式。 */
Ptr fromPolynomial(const UnivariatePoly &polynomial);

/* ------------------------------------------------------------------ 展开 */

/**
 * @brief 展开乘积与整数幂：`(x+1)^3` → `x^3+3x^2+3x+1`。
 * @param termLimit 展开后项数上限，超过则放弃展开原样返回，防止组合爆炸。
 */
Ptr expand(const Ptr &expression, std::size_t termLimit = 100000);

/* -------------------------------------------------------- 多项式运算 */

/** @brief 多项式加法。 */
Ptr polyAdd(const UnivariatePoly &left, const UnivariatePoly &right);
/** @brief 多项式乘法。 */
Ptr polyMultiply(const UnivariatePoly &left, const UnivariatePoly &right);
/** @brief 多项式带余除法：返回商与余式。 */
bool polyDivide(const UnivariatePoly &dividend, const UnivariatePoly &divisor, UnivariatePoly *quotient,
                UnivariatePoly *remainder);
/** @brief 多项式 GCD（首一化）。 */
Ptr polyGcd(const NodePtr &left, const NodePtr &right, std::uint32_t variable);
/** @brief 多项式 Resultant（用于消元与判别式）。 */
Ptr polyResultant(const NodePtr &left, const NodePtr &right, std::uint32_t variable);

/* ------------------------------------------------------------ 因式分解 */

/**
 * @brief 因式分解（有理数域上的单变量情形）。
 *
 * 覆盖：公因子提取、平方差/完全平方识别、有理根定理求一次因子、
 * 二次式判别式分解。无法完全分解时返回"已分解部分 × 剩余部分"。
 */
Ptr factor(const Ptr &expression);

/** @brief 通分合并各项，返回单一分式。 */
Ptr together(const Ptr &expression);

/** @brief 有理函数约分：分子分母取多项式 GCD 后相约。 */
Ptr cancelCommon(const Ptr &numerator, const Ptr &denominator, std::uint32_t variable);

/** @brief 对已展开多项式按 variable 整理成降幂形式。 */
Ptr collect(const Ptr &expression, const Memory::Vector<std::uint32_t> &variables);

} // namespace Cas::Detail

#endif // CAS_DETAIL_POLYNOMIAL_H

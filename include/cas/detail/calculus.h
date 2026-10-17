#ifndef CAS_DETAIL_CALCULUS_H
#define CAS_DETAIL_CALCULUS_H

/**
 * @file calculus.h
 * @brief 求导、不定积分、定积分、极限与 Taylor 展开。
 *
 * 各模块的能力边界写得很明确：做不了的返回"原样保留"或 null，绝不猜测。
 * 对 CAS 来说，返回一个错误答案比返回"未求值形式"危害大得多。
 */

#include "cas/detail/node.h"
#include "cas/detail/simplify.h"

namespace Cas::Detail {

/* -------------------------------------------------------------- 求导 */

/** @brief 对 variable 求导（符号求导，含链式法则与常见函数）。 */
Ptr differentiate(const Ptr &expression, std::uint32_t variable);
/** @brief 求 n 阶导。 */
Ptr differentiateN(const Ptr &expression, std::uint32_t variable, int order);

/* -------------------------------------------------------------- 积分 */

/**
 * @brief 求不定积分。
 *
 * 覆盖：幂函数、指数/对数、三角与反三角基本式、线性换元、
 * 乘积分部、部分分式（分母为一次或二次因式）、多项式除法。
 * 无法积分时返回 null，由调用方呈现为未求值积分。
 */
Ptr integrate(const Ptr &expression, std::uint32_t variable);
/** @brief 定积分：先求原函数再代上下限；上下限为无穷时走极限路径。 */
Ptr definiteIntegral(const Ptr &expression, std::uint32_t variable, const Ptr &lower, const Ptr &upper);

/* -------------------------------------------------------------- 极限 */

/**
 * @brief 求 x → point 的极限。
 *
 * 策略依次为：直接代入 → 有理式约分 → 洛必达（0/0 与 ∞/∞）→
 * 主导项分析 → Taylor 展开比较阶。
 */
Ptr limit(const Ptr &expression, std::uint32_t variable, const Ptr &point);
/** @brief x → ±∞ 的极限。 */
Ptr limitAtInfinity(const Ptr &expression, std::uint32_t variable, bool positive);

/* ---------------------------------------------------------- Taylor */

/** @brief 在 point 处展开到 order 阶，返回多项式。 */
Ptr taylorExpand(const Ptr &expression, std::uint32_t variable, const Ptr &point, int order);

/* ------------------------------------------------------------ 辅助 */

/** @brief 把变量代入表达式的值。 */
Ptr substitute(const Ptr &expression, std::uint32_t variable, const Ptr &value);
/** @brief 表达式是否含指定变量的函数依赖（含幂、函数实参等）。 */
bool dependsOn(const NodePtr &expression, std::uint32_t variable);
/** @brief 是否为常数（不含任何符号）。 */
bool isConstant(const NodePtr &expression);

} // namespace Cas::Detail

#endif // CAS_DETAIL_CALCULUS_H

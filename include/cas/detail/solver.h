#ifndef CAS_DETAIL_SOLVER_H
#define CAS_DETAIL_SOLVER_H

/**
 * @file solver.h
 * @brief 方程求解：把 f(x) = 0 的求解分解为若干可判定情形。
 *
 * 支持范围：
 *  - 多项式：一次、二次（含复根时不返回）、三次及以上走有理根 + 降次；
 *  - 可整体换元：`x^4 - 5x^2 + 4 = 0`（把 x^2 当新变量）；
 *  - 简单超越方程：`e^x = a`、`ln x = a`、`sin x = a`（返回主值分支）。
 *
 * 无解析解或超出范围时返回空列表，调用方据此显示"未求解"。
 */

#include "cas/detail/node.h"
#include "cas/detail/simplify.h"

namespace Cas::Detail {

/** @brief 一个方程的解。 */
struct Solution {
    Ptr value;     ///< 解的值（表达式）
    bool verified = false; ///< 是否已代回验证通过
};

/**
 * @brief 解方程 left = right 关于 variable。
 * @param solutions 输出解列表，可能为空（无解或无法求解）。
 * @return 成功进入求解流程返回 true；表达式无法识别为可解类型时返回 false。
 */
bool solve(const NodePtr &left, const NodePtr &right, std::uint32_t variable,
           Memory::Vector<Solution> *solutions);

/** @brief 解不等式（返回区间描述文本，暂只支持线性与二次）。 */
bool solveInequality(const NodePtr &left, const NodePtr &right, const std::uint32_t variable,
                     const char *relation, std::string *description);

/** @brief 解线性方程组（高斯消元，系数必须为数值）。 */
bool solveLinearSystem(const Memory::Vector<Memory::Vector<Ptr>> &equations,
                       const Memory::Vector<std::uint32_t> &variables,
                       Memory::Vector<Ptr> *solution);

} // namespace Cas::Detail

#endif // CAS_DETAIL_SOLVER_H

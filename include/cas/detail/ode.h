// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#ifndef CAS_DETAIL_ODE_H
#define CAS_DETAIL_ODE_H

/**
 * @file ode.h
 * @brief 常微分方程（ODE）的符号求解。
 *
 * 表达式约定：
 * - 未知函数写成 `y`（符号），自变量写成 `x`；先支持一阶 `dy/dx` 与二阶 `d²y/dx²`。
 * - 导数节点写作 `Derivative(y, x)` 与 `Derivative(y, x, 2)`（函数名 `Derivative`）。
 * - 积分常数写作 `C1`、`C2`（符号）。
 *
 * 能力边界与 sympy 的 `classify_ode` 对齐，但只实现"能给出通解"的那几类：
 *
 * 一阶：
 * - 可分离变量      `y' = f(x)·g(y)`
 * - 一阶线性        `y' + P(x)·y = Q(x)`，用积分因子
 * - 齐次型          `y' = F(y/x)`，换元 `u = y/x`
 * - 伯努利          `y' + P(x)·y = Q(x)·y^n`
 * - 恰当方程        `M dx + N dy = 0` 且 `M_y = N_x`
 *
 * 二阶常系数（`a·y'' + b·y' + c·y = f(x)`）：
 * - 齐次：特征方程三个情形（两实根 / 重根 / 共轭复根）
 * - 非齐次：待定系数法（`f` 为多项式、指数、正弦余弦、及它们的乘积）
 *
 * 做不了的一律返回 `Nullable` 并能给出"为什么不行"的诊断，不返回猜测的答案。
 */

#include "cas/detail/node.h"
#include "cas/detail/simplify.h"

#include <string>

namespace Cas::Detail {

/** @brief 构造导数节点 `Derivative(function, variable [, order])`。 */
Ptr makeDerivative(const Ptr &function, std::uint32_t variable, int order = 1);

/**
 * @brief 识别导数节点。
 * @param order 输出阶数（缺省 1）；不需要时传 nullptr。
 * @return 是导数节点时返回被求导的函数，否则返回空。
 */
Ptr asDerivative(const Ptr &node, std::uint32_t *variable, int *order);

/** @brief 是否含未知函数 `y` 或其导数。 */
bool containsUnknown(const Ptr &expression, std::uint32_t unknown);

/** @brief 把 `y'`、`y''` 收集出来，判断最高阶；返回 -1 表示不是 ODE。 */
int derivativeOrder(const Ptr &equation, std::uint32_t unknown);

/**
 * @brief 求解 ODE。
 *
 * @param equation 方程，形如 `y' = x*y`（`=` 已由调用方拆成等号两边的差）。
 * @param unknown  未知函数的符号编号（`y`）。
 * @param variable 自变量符号编号（`x`）。
 * @param description 成功时写入所用方法的说明（如"可分离变量"）。
 * @param message  失败时写入原因。
 * @return 通解（`y = ...` 的右端），失败返回空。
 */
Ptr solveOde(const Ptr &equation, std::uint32_t unknown, std::uint32_t variable,
             std::string *description, std::string *message);

/**
 * @brief 带初值的特解。
 *
 * 初值形如 `y(0) = 1` 与 `y'(0) = 0`（`y` 的初值 / 一阶导数初值）。
 * 目前支持把这些条件代入通解、解出 `C1`/`C2`。
 */
Ptr solveOdeWithInitialValues(const Ptr &equation, std::uint32_t unknown, std::uint32_t variable,
                              const Memory::Vector<Ptr> &conditions, std::string *description,
                              std::string *message);

} // namespace Cas::Detail

#endif // CAS_DETAIL_ODE_H

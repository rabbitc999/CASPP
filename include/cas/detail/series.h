// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#ifndef CAS_DETAIL_SERIES_H
#define CAS_DETAIL_SERIES_H

/**
 * @file series.h
 * @brief 级数求和、无穷级数判别与渐近展开。
 *
 * 求和策略对标 sympy 的 `concrete.summations.eval_sum`，但按本工程的表达力裁剪成
 * 一条固定流水线，每步失配就落到下一步：
 *
 * 1. **平凡情形**：被加项不含求和变量 → `f·(b-a+1)`；上下限相等 → 直接代入。
 * 2. **线性拆分**：和式按 Add 拆开，按 Mul 提出与变量无关的因子。
 * 3. **多项式（Faulhaber）**：`Σ i^k` 用伯努利数闭式，`Σ 1/i^r` 用广义谐波数。
 * 4. **几何级数**：`Σ q^i` 闭式，含 `|q|<1` 的无穷和。
 * 5. **裂项（telescoping）**：`Σ (a_i - a_{i-1})` 只留两端。
 * 6. **p 级数**：`Σ 1/i^p` 在 `p>1` 时收敛，给出 ζ(p) 的精确值（p 为偶数）。
 * 7. **直接展开**：上限为具体整数且项数不多时逐项相加。
 *
 * 全部失败时返回 null，由调用方呈现为未求值的 `sum(...)`，绝不猜答案。
 */

#include "cas/detail/node.h"
#include "cas/detail/simplify.h"

namespace Cas::Detail {

/** @brief 求和上下限；`upper` 为 Infinity 表示无穷和。 */
struct SumLimits {
    std::uint32_t variable = 0;
    Ptr lower;
    Ptr upper;
};

/**
 * @brief 符号求和 `Σ_{i=lower}^{upper} term`。
 * @return 闭式结果；无法求和时返回 `Ptr()`。
 */
Ptr summation(const Ptr &term, const SumLimits &limits);

/**
 * @brief 无穷级数敛散性判别。
 * @param verdict 输出判定结果（Convergent / Divergent / Unknown）
 * @param value 收敛时输出和（可能为 null，表示只知道收敛）
 * @return 是否做出了判定。
 */
enum class Convergence { Convergent, Divergent, Unknown };
bool seriesConvergence(const Ptr &term, const SumLimits &limits, Convergence *verdict, Ptr *value);

/**
 * @brief 在 point 处做渐近/幂级数展开，直到出现 order 阶以上的项。
 *
 * 与 `taylorExpand` 的区别：内部按"逐项求导 + 取极限"生成系数，因此能处理
 * `exp(-1/x)` 这类在 0 处没有 Taylor 级数的函数（此时给出 0 阶常数项）。
 */
Ptr seriesExpand(const Ptr &expression, std::uint32_t variable, const Ptr &point, int order);

/** @brief 伯努利数 B_n（约定 B_1 = -1/2）。 */
Ptr bernoulliNumber(long long n);

/**
 * @brief 识别简单序列的通项（前若干项的差分/比值模式）。
 *
 * 仅用于 `guess` 命令：给出前几项，尝试识别等差、等比、多项式与阶乘型序列。
 * @return 通项表达式（以 index 为变量）；无法识别返回 `Ptr()`。
 */
Ptr guessSequence(const Memory::Vector<Ptr> &terms, std::uint32_t index);

} // namespace Cas::Detail

#endif // CAS_DETAIL_SERIES_H

// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#ifndef CAS_DETAIL_TRIGONOMETRY_H
#define CAS_DETAIL_TRIGONOMETRY_H

/**
 * @file trigonometry.h
 * @brief 三角与超越函数的化简与恒等式重写。
 *
 * 对标 SymPy 的 `trigsimp` / `expand_trig` / `fu` 三条思路，本模块只做**有向**重写：
 * 每个规则都把表达式推向"更短或更规范"的方向，因此不需要搜索，也不会来回打转。
 *
 * 规则集：
 *  1. 奇偶性：`sin(-x) → -sin(x)`、`cos(-x) → cos(x)`、`tan(-x) → -tan(x)`；
 *  2. 勾股恒等式：`sin(x)^2 + cos(x)^2 → 1`（含其倍数与移位形式）；
 *  3. 特殊角：`sin(0)=0`、`sin(pi/6)=1/2`、`sin(pi/4)=sqrt(2)/2`、`cos(pi/3)=1/2` …；
 *  4. 反函数：`sin(asin(x)) → x`、`asin(sin(x)) → x`（主值区间内成立时）；
 *  5. 对数/指数：`exp(0)=1`、`log(1)=0`、`log(e)=1`、`log(a^b) → b·log(a)`、
 *     `exp(a)·exp(b) → exp(a+b)`、`a^x · a^y → a^(x+y)`（已在 Mul 的指数合并里处理）。
 */

#include "cas/detail/node.h"

namespace Cas::Detail {

/** @brief 三角化简主入口：反复应用规则直到不动点（最多 8 轮，防振荡）。 */
Ptr trigSimplify(const Ptr &expression);

/**
 * @brief 展开三角加法公式：`sin(a+b) → sin(a)cos(b) + cos(a)sin(b)`。
 *
 * 仅对二倍角与和角展开，避免对 n 项和做全展开导致的项数膨胀。
 */
Ptr expandTrig(const Ptr &expression);

/**
 * @brief 把 `sin(x)^2` 型项按 `1 - cos(2x))/2` 降幂，用于积分与线性化。
 */
Ptr reduceTrigPower(const Ptr &expression);

} // namespace Cas::Detail

#endif // CAS_DETAIL_TRIGONOMETRY_H

// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#ifndef CAS_DETAIL_SIMPLIFY_H
#define CAS_DETAIL_SIMPLIFY_H

/**
 * @file simplify.h
 * @brief 表达式自动化简：Add/Mul/Pow 的规范化、同类项合并、数值折叠、恒等式消去。
 *
 * 设计原则（对标 SymPy 的 `Add.flatten` / `Mul.flatten` / `Pow.__new__`）：
 *  - **扁平化**：`(a+b)+c` 与 `a+(b+c)` 都变成三元 Add，括号差异不影响结果；
 *  - **规范排序**：操作数按 `Node::compareTo` 排序，结构相等即可用 `equals` 判定；
 *  - **同类项合并**：Add 收集系数，Mul 收集指数，`2x + 3x → 5x`、`x·x → x²`；
 *  - **数值优先**：所有纯数值运算在化简期完成，结果保持精确有理数；
 *  - **不引入新的未定义行为**：`0^0` 保持未求值，不擅自约定为 1。
 */

#include "cas/detail/node.h"

namespace Cas::Detail {

/** @brief 创建并化简 n 元加法。 */
Ptr makeAdd(Ptr left, Ptr right);
/** @brief 创建并化简 n 元加法（多操作数）。 */
Ptr makeAddMany(const Memory::Vector<NodePtr> &terms);
/** @brief 创建并化简 n 元乘法。 */
Ptr makeMul(Ptr left, Ptr right);
/** @brief 创建并化简 n 元乘法（多操作数）。 */
Ptr makeMulMany(const Memory::Vector<NodePtr> &factors);
/** @brief 创建并化简幂 base^exponent。 */
Ptr makePow(Ptr base, Ptr exponent);
/** @brief 创建并化简相反数 -value。 */
Ptr makeNegate(Ptr value);
/** @brief 创建并化简除法 left / right。 */
Ptr makeDivide(Ptr numerator, Ptr denominator);

/** @brief 数值：把两个 Node 转成有理数（要求 isNumber）。 */
bool numberValue(const NodePtr &node, BigInt *numerator, BigInt *denominator);
/** @brief 由分子分母构造约分后的数值节点。 */
Ptr makeNumber(const BigInt &numerator, const BigInt &denominator);

/** @brief 递归化简（幂等）：先化简子式，再按规则折叠。 */
Ptr simplify(const Ptr &expression);

/** @brief 表达式里出现的全部符号（去重，按编号排序）。 */
Memory::Vector<std::uint32_t> freeSymbols(const NodePtr &expression);
/** @brief 表达式是否含有指定符号。 */
bool containsSymbol(const NodePtr &expression, std::uint32_t symbolId);
/** @brief 表达式是否含有 Infinity。 */
bool containsInfinity(const NodePtr &expression);

} // namespace Cas::Detail

#endif // CAS_DETAIL_SIMPLIFY_H

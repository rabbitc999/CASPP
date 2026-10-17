// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#ifndef CAS_DETAIL_PARSER_H
#define CAS_DETAIL_PARSER_H

/**
 * @file parser.h
 * @brief 表达式解析器：把文本转成表达式树。
 *
 * 手写递归下降，优先级自低向高：
 * ```
 *   expr    := assign
 *   assign  := add (('=' | '==') add)?        // 等式，产出 Equal 函数节点
 *   add     := mul (('+'|'-') mul)*
 *   mul     := unary (('*'|'/'|implicit) unary)*
 *   unary   := ('-'|'+') unary | power
 *   power   := postfix ('^' unary)?           // 右结合，且 -2^2 = -(2^2)
 *   postfix := atom ('!' | "'")*
 *   atom    := number | symbol | '(' expr ')' | '|' expr '|' | func '(' args ')'
 * ```
 *
 * 隐式乘法（`2x`、`3(x+1)`、`(x+1)(x-1)`、`2sin(x)`）在 `mul` 层识别：
 * 只要下一个 token 能开启一个 atom，就当作乘法，这样手写输入不需要到处加 `*`。
 *
 * 语法糖：
 *  - `x^2` 与 `x**2` 等价，`^` 右结合；
 *  - `sqrt(x)`、`abs(x)`、`cbrt(x)` 是函数；`√x` 也可用；
 *  - `pi`、`e`、`oo`、`infinity` 是常数；`E` 视为自然对数底；
 *  - 赋值 `f := x^2+1` 只记录表达式，不做符号替换（v1 不实现函数定义）。
 */

#include "cas/detail/simplify.h"

#include <functional>

namespace Cas::Detail {

/** @brief 解析结果与诊断信息。 */
struct ParseResult {
    Ptr expression;               ///< 成功时的表达式
    bool ok = false;              ///< 是否成功
    std::string message;          ///< 失败原因（含位置）
    std::size_t position = 0;     ///< 出错位置（从 0 起的字符下标）
    std::size_t length = 0;       ///< 表达式消耗的字符数
};

/**
 * @brief 解析选项。
 *
 * 缺省行为是「只有内置名字后跟括号才算函数调用，其余按隐式乘法」，
 * 这样 `f(x)` 仍读作 `f` 乘 `x`。接入用户定义后需要用 `isUserFunction`
 * 把已登记的名字标出来，它们在 `name(...)` 位置才按调用解析。
 */
struct ParseOptions {
    /** @brief 判断某个名字是否为用户定义的函数（返回 true 时 `name(...)` 按调用解析）。 */
    std::function<bool(const std::string &)> isUserFunction;
    /** @brief 是否把 `define`/`let` 的第一个参数整体当作函数签名解析。 */
    bool allowSignature = false;
    /**
     * @brief 允许任意裸符号后跟括号按调用解析。
     *
     * `define` 解析签名时用：此时还不存在定义，无法靠 `isUserFunction` 判断。
     */
    bool allowAnyCall = false;
};

/** @brief 解析单个表达式。 */
ParseResult parseExpression(const std::string &text);
/** @brief 带选项解析单个表达式。 */
ParseResult parseExpression(const std::string &text, const ParseOptions &options);

/**
 * @brief 解析一条初值条件，如 `y(0) = 1`、`y'(0) = 0`。
 *
 * 与 `parseExpression` 的区别只在于：裸符号后跟括号按**函数调用**处理，
 * 这样 `y(0)` 才不会被读成 `y` 乘 `0`。
 */
ParseResult parseInitialConditionText(const std::string &text);

/** @brief 把常数的名字登记到符号表（`pi`、`e`、`infinity` 等按常数处理，其余为符号）。 */
bool isReservedConstantName(const std::string &name);

} // namespace Cas::Detail

#endif // CAS_DETAIL_PARSER_H

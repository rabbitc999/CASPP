#ifndef CAS_DETAIL_PRINTER_H
#define CAS_DETAIL_PRINTER_H

/**
 * @file printer.h
 * @brief 表达式输出：纯文本（前缀/中缀）与 LaTeX。
 *
 * 排版规则参考常规数学书写习惯：
 *  - 乘法在纯文本中输出 `2*x`，在 LaTeX 中输出 `2 x`（省略乘号）；
 *  - 负数加法输出 `x^2 - 2x + 1` 而不是 `x^2 + (-2x) + 1`；
 *  - 分式在 LaTeX 中输出 `\frac{p}{q}`，纯文本输出 `p/q`；
 *  - 幂输出 `x^{2}`，三角函数输出 `\sin`、`\cos` 等。
 */

#include "cas/detail/node.h"

#include <string>

namespace Cas::Detail {

/** @brief 输出选项。 */
struct PrintOptions {
    bool useUnicode = false; ///< 是否用 `π`、`√`、`·`、`≠` 等 Unicode 字符
};

/** @brief 输出为纯文本（适合命令行 REPL）。 */
std::string toText(const NodePtr &expression, const PrintOptions &options = PrintOptions());

/** @brief 输出为 LaTeX 公式。 */
std::string toLatex(const NodePtr &expression, const PrintOptions &options = PrintOptions());

/** @brief 把有理数格式化为文本：整数直接输出，分数输出 `p/q`。 */
std::string numberToText(const NodePtr &number);

} // namespace Cas::Detail

#endif // CAS_DETAIL_PRINTER_H

// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#ifndef CAS_ENGINE_H
#define CAS_ENGINE_H

/**
 * @file engine.h
 * @brief 对外门面：`Cas::Engine`。
 *
 * 命名沿用参考工程 cas 的 Qt 风格驼峰约定（`Cas::Engine`、`calculate()`、`memoryStats()`）。
 * 内部实现全部在 `Cas::Detail` 命名空间里，外部只需要这一个头文件。
 *
 * 典型用法：
 * @code
 * Cas::Engine engine;
 * const Cas::Result r = engine.calculate("diff(x^3 + 2x, x)");
 * if (r.ok) std::cout << r.text << "\n";   // 3*x^2 + 2
 * @endcode
 *
 * 线程安全：`Engine` 不做内部加锁。内存池是全局共享的，因此**同一时刻只应有一个线程**
 * 通过 Engine 修改表达式；只读输出（`toText`/`toLatex`）在多线程下是安全的。
 */

#include "cas/allocator.h"
#include "cas/detail/parser.h"
#include "cas/memory.h"

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace Cas {

namespace Detail {
class Node;
class Ptr;
using NodePtr = Node *;
} // namespace Detail

/** @brief 一次计算的返回值。 */
struct Result {
    bool ok = false;        ///< 是否成功
    std::string text;       ///< 纯文本结果
    std::string latex;      ///< LaTeX 结果
    std::string message;    ///< 失败原因
    std::size_t position = 0; ///< 语法错误位置（0 起）
    double elapsedMs = 0.0; ///< 耗时（毫秒）
};

/** @brief 代数符号计算引擎。 */
class Engine {
public:
    Engine();
    ~Engine();
    Engine(const Engine &) = delete;
    Engine &operator=(const Engine &) = delete;

    /* ---------------------------------------------------- 求值入口 */

    /** @brief 按命令语法求值：`diff(f, x)`、`integrate(f, x)`、`solve(f = 0, x)` … 或纯表达式。 */
    Result evaluate(const std::string &input);
    /** @brief `evaluate` 的别名，便于沿用参考工程的调用习惯。 */
    Result calculate(const std::string &input) { return evaluate(input); }

    /* ------------------------------------------- 逐项 API（无需写命令） */

    Result simplify(const std::string &expression);
    Result expand(const std::string &expression);
    Result factor(const std::string &expression);
    Result together(const std::string &expression);
    Result trigSimplify(const std::string &expression);
    Result differentiate(const std::string &expression, const std::string &variable = "x");
    Result integrate(const std::string &expression, const std::string &variable = "x");
    /** @brief 定积分；上下限可为 `0`、`pi`、`oo` 等表达式。 */
    Result integrateDefinite(const std::string &expression, const std::string &lower, const std::string &upper,
                             const std::string &variable = "x");
    /**
     * @brief 变量替换：`subs(x^2+1, x, 3)`。
     *
     * @param expression 待替换的表达式。
     * @param variable   被替换的变量名（也可传 `[x, y]` 这样的列表，配合 `replacements` 列表）。
     * @param value      代入值（或与变量列表等长的值列表）。
     * @return 替换并化简后的结果。
     */
    Result substitute(const std::string &expression, const std::string &variable, const std::string &value);

    /**
     * @brief 定义函数或常量。
     *
     * - `define(f(x), x^2+1)` 定义单参函数，之后 `f(2)` 求值得 `5`；
     * - `define(f(x, y), x+y)` 定义多参函数；
     * - `define(a, 3)` 定义符号常量，之后 `a+1` 得 `4`；
     * - 同名重复定义会覆盖，`undef(f)` 删除。
     *
     * @param signature 函数签名（`f(x)`、`f(x, y)`）或常量名（`a`）。
     * @param body      函数体 / 常量值。
     */
    Result define(const std::string &signature, const std::string &body);
    /** @brief 删除定义；返回是否确实存在。 */
    Result undefine(const std::string &name);
    /** @brief 列出全部定义（`defs` 命令）。 */
    Result definitions() const;
    Result limit(const std::string &expression, const std::string &point, const std::string &variable = "x");
    Result taylor(const std::string &expression, const std::string &point, int order,
                  const std::string &variable = "x");
    /**
     * @brief 求驻点并判定极值。
     *
     * 解 `f'(x) = 0` 得驻点，再用 `f''` 的符号区分极小、极大与拐点。
     * 输出形如 `x = 1（极小值，f = -1）; x = -1（极大值，f = 3）`。
     */
    Result criticalPoints(const std::string &expression, const std::string &variable = "x");
    /** @brief 解方程；`equation` 可写成 `x^2-1=0` 或 `x^2-1`。 */
    Result solve(const std::string &equation, const std::string &variable = "x");
    /** @brief 解方程组：`{"x+y=3", "x-y=1"}`。 */
    Result solveSystem(const std::vector<std::string> &equations, const std::vector<std::string> &variables);
    /** @brief 解不等式：`x^2-1>0`。 */
    Result solveInequality(const std::string &inequality, const std::string &variable = "x");
    Result gcd(const std::string &left, const std::string &right, const std::string &variable = "x");
    Result resultant(const std::string &left, const std::string &right, const std::string &variable = "x");
    Result groebner(const std::vector<std::string> &polynomials, const std::vector<std::string> &variables);
    /**
     * @brief 解常微分方程。
     *
     * @param equation   方程，如 `y' = x*y`、`y'' + y = 0`（`=` 可省略，缺省视为 = 0）。
     * @param unknown    未知函数名，缺省 `y`。
     * @param variable   自变量名，缺省 `x`。
     * @param conditions 初值条件，如 `y(0) = 1`、`y'(0) = 0`；为空时给通解。
     *
     * 初值条件只对能给出显式解的类型（一阶线性、二阶常系数）生效；隐式解会明确提示。
     */
    Result dsolve(const std::string &equation, const std::string &unknown = "y",
                  const std::string &variable = "x",
                  const std::vector<std::string> &conditions = std::vector<std::string>());

    /* ------------------------------------------------------ 工具 */

    /** @brief 文本转 LaTeX（输入为任意表达式）。 */
    Result toLatexResult(const std::string &expression);
    /** @brief 当前内存池统计。 */
    MemoryStats memoryStats() const;
    /** @brief 把池内所有可回收块归还（表达式都已析构时调用）。 */
    void releaseMemory();
    /** @brief 最近一次成功的表达式（供 REPL 的 `%` 引用）。 */
    std::string lastResult() const { return m_lastResult; }

    /** @brief 输出用选项：true 时使用 Unicode 字符（π、√、·）。 */
    void setUnicodeOutput(bool enabled) { m_unicodeOutput = enabled; }
    bool unicodeOutput() const { return m_unicodeOutput; }

    /** @brief 帮助文本（REPL 的 `help`）。 */
    static std::string helpText();

private:
    /** @brief 按命令名分发到具体方法；不认识返回 message 为空的 Result。 */
    Result dispatch(const std::string &name, const Detail::Node *call, const std::string &source);

    /**
     * @brief 把表达式里内嵌的矩阵命令就地求值，输出改写后的文本。
     *
     * `det(inv(A))`、`matmul(A, B)` 这类写法里，内层命令在简化流程中只是普通的
     * 未定义函数。这里自底向上（后序遍历）把它们替换成结果文本，再由调用方重新解析，
     * 从而不必在引擎内部手工维护表达式树的引用计数。
     * @return 是否有任意一层被替换。
     */
    bool resolveCommands(const Detail::Ptr &node, std::string *out);

    /* 节点直传版本：避免先把参数渲染成文本再重新解析，从而丢失结构 */
    Result derivativeOfNode(const Detail::Ptr &expression, const std::string &variable);
    Result limitOfNode(const Detail::Ptr &expression, const Detail::Ptr &point, const std::string &variable);
    Result taylorOfNode(const Detail::Ptr &expression, const Detail::Ptr &point, int order,
                        const std::string &variable);

    /**
     * @brief `sum(term, i=a..b)` 的实现。
     *
     * `term` 是已经解析好的表达式（未求值，保留对求和变量的依赖），
     * `range` 是 `i=a..b` 这样的区间写法，`variable` 是可选的显式变量名。
     */
    Result sumCommand(const Detail::Ptr &term, const Detail::Ptr &range, const Detail::Ptr &variable);
    /** @brief `converge(term, i=a..oo)` 的实现：判定级数敛散性。 */
    Result convergeCommand(const Detail::Ptr &term, const Detail::Ptr &range);

    /** @brief 用户定义：签名 `f(x, y)` 或常量 `a` → 定义体文本。 */
    struct Definition {
        std::string name;                  ///< 函数名 / 常量名
        std::vector<std::string> parameters; ///< 形参；常量定义时为空
        std::string body;                  ///< 定义体文本
    };

    /** @brief 查表并展开定义：把 `f(2)` 替换成 `2^2+1` 再化简。 */
    Result expandDefinition(const std::string &name, const Detail::Ptr &call);

    /**
     * @brief 把表达式树里的用户定义调用就地替换成展开后的文本。
     *
     * 与 `resolveCommands` 同一思路：自底向上替换成文本，再由调用方重新解析，
     * 从而让化简/求导/积分等流程完全不需要感知定义表。
     * @return 是否有任意一层被替换。
     */
    bool foldDefinitions(const Detail::Ptr &node, std::string *out);

    bool m_unicodeOutput = false;
    std::string m_lastResult;
    /** @brief 用户定义表（键为名字，按插入顺序输出）。 */
    std::vector<Definition> m_definitions;
    /** @brief 定义折叠的递归深度（防止自引用定义无限展开）。 */
    int m_foldDepth = 0;

    /** @brief 按当前定义表构造解析选项（让 `f(2)` 走函数调用而非乘法）。 */
    Detail::ParseOptions parseOptions() const;
    /** @brief 使用当前定义表解析表达式。 */
    Detail::ParseResult parse(const std::string &text) const;
};

} // namespace Cas

#endif // CAS_ENGINE_H

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

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace Cas {

namespace Detail {
class Node;
class Ptr;
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
    Result limit(const std::string &expression, const std::string &point, const std::string &variable = "x");
    Result taylor(const std::string &expression, const std::string &point, int order,
                  const std::string &variable = "x");
    /** @brief 解方程；`equation` 可写成 `x^2-1=0` 或 `x^2-1`。 */
    Result solve(const std::string &equation, const std::string &variable = "x");
    /** @brief 解方程组：`{"x+y=3", "x-y=1"}`。 */
    Result solveSystem(const std::vector<std::string> &equations, const std::vector<std::string> &variables);
    /** @brief 解不等式：`x^2-1>0`。 */
    Result solveInequality(const std::string &inequality, const std::string &variable = "x");
    Result gcd(const std::string &left, const std::string &right, const std::string &variable = "x");
    Result resultant(const std::string &left, const std::string &right, const std::string &variable = "x");
    Result groebner(const std::vector<std::string> &polynomials, const std::vector<std::string> &variables);

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

    /* 节点直传版本：避免先把参数渲染成文本再重新解析，从而丢失结构 */
    Result derivativeOfNode(const Detail::Ptr &expression, const std::string &variable);
    Result limitOfNode(const Detail::Ptr &expression, const Detail::Ptr &point, const std::string &variable);
    Result taylorOfNode(const Detail::Ptr &expression, const Detail::Ptr &point, int order,
                        const std::string &variable);

    bool m_unicodeOutput = false;
    std::string m_lastResult;
};

} // namespace Cas

#endif // CAS_ENGINE_H

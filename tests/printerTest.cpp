// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

// 打印器自检：内部节点（equal / list）不能泄露到用户界面，
// 且 solveSystem 这类模块入口不能被「未知函数」诊断误伤。
#include "cas/allocator.h"
#include "cas/engine.h"

#include <cstdio>
#include <string>

namespace {

int failures = 0;

/** @brief 断言一次求值的文本结果。 */
void expectText(Cas::Engine *engine, const char *input, const char *expected) {
    const Cas::Result result = engine->evaluate(input);
    if (!result.ok) {
        std::printf("FAIL %-42s ← 错误：%s\n", input, result.message.c_str());
        ++failures;
        return;
    }
    if (result.text != expected) {
        std::printf("FAIL %-42s ← 得到 %s，期望 %s\n", input, result.text.c_str(), expected);
        ++failures;
        return;
    }
    std::printf("  %-40s = %s\n", input, result.text.c_str());
}

/** @brief 断言一次求值的 LaTeX 结果。 */
void expectLatex(Cas::Engine *engine, const char *input, const char *expected) {
    const Cas::Result result = engine->evaluate(input);
    if (!result.ok) {
        std::printf("FAIL %-42s ← 错误：%s\n", input, result.message.c_str());
        ++failures;
        return;
    }
    if (result.latex != expected) {
        std::printf("FAIL %-42s ← LaTeX 得到 %s，期望 %s\n", input, result.latex.c_str(), expected);
        ++failures;
        return;
    }
    std::printf("  %-40s ⇒ %s\n", input, result.latex.c_str());
}

/** @brief 断言一次求值应当失败。 */
void expectError(Cas::Engine *engine, const char *input) {
    const Cas::Result result = engine->evaluate(input);
    if (result.ok) {
        std::printf("FAIL %-42s ← 本应报错，却得到 %s\n", input, result.text.c_str());
        ++failures;
        return;
    }
    std::printf("  %-40s ↯ %s\n", input, result.message.c_str());
}

/* --------------------------------------------------- 等式不泄露 equal() */

void testEqualPrinting() {
    Cas::Engine engine;

    expectText(&engine, "y=x", "y = x");
    expectText(&engine, "x=y", "x = y");
    expectText(&engine, "y=x+1", "y = x + 1");
    expectText(&engine, "y=x^2", "y = x^2");
    expectText(&engine, "y=2*x+3", "y = 2*x + 3");

    expectLatex(&engine, "y=x", "y = x");
    expectLatex(&engine, "x=y", "x = y");
    expectLatex(&engine, "y=2*x+3", "y = 2 x + 3");

    // latex(...) 命令走的是同一套 LaTeX 输出
    expectText(&engine, "latex(y=x)", "y = x");
    expectText(&engine, "latex(x=y^2)", "x = y^2");
}

/* ------------------------------------------------- 列表不泄露 list(...) */

void testListPrinting() {
    Cas::Engine engine;

    expectText(&engine, "[1,2,3]", "[1, 2, 3]");
    expectText(&engine, "list(1,2,3)", "[1, 2, 3]");
    expectText(&engine, "[x, y]", "[x, y]");
    expectText(&engine, "[]", "[]");
    expectText(&engine, "[x+y=3, x-y=1]", "[x + y = 3, x - y = 1]");

    expectLatex(&engine, "[1,2,3]", "\\left[1, 2, 3\\right]");
    expectLatex(&engine, "[x+y=3]", "\\left[x + y = 3\\right]");

    // 嵌套列表
    expectText(&engine, "[[1,2],[3,4]]", "[[1, 2], [3, 4]]");
}

/* ------------------------------------- 模块入口不被「未知函数」诊断误伤 */

void testModuleEntries() {
    Cas::Engine engine;

    expectText(&engine, "solveSystem([x+y=3, x-y=1], [x, y])", "x = 2, y = 1");
    expectText(&engine, "solve([x+y=3, x-y=1], [x, y])", "x = 2, y = 1");
    expectText(&engine, "groebner([x^2+y^2-1, x-y], [x, y])", "{x - y, y^2 - 1/2}");

    // 大小写变体与别名
    expectText(&engine, "solvesystem([x+y=3, x-y=1], [x, y])", "x = 2, y = 1");

    // 参数不足仍是清晰报错，而不是被当成乘法
    expectError(&engine, "solveSystem([x+y=3])");
}

/* -------------------------------- 未知函数诊断本身仍然生效（没有被放松） */

void testUnknownCallStillDiagnosed() {
    Cas::Engine engine;

    // 含逗号的括号 → 用户显然想写调用，但既不是内置函数，也没有定义
    expectError(&engine, "g(x, t)");

    // 定义之后就能调用
    expectText(&engine, "define(g(x, y), x*y+1)", "已定义 g(x, y) = x*y + 1");
    expectText(&engine, "g(2, 3)", "7");

    // 参数个数不匹配
    expectError(&engine, "g(1)");
    expectError(&engine, "g(1, 2, 3)");
}

/* ---------------------------------------- 其他命令的输出没有受到连带影响 */

void testOtherPrintsUnchanged() {
    Cas::Engine engine;

    expectText(&engine, "diff(x^2, x)", "2*x");
    expectText(&engine, "factor(x^2-1)", "(x - 1)*(x + 1)");
    expectText(&engine, "det([[1,2],[3,4]])", "-2");
    expectText(&engine, "sin(0)", "0");
    expectLatex(&engine, "sin(x)", "\\sin\\left(x\\right)");
}

} // namespace

int main() {
    std::printf("打印器与内部表示自检\n");
    Cas::Allocator::initialize();

    testEqualPrinting();            std::printf("equal printing ok\n");
    testListPrinting();             std::printf("list printing ok\n");
    testModuleEntries();            std::printf("module entries ok\n");
    testUnknownCallStillDiagnosed(); std::printf("unknown call diagnosis ok\n");
    testOtherPrintsUnchanged();     std::printf("other prints ok\n");

    if (failures == 0) {
        std::printf("全部通过\n");
        return 0;
    }
    std::printf("%d 项失败\n", failures);
    return 1;
}

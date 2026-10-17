// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

// 变量替换与用户定义的自检：subs 单/多变量、define 函数与常量、嵌套定义、错误处理。
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
        std::printf("FAIL %-40s ← 错误：%s\n", input, result.message.c_str());
        ++failures;
        return;
    }
    if (result.text != expected) {
        std::printf("FAIL %-40s ← 得到 %s，期望 %s\n", input, result.text.c_str(), expected);
        ++failures;
        return;
    }
    std::printf("  %-38s = %s\n", input, result.text.c_str());
}

/** @brief 断言一次求值应当失败。 */
void expectError(Cas::Engine *engine, const char *input) {
    const Cas::Result result = engine->evaluate(input);
    if (result.ok) {
        std::printf("FAIL %-40s ← 本应报错，却得到 %s\n", input, result.text.c_str());
        ++failures;
        return;
    }
    std::printf("  %-38s ↯ %s\n", input, result.message.c_str());
}

/* ------------------------------------------------------------ 变量替换 */

void testSubstitute() {
    Cas::Engine engine;

    expectText(&engine, "subs(x^2+1, x, 3)", "10");
    expectText(&engine, "subs(x^2+y^2, x, 1)", "y^2 + 1");
    expectText(&engine, "subs(x^2+y^2, x, a)", "a^2 + y^2");
    expectText(&engine, "subs(sin(x), x, 0)", "0");
    expectText(&engine, "subs(x^2+1, y, 3)", "x^2 + 1"); // 不含 y：原样返回
    expectText(&engine, "subs(x*y, [x, y], [2, 3])", "6");
    expectText(&engine, "subs(x^2+y^2, [x, y], [3, 4])", "25");
    // 代入值本身可以是表达式
    expectText(&engine, "subs(x^2, x, y+1)", "(y + 1)^2");

    // 错误路径
    expectError(&engine, "subs(x^2, 3, 1)");                        // 目标不是变量
    expectError(&engine, "subs(x^2, [x, y], [1])");                 // 长度不一致
    expectError(&engine, "subs(x^2, [x, y], 1)");                   // 变量是列表、值不是
    expectError(&engine, "subs(x^2, x)");                           // 参数不足
}

/* ------------------------------------------------------------ 函数定义 */

void testDefineFunction() {
    Cas::Engine engine;

    expectText(&engine, "define(f(x), x^2+1)", "已定义 f(x) = x^2 + 1");
    expectText(&engine, "f(2)", "5");
    expectText(&engine, "f(0)", "1");
    expectText(&engine, "f(x)", "x^2 + 1");     // 形参代入自身：还原定义体
    expectText(&engine, "f(a)", "a^2 + 1");     // 参数可以是符号
    expectText(&engine, "f(y+1)", "(y + 1)^2 + 1");

    // 多参
    expectText(&engine, "define(g(x, y), x^2+y^2)", "已定义 g(x, y) = x^2 + y^2");
    expectText(&engine, "g(3, 4)", "25");
    expectText(&engine, "g(x, 1)", "x^2 + 1");

    // 与命令组合：定义先折叠成普通表达式，再走原有流程
    expectText(&engine, "diff(f(x), x)", "2*x");
    expectText(&engine, "integrate(f(x), x)", "1/3*x^3 + x");
    expectText(&engine, "limit(f(x), x, 1)", "2");
    expectText(&engine, "f(2)+g(3, 4)", "30");
    expectText(&engine, "expand(f(x))", "x^2 + 1");
    expectText(&engine, "f(1)*g(1, 2)", "10");

    // 嵌套定义：k 的定义体调用 h
    expectText(&engine, "define(h(x), x^2)", "已定义 h(x) = x^2");
    expectText(&engine, "define(k(x), h(x)+1)", "已定义 k(x) = h(x) + 1");
    expectText(&engine, "k(2)", "5");

    // 覆盖定义
    expectText(&engine, "define(f(x), x^3)", "已定义 f(x) = x^3");
    expectText(&engine, "f(2)", "8");
}

/* ------------------------------------------------------------ 常量定义 */

void testDefineConstant() {
    Cas::Engine engine;

    expectText(&engine, "define(a, 3)", "已定义 a = 3");
    expectText(&engine, "a+1", "4");
    expectText(&engine, "a^2", "9");
    expectText(&engine, "2*a", "6");

    // 常量可以引用另一个常量
    expectText(&engine, "define(b, a*2)", "已定义 b = 2*a");
    expectText(&engine, "b+1", "7");

    // 未定义的名字仍是自由符号
    expectText(&engine, "zzz+1", "zzz + 1");
}

/* ------------------------------------------------------------ 定义管理 */

void testDefinitions() {
    Cas::Engine engine;

    expectText(&engine, "defs", "（没有定义）");
    expectText(&engine, "define(f(x), x+1)", "已定义 f(x) = x + 1");
    expectText(&engine, "define(a, 2)", "已定义 a = 2");
    expectText(&engine, "defs", "f(x) = x + 1\na = 2");

    expectText(&engine, "undef(f)", "已删除定义 f");
    // 删除后 `f(1)` 退回隐式乘法 `f*1 = f`（定义表里没有 f，不再当调用）
    expectText(&engine, "f(1)", "f");
    expectText(&engine, "a+1", "3");   // 常量仍在

    expectText(&engine, "undef(a)", "已删除定义 a");
    expectError(&engine, "undef(a)");   // 重复删除
}

/* ------------------------------------------------------------ 错误处理 */

void testDefineErrors() {
    Cas::Engine engine;

    expectError(&engine, "define(sin, 1)");          // 内置函数名
    expectError(&engine, "define(diff, 1)");         // 内置命令名
    expectError(&engine, "define(f(x, x), x)");      // 形参重名
    expectError(&engine, "define(f(x), )");          // 定义体为空
    expectError(&engine, "define(c, c+1)");          // 自引用
    expectError(&engine, "define(f(x), x");
    expectError(&engine, "define(f(x))");            // 参数不足

    // 参数个数不匹配
    expectText(&engine, "define(g(x, y), x+y)", "已定义 g(x, y) = x + y");
    expectError(&engine, "g(1)");
    expectError(&engine, "g(1, 2, 3)");
}

/* ---------------------------------------------- 既有隐式乘法语义不受影响 */

void testNoRegression() {
    Cas::Engine engine;

    expectText(&engine, "f*x", "f*x");       // 未定义时仍是乘法
    expectText(&engine, "2x", "2*x");
    expectText(&engine, "a*b", "a*b");
    expectText(&engine, "x(1)", "x");
    expectText(&engine, "sin(x)", "sin(x)");

    // 定义 f 之后，`f*x` 仍是乘法（只有 `f(...)` 才是调用）
    expectText(&engine, "define(f(x), x^2+1)", "已定义 f(x) = x^2 + 1");
    expectText(&engine, "f*x", "f*x");
    expectText(&engine, "f(2)", "5");
}

} // namespace

int main() {
    std::printf("变量替换与用户定义自检\n");
    Cas::Allocator::initialize();

    testSubstitute();    std::printf("subs ok\n");
    testDefineFunction(); std::printf("define function ok\n");
    testDefineConstant(); std::printf("define constant ok\n");
    testDefinitions();   std::printf("definitions ok\n");
    testDefineErrors();  std::printf("define errors ok\n");
    testNoRegression();  std::printf("no regression ok\n");

    if (failures == 0) {
        std::printf("全部通过\n");
        return 0;
    }
    std::printf("%d 项失败\n", failures);
    return 1;
}

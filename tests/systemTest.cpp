// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

// 方程组求解自检：线性走消元、非线性走三角消元（Groebner 基 + 代回剪枝），
// 以及欠定方程组的参数解。
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
        std::printf("FAIL %-52s ← 错误：%s\n", input, result.message.c_str());
        ++failures;
        return;
    }
    if (result.text != expected) {
        std::printf("FAIL %-52s ← 得到 %s，期望 %s\n", input, result.text.c_str(), expected);
        ++failures;
        return;
    }
    std::printf("  %-50s = %s\n", input, result.text.c_str());
}

/**
 * @brief 断言结果里含有某个子串（用于解的顺序不稳定、但集合确定的场合）。
 */
void expectContains(Cas::Engine *engine, const char *input, const char *needle) {
    const Cas::Result result = engine->evaluate(input);
    if (!result.ok) {
        std::printf("FAIL %-52s ← 错误：%s\n", input, result.message.c_str());
        ++failures;
        return;
    }
    if (result.text.find(needle) == std::string::npos) {
        std::printf("FAIL %-52s ← 得到 %s，其中没有 %s\n", input, result.text.c_str(), needle);
        ++failures;
        return;
    }
    std::printf("  %-50s = %s\n", input, result.text.c_str());
}

/** @brief 断言一次求值应当失败。 */
void expectError(Cas::Engine *engine, const char *input) {
    const Cas::Result result = engine->evaluate(input);
    if (result.ok) {
        std::printf("FAIL %-52s ← 本应报错，却得到 %s\n", input, result.text.c_str());
        ++failures;
        return;
    }
    std::printf("  %-50s ↯ %s\n", input, result.message.c_str());
}

/* ------------------------------------------------------------ 线性方程组 */

void testLinearSystems() {
    Cas::Engine engine;

    expectText(&engine, "solveSystem([x+y=3, x-y=1], [x, y])", "x = 2, y = 1");
    expectText(&engine, "solveSystem([x=1, y=2], [x, y])", "x = 1, y = 2");
    expectText(&engine, "solveSystem([2*x+3*y=7, x-y=1], [x, y])", "x = 2, y = 1");

    // 三元的线性方程组
    expectText(&engine, "solveSystem([x+y+z=6, x-y=0, z=2], [x, y, z])", "x = 2, y = 2, z = 2");

    // 不相容：无解，而不是报错
    expectText(&engine, "solveSystem([x+y=1, x+y=2], [x, y])", "无解");
}

/* ---------------------------------------------------------- 非线性：两元 */

void testNonlinearTwo() {
    Cas::Engine engine;

    // 圆与双曲线的交点
    expectContains(&engine, "solveSystem([x^2+y^2=25, x*y=12], [x, y])", "x = 3, y = 4");
    expectContains(&engine, "solveSystem([x^2+y^2=25, x*y=12], [x, y])", "x = -3, y = -4");
    expectContains(&engine, "solveSystem([x^2+y^2=25, x*y=12], [x, y])", "x = 4, y = 3");
    expectContains(&engine, "solveSystem([x^2+y^2=25, x*y=12], [x, y])", "x = -4, y = -3");

    // 两个二次式相减（首项相同的坑）
    expectContains(&engine, "solveSystem([x^2+y^2=10, x^2-y^2=8], [x, y])", "x = 3, y = 1");
    expectContains(&engine, "solveSystem([x^2+y^2=10, x^2-y^2=8], [x, y])", "x = -3, y = -1");

    // 二次与一次混合
    expectContains(&engine, "solveSystem([x^2+y^2=5, x-y=1], [x, y])", "x = 2, y = 1");
    expectContains(&engine, "solveSystem([x^2+y^2=5, x-y=1], [x, y])", "x = -1, y = -2");

    // 乘积型：线性路径会给出 x = 6/y 这种伪解，必须走非线性
    expectContains(&engine, "solveSystem([x*y=6, x+y=5], [x, y])", "x = 2, y = 3");
    expectContains(&engine, "solveSystem([x*y=6, x+y=5], [x, y])", "x = 3, y = 2");

    // 只有复根 → 无实解
    expectText(&engine, "solveSystem([x^2+y^2=1, x*y=1], [x, y])", "无解");
}

/* ---------------------------------------------------------- 非线性：三元 */

void testNonlinearThree() {
    Cas::Engine engine;

    // 对称方程组：三数两两组合，共 6 个排列
    expectContains(&engine, "solveSystem([x^2+y^2+z^2=14, x+y+z=6, x*y*z=6], [x, y, z])", "x = 1, y = 2, z = 3");
    expectContains(&engine, "solveSystem([x^2+y^2+z^2=14, x+y+z=6, x*y*z=6], [x, y, z])", "x = 3, y = 2, z = 1");
    expectContains(&engine, "solveSystem([x^2+y^2+z^2=14, x+y+z=6, x*y*z=6], [x, y, z])", "x = 2, y = 3, z = 1");

    // 圆与坐标轴的交点
    expectContains(&engine, "solveSystem([x^2+y^2=1, x*y=0], [x, y])", "x = 1, y = 0");
    expectContains(&engine, "solveSystem([x^2+y^2=1, x*y=0], [x, y])", "x = 0, y = 1");
}

/* ------------------------------------------------------------ 单变量方程 */

void testSingleVariable() {
    Cas::Engine engine;

    expectText(&engine, "solveSystem([x^2=1], [x])", "x = 1; x = -1");
    expectContains(&engine, "solveSystem([x^3-2*x^2-x+2=0], [x])", "x = 2");
    expectText(&engine, "solveSystem([x^2=-1], [x])", "无解");
}

/* -------------------------------------------------------------- 欠定情形 */

void testUnderdetermined() {
    Cas::Engine engine;

    // 一条直线：应给出参数解而不是「无解」
    expectContains(&engine, "solveSystem([x+y=1], [x, y])", "存在自由变量");
    expectContains(&engine, "solveSystem([x+y=1], [x, y])", "x =");

    // 两条成比例的直线（同一理想）
    expectContains(&engine, "solveSystem([x+y=1, 2*x+2*y=2], [x, y])", "存在自由变量");

    // 一条曲线：解不唯一，也不该报「无解」
    expectContains(&engine, "solveSystem([x^2+y^2=1], [x, y])", "存在自由变量");

    // 三元自由度 1
    expectContains(&engine, "solveSystem([x-y=0, y-z=0], [x, y, z])", "存在自由变量");
}

/* ---------------------------------------------------------------- 边界 */

void testBoundaries() {
    Cas::Engine engine;

    // 非多项式方程超出求解范围，要报错而不是给出错误答案
    expectError(&engine, "solveSystem([sin(x)=0], [x])");

    // 变量表为空
    expectError(&engine, "solveSystem([x=1], [])");
}

} // namespace

int main() {
    std::printf("方程组求解自检\n");
    Cas::Allocator::initialize();

    testLinearSystems();       std::printf("linear systems ok\n");
    testNonlinearTwo();        std::printf("nonlinear two vars ok\n");
    testNonlinearThree();      std::printf("nonlinear three vars ok\n");
    testSingleVariable();      std::printf("single variable ok\n");
    testUnderdetermined();     std::printf("underdetermined ok\n");
    testBoundaries();          std::printf("boundaries ok\n");

    if (failures == 0) {
        std::printf("全部通过\n");
        return 0;
    }
    std::printf("%d 项失败\n", failures);
    return 1;
}

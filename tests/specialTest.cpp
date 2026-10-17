// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#include "cas/engine.h"

#include <cstdio>
#include <string>

/**
 * @file specialTest.cpp
 * @brief 特殊函数、级数求和与敛散性判定的回归测试。
 *
 * 覆盖范围：
 * - 特殊函数：Gamma（整数 / 半整数 / 极点）、阶乘与双阶乘、Beta、组合数、谐波数、
 *   erf / erfc、Bessel 各类、伯努利数，以及 `specialsimplify` 的折叠行为。
 * - 级数求和：多项式（Faulhaber）、p 级数（含 ζ(2n) 的 π 闭式）、几何级数、
 *   算术-几何级数、裂项级数，有限与无穷两种上限。
 * - 敛散性：`converge` 的收敛求和与发散判定。
 * - 序列识别：`guess` 的等差 / 等比 / 多项式。
 *
 * 输出约定与其它测试一致：每行 `ok   <描述>` 或 `FAIL <描述> ← 实际值`。
 */

namespace {

int g_failures = 0;
int g_total = 0;

/** @brief 断言文本结果等于期望值。 */
void expect(Cas::Engine &engine, const char *input, const char *expected) {
    ++g_total;
    const Cas::Result result = engine.evaluate(input);
    if (!result.ok) {
        std::printf("FAIL %-44s ← 错误：%s\n", input, result.message.c_str());
        ++g_failures;
        return;
    }
    if (result.text == expected) {
        std::printf("ok   %-44s = %s\n", input, result.text.c_str());
        return;
    }
    std::printf("FAIL %-44s ← 得到 %s，期望 %s\n", input, result.text.c_str(), expected);
    ++g_failures;
}

/** @brief 断言求值失败（用于极点、发散等应当拒绝的输入）。 */
void expectError(Cas::Engine &engine, const char *input) {
    ++g_total;
    const Cas::Result result = engine.evaluate(input);
    if (result.ok) {
        std::printf("FAIL %-44s ← 本该报错，却得到 %s\n", input, result.text.c_str());
        ++g_failures;
        return;
    }
    std::printf("ok   %-44s （按预期报错）\n", input);
}

/** @brief 断言结果文本包含指定片段（用于"收敛，和为 …"这类带说明的输出）。 */
void expectContains(Cas::Engine &engine, const char *input, const char *needle) {
    ++g_total;
    const Cas::Result result = engine.evaluate(input);
    if (!result.ok) {
        std::printf("FAIL %-44s ← 错误：%s\n", input, result.message.c_str());
        ++g_failures;
        return;
    }
    if (result.text.find(needle) != std::string::npos) {
        std::printf("ok   %-44s = %s\n", input, result.text.c_str());
        return;
    }
    std::printf("FAIL %-44s ← 得到 %s，未包含 %s\n", input, result.text.c_str(), needle);
    ++g_failures;
}

/* ------------------------------------------------------------ Gamma */

void testGamma(Cas::Engine &engine) {
    std::printf("\n== Gamma 函数 ==\n");
    expect(engine, "gamma(1)", "1");
    expect(engine, "gamma(2)", "1");
    expect(engine, "gamma(5)", "24");
    expect(engine, "gamma(10)", "362880");
    // Γ(1/2) = √π
    expect(engine, "gamma(1/2)", "sqrt(pi)");
    // Γ(1/2 + n) = (2n)!/(4^n n!) √π
    expect(engine, "gamma(3/2)", "1/2*sqrt(pi)");
    expect(engine, "gamma(5/2)", "3/4*sqrt(pi)");
    expect(engine, "gamma(7/2)", "15/8*sqrt(pi)");
    // Γ(1/2 - n) = (-4)^n n!/(2n)! √π，注意 n = 0 的情形已由 Γ(1/2) 覆盖
    expect(engine, "gamma(-1/2)", "-2*sqrt(pi)");
    expect(engine, "gamma(-3/2)", "4/3*sqrt(pi)");
    expect(engine, "gamma(-5/2)", "-8/15*sqrt(pi)");
    // 非极点、非半整数的值应保持未求值
    expect(engine, "gamma(0)", "gamma(0)");
    expect(engine, "gamma(1/3)", "gamma(1/3)");
    // 极点：0 与负整数
    expect(engine, "gamma(-1)", "gamma(-1)");
    expect(engine, "gamma(-7)", "gamma(-7)");
}

/* ------------------------------------------- 阶乘 / 双阶乘 / 组合数 */

void testFactorialFamily(Cas::Engine &engine) {
    std::printf("\n== 阶乘与组合数 ==\n");
    expect(engine, "0!", "1");
    expect(engine, "5!", "120");
    expect(engine, "10!", "3628800");
    expect(engine, "20!", "2432902008176640000");
    // 双阶乘：n!! = n(n-2)(n-4)…
    expect(engine, "7!!", "105");
    expect(engine, "6!!", "48");
    expect(engine, "9!!", "945");
    expect(engine, "10!!", "3840");
    expect(engine, "1!!", "1");
    expect(engine, "0!!", "1");
    // (-1)!! = 1（空积约定）
    expect(engine, "(-1)!!", "1");

    expect(engine, "binomial(5,2)", "10");
    expect(engine, "binomial(10,3)", "120");
    expect(engine, "binomial(100,0)", "1");
    expect(engine, "binomial(100,100)", "1");
    expect(engine, "binomial(20,10)", "184756");
    // k > n 时按组合数约定取 0（与 sympy 一致）
    expect(engine, "binomial(3,5)", "0");
}

/* ----------------------------------------- Beta / erf / Bessel / 谐波 */

void testBetaAndErrorFunctions(Cas::Engine &engine) {
    std::printf("\n== Beta / erf / 谐波数 ==\n");
    expect(engine, "beta(2,3)", "1/12");
    expect(engine, "beta(1,1)", "1");
    expect(engine, "beta(3,4)", "1/60");

    expect(engine, "erf(0)", "0");
    expect(engine, "erfc(0)", "1");
    expect(engine, "erfi(0)", "0");

    expect(engine, "harmonic(1)", "1");
    expect(engine, "harmonic(5)", "137/60");
    expect(engine, "harmonic(10)", "7381/2520");
    expect(engine, "harmonic(0)", "0");
    // 二阶谐波数 H(n, 2) = Σ 1/k²
    expect(engine, "harmonic(3,2)", "49/36");
}

void testBessel(Cas::Engine &engine) {
    std::printf("\n== Bessel 函数 ==\n");
    // J_0(0) = 1，J_n(0) = 0（n ≥ 1）
    expect(engine, "besselj(0,0)", "1");
    expect(engine, "besselj(1,0)", "0");
    expect(engine, "besselj(3,0)", "0");
    // Y_n(0) 在 0 处发散，保持未求值
    expect(engine, "bessely(1,0)", "bessely(1, 0)");
    // 非零点保持符号形式
    expect(engine, "besselj(0,1)", "besselj(0, 1)");
}

/* --------------------------------------------------------- 伯努利数 */

void testBernoulli(Cas::Engine &engine) {
    std::printf("\n== 伯努利数 ==\n");
    expect(engine, "bernoulli(0)", "1");
    expect(engine, "bernoulli(1)", "-1/2");
    expect(engine, "bernoulli(2)", "1/6");
    expect(engine, "bernoulli(3)", "0");
    expect(engine, "bernoulli(4)", "-1/30");
    expect(engine, "bernoulli(5)", "0");
    expect(engine, "bernoulli(6)", "1/42");
    expect(engine, "bernoulli(8)", "-1/30");
    expect(engine, "bernoulli(10)", "5/66");
    expect(engine, "bernoulli(12)", "-691/2730");
    // 奇数下标（> 1）恒为 0
    expect(engine, "bernoulli(7)", "0");
    expect(engine, "bernoulli(9)", "0");
}

/* -------------------------------------------- specialsimplify 折叠 */

void testSpecialSimplify(Cas::Engine &engine) {
    std::printf("\n== specialsimplify ==\n");
    expect(engine, "specialsimplify(gamma(6))", "120");
    expect(engine, "specialsimp(gamma(1/2))", "sqrt(pi)");
    expect(engine, "spsimp(binomial(7,2))", "21");
    // 复合表达式里的调用也会被折叠
    expect(engine, "specialsimplify(gamma(5) + binomial(4,2))", "30");
}

/* ------------------------------------------------------- 多项式求和 */

void testPolynomialSum(Cas::Engine &engine) {
    std::printf("\n== 幂和（Faulhaber）==\n");
    expect(engine, "sum(1, i=1..n)", "n");
    expect(engine, "sum(i, i=1..n)", "1/2*n^2 + 1/2*n");
    expect(engine, "sum(i^2, i=1..n)", "1/3*n^3 + 1/2*n^2 + 1/6*n");
    expect(engine, "sum(i^3, i=1..n)", "1/4*n^4 + 1/2*n^3 + 1/4*n^2");
    expect(engine, "sum(i^4, i=1..n)", "1/5*n^5 + 1/2*n^4 + 1/3*n^3 - 1/30*n");
    // 带系数的线性组合：必须整体走多项式分支，而不是逐项拆分后各自打包
    expect(engine, "sum(3i^2-2i+1, i=1..n)", "n^3 + 1/2*n^2 + 1/2*n");
    // 有限上限直接算出数值
    expect(engine, "sum(i, i=1..10)", "55");
    expect(engine, "sum(i, i=1..100)", "5050");
    expect(engine, "sum(i^2, i=1..10)", "385");
    expect(engine, "sum(i^2, i=1..100)", "338350");
}

/* --------------------------------------------------------- p 级数 */

void testPSeriesSum(Cas::Engine &engine) {
    std::printf("\n== p 级数与调和和 ==\n");
    expect(engine, "sum(1/i, i=1..5)", "137/60");
    expect(engine, "sum(1/i, i=1..10)", "7381/2520");
    expect(engine, "sum(1/i^2, i=1..10)", "1968329/1270080");
    // ζ(2n) = (-1)^(n+1) B_2n (2π)^(2n) / (2 (2n)!)
    expect(engine, "sum(1/i^2, i=1..oo)", "1/6*pi^2");
    expect(engine, "sum(1/i^4, i=1..oo)", "1/90*pi^4");
    expect(engine, "sum(1/i^6, i=1..oo)", "1/945*pi^6");
    // 从 2 开始的 p 级数要扣掉首项
    expect(engine, "sum(1/i^2, i=2..oo)", "1/6*pi^2 - 1");
}

/* ------------------------------------------------------- 几何级数 */

void testGeometricSum(Cas::Engine &engine) {
    std::printf("\n== 几何级数 ==\n");
    expect(engine, "sum((1/2)^i, i=0..oo)", "2");
    expect(engine, "sum(2^(-i), i=0..oo)", "2");
    expect(engine, "sum((1/3)^i, i=0..oo)", "3/2");
    expect(engine, "sum(x^i, i=0..3)", "x^3 + x^2 + x + 1");
    // 与 (1 - x^(n+1))/(1 - x) 等价，化简层未做分母符号归并
    expect(engine, "sum(x^i, i=0..n)", "(-x^(n + 1) + 1)/(-x + 1)");
    // 公比绝对值 ≥ 1 时无穷和不收敛
    expectError(engine, "sum(2^i, i=0..oo)");
}

/* ---------------------------------------------- 算术-几何级数 */

void testArithmeticGeometricSum(Cas::Engine &engine) {
    std::printf("\n== 算术-几何级数 ==\n");
    // Σ i q^i = q/(1-q)²
    expect(engine, "sum(i*(1/2)^i, i=0..oo)", "2");
    expect(engine, "sum(i*2^(-i), i=0..oo)", "2");
    expect(engine, "sum(i*(1/3)^i, i=1..oo)", "3/4");
    // Σ (i+1) q^i = 1/(1-q)²
    expect(engine, "sum((i+1)*(1/2)^i, i=0..oo)", "4");
}

/* --------------------------------------------------------- 裂项级数 */

void testTelescopingSum(Cas::Engine &engine) {
    std::printf("\n== 裂项级数 ==\n");
    // 1/(i(i+1)) = 1/i - 1/(i+1)，无穷和只剩首项
    expect(engine, "sum(1/(i*(i+1)), i=1..oo)", "1");
    expect(engine, "sum(1/(i^2+i), i=1..oo)", "1");
    // 1/(i(i+2)) = 1/2·(1/i - 1/(i+2))，剩 1/2·(1+1/2)
    expect(engine, "sum(1/(i*(i+2)), i=1..oo)", "3/4");
    // 偏移不为零的裂项：1/((i+2)(i+5)) = 1/3·(1/(i+2) - 1/(i+5))
    expect(engine, "sum(1/((i+2)*(i+5)), i=0..oo)", "13/36");
    // 1/(i²-1) 在 i ≥ 2 上裂项
    expect(engine, "sum(1/(i^2-1), i=2..oo)", "3/4");
    // 有限上限
    expect(engine, "sum(1/(i*(i+1)), i=1..5)", "5/6");
    expect(engine, "sum(1/(i*(i+1)), i=1..n)", "harmonic(n) - harmonic(0) - (harmonic(n + 1) - harmonic(1))");
}

/* --------------------------------------------------------- 敛散性 */

void testConvergence(Cas::Engine &engine) {
    std::printf("\n== 敛散性 ==\n");
    expectContains(engine, "converge(1/i^2, i=1..oo)", "收敛");
    expectContains(engine, "converge(1/i^2, i=1..oo)", "1/6*pi^2");
    expectContains(engine, "converge(1/i, i=1..oo)", "发散");
    // 等比级数：|q| < 1 收敛，否则发散
    expectContains(engine, "converge((1/2)^i, i=0..oo)", "收敛");
    expectContains(engine, "converge(2^i, i=0..oo)", "发散");
}

/* ----------------------------------------------------- 序列识别 */

void testGuess(Cas::Engine &engine) {
    std::printf("\n== 序列识别 ==\n");
    expect(engine, "guess([1,3,5,7,9])", "2*n + 1");
    expect(engine, "guess([1,4,9,16,25])", "n^2 + 2*n + 1");
    expect(engine, "guess([5,5,5,5])", "5");
    // 等比
    expect(engine, "guess([2,4,8,16,32])", "2*2^n");
    // 参数不是列表时报错
    expectError(engine, "guess(1,4,9,16)");
}

/* ------------------------------------------------------- 综合 / 交互 */

void testCrossFeature(Cas::Engine &engine) {
    std::printf("\n== 跨模块组合 ==\n");
    // 求和结果可以直接进入化简 / 求值链
    expect(engine, "sum(i, i=1..10) + sum(i^2, i=1..10)", "440");
    // 特殊函数参与代数运算
    expect(engine, "gamma(5)/factorial(4)", "1");
    expect(engine, "binomial(6,3) + harmonic(3)", "131/6");
}

} // namespace

int main() {
    Cas::Engine engine;

    testGamma(engine);
    testFactorialFamily(engine);
    testBetaAndErrorFunctions(engine);
    testBessel(engine);
    testBernoulli(engine);
    testSpecialSimplify(engine);
    testPolynomialSum(engine);
    testPSeriesSum(engine);
    testGeometricSum(engine);
    testArithmeticGeometricSum(engine);
    testTelescopingSum(engine);
    testConvergence(engine);
    testGuess(engine);
    testCrossFeature(engine);

    std::printf("\n------------------------------------\n");
    std::printf("共 %d 项，失败 %d 项\n", g_total, g_failures);
    return g_failures == 0 ? 0 : 1;
}

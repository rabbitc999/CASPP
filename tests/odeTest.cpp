#include "cas/engine.h"

#include <cstdio>
#include <string>

/**
 * @file odeTest.cpp
 * @brief 常微分方程求解的回归测试。
 *
 * 覆盖范围：
 * - 一阶：可分离变量、一阶线性（含积分因子与 `e^{c·log u} = u^c` 化简）、
 *   齐次型（换元 u = y/x）、伯努利方程。
 * - 二阶常系数线性齐次：两相异实根、重根、共轭复根；以及缺 y' 项、缺常数项、
 *   同时缺 y' 与 y 的退化情形。
 * - 二阶常系数线性非齐次：右端为多项式、`e^{kx}`、`cos(kx)`/`sin(kx)`，
 *   以及右端恰为齐次解时的共振（试探式乘 x）。
 * - 边界：三阶以上拒绝、方程里没有未知函数时报错。
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

/** @brief 断言求值失败。 */
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

/** @brief 断言结果文本包含指定片段。 */
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

/* ------------------------------------------------------ 一阶：线性 */

void testLinearFirstOrder(Cas::Engine &engine) {
    std::printf("\n== 一阶线性方程 ==\n");

    // y' = y：最简线性，积分因子退化为 1
    expect(engine, "dsolve(y' = y)", "y = C1*e^x");
    // y' = x·y：P(x) = x，指数上出现 x²/2
    expect(engine, "dsolve(y' = x*y)", "y = C1*e^(1/2*x^2)");
    // y'/y = 1/x：P(x) = 1/x，要验证 log 型积分因子
    expect(engine, "dsolve(y' = y/x)", "y = C1*x");

    // y' + y/x = x²：Q(x) = x²，积分因子为 x，通解要通分
    expect(engine, "dsolve(y' + y/x = x^2)", "y = (1/4*x^4 + C1)/x");
    // 同上但符号相反，考验 e^{-log x} = 1/x 的化简
    expect(engine, "dsolve(y' - y/x = x^2)", "y = x*(1/2*x^2 + C1)");

    // y' + 2xy = x：特解含 e^{-x²}·e^{x²} 的配对，考验一般换元积分 ∫x·e^{x²}dx
    expect(engine, "dsolve(y' + 2*x*y = x)", "y = e^(-x^2)*(C1 + 1/2*e^(x^2))");
    // y' + y = x：积分因子 e^x，答案含 x·e^x 的展开
    expect(engine, "dsolve(y' + y = x)", "y = e^(-x)*(x*e^x + C1 - e^x)");
}

/* -------------------------------------------------- 一阶：可分离 */

void testSeparable(Cas::Engine &engine) {
    std::printf("\n== 可分离变量方程 ==\n");

    // 右端只含 y：∫dy/y² = ∫dx
    expect(engine, "dsolve(y' = y^2)", "隐式解（可分离变量）：-C1 - x - 1/y = 0");
    // 右端为 x·y²：∫dy/y² = ∫x dx
    expect(engine, "dsolve(y' = x*y^2)", "隐式解（可分离变量）：-1/2*x^2 - C1 - 1/y = 0");
    // 右端含 y³
    expect(engine, "dsolve(y' = y^2*x^3)", "隐式解（可分离变量）：-1/4*x^4 - C1 - 1/y = 0");

    // 右端不含 y：退化为纯积分，可直接给出显式解
    expect(engine, "dsolve(y' = 3*x^2)", "y = x^3 + C1");
    expect(engine, "dsolve(y' = e^x)", "y = C1 + e^x");
}

/* ---------------------------------------------------- 一阶：齐次型 */

void testHomogeneous(Cas::Engine &engine) {
    std::printf("\n== 齐次型方程 ==\n");

    // y' = y/x + 1：F(u) = u + 1，换元后立即退化为线性
    expect(engine, "dsolve(y' = y/x + 1)", "y = x*(C1 + log(x))");
    // y' = (x²+y²)/(xy)：F(u) = (1+u²)/u，F(u)−u = 1/u，积出 u²/2 − log x
    expect(engine, "dsolve(y' = (x^2+y^2)/(x*y))",
           "隐式解（齐次型（换元 u = y/x））：y^2/(2*x^2) + C1 - log(x) = 0");
    // y' = (y/x)²：F(u) − u = u(u−1)/u … 积出对数组合
    expect(engine, "dsolve(y' = (y/x)^2)",
           "隐式解（齐次型（换元 u = y/x））：C1 + log(y/x - 1) - log(x) - log(y/x) = 0");
}

/* ------------------------------------------------ 一阶：伯努利方程 */

void testBernoulli(Cas::Engine &engine) {
    std::printf("\n== 伯努利方程 ==\n");

    // y' − y = x·y²：n = 2，换元 v = 1/y 得 v' + v = −x
    expect(engine, "dsolve(y' - y = x*y^2)",
           "隐式解（伯努利方程（换元 v = y^(1-n)））：e^x/(-x*e^x + C1 + e^x) - y = 0");
    // y' + y/x = y²：v' − v/x = −1，积分因子 1/x
    expect(engine, "dsolve(y' + y/x = y^2)", "隐式解（伯努利方程（换元 v = y^(1-n)））：1/(x*(C1 - log(x))) - y = 0");
    // y' − 2x·y = x·y²：v' + 2x·v = −x
    expect(engine, "dsolve(y' - 2*x*y = x*y^2)",
           "隐式解（伯努利方程（换元 v = y^(1-n)））：e^(x^2)/(C1 - 1/2*e^(x^2)) - y = 0");
}

/* ---------------------------------------------- 二阶：齐次三情形 */

void testSecondOrderHomogeneous(Cas::Engine &engine) {
    std::printf("\n== 二阶常系数齐次 ==\n");

    // 两相异实根 r = 1, −1
    expect(engine, "dsolve(y'' = y)", "y = C1*e^x + C2*e^(-x)");
    // 两相异实根 r = 2, 1
    expect(engine, "dsolve(y'' - 3*y' + 2*y = 0)", "y = C1*e^(2*x) + C2*e^x");
    // 共轭复根 r = ±i
    expect(engine, "dsolve(y'' + y = 0)", "y = C1*cos(x) + C2*sin(x)");
    // 共轭复根 r = ±2i
    expect(engine, "dsolve(y'' + 4*y = 0)", "y = C1*cos(2*x) + C2*sin(2*x)");
    // 重根 r = −1：解要乘 x
    expect(engine, "dsolve(y'' + 2*y' + y = 0)", "y = e^(-x)*(C2*x + C1)");

    // 缺 y' 与常数项（只有 y'' 与 y）
    expect(engine, "dsolve(y'' = y)", "y = C1*e^x + C2*e^(-x)");
    // 完全退化：y'' = 0，解是一次多项式
    expect(engine, "dsolve(y'' = 0)", "y = C2*x + C1");
}

/* -------------------------------------------- 二阶：非齐次特解 */

void testSecondOrderParticular(Cas::Engine &engine) {
    std::printf("\n== 二阶非齐次（待定系数法）==\n");

    // 常数右端
    expect(engine, "dsolve(y'' + y = 1)", "y = C1*cos(x) + C2*sin(x) + 1");
    // 一次多项式
    expect(engine, "dsolve(y'' - 3*y' + 2*y = x)", "y = C1*e^(2*x) + C2*e^x + 1/2*x + 3/4");
    // 二次多项式（右端缺一次项，同样要补全试探式）
    expect(engine, "dsolve(y'' - y = x^2)", "y = -x^2 + C1*e^x + C2*e^(-x) - 2");
    expect(engine, "dsolve(y'' + 2*y' + y = x^2)", "y = x^2 + e^(-x)*(C2*x + C1) - 4*x + 6");

    // 指数右端（非共振）
    expect(engine, "dsolve(y'' + y = e^x)", "y = C1*cos(x) + C2*sin(x) + 1/2*e^x");
    // 三角右端（非共振）
    expect(engine, "dsolve(y'' + y = cos(2*x))", "y = C1*cos(x) + C2*sin(x) - 1/3*cos(2*x)");
}

/* ------------------------------------------------- 二阶：共振情形 */

void testSecondOrderResonance(Cas::Engine &engine) {
    std::printf("\n== 二阶共振（试探式乘 x）==\n");

    // y'' + y = sin x：sin x 本身是齐次解，特解要乘 x
    expect(engine, "dsolve(y'' + y = sin(x))", "y = -1/2*x*cos(x) + C1*cos(x) + C2*sin(x)");
    // y'' + 4y = sin 2x：同样是共振
    expect(engine, "dsolve(y'' + 4*y = sin(2*x))", "y = -1/4*x*cos(2*x) + C1*cos(2*x) + C2*sin(2*x)");
    // y'' − 3y' + 2y = e^x：e^x 对应单根 r = 1，特解乘一次 x
    expect(engine, "dsolve(y'' - 3*y' + 2*y = e^x)", "y = -x*e^x + C1*e^(2*x) + C2*e^x");
}

/* --------------------------------------------------------- 边界 */

void testBoundary(Cas::Engine &engine) {
    std::printf("\n== 边界与报错 ==\n");

    // 三阶以上目前不支持
    expectError(engine, "dsolve(y''' = y)");
    // 方程里没有 y：不是微分方程
    expectError(engine, "dsolve(x = 1)");
    // 指定了不存在的未知函数
    expectError(engine, "dsolve(y' = x, z)");

    // 隐式解要带积分常数，否则答案不完整
    expectContains(engine, "dsolve(y' = y^2)", "C1");
}

} // namespace

int main() {
    Cas::Engine engine;

    testLinearFirstOrder(engine);
    testSeparable(engine);
    testHomogeneous(engine);
    testBernoulli(engine);
    testSecondOrderHomogeneous(engine);
    testSecondOrderParticular(engine);
    testSecondOrderResonance(engine);
    testBoundary(engine);

    std::printf("\n------------------------------------\n");
    std::printf("共 %d 项，失败 %d 项\n", g_total, g_failures);
    return g_failures == 0 ? 0 : 1;
}

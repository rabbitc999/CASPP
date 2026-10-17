#include "cas/engine.h"

#include <cstdio>
#include <string>
#include <vector>

/**
 * @file calculusTest.cpp
 * @brief 微积分 / 求解 / 解析 / 输出 的回归测试。
 *
 * 输出约定：每行 `ok   <描述>` 或 `FAIL <描述> ← 实际值`。
 * 任一失败都会让进程返回 1，便于 make test 汇总。
 */

namespace {

int g_failures = 0;
int g_total = 0;

/** @brief 断言文本结果等于期望值。 */
void expect(Cas::Engine &engine, const char *input, const char *expected) {
    ++g_total;
    const Cas::Result result = engine.evaluate(input);
    if (!result.ok) {
        std::printf("FAIL %-42s ← 错误：%s\n", input, result.message.c_str());
        ++g_failures;
        return;
    }
    if (result.text == expected) {
        std::printf("ok   %-42s = %s\n", input, result.text.c_str());
        return;
    }
    std::printf("FAIL %-42s ← 得到 %s，期望 %s\n", input, result.text.c_str(), expected);
    ++g_failures;
}

/** @brief 断言求值失败（用于"应当报错"的用例）。 */
void expectError(Cas::Engine &engine, const char *input) {
    ++g_total;
    const Cas::Result result = engine.evaluate(input);
    if (!result.ok) {
        std::printf("ok   %-42s 正确报错\n", input);
        return;
    }
    std::printf("FAIL %-42s ← 本应报错，却得到 %s\n", input, result.text.c_str());
    ++g_failures;
}

/** @brief 只校验解析成功（表达式可解析即可）。 */
void expectParse(Cas::Engine &engine, const char *input) {
    ++g_total;
    const Cas::Result result = engine.evaluate(input);
    if (result.ok) {
        std::printf("ok   %-42s → %s\n", input, result.text.c_str());
        return;
    }
    std::printf("FAIL %-42s ← 解析失败：%s (位置 %zu)\n", input, result.message.c_str(), result.position);
    ++g_failures;
}

/* ------------------------------------------------------------- 用例组 */

void testParsing(Cas::Engine &engine) {
    std::printf("\n== 解析 ==\n");
    expect(engine, "2x + 3x", "5*x");
    expect(engine, "2(x+1)", "2*(x + 1)"); // 隐式乘法，默认不展开
    expect(engine, "(x+1)(x-1)", "(x - 1)*(x + 1)"); // 隐式乘法，默认不展开
    expect(engine, "x^2^3", "x^8");
    expect(engine, "-x^2", "-x^2");
    expect(engine, "1/2 + 1/3", "5/6");
    expect(engine, "0.5", "1/2");
    expect(engine, "2.5e-2", "1/40");
    expect(engine, "3!", "6");
    expect(engine, "sin(0)", "0");
    expect(engine, "cos(0)", "1");
    expect(engine, "sqrt(4)", "2");
    expect(engine, "sqrt(2)^2", "2");
    expect(engine, "exp(0)", "1");
    expect(engine, "log(1)", "0");
    expectError(engine, "1 +");
    expectError(engine, "(x+1");
    expectError(engine, "x + @");
}

void testSimplification(Cas::Engine &engine) {
    std::printf("\n== 化简 ==\n");
    expect(engine, "simplify(x + x + x)", "3*x");
    expect(engine, "simplify(x*x*x)", "x^3");
    expect(engine, "simplify((x^2)^3)", "x^6");
    expect(engine, "simplify(2x/2)", "x");
    expect(engine, "simplify(x - x)", "0");
    expect(engine, "simplify(x^0)", "1");
    expect(engine, "simplify(1 + 1/2)", "3/2");
    expect(engine, "simplify(2^10)", "1024");
}

void testExpansion(Cas::Engine &engine) {
    std::printf("\n== 展开 ==\n");
    expect(engine, "expand((x+1)^2)", "x^2 + 2*x + 1");
    expect(engine, "expand((x+1)^3)", "x^3 + 3*x^2 + 3*x + 1");
    expect(engine, "expand((x-1)(x+1))", "x^2 - 1");
    expect(engine, "expand((x+y)^2)", "x^2 + y^2 + 2*x*y");
    // 三元及以上：验证卷积不会残留嵌套括号
    expect(engine, "expand((x+y+z)^3)",
           "x^3 + y^3 + z^3 + 3*x*y^2 + 3*x*z^2 + 3*y*x^2 + 3*y*z^2 + 3*z*x^2 + 3*z*y^2 + 6*x*y*z");
    expect(engine, "expand((a+b+c+d)^2)",
           "a^2 + b^2 + c^2 + d^2 + 2*a*b + 2*a*c + 2*a*d + 2*b*c + 2*b*d + 2*c*d");
    // 二项高次：系数应为组合数
    expect(engine, "expand((x+y)^8)", "x^8 + y^8 + 8*x*y^7 + 8*y*x^7 + 28*x^2*y^6 + 28*x^6*y^2 + "
                                      "56*x^3*y^5 + 56*x^5*y^3 + 70*x^4*y^4");
}

void testFactoring(Cas::Engine &engine) {
    std::printf("\n== 因式分解 ==\n");
    expect(engine, "factor(x^2-1)", "(x - 1)*(x + 1)");
    expect(engine, "factor(x^2-9)", "(x - 3)*(x + 3)");
    expect(engine, "factor(x^2+2x+1)", "(x + 1)^2");
    expect(engine, "factor(2x^2+7x+3)", "(x + 3)*(2*x + 1)");
    expect(engine, "factor(x^2+1)", "x^2 + 1"); // 有理数域上不可约
}

void testDerivative(Cas::Engine &engine) {
    std::printf("\n== 求导 ==\n");
    expect(engine, "diff(x^2, x)", "2*x");
    expect(engine, "diff(x^3+2x, x)", "3*x^2 + 2");
    expect(engine, "diff(sin(x), x)", "cos(x)");
    expect(engine, "diff(cos(x), x)", "-sin(x)");
    expect(engine, "diff(exp(x), x)", "exp(x)");
    expect(engine, "diff(log(x), x)", "1/x");
    expect(engine, "diff(sin(x^2), x)", "2*x*cos(x^2)");
    expect(engine, "diff(x*exp(x), x)", "x*exp(x) + exp(x)");
    expect(engine, "diff(1/x, x)", "-1/x^2");
    expect(engine, "diff(sqrt(x), x)", "1/(2*sqrt(x))");
}

void testIntegration(Cas::Engine &engine) {
    std::printf("\n== 积分 ==\n");
    expect(engine, "integrate(x^2, x)", "1/3*x^3");
    expect(engine, "integrate(x, x)", "1/2*x^2");
    expect(engine, "integrate(1/x, x)", "log(x)");
    expect(engine, "integrate(sin(x), x)", "-cos(x)");
    expect(engine, "integrate(cos(x), x)", "sin(x)");
    expect(engine, "integrate(exp(x), x)", "exp(x)");
    expect(engine, "integrate(x*exp(x), x)", "x*exp(x) - exp(x)");
    expect(engine, "integrate(x*sin(x), x)", "-x*cos(x) + sin(x)");
    expect(engine, "integrate(3x^2 + 2x, x)", "x^3 + x^2");
    // 分部积分
    expect(engine, "integrate(x^2*exp(x), x)", "x^2*exp(x) - 2*x*exp(x) + 2*exp(x)");
    expect(engine, "integrate(x^3*exp(x), x)", "x^3*exp(x) - 3*x^2*exp(x) + 6*x*exp(x) - 6*exp(x)");
    // 定积分
    expect(engine, "integrate(x^2, 0, 1)", "1/3");
    expect(engine, "integrate(2x, 0, 2)", "4");
    expect(engine, "integrate(sin(x), 0, pi)", "2");
    // 变量写在第二位的写法，以及上下限代入后的特殊角折叠
    expect(engine, "integrate(sin(x), x, 0, pi)", "2");
    expect(engine, "integrate(cos(x), 0, pi/2)", "1");
    expect(engine, "integrate(1/x, x, 1, 2)", "log(2)");
}

void testRationalIntegration(Cas::Engine &engine) {
    std::printf("\n== 有理函数积分 ==\n");
    // 分母一次：多项式除法 + 对数
    expect(engine, "integrate(1/x, x)", "log(x)");
    expect(engine, "integrate(1/(x+1), x)", "log(x + 1)");
    // 分母二次、判别式 > 0：部分分式的两个对数
    expect(engine, "integrate(1/(x^2-1), x)", "-1/2*log(x + 1) + 1/2*log(x - 1)");
    expect(engine, "integrate(1/(x^2-4), x)", "-1/4*log(x + 2) + 1/4*log(x - 2)");
    // 分母二次、判别式 < 0：对数 + atan
    expect(engine, "integrate(1/(x^2+1), x)", "atan(x)");
    expect(engine, "integrate(1/(x^2+4), x)", "1/2*atan(1/2*x)");
    // 定积分得到精确的 π 有理倍数
    expect(engine, "integrate(1/(x^2+1), x, 0, 1)", "1/4*pi");
    expect(engine, "integrate(1/(x^2+4), x, 0, 2)", "1/8*pi");
}

void testLimit(Cas::Engine &engine) {
    std::printf("\n== 极限 ==\n");
    expect(engine, "limit(sin(x)/x, x, 0)", "1");
    expect(engine, "limit((1-cos(x))/x^2, x, 0)", "1/2");
    expect(engine, "limit(tan(x)/x, x, 0)", "1");
    expect(engine, "limit((exp(x)-1)/x, x, 0)", "1");
    expect(engine, "limit((x^2-1)/(x-1), x, 1)", "2");
    expect(engine, "limit(1/x, x, oo)", "0");
    expect(engine, "limit((3x^2+1)/(x^2-2), x, oo)", "3");
    expect(engine, "limit((2x+1)/(x-3), x, oo)", "2");
    expect(engine, "limit(x^2, x, -oo)", "inf");
    expect(engine, "limit(x^3, x, -oo)", "-inf");
}

void testTaylor(Cas::Engine &engine) {
    std::printf("\n== Taylor 展开 ==\n");
    expect(engine, "taylor(exp(x), 0, 4)", "1/24*x^4 + 1/6*x^3 + 1/2*x^2 + x + 1");
    expect(engine, "taylor(sin(x), 0, 5)", "1/120*x^5 - 1/6*x^3 + x");
    expect(engine, "taylor(cos(x), 0, 4)", "1/24*x^4 - 1/2*x^2 + 1");
    expect(engine, "taylor(1/(1-x), 0, 3)", "x^3 + x^2 + x + 1");
}

void testSolve(Cas::Engine &engine) {
    std::printf("\n== 方程求解 ==\n");
    expect(engine, "solve(x^2-5x+6=0, x)", "x = 2, x = 3");
    expect(engine, "solve(x^2-1=0, x)", "x = -1, x = 1");
    expect(engine, "solve(2x+4=0, x)", "x = -2");
    expect(engine, "solve(x^2+x+1=0, x)", "无解（或在当前支持的范围内无法求出解析解）");
    expect(engine, "solve([x+y=3, x-y=1], [x, y])", "x = 2, y = 1");
    expect(engine, "solve([2x+3y=7, x-y=1], [x, y])", "x = 2, y = 1");
    expect(engine, "solve([x+y+z=6, x-y=0, z=2], [x,y,z])", "x = 2, y = 2, z = 2");
}

void testInequality(Cas::Engine &engine) {
    std::printf("\n== 不等式 ==\n");
    expectParse(engine, "x^2-1>0");
    expectParse(engine, "2x+1<3");
}

void testTrig(Cas::Engine &engine) {
    std::printf("\n== 三角化简 ==\n");
    expect(engine, "trigsimp(sin(x)^2+cos(x)^2)", "1");
    expect(engine, "trigsimp(2sin(x)^2+2cos(x)^2)", "2");
    expect(engine, "sin(pi/6)", "1/2");
    expect(engine, "sin(pi/4)", "1/2*sqrt(2)");
    expect(engine, "sin(pi/3)", "1/2*sqrt(3)");
    expect(engine, "cos(pi/3)", "1/2");
    expect(engine, "cos(pi/4)", "1/2*sqrt(2)");
    expect(engine, "sin(pi/2)", "1");
    expect(engine, "cos(pi)", "-1");
    expect(engine, "sin(pi)", "0");
    expect(engine, "sin(-pi/6)", "-1/2");
    // 反三角特殊值
    expect(engine, "atan(1)", "1/4*pi");
    expect(engine, "atan(0)", "0");
    expect(engine, "asin(1)", "1/2*pi");
    expect(engine, "acos(0)", "1/2*pi");
    expect(engine, "acos(-1)", "pi");
    // 根式精确开方
    expect(engine, "sqrt(4)", "2");
    expect(engine, "sqrt(16)", "4");
    expect(engine, "8^(1/3)", "2");
    expect(engine, "sqrt(2)", "sqrt(2)"); // 无理数保持根式
}

void testGroebner(Cas::Engine &engine) {
    std::printf("\n== Groebner 基 ==\n");
    expectParse(engine, "groebner([x^2+y^2-1, x-y], [x, y])");
    expectParse(engine, "groebner([x*y-1, x+y-3], [x, y])");
}

void testLatex(Cas::Engine &engine) {
    std::printf("\n== LaTeX 输出 ==\n");
    ++g_total;
    const Cas::Result result = engine.toLatexResult("1/x + x^2");
    const bool good = result.ok && result.latex.find("\\frac") != std::string::npos &&
                      result.latex.find("x^{2}") != std::string::npos;
    if (good) {
        std::printf("ok   %-42s → %s\n", "latex(1/x + x^2)", result.latex.c_str());
    } else {
        std::printf("FAIL %-42s ← %s\n", "latex(1/x + x^2)", result.latex.c_str());
        ++g_failures;
    }
}

void testMemory(Cas::Engine &engine) {
    std::printf("\n== 内存池 ==\n");
    ++g_total;
    const Cas::MemoryStats stats = engine.memoryStats();
    if (stats.capacity > 0) {
        std::printf("ok   池容量 %zu 字节，活动 %zu 字节 / %zu 块，峰值 %zu 字节\n", stats.capacity,
                    stats.liveBytes, stats.liveBlocks, stats.peakBytes);
    } else {
        std::printf("FAIL 内存池未初始化\n");
        ++g_failures;
    }

    // 反复计算后活动块应能回落到较低水平（表达式随作用域释放）
    for (int i = 0; i < 200; ++i) engine.evaluate("expand((x+1)^6)");
    ++g_total;
    const Cas::MemoryStats after = engine.memoryStats();
    if (after.liveBytes < stats.capacity / 2) {
        std::printf("ok   200 次展开后活动 %zu 字节（未泄漏）\n", after.liveBytes);
    } else {
        std::printf("FAIL 200 次展开后活动 %zu 字节，疑似泄漏\n", after.liveBytes);
        ++g_failures;
    }
}

/**
 * @brief 微积分加深：三角有理式、根式代换、部分分式、幂指极限、Laurent、驻点。
 *
 * 这一组覆盖的规则都是"旧版返回原式/报错"的场景，因此期望值都写死成闭式，
 * 一旦某条规则退化就能直接看出来。
 */
void testCalculusExtensions(Cas::Engine &engine) {
    // 三角有理式：万能代换 t = tan(u/2)
    expect(engine, "integrate(1/sin(x), x)", "log(tan(1/2*x))");
    expect(engine, "integrate(1/(1+cos(x)), x)", "tan(1/2*x)");

    // 根式三角代换
    expect(engine, "integrate(1/sqrt(1-x^2), x)", "asin(x)");
    expect(engine, "integrate(x/sqrt(1-x^2), x)", "-sqrt((-x^2 + 1))");

    // 部分分式分解（含重数与不可约二次因子）
    expect(engine, "integrate(1/(t^3-1), t)",
           "-1/2*atan((t + 1/2)/sqrt(3/4))*1/sqrt(3/4) - 1/6*log(t^2 + t + 1) + 1/3*log(t - 1)");

    // sin/cos 幂与 log 幂
    expect(engine, "integrate(sin(x)^2, x)", "-1/4*sin(2*x) + 1/2*x");
    expect(engine, "integrate(cos(x)^3, x)", "-1/3*sin(x)^3 + sin(x)");
    expect(engine, "integrate(log(x), x)", "x*log(x) - x");

    // 幂指型未定式：F^G = e^(G·ln F)
    expect(engine, "limit((1+x)^(1/x), x, 0)", "e");
    expect(engine, "limit(x^x, x, 0)", "1");
    // x→∞ 换元兜底：x = 1/t 转 t→0
    expect(engine, "limit(x*log(1+1/x), x, oo)", "1");
    expect(engine, "limit(x*sin(1/x), x, oo)", "1");

    // Laurent 展开：奇点处出现负幂项
    expect(engine, "taylor(exp(x)/x, 0, 3)", "1/6*x^2 + 1/x + 1/2*x + 1");
    expect(engine, "taylor(cos(x)/x^2, 0, 3)", "1/x^2 - 1/2");

    // 驻点与极值判定
    expect(engine, "critical(x^2-2*x)", "x = 1（极小值，f = -1）");
    expect(engine, "critical(x^3-3*x)", "x = -1（极大值，f = 2）; x = 1（极小值，f = -2）");
    expect(engine, "critical(x^3)", "x = 0（无法判定）");
    expect(engine, "critical(sin(x), x)", "x = 1/2*pi（极大值，f = 1）");
    expect(engine, "critical(cos(x), x)", "x = 0（极大值，f = 1）");
    // 严格单调：没有驻点，应明确报错而不是返回空串
    expectError(engine, "critical(2*x+1)");
}

} // namespace

int main() {
    Cas::Engine engine;

    testParsing(engine);
    testSimplification(engine);
    testExpansion(engine);
    testFactoring(engine);
    testDerivative(engine);
    testIntegration(engine);
    testRationalIntegration(engine);
    testCalculusExtensions(engine);
    testLimit(engine);
    testTaylor(engine);
    testSolve(engine);
    testInequality(engine);
    testTrig(engine);
    testGroebner(engine);
    testLatex(engine);
    testMemory(engine);

    std::printf("\n------------------------------------\n");
    std::printf("共 %d 项，失败 %d 项\n", g_total, g_failures);
    return g_failures == 0 ? 0 : 1;
}

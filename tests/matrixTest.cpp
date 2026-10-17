// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#include "cas/engine.h"

#include <cstdio>
#include <string>

/**
 * @file matrixTest.cpp
 * @brief 线性代数模块的回归测试。
 *
 * 引擎把矩阵结果渲染成 `[[a, b], [c, d]]`，并在文本开头附一个内部标记
 * （供 `det(inv(A))` 这类嵌套调用识别"上一步的结果是矩阵"）。测试比较前统一剥掉，
 * 只比较用户可见的部分。
 *
 * 输出约定与其它测试一致：每行 `ok   <描述>` 或 `FAIL <描述> ← 实际值`。
 */

namespace {

int g_failures = 0;
int g_total = 0;

/** @brief 去掉引擎内部使用的矩阵标记。 */
std::string display(const std::string &text) {
    std::size_t begin = 0;
    while (begin < text.size() && text[begin] == '\x1F') ++begin;
    return text.substr(begin);
}

/** @brief 断言文本结果等于期望值。 */
void expect(Cas::Engine &engine, const char *input, const char *expected) {
    ++g_total;
    const Cas::Result result = engine.evaluate(input);
    if (!result.ok) {
        std::printf("FAIL %-46s ← 错误：%s\n", input, result.message.c_str());
        ++g_failures;
        return;
    }
    const std::string actual = display(result.text);
    if (actual == expected) {
        std::printf("ok   %-46s = %s\n", input, actual.c_str());
        return;
    }
    std::printf("FAIL %-46s ← 得到 %s，期望 %s\n", input, actual.c_str(), expected);
    ++g_failures;
}

/** @brief 断言求值失败（用于"应当报错"的用例）。 */
void expectError(Cas::Engine &engine, const char *input) {
    ++g_total;
    const Cas::Result result = engine.evaluate(input);
    if (!result.ok) {
        std::printf("ok   %-46s 正确报错\n", input);
        return;
    }
    std::printf("FAIL %-46s ← 本应报错，却得到 %s\n", input, display(result.text).c_str());
    ++g_failures;
}

/** @brief 断言 LaTeX 输出里包含指定片段。 */
void expectLatexContains(Cas::Engine &engine, const char *input, const char *needle) {
    ++g_total;
    const Cas::Result result = engine.evaluate(input);
    if (result.ok && result.latex.find(needle) != std::string::npos) {
        std::printf("ok   %-46s LaTeX 含 %s\n", input, needle);
        return;
    }
    std::printf("FAIL %-46s ← LaTeX 为 %s\n", input, result.latex.c_str());
    ++g_failures;
}

/* ------------------------------------------------------------- 用例组 */

void testDeterminant(Cas::Engine &engine) {
    std::printf("\n== 行列式 ==\n");
    expect(engine, "det([[1,2],[3,4]])", "-2");
    expect(engine, "det([[1,2,3],[4,5,6],[7,8,10]])", "-3");
    // 两行成比例 → 0
    expect(engine, "det([[1,2,3],[4,5,6],[7,8,9]])", "0");
    expect(engine, "det([[5]])", "5");
    expect(engine, "determinant([[2,0],[0,3]])", "6");
    // 符号矩阵：ad - bc
    expect(engine, "det([[a,b],[c,d]])", "-b*c + a*d");
    // 三角矩阵的行列式 = 对角元乘积
    expect(engine, "det([[1,2,3],[0,5,6],[0,0,9]])", "45");
}

void testInverse(Cas::Engine &engine) {
    std::printf("\n== 逆矩阵 ==\n");
    expect(engine, "inv([[1,2],[3,4]])", "[[-2, 1], [3/2, -1/2]]");
    expect(engine, "inv([[2,0],[0,4]])", "[[1/2, 0], [0, 1/4]]");
    expect(engine, "inverse([[1,0],[0,1]])", "[[1, 0], [0, 1]]");
    expect(engine, "invert([[2,1],[1,2]])", "[[2/3, -1/3], [-1/3, 2/3]]");
    // 奇异矩阵必须报错，而不是给出伪逆
    expectError(engine, "inv([[1,2],[2,4]])");
    expectError(engine, "inv([[1,2,3],[4,5,6],[7,8,9]])");
}

void testTransposeAndTrace(Cas::Engine &engine) {
    std::printf("\n== 转置与迹 ==\n");
    expect(engine, "transpose([[1,2,3],[4,5,6]])", "[[1, 4], [2, 5], [3, 6]]");
    expect(engine, "trans([[1,2],[3,4]])", "[[1, 3], [2, 4]]");
    expect(engine, "transpose([[7]])", "[[7]]");
    expect(engine, "trace([[1,2],[3,4]])", "5");
    expect(engine, "tr([[1,2],[3,4]])", "5");
    expect(engine, "trace([[a,b],[c,d]])", "a + d");
    expect(engine, "trace(identity(4))", "4");
    // 迹只对方阵有定义
    expectError(engine, "trace([[1,2,3],[4,5,6]])");
}

void testMultiplication(Cas::Engine &engine) {
    std::printf("\n== 矩阵乘法 ==\n");
    expect(engine, "matmul([[1,2],[3,4]],[[5,6],[7,8]])", "[[19, 22], [43, 50]]");
    // 与单位阵相乘不变
    expect(engine, "matmul([[1,2],[3,4]],identity(2))", "[[1, 2], [3, 4]]");
    expect(engine, "matmul(identity(2),[[1,2],[3,4]])", "[[1, 2], [3, 4]]");
    // 非方阵：(2×3)·(3×2) → 2×2
    expect(engine, "matmul([[1,2,3],[4,5,6]],[[1,0],[0,1],[1,1]])", "[[4, 5], [10, 11]]");
    // 维数不匹配必须报错
    expectError(engine, "matmul([[1,2],[3,4]],[[1,2,3]])");
}

void testReducedRowEchelon(Cas::Engine &engine) {
    std::printf("\n== 行最简形与秩 ==\n");
    expect(engine, "rref([[1,2,3],[4,5,6]])", "[[1, 0, -1], [0, 1, 2]]");
    expect(engine, "rref([[1,1,1],[2,2,2]])", "[[1, 1, 1], [0, 0, 0]]");
    expect(engine, "rref([[0,1],[1,0]])", "[[1, 0], [0, 1]]");
    expect(engine, "rank([[1,2],[2,4]])", "1");
    expect(engine, "rank([[1,2],[3,4]])", "2");
    expect(engine, "rank([[1,2,3],[4,5,6],[7,8,9]])", "2");
    expect(engine, "rank(identity(3))", "3");
}

void testCharacteristicPolynomial(Cas::Engine &engine) {
    std::printf("\n== 特征多项式 ==\n");
    expect(engine, "charpoly([[2,0],[0,1]])", "x^2 - 3*x + 2");
    expect(engine, "charpoly([[1,2],[3,4]])", "x^2 - 5*x - 2");
    // 三角矩阵的特征多项式 = (x-a11)(x-a22)
    expect(engine, "charpoly([[3,1],[0,5]])", "x^2 - 8*x + 15");
    expect(engine, "charpoly(identity(2))", "x^2 - 2*x + 1");
}

void testEigenvalues(Cas::Engine &engine) {
    std::printf("\n== 特征值 ==\n");
    expect(engine, "eigenvals([[2,0],[0,1]])", "{1, 2}");
    expect(engine, "eigenvalues([[1,2],[3,4]])", "{1/2*(sqrt(33) + 5), 1/2*(-sqrt(33) + 5)}");
    expect(engine, "eigenvals([[1,0],[0,1]])", "{1}");
    // 无实重根：判别式 < 0，本实现不给复数解
    expectError(engine, "eigenvals([[0,1],[-1,0]])");
}

void testIdentity(Cas::Engine &engine) {
    std::printf("\n== 单位阵 ==\n");
    expect(engine, "identity(1)", "[[1]]");
    expect(engine, "identity(2)", "[[1, 0], [0, 1]]");
    expect(engine, "identity(3)", "[[1, 0, 0], [0, 1, 0], [0, 0, 1]]");
    expect(engine, "eye(2)", "[[1, 0], [0, 1]]");
}

void testNestedCommands(Cas::Engine &engine) {
    std::printf("\n== 嵌套调用 ==\n");
    // 矩阵作为中间结果参与外层标量计算
    expect(engine, "det(inv([[1,2],[3,4]]))", "-1/2");
    expect(engine, "trace(matmul([[1,2],[3,4]],inv([[1,2],[3,4]])))", "2");
    expect(engine, "det(matmul(inv([[1,2],[3,4]]),[[1,2],[3,4]]))", "1");
    expect(engine, "det(matmul([[1,2],[3,4]],[[1,2],[3,4]]))", "4");
    // 矩阵作为最终结果，保持矩阵渲染
    expect(engine, "transpose(inv([[1,2],[3,4]]))", "[[-2, 3/2], [1, -1/2]]");
    expect(engine, "inv(transpose([[1,2],[3,4]]))", "[[-2, 3/2], [1, -1/2]]");
    // 内层失败要有明确诊断
    expectError(engine, "det(inv([[1,2],[2,4]]))");
}

void testIdentityProperties(Cas::Engine &engine) {
    std::printf("\n== 代数恒等式 ==\n");
    // det(AB) = det(A)·det(B)
    expect(engine, "det(matmul([[1,2],[3,4]],[[2,0],[1,3]]))", "-12");
    // det(A^-1) = 1/det(A)
    expect(engine, "inv([[2,0],[0,5]])", "[[1/2, 0], [0, 1/5]]");
    // (A^-1)^-1 = A
    expect(engine, "inv(inv([[1,2],[3,4]]))", "[[1, 2], [3, 4]]");
    // det(A^T) = det(A)
    expect(engine, "det(transpose([[1,2],[3,4]]))", "-2");
}

void testLatex(Cas::Engine &engine) {
    std::printf("\n== 矩阵 LaTeX ==\n");
    expectLatexContains(engine, "inv([[1,2],[3,4]])", "\\begin{pmatrix}");
    expectLatexContains(engine, "inv([[1,2],[3,4]])", "\\frac{3}{2}");
    expectLatexContains(engine, "transpose([[1,2,3],[4,5,6]])", "\\begin{pmatrix}1 & 4");
}

void testMalformedInput(Cas::Engine &engine) {
    std::printf("\n== 非法输入 ==\n");
    // 行列长度不一致
    expectError(engine, "det([[1,2],[3]])");
    // 非方阵做行列式
    expectError(engine, "det([[1,2,3],[4,5,6]])");
    // 非方阵求逆
    expectError(engine, "inv([[1,2,3],[4,5,6]])");
}

} // namespace

int main() {
    Cas::Engine engine;

    testDeterminant(engine);
    testInverse(engine);
    testTransposeAndTrace(engine);
    testMultiplication(engine);
    testReducedRowEchelon(engine);
    testCharacteristicPolynomial(engine);
    testEigenvalues(engine);
    testIdentity(engine);
    testNestedCommands(engine);
    testIdentityProperties(engine);
    testLatex(engine);
    testMalformedInput(engine);

    std::printf("\n------------------------------------\n");
    std::printf("共 %d 项，失败 %d 项\n", g_total, g_failures);
    return g_failures == 0 ? 0 : 1;
}

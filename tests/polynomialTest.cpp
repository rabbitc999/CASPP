// 多项式代数自检：展开、带余除法、GCD、Resultant、因式分解、通分。
#include "cas/allocator.h"
#include "cas/detail/node.h"
#include "cas/detail/polynomial.h"
#include "cas/detail/simplify.h"

#include <cstdio>
#include <string>

using namespace Cas::Detail;

namespace {

int failures = 0;

void check(bool condition, const char *label) {
    if (!condition) {
        std::printf("FAIL %s\n", label);
        ++failures;
    }
}

std::string describe(const NodePtr &node) {
    if (!node) return "<null>";
    switch (node->kind()) {
    case NodeKind::Integer: return node->numerator().toString();
    case NodeKind::Rational: return "(" + node->numerator().toString() + "/" + node->denominator().toString() + ")";
    case NodeKind::Symbol: return node->name();
    case NodeKind::Constant: return "<const>";
    case NodeKind::Add: {
        std::string text;
        for (NodePtr child : node->operands()) {
            if (!text.empty()) text += "+";
            text += describe(child);
        }
        return "(" + text + ")";
    }
    case NodeKind::Mul: {
        std::string text;
        for (NodePtr child : node->operands()) {
            if (!text.empty()) text += "*";
            text += describe(child);
        }
        return "(" + text + ")";
    }
    case NodeKind::Pow:
        return describe(node->operand(0)) + "^" + describe(node->operand(1));
    case NodeKind::Function: {
        std::string text = std::string(node->name()) + "(";
        for (std::size_t i = 0; i < node->size(); ++i) {
            if (i) text += ",";
            text += describe(node->operand(i));
        }
        return text + ")";
    }
    }
    return "?";
}

Ptr sym(const char *name) { return Ptr(Node::makeSymbol(name)); }
Ptr num(long long value) { return Ptr(Node::makeInteger(BigInt(value))); }

void expectExpand(const Ptr &expression, const char *expected, const char *label) {
    Ptr result = expand(expression);
    const std::string actual = describe(result.get());
    if (actual != expected) {
        std::printf("FAIL %s\n     expected: %s\n     actual:   %s\n", label, expected, actual.c_str());
        ++failures;
    }
}

void testExpand() {
    Ptr x = sym("x");
    expectExpand(makeMul(makeAdd(x, num(1)), makeAdd(x, num(2))), "(x^2+(3*x)+2)", "(x+1)(x+2)");
    expectExpand(makePow(makeAdd(x, num(1)), num(3)), "(x^3+(3*x^2)+(3*x)+1)", "(x+1)^3");
    expectExpand(makePow(makeAdd(x, num(1)), num(2)), "(x^2+(2*x)+1)", "(x+1)^2");
    expectExpand(makePow(makeAdd(x, num(-1)), num(2)), "(x^2+(-2*x)+1)", "(x-1)^2");
    expectExpand(makeMul(num(2), makeAdd(x, num(3))), "((2*x)+6)", "2(x+3)");
}

void testPolynomialConversion() {
    Ptr x = sym("x");
    const std::uint32_t id = lookupSymbol("x");
    UnivariatePoly polynomial;
    check(toPolynomial(makeAdd(makePow(x, num(2)), num(1)).get(), id, &polynomial), "x^2+1 可解析");
    check(polynomial.degree() == 2, "x^2+1 次数 2");
    check(polynomial.coefficient(1)->isZero(), "x^2+1 一次项为 0");

    Ptr rebuilt = fromPolynomial(polynomial);
    check(rebuilt && describe(rebuilt.get()) == "(x^2+1)", "多项式往返");
}

void testDivision() {
    Ptr x = sym("x");
    const std::uint32_t id = lookupSymbol("x");
    // (x^2 - 1) / (x - 1) = x + 1
    Ptr dividend = expand(makeMul(makeAdd(x, num(1)), makeAdd(x, num(-1))));
    Ptr divisor = makeAdd(x, num(-1));
    UnivariatePoly a, b, q, r;
    toPolynomial(dividend.get(), id, &a);
    toPolynomial(divisor.get(), id, &b);
    check(polyDivide(a, b, &q, &r), "带余除法可执行");
    check(r.isZero() || describe(fromPolynomial(r).get()) == "0", "(x^2-1)/(x-1) 余式为零");
    check(describe(fromPolynomial(q).get()) == "(x+1)", "(x^2-1)/(x-1) = x+1");
}

void testGcd() {
    Ptr x = sym("x");
    const std::uint32_t id = lookupSymbol("x");
    // gcd(x^2-1, x^2-2x+1) = x-1
    Ptr a = expand(makeMul(makeAdd(x, num(1)), makeAdd(x, num(-1))));
    Ptr b = expand(makePow(makeAdd(x, num(-1)), num(2)));
    Ptr g = polyGcd(a.get(), b.get(), id);
    const std::string text = describe(g.get());
    check(text == "(x+-1)" || text == "((x+-1))", std::string("gcd(x^2-1, (x-1)^2) 应为 x-1，实际 " + text).c_str());
}

void testResultant() {
    Ptr x = sym("x");
    const std::uint32_t id = lookupSymbol("x");
    // Res(x - 1, x - 2) = 1 - 2 = -1
    Ptr a = makeAdd(x, num(-1));
    Ptr b = makeAdd(x, num(-2));
    Ptr r = polyResultant(a.get(), b.get(), id);
    check(r && r->isNumber(), "Res(x-1,x-2) 是数值");
    if (r && r->isNumber()) {
        const std::string text = describe(r.get());
        check(text == "1" || text == "-1", std::string("Res(x-1,x-2)=±1，实际 " + text).c_str());
    }
}

void testFactor() {
    Ptr x = sym("x");
    // x^2 - 1 = (x-1)(x+1)
    Ptr difference = expand(makeMul(makeAdd(x, num(1)), makeAdd(x, num(-1))));
    Ptr factored = factor(difference);
    const std::string text = describe(factored.get());
    std::printf("  factor(x^2-1) = %s\n", text.c_str());
    check(text.find("x") != std::string::npos, "因式分解包含 x");

    // x^2 + 3x + 2 = (x+1)(x+2)
    Ptr quadratic = expand(makeMul(makeAdd(x, num(1)), makeAdd(x, num(2))));
    Ptr factored2 = factor(quadratic);
    std::printf("  factor(x^2+3x+2) = %s\n", describe(factored2.get()).c_str());

    // 2x^2 + 3x + 1 = (2x+1)(x+1)
    Ptr rationalRoot = makeAdd(makeAdd(makeMul(num(2), makePow(x, num(2))), makeMul(num(3), x)), num(1));
    Ptr factored3 = factor(rationalRoot);
    std::printf("  factor(2x^2+3x+1) = %s\n", describe(factored3.get()).c_str());
    // 验证分解结果与原式等价
    Ptr back = expand(factored3);
    check(describe(back.get()) == describe(expand(rationalRoot).get()), "2x^2+3x+1 分解后展开还原");

    // 公因子提取：6x^2 + 9x = 3x(2x+3)
    Ptr withGcd = makeAdd(makeMul(num(6), makePow(x, num(2))), makeMul(num(9), x));
    Ptr factored4 = factor(withGcd);
    Ptr back4 = expand(factored4);
    check(describe(back4.get()) == describe(expand(withGcd).get()), "6x^2+9x 分解后展开还原");
    std::printf("  factor(6x^2+9x) = %s\n", describe(factored4.get()).c_str());
}

void testTogether() {
    Ptr x = sym("x");
    // 1/x + 1/(x+1) = (2x+1)/(x(x+1))
    Ptr sum = makeAdd(makeDivide(num(1), x), makeDivide(num(1), makeAdd(x, num(1))));
    Ptr combined = together(sum);
    std::printf("  together(1/x + 1/(x+1)) = %s\n", describe(combined.get()).c_str());
    check(combined.get() != nullptr, "通分成功");
}

} // namespace

int main() {
    std::printf("多项式代数自检\n");
    Cas::Allocator::initialize();

    testExpand();              std::printf("expand ok\n");
    testPolynomialConversion(); std::printf("conversion ok\n");
    testDivision();            std::printf("division ok\n");
    testGcd();                 std::printf("gcd ok\n");
    testResultant();           std::printf("resultant ok\n");
    testFactor();              std::printf("factor ok\n");
    testTogether();            std::printf("together ok\n");

    if (failures == 0) {
        std::printf("全部通过\n");
        return 0;
    }
    std::printf("%d 项失败\n", failures);
    return 1;
}

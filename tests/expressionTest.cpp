// 表达式树与自动化简的自检：验证扁平化、同类项合并、数值折叠、恒等式消去。
#include "cas/allocator.h"
#include "cas/detail/node.h"
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

/** @brief 极简结构打印：便于失败时肉眼比对。 */
std::string describe(const NodePtr &node) {
    if (!node) return "<null>";
    switch (node->kind()) {
    case NodeKind::Integer:
        return node->numerator().toString();
    case NodeKind::Rational:
        return "(" + node->numerator().toString() + "/" + node->denominator().toString() + ")";
    case NodeKind::Symbol:
        return node->name();
    case NodeKind::Constant:
        switch (node->constant()) {
        case ConstantId::Pi: return "pi";
        case ConstantId::E: return "e";
        case ConstantId::EulerGamma: return "gamma";
        case ConstantId::Infinity: return "inf";
        case ConstantId::NegativeInfinity: return "-inf";
        }
        return "?";
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

/** @brief 断言化简结果的结构文本等于期望。 */
void expect(const Ptr &expression, const char *expected, const char *label) {
    Ptr simplified = simplify(expression);
    const std::string actual = describe(simplified.get());
    if (actual != expected) {
        std::printf("FAIL %s\n     expected: %s\n     actual:   %s\n", label, expected, actual.c_str());
        ++failures;
    }
}

Ptr sym(const char *name) { return Ptr(Node::makeSymbol(name)); }
Ptr num(long long value) { return Ptr(Node::makeInteger(BigInt(value))); }

void testNumbers() {
    expect(makeAdd(num(2), num(3)), "5", "2+3");
    expect(makeAdd(num(1), makeAdd(num(2), num(3))), "6", "1+(2+3) 展平");
    expect(makeMul(num(2), num(3)), "6", "2*3");
    expect(makeMul(num(2), makeMul(num(3), num(4))), "24", "2*(3*4) 展平");
    expect(makePow(num(2), num(10)), "1024", "2^10");
    expect(makePow(num(2), num(-3)), "(1/8)", "2^-3");
    expect(makePow(num(2), num(0)), "1", "2^0");
    expect(makePow(num(0), num(5)), "0", "0^5");
    expect(makeAdd(num(1), num(-1)), "0", "1+(-1)");
    expect(makeMul(num(0), sym("x")), "0", "0*x");
}

void testLikeTerms() {
    Ptr x = sym("x"), y = sym("y");
    expect(makeAdd(makeMul(num(2), x), makeMul(num(3), x)), "(5*x)", "2x+3x");
    expect(makeAdd(x, x), "(2*x)", "x+x");
    expect(makeAdd(x, makeMul(num(-1), x)), "0", "x-x");
    expect(makeMul(x, x), "x^2", "x*x");
    expect(makeMul(makePow(x, num(2)), x), "x^3", "x^2*x");
    expect(makeMul(makePow(x, num(2)), makePow(x, num(3))), "x^5", "x^2*x^3");
    expect(makeMul(num(2), makeMul(num(3), x)), "(6*x)", "2*(3x)");
}

void testIdentities() {
    Ptr x = sym("x");
    expect(makeMul(num(1), x), "x", "1*x");
    expect(makeAdd(num(0), x), "x", "x+0");
    expect(makePow(x, num(1)), "x", "x^1");
    expect(makePow(makePow(x, num(2)), num(3)), "x^6", "(x^2)^3");
    expect(makePow(makeMul(sym("a"), sym("b")), num(2)), "(a^2*b^2)", "(a*b)^2");
    expect(makeDivide(makeMul(num(4), x), num(2)), "(2*x)", "4x/2");
    expect(makeNegate(makeNegate(x)), "x", "-(-x)");
}

void testOrdering() {
    Ptr a = sym("a"), b = sym("b");
    Ptr left = makeAdd(a, b), right = makeAdd(b, a);
    check(left->equals(right.get()), "a+b 与 b+a 结构相同");
    Ptr mleft = makeMul(a, b), mright = makeMul(b, a);
    check(mleft->equals(mright.get()), "a*b 与 b*a 结构相同");
}

void testNesting() {
    Ptr x = sym("x");
    Ptr squared = makePow(makeAdd(x, num(1)), num(2));
    Ptr again = simplify(squared);
    check(squared->equals(again.get()), "化简幂等");
    Ptr twice = simplify(again);
    check(again->equals(twice.get()), "二次化简稳定");
}

void testQuery() {
    Ptr x = sym("x"), y = sym("y");
    Ptr expression = makeAdd(makeMul(num(2), x), makeMul(num(3), y));
    const auto symbols = freeSymbols(expression.get());
    check(symbols.size() == 2, "自由符号 2 个");
    check(containsSymbol(expression.get(), lookupSymbol("x")), "含 x");
    check(!containsSymbol(expression.get(), lookupSymbol("z")), "不含 z");
    check(!containsInfinity(expression.get()), "不含无穷");
    Ptr infinity(Node::makeConstant(ConstantId::Infinity));
    check(containsInfinity(infinity.get()), "含无穷");
}

void testPoolAccounting() {
    Cas::Allocator::initialize();
    {
        Ptr x = sym("x");
        Ptr expression = makePow(makeAdd(x, num(1)), num(2));
        Ptr simplified = simplify(expression);
        check(simplified.get() != nullptr, "复杂表达式化简成功");
    }
    const auto stats = Cas::Allocator::statistics();
    check(stats.liveBytes < 4 * 1024 * 1024, "无大额内存泄漏");
}

} // namespace

int main() {
    std::printf("表达式树自检\n");
    Cas::Allocator::initialize();

    testNumbers();   std::printf("numbers ok\n");
    testLikeTerms(); std::printf("liketerms ok\n");
    testIdentities(); std::printf("identities ok\n");
    testOrdering();  std::printf("ordering ok\n");
    testNesting();   std::printf("nesting ok\n");
    testQuery();     std::printf("query ok\n");
    testPoolAccounting(); std::printf("pool ok\n");

    if (failures == 0) {
        std::printf("全部通过\n");
        return 0;
    }
    std::printf("%d 项失败\n", failures);
    return 1;
}

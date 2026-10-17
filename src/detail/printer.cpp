#include "cas/detail/printer.h"
#include "cas/detail/simplify.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace Cas::Detail {
namespace {

/** @brief 运算符优先级：用于决定是否补括号。 */
enum class Precedence { Add = 1, Mul = 2, Pow = 3, Atom = 4 };

/** @brief 取节点在书写时的优先级。 */
Precedence precedenceOf(const NodePtr &node) {
    if (!node) return Precedence::Atom;
    switch (node->kind()) {
    case NodeKind::Add: return Precedence::Add;
    case NodeKind::Mul: return Precedence::Mul;
    case NodeKind::Pow: return Precedence::Pow;
    default: return Precedence::Atom;
    }
}

/** @brief 是否为负数值。 */
bool isNegativeNumber(const NodePtr &node) { return node && node->isNumber() && node->isNegative(); }

/** @brief 取数值的相反数文本（去掉开头的负号，用作减法右项）。 */
Ptr absoluteNumber(const NodePtr &node) {
    if (!node || !node->isNumber()) return Ptr();
    BigInt numerator = node->numerator();
    if (numerator.isNegative()) numerator = numerator.negated();
    return makeNumber(numerator, node->denominator());
}

/* ---------------------------------------------------------- 文本输出 */

std::string textOf(const NodePtr &node, const PrintOptions &options, Precedence parent);

/** @brief 幂指数加括号的条件：指数不是单个原子。 */
std::string exponentText(const Ptr &exponent, const PrintOptions &options) {
    const std::string inner = textOf(exponent.get(), options, Precedence::Atom);
    if (exponent->isNumber() && exponent->isInteger()) return inner;
    if (exponent->kind() == NodeKind::Symbol || exponent->kind() == NodeKind::Constant) return inner;
    return "(" + inner + ")";
}

/** @brief 文本是否已被整对圆括号包住。 */
bool wrapped(const std::string &text) {
    if (text.size() < 2 || text.front() != '(' || text.back() != ')') return false;
    int depth = 0;
    for (std::size_t i = 0; i < text.size(); ++i) {
        if (text[i] == '(') ++depth;
        else if (text[i] == ')') {
            --depth;
            if (depth == 0) return i == text.size() - 1; // 第一层就闭合，说明整串被包住
        }
    }
    return false;
}

/** @brief 分式文本：分子/分母。 */
std::string fractionText(const Ptr &numerator, const Ptr &denominator, const PrintOptions &options) {
    const std::string top = textOf(numerator.get(), options, Precedence::Pow);
    const std::string bottom = textOf(denominator.get(), options, Precedence::Pow);
    const bool simpleTop = numerator->kind() != NodeKind::Add && numerator->kind() != NodeKind::Mul;
    const bool simpleBottom = denominator->kind() != NodeKind::Add && denominator->kind() != NodeKind::Mul;
    std::string result;
    result += (simpleTop || wrapped(top)) ? top : "(" + top + ")";
    result += "/";
    result += (simpleBottom || wrapped(bottom)) ? bottom : "(" + bottom + ")";
    return result;
}

std::string textOf(const NodePtr &node, const PrintOptions &options, Precedence parent) {
    if (!node) return "<null>";

    switch (node->kind()) {
    case NodeKind::Integer:
        return node->numerator().toString();

    case NodeKind::Rational: {
        // 直接输出 p/q，不再走递归（分子是自身会无限递归）
        const Ptr numerator(makeNumber(node->numerator(), BigInt(1)));
        const Ptr denominator(makeNumber(node->denominator(), BigInt(1)));
        return fractionText(numerator, denominator, options);
    }

    case NodeKind::Symbol:
        return node->name();

    case NodeKind::Constant:
        switch (node->constant()) {
        case ConstantId::Pi: return options.useUnicode ? "\u03c0" : "pi";
        case ConstantId::E: return "e";
        case ConstantId::EulerGamma: return options.useUnicode ? "\u03b3" : "gamma";
        case ConstantId::Infinity: return options.useUnicode ? "\u221e" : "inf";
        case ConstantId::NegativeInfinity: return options.useUnicode ? "-\u221e" : "-inf";
        }
        return "?";

    case NodeKind::Add: {
        std::string text;
        for (std::size_t i = 0; i < node->size(); ++i) {
            NodePtr child = node->operand(i);
            if (i == 0) {
                text += textOf(child, options, Precedence::Add);
                continue;
            }
            if (isNegativeNumber(child)) {
                // a + (-b) 写成 a - b
                text += options.useUnicode ? " - " : " - ";
                text += textOf(absoluteNumber(child).get(), options, Precedence::Mul);
                continue;
            }
            if (child->kind() == NodeKind::Mul && child->size() >= 2 && isNegativeNumber(child->operand(0))) {
                // a + (-2*x) 写成 a - 2*x；系数为 -1 时丢弃它，得到 a - x
                Ptr positive = absoluteNumber(child->operand(0));
                std::string body;
                if (!positive->isOne()) body += textOf(positive.get(), options, Precedence::Mul) + "*";
                for (std::size_t k = 1; k < child->size(); ++k) {
                    NodePtr factor = child->operand(k);
                    if (k > 1) body += "*";
                    if (precedenceOf(factor) == Precedence::Add) {
                        body += "(" + textOf(factor, options, Precedence::Add) + ")";
                    } else {
                        body += textOf(factor, options, Precedence::Mul);
                    }
                }
                text += " - ";
                text += body;
                continue;
            }
            text += " + ";
            text += textOf(child, options, Precedence::Mul);
        }
        if (parent == Precedence::Mul || parent == Precedence::Pow) return "(" + text + ")";
        return text;
    }

    case NodeKind::Mul: {
        // 抽成单一分式：形如 (1/2)*x^(-1/2) 输出 1/(2*sqrt(x))，比 1/2*1/sqrt(x) 自然
        {
            // 有理系数拆成 分子/分母；带负指数的因子把底数（正指数）搬进分母
            bool sign = false;
            BigInt coefNumerator(1);
            BigInt coefDenominator(1);
            std::vector<std::string> topFactors;
            std::vector<std::string> bottomFactors;
            bool hasNegativePower = false;
            bool okay = true;

            for (std::size_t i = 0; i < node->size() && okay; ++i) {
                NodePtr child = node->operand(i);
                if (child->isNumber()) {
                    coefNumerator = coefNumerator * child->numerator();
                    coefDenominator = coefDenominator * child->denominator();
                    if (coefNumerator.isNegative()) {
                        coefNumerator = coefNumerator.negated();
                        sign = !sign;
                    }
                    continue;
                }
                if (child->kind() == NodeKind::Pow && child->operand(1)->isNumber() &&
                    child->operand(1)->isNegative()) {
                    Ptr positive = makeNegate(Ptr::borrow(child->operand(1)));
                    Ptr flipped = makePow(Ptr::borrow(child->operand(0)), positive);
                    bottomFactors.push_back(textOf(flipped.get(), options, Precedence::Mul));
                    hasNegativePower = true;
                    continue;
                }
                // 加法因子（如分子 2x+1）作为整体进分子，渲染时自带括号
                if (precedenceOf(child) == Precedence::Add) {
                    topFactors.push_back(textOf(child, options, Precedence::Mul));
                    continue;
                }
                topFactors.push_back(textOf(child, options, Precedence::Mul));
            }

            // 只有存在负指数因子时才聚合成分式；`1/2*x` 之类保持原样
            if (okay && hasNegativePower) {
                std::string result = sign ? "-" : "";
                std::string top;
                for (const std::string &piece : topFactors) {
                    if (!top.empty()) top += "*";
                    top += piece;
                }
                std::string bottom;
                for (const std::string &piece : bottomFactors) {
                    if (!bottom.empty()) bottom += "*";
                    bottom += piece;
                }

                // 分子 = 系数分子 × 各正指数因子
                std::string topText;
                if (!coefNumerator.isOne()) topText += coefNumerator.toString();
                if (!top.empty()) {
                    if (!topText.empty()) topText += "*";
                    topText += top;
                }
                if (topText.empty()) topText = "1";

                // 分母 = 系数分母 × 各负指数因子
                std::string bottomText;
                if (!coefDenominator.isOne()) bottomText += coefDenominator.toString();
                if (!bottom.empty()) {
                    if (!bottomText.empty()) bottomText += "*";
                    bottomText += bottom;
                }

                const bool wrapTop =
                    !wrapped(topText) && (topText.find('*') != std::string::npos || topText.find(' ') != std::string::npos);
                const bool wrapBottom = !wrapped(bottomText) && (bottomText.find('*') != std::string::npos ||
                                                                bottomText.find(' ') != std::string::npos);
                result += wrapTop ? "(" + topText + ")" : topText;
                if (!bottomText.empty()) {
                    result += "/";
                    result += wrapBottom ? "(" + bottomText + ")" : bottomText;
                }
                if (parent == Precedence::Pow) return "(" + result + ")";
                return result;
            }
        }

        // 首元素是负系数时提出负号，写成 `-2*x` 而不是 `-1*2*x`
        std::string text;
        std::size_t first = 0;
        if (node->size() >= 1) {
            NodePtr head = node->operand(0);
            if (head->isNumber() && head->numerator().isNegative()) {
                if (head->isInteger() && head->numerator().abs().isOne()) {
                    text += "-";
                    first = 1; // 单位负号：直接丢弃 -1
                } else {
                    text += "-";
                    Ptr magnitude = makeNegate(Ptr::borrow(head));
                    if (node->size() == 1) return text + textOf(magnitude.get(), options, parent);
                    text += textOf(magnitude.get(), options, Precedence::Mul) + "*";
                    first = 1;
                }
            }
        }
        for (std::size_t i = first; i < node->size(); ++i) {
            NodePtr child = node->operand(i);
            const bool needsParen = precedenceOf(child) == Precedence::Add;
            if (i > first) text += "*";
            if (needsParen) {
                text += "(" + textOf(child, options, Precedence::Add) + ")";
            } else {
                text += textOf(child, options, Precedence::Mul);
            }
        }
        if (first == node->size()) return text.empty() ? "1" : "-1"; // 全部被提走
        if (parent == Precedence::Pow) return "(" + text + ")";
        return text;
    }

    case NodeKind::Pow: {
        const std::string base = textOf(node->operand(0), options, Precedence::Pow);
        Ptr exponent(Ptr::borrow(node->operand(1)).share());
        // 负指数写成 1/x^n 更贴近习惯（含负有理指数，如 x^(-1/2) → 1/sqrt(x)）
        if (exponent->isNumber() && exponent->isNegative()) {
            Ptr positive = makeNegate(exponent);
            Ptr reciprocal = makePow(Ptr::borrow(node->operand(0)), positive);
            return fractionText(Ptr(makeNumber(BigInt(1), BigInt(1))), reciprocal, options);
        }
        // 开方：x^(1/2) → sqrt(x)
        if (exponent->isNumber() && !exponent->isInteger() && exponent->numerator().isOne() &&
            exponent->denominator() == BigInt(2)) {
            // `sqrt(...)` 自带一层括号，而 base 已按 Pow 优先级求值
            // （加法子式在那里就加过括号了），这里不能再补一层，否则
            // 会打印成 sqrt(((-x^2 + 1)))。
            return "sqrt(" + base + ")";
        }
        return base + "^" + exponentText(exponent, options);
    }

    case NodeKind::Function: {
        std::string text = std::string(node->name()) + "(";
        for (std::size_t i = 0; i < node->size(); ++i) {
            if (i) text += ", ";
            text += textOf(node->operand(i), options, Precedence::Add);
        }
        return text + ")";
    }
    }
    return "?";
}

/* ---------------------------------------------------------- LaTeX 输出 */

std::string latexOf(const NodePtr &node, const PrintOptions &options, Precedence parent);

/** @brief LaTeX 分式。 */
std::string latexFraction(const Ptr &numerator, const Ptr &denominator, const PrintOptions &options) {
    return "\\frac{" + latexOf(numerator.get(), options, Precedence::Add) + "}{" +
           latexOf(denominator.get(), options, Precedence::Add) + "}";
}

/** @brief LaTeX 幂。 */
std::string latexPower(const Ptr &base, const Ptr &exponent, const PrintOptions &options) {
    const std::string baseText = latexOf(base.get(), options, Precedence::Pow);
    // 负整数指数写成 \frac{1}{x^n}
    if (exponent->isNumber() && exponent->isInteger() && exponent->isNegative()) {
        Ptr positive = makeNegate(exponent);
        return "\\frac{1}{" + latexPower(base, positive, options) + "}";
    }
    // 分数指数写成 \sqrt[n]{x} 或 \sqrt{x}
    return baseText + "^{" + latexOf(exponent.get(), options, Precedence::Atom) + "}";
}

std::string latexOf(const NodePtr &node, const PrintOptions &options, Precedence parent) {
    if (!node) return "<null>";

    switch (node->kind()) {
    case NodeKind::Integer:
        return node->numerator().toString();

    case NodeKind::Rational: {
        const Ptr numerator(makeNumber(node->numerator(), BigInt(1)));
        const Ptr denominator(makeNumber(node->denominator(), BigInt(1)));
        return latexFraction(numerator, denominator, options);
    }

    case NodeKind::Symbol: {
        // 多字符变量名用 \mathrm 保持正体，避免被当成多个变量相乘
        const std::string name = node->name();
        if (name.size() > 1) return "\\mathrm{" + name + "}";
        return name;
    }

    case NodeKind::Constant:
        switch (node->constant()) {
        case ConstantId::Pi: return "\\pi";
        case ConstantId::E: return "e";
        case ConstantId::EulerGamma: return "\\gamma";
        case ConstantId::Infinity: return "\\infty";
        case ConstantId::NegativeInfinity: return "-\\infty";
        }
        return "?";

    case NodeKind::Add: {
        std::string text;
        for (std::size_t i = 0; i < node->size(); ++i) {
            NodePtr child = node->operand(i);
            if (i == 0) {
                text += latexOf(child, options, Precedence::Add);
                continue;
            }
            if (isNegativeNumber(child)) {
                text += " - " + latexOf(absoluteNumber(child).get(), options, Precedence::Mul);
                continue;
            }
            if (child->kind() == NodeKind::Mul && child->size() >= 2 &&
                isNegativeNumber(child->operand(0))) {
                // a + (-2*x) 写成 a - 2 x；系数为 -1 时丢弃它
                Ptr positive = absoluteNumber(child->operand(0));
                std::string body;
                if (!positive->isOne()) body += latexOf(positive.get(), options, Precedence::Mul) + " ";
                for (std::size_t k = 1; k < child->size(); ++k) {
                    NodePtr factor = child->operand(k);
                    if (k > 1) {
                        const NodePtr previous = child->operand(k - 1);
                        body += (previous->isNumber() && factor->isNumber()) ? " \\cdot " : " ";
                    }
                    if (precedenceOf(factor) == Precedence::Add) {
                        body += "\\left(" + latexOf(factor, options, Precedence::Add) + "\\right)";
                    } else {
                        body += latexOf(factor, options, Precedence::Mul);
                    }
                }
                text += " - " + body;
                continue;
            }
            text += " + " + latexOf(child, options, Precedence::Mul);
        }
        if (parent == Precedence::Mul || parent == Precedence::Pow) return "\\left(" + text + "\\right)";
        return text;
    }

    case NodeKind::Mul: {
        std::string text;
        std::size_t first = 0;
        if (node->size() >= 1) {
            NodePtr head = node->operand(0);
            if (head->isNumber() && head->numerator().isNegative()) {
                text += "-";
                if (head->isInteger() && head->numerator().abs().isOne()) {
                    first = 1;
                } else {
                    Ptr magnitude = makeNegate(Ptr::borrow(head));
                    if (node->size() == 1) return text + latexOf(magnitude.get(), options, parent);
                    text += latexOf(magnitude.get(), options, Precedence::Mul) + " ";
                    first = 1;
                }
            }
        }
        for (std::size_t i = first; i < node->size(); ++i) {
            NodePtr child = node->operand(i);
            if (i > first) {
                // 数字与字母之间省略乘号；两个数字之间用 \cdot
                const NodePtr previous = node->operand(i - 1);
                const bool bothNumbers = previous->isNumber() && child->isNumber();
                text += bothNumbers ? " \\cdot " : " ";
            }
            if (precedenceOf(child) == Precedence::Add) {
                text += "\\left(" + latexOf(child, options, Precedence::Add) + "\\right)";
            } else {
                text += latexOf(child, options, Precedence::Mul);
            }
        }
        if (first == node->size()) return text.empty() ? "1" : "-1";
        if (parent == Precedence::Pow) return "\\left(" + text + "\\right)";
        return text;
    }

    case NodeKind::Pow: {
        Ptr base = Ptr::borrow(node->operand(0));
        Ptr exponent = Ptr::borrow(node->operand(1));
        // 负指数：x^(-n) → \frac{1}{x^{n}}，x^(-1/2) → \frac{1}{\sqrt{x}}
        if (exponent->isNumber() && exponent->isNegative()) {
            Ptr positive = makeNegate(exponent);
            return latexFraction(Ptr(makeNumber(BigInt(1), BigInt(1))),
                                 makePow(base, positive), options);
        }
        // 开方：x^(1/2) → \sqrt{x}
        if (exponent->isNumber() && !exponent->isInteger()) {
            const BigInt denominator = exponent->denominator();
            if (exponent->numerator().isOne()) {
                if (denominator == BigInt(2)) {
                    return "\\sqrt{" + latexOf(base.get(), options, Precedence::Add) + "}";
                }
                return "\\sqrt[" + denominator.toString() + "]{" +
                       latexOf(base.get(), options, Precedence::Add) + "}";
            }
        }
        return latexPower(base, exponent, options);
    }

    case NodeKind::Function: {
        const std::string name = node->name();
        static const char *kKnown[] = {"sin", "cos", "tan", "cot", "sec", "csc", "arcsin", "arccos",
                                       "arctan", "sinh", "cosh", "tanh", "log", "ln", "exp", "abs"};
        bool known = false;
        for (const char *candidate : kKnown) {
            if (name == candidate) {
                known = true;
                break;
            }
        }
        std::string text;
        if (known) {
            text = "\\" + name;
            // 三角函数直接接括号；log/exp 等用 \left(\right)
            text += "\\left(";
        } else {
            text = "\\operatorname{" + name + "}\\left(";
        }
        for (std::size_t i = 0; i < node->size(); ++i) {
            if (i) text += ", ";
            text += latexOf(node->operand(i), options, Precedence::Add);
        }
        return text + "\\right)";
    }
    }
    return "?";
}

} // namespace

std::string numberToText(const NodePtr &number) {
    if (!number || !number->isNumber()) return "<not a number>";
    if (number->denominator().isOne()) return number->numerator().toString();
    return number->numerator().toString() + "/" + number->denominator().toString();
}

std::string toText(const NodePtr &expression, const PrintOptions &options) {
    return textOf(expression, options, Precedence::Add);
}

std::string toLatex(const NodePtr &expression, const PrintOptions &options) {
    return latexOf(expression, options, Precedence::Add);
}

} // namespace Cas::Detail

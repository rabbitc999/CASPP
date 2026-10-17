// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#include "cas/detail/parser.h"
#include "cas/detail/trigonometry.h"

#include <cctype>
#include <cmath>
#include <cstring>

namespace Cas::Detail {
namespace {

/** @brief 快捷构造整数。 */
Ptr integer(long long value) { return Ptr(Node::makeInteger(BigInt(value))); }

/** @brief 词法单元种类。 */
enum class TokenKind {
    End,
    Number,
    Identifier,
    Operator, // + - * / ^ ! ' = ,
    LeftParen,
    RightParen,
    LeftBracket,
    RightBracket,
    Bar,
};

/** @brief 词法单元。 */
struct Token {
    TokenKind kind = TokenKind::End;
    std::string text;   ///< 标识符名字/运算符文本/数字原文
    std::size_t start = 0;
    std::size_t end = 0;
    bool isInteger = true; ///< 数字是否为整数
    BigInt integerValue;   ///< 整数值
    BigInt floatNumerator;  ///< 小数转成的有理数分子
    BigInt floatDenominator;///< 小数转成的有理数分母
};

/** @brief 递归下降解析器。 */
class Parser {
public:
    Parser(const std::string &text, bool allowUnknownCalls = false)
        : m_text(text), m_allowUnknownCalls(allowUnknownCalls) {
        m_position = 0;
    }

    Parser(const std::string &text, const ParseOptions &options, bool allowUnknownCalls = false)
        : m_text(text), m_allowUnknownCalls(allowUnknownCalls), m_options(options) {
        m_position = 0;
    }

    /** @brief 解析整个表达式。 */
    ParseResult run() {
        ParseResult result;
        advance();
        if (m_token.kind == TokenKind::End) {
            result.message = "表达式为空";
            return result;
        }

        Ptr expression = parseAssignment();
        if (!expression) {
            result.message = m_error.empty() ? "无法解析表达式" : m_error;
            result.position = m_errorPosition;
            return result;
        }
        if (m_token.kind != TokenKind::End) {
            result.message = "多余的输入：'" + m_token.text + "'";
            result.position = m_token.start;
            return result;
        }
        result.expression = expression;
        result.ok = true;
        result.length = m_text.size();
        return result;
    }

private:
    /* ------------------------------------------------------ 词法 */

    void skipSpace() {
        while (m_position < m_text.size() && std::isspace(static_cast<unsigned char>(m_text[m_position]))) {
            ++m_position;
        }
    }

    void advance() {
        skipSpace();
        m_token = Token();
        m_token.start = m_position;
        if (m_position >= m_text.size()) {
            m_token.kind = TokenKind::End;
            m_token.end = m_position;
            return;
        }

        const char c = m_text[m_position];

        // 区间运算符 `..`（sum(i^2, i=1..n) 的上下限分隔）要先于数字识别，
        // 否则 `.` 会被当成小数点送进 readNumber，遇到 `..` 又没有字符可读，
        // 结果既不产出 token 也不推进位置——死循环。
        if (c == '.' && m_position + 1 < m_text.size() && m_text[m_position + 1] == '.') {
            m_token.kind = TokenKind::Operator;
            m_token.text = "..";
            m_position += 2;
            m_token.end = m_position;
            return;
        }
        if (std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
            readNumber();
            return;
        }
        if (std::isalpha(static_cast<unsigned char>(c))) {
            readIdentifier();
            return;
        }
        // UTF-8 多字节：√ 与 π 需要整体识别
        if (static_cast<unsigned char>(c) >= 0x80) {
            readUtf8Symbol();
            return;
        }

        switch (c) {
        case '(': m_token.kind = TokenKind::LeftParen; m_token.text = "("; m_position += 1; break;
        case ')': m_token.kind = TokenKind::RightParen; m_token.text = ")"; m_position += 1; break;
        case '[': m_token.kind = TokenKind::LeftBracket; m_token.text = "["; m_position += 1; break;
        case ']': m_token.kind = TokenKind::RightBracket; m_token.text = "]"; m_position += 1; break;
        case '|': m_token.kind = TokenKind::Bar; m_token.text = "|"; m_position += 1; break;
        case '+': case '-': case '*': case '/': case '^': case '!': case '\'': case ',':
            m_token.kind = TokenKind::Operator;
            m_token.text = std::string(1, c);
            m_position += 1;
            // ** 与 := 当作单 token
            if (c == '*' && m_position < m_text.size() && m_text[m_position] == '*') {
                m_token.text = "^";
                m_position += 1;
            } else if (c == ':' && m_position < m_text.size() && m_text[m_position] == '=') {
                m_token.text = "=";
                m_position += 1;
            }
            break;
        case '=':
            m_token.kind = TokenKind::Operator;
            m_token.text = "=";
            m_position += 1;
            if (m_position < m_text.size() && m_text[m_position] == '=') ++m_position;
            break;
        case ':':
            // 单独的 ':' 出现在 := 之外时按错误处理
            m_token.kind = TokenKind::Operator;
            m_token.text = ":";
            m_position += 1;
            break;
        case '.':
            // 走到这里说明是孤立的 `.`（`..` 已在函数开头处理）
            m_token.kind = TokenKind::Operator;
            m_token.text = ".";
            m_position += 1;
            fail("无法识别的字符 '.'");
            break;
        default:
            m_token.kind = TokenKind::Operator;
            m_token.text = std::string(1, c);
            m_position += 1;
            fail("无法识别的字符");
            break;
        }
        m_token.end = m_position;
    }

    void readNumber() {
        const std::size_t begin = m_position;
        bool seenDot = false, seenExponent = false;
        while (m_position < m_text.size()) {
            const char c = m_text[m_position];
            if (std::isdigit(static_cast<unsigned char>(c))) { ++m_position; continue; }
            // `..` 是区间运算符（`1..n`），不是小数点，禁止在这里吞掉它
            if (c == '.' && !seenDot && !seenExponent) {
                if (m_position + 1 < m_text.size() && m_text[m_position + 1] == '.') break;
                seenDot = true;
                ++m_position;
                continue;
            }
            if ((c == 'e' || c == 'E') && !seenExponent &&
                m_position + 1 < m_text.size() &&
                (std::isdigit(static_cast<unsigned char>(m_text[m_position + 1])) ||
                 ((m_text[m_position + 1] == '+' || m_text[m_position + 1] == '-') &&
                  m_position + 2 < m_text.size() &&
                  std::isdigit(static_cast<unsigned char>(m_text[m_position + 2]))))) {
                seenExponent = true;
                m_position += 1;
                if (m_text[m_position] == '+' || m_text[m_position] == '-') ++m_position;
                continue;
            }
            break;
        }
        m_token.kind = TokenKind::Number;
        m_token.text = m_text.substr(begin, m_position - begin);
        m_token.end = m_position;
        m_token.isInteger = !seenDot && !seenExponent;

        if (m_token.isInteger) {
            m_token.integerValue = BigInt(m_token.text, 10);
            return;
        }

        // 小数按精确十进制转有理数：0.125 → 1/8，1.2e-3 → 3/2500
        std::string digits;
        int exponentShift = 0;
        int fractionalDigits = 0;
        bool afterDot = false;
        std::size_t cursor = 0;
        for (; cursor < m_token.text.size(); ++cursor) {
            const char c = m_token.text[cursor];
            if (std::isdigit(static_cast<unsigned char>(c))) {
                digits.push_back(c);
                if (afterDot) ++fractionalDigits;
            } else if (c == '.') {
                afterDot = true;
            } else if (c == 'e' || c == 'E') {
                exponentShift = std::stoi(m_token.text.substr(cursor + 1));
                break;
            }
        }
        if (digits.empty()) digits = "0";
        BigInt numerator(digits, 10);
        BigInt denominator(1);
        const long long scale = static_cast<long long>(fractionalDigits) - exponentShift;
        if (scale > 0) denominator = BigInt::power(BigInt(10), static_cast<unsigned long long>(scale));
        else if (scale < 0) numerator = numerator * BigInt::power(BigInt(10), static_cast<unsigned long long>(-scale));
        m_token.floatNumerator = numerator;
        m_token.floatDenominator = denominator;
    }

    void readIdentifier() {
        const std::size_t begin = m_position;
        while (m_position < m_text.size()) {
            const unsigned char c = static_cast<unsigned char>(m_text[m_position]);
            if (std::isalnum(c) || m_text[m_position] == '_' || c >= 0x80) {
                ++m_position;
            } else if (m_text[m_position] == '\'' && m_position == begin) {
                break;
            } else {
                break;
            }
        }
        m_token.kind = TokenKind::Identifier;
        m_token.text = m_text.substr(begin, m_position - begin);
        m_token.end = m_position;
    }

    void readUtf8Symbol() {
        const std::size_t begin = m_position;
        // 识别常见的单字符 UTF-8 记号
        static const struct {
            const char *utf8;
            std::size_t length;
            const char *ascii;
        } kSymbols[] = {
            {"\xE2\x88\x9A", 3, "sqrt"}, // √
            {"\xCF\x80", 2, "pi"},       // π
            {"\xE2\x88\x9E", 3, "infinity"}, // ∞
            {"\xC3\xB7", 2, "/"},        // ÷
            {"\xC3\x97", 2, "*"},        // ×
            {"\xE2\x88\x92", 3, "-"},    // −
        };
        for (const auto &symbol : kSymbols) {
            if (m_text.compare(begin, symbol.length, symbol.utf8) == 0) {
                m_position += symbol.length;
                m_token.kind = TokenKind::Identifier;
                m_token.text = symbol.ascii;
                m_token.end = m_position;
                return;
            }
        }
        m_position += 1;
        m_token.kind = TokenKind::Identifier;
        m_token.text = m_text.substr(begin, m_position - begin);
        m_token.end = m_position;
    }

    /* ------------------------------------------------------ 语法 */

    void fail(const std::string &message) {
        if (m_error.empty()) {
            m_error = message;
            m_errorPosition = m_token.start;
        }
    }

    bool isOperator(const char *text) const {
        return m_token.kind == TokenKind::Operator && m_token.text == text;
    }

    /** @brief 下一个 token 是否能开启一个原子（用于隐式乘法判定）。
     *
     * `,` `)` `]` `|` 都不会开启原子，因此 `diff(f, x)` 里的逗号不会
     * 被误当成隐式乘法——这是手写解析器最容易踩的坑。
     */
    bool startsAtom() const {
        switch (m_token.kind) {
        case TokenKind::Number:
        case TokenKind::LeftParen:
        case TokenKind::LeftBracket:
        case TokenKind::Bar:
        case TokenKind::Identifier: return true;
        default: return false;
        }
    }

    Ptr parseAssignment() {
        Ptr left = parseAdditive();
        if (!left) return Ptr();
        if (isOperator("=")) {
            advance();
            Ptr right = parseAdditive();
            if (!right) return Ptr();
            // 区间写法 `i=1..n`（求和/敛散性判别的上下限）→ range(i, 1, n)
            if (isOperator("..")) {
                advance();
                Ptr upper = parseAdditive();
                if (!upper) return Ptr();
                Memory::Vector<NodePtr> arguments{left.get(), right.get(), upper.get()};
                return Ptr(Node::makeFunction("range", arguments));
            }
            // 等式表示为 equal(left, right)，供 solve 直接使用
            Memory::Vector<NodePtr> arguments{left.get(), right.get()};
            return Ptr(Node::makeFunction("equal", arguments));
        }
        return left;
    }

    Ptr parseAdditive() {
        Ptr left = parseMultiplicative();
        if (!left) return Ptr();
        while (m_token.kind == TokenKind::Operator && (m_token.text == "+" || m_token.text == "-")) {
            const bool minus = m_token.text == "-";
            advance();
            Ptr right = parseMultiplicative();
            if (!right) return Ptr();
            left = minus ? makeAdd(left, makeNegate(right)) : makeAdd(left, right);
            if (!left) return Ptr();
        }
        return left;
    }

    Ptr parseMultiplicative() {
        Ptr left = parseUnary();
        if (!left) return Ptr();
        for (;;) {
            if (m_token.kind == TokenKind::Operator && (m_token.text == "*" || m_token.text == "/")) {
                const bool divide = m_token.text == "/";
                advance();
                Ptr right = parseUnary();
                if (!right) return Ptr();
                left = divide ? makeDivide(left, right) : makeMul(left, right);
                if (!left) return Ptr();
                continue;
            }
            if (startsAtom()) {
                Ptr right = parseUnary();
                if (!right) return Ptr();
                left = makeMul(left, right);
                if (!left) return Ptr();
                continue;
            }
            break;
        }
        return left;
    }

    Ptr parseUnary() {
        if (m_token.kind == TokenKind::Operator && (m_token.text == "-" || m_token.text == "+")) {
            const bool minus = m_token.text == "-";
            advance();
            Ptr operand = parseUnary();
            if (!operand) return Ptr();
            return minus ? makeNegate(operand) : operand;
        }
        return parsePower();
    }

    Ptr parsePower() {
        Ptr base = parsePostfix();
        if (!base) return Ptr();
        if (isOperator("^")) {
            advance();
            // 右结合，且允许 2^-1
            Ptr exponent = parseUnary();
            if (!exponent) return Ptr();
            return makePow(base, exponent);
        }
        return base;
    }

    /** @brief 把 `callee(...)` 的括号部分接到 callee 上，生成函数调用节点。 */
    Ptr parseCallOn(const Ptr &callee) {
        if (!callee || m_token.kind != TokenKind::LeftParen) return Ptr();
        const char *name = callee->name();
        if (callee->kind() != NodeKind::Symbol || !name) return Ptr();
        advance();
        Memory::Vector<Ptr> arguments;
        if (m_token.kind != TokenKind::RightParen) {
            for (;;) {
                Ptr argument = parseAssignment();
                if (!argument) return Ptr();
                arguments.push_back(argument);
                if (m_token.kind == TokenKind::Operator && m_token.text == ",") {
                    advance();
                    continue;
                }
                break;
            }
        }
        if (m_token.kind != TokenKind::RightParen) {
            fail("函数调用缺少右括号");
            return Ptr();
        }
        advance();
        Memory::Vector<NodePtr> raw;
        raw.reserve(arguments.size());
        for (const Ptr &argument : arguments) raw.push_back(argument.get());
        return Ptr(Node::makeFunction(name, raw));
    }

    Ptr parsePostfix() {
        Ptr value = parseAtom();
        if (!value) return Ptr();
        for (;;) {
            if (isOperator("!")) {
                advance();
                // `!!` 是双阶乘（一次消费两个 `!`），不能当成两次阶乘
                const bool doubleFactorial = isOperator("!");
                if (doubleFactorial) advance();
                // 阶乘只对非负整数求值
                long long number = 0;
                if (doubleFactorial) {
                    Memory::Vector<NodePtr> arguments{value.get()};
                    value = Ptr(Node::makeFunction("factorial2", arguments));
                    continue;
                }
                if (value->isSmallInteger(&number) && number >= 0 && number <= 20) {
                    long long product = 1;
                    for (long long k = 2; k <= number; ++k) product *= k;
                    value = integer(product);
                } else {
                    Memory::Vector<NodePtr> arguments{value.get()};
                    value = Ptr(Node::makeFunction("factorial", arguments));
                }
                continue;
            }
            if (isOperator("'")) {
                // 连续撇号 `y''` 表示二阶导，整段一次消费
                int order = 0;
                while (isOperator("'")) {
                    advance();
                    ++order;
                }
                // 简单符号的撇号按 ODE 记法解释：y' ⇒ Derivative(y, x)
                // 其它形态（如 (x^2+1)'）沿用 diff(expr) 的写法。
                //
                // `y'(0)` 这种"导数在某点取值"：撇号后面紧跟括号时，括号归属于
                // 被求导的函数本身，于是生成 Derivative(y(0), x, n)——初值条件
                // 的解析就靠这个形状识别。
                if (value->kind() == NodeKind::Symbol && order <= 3) {
                    Ptr target = value;
                    if (m_allowUnknownCalls && m_token.kind == TokenKind::LeftParen) {
                        Ptr call = parseCallOn(value);
                        if (!call) return Ptr();
                        target = call;
                    }
                    const std::string name = value->name() ? value->name() : std::string();
                    const std::uint32_t variableId = internSymbol("x");
                    Ptr variablePtr = Ptr(Node::makeSymbolById(variableId));
                    Ptr orderPtr = integer(order);
                    Memory::Vector<NodePtr> arguments{target.get(), variablePtr.get(), orderPtr.get()};
                    value = Ptr(Node::makeFunction("Derivative", arguments));
                    continue;
                }
                for (int k = 0; k < order; ++k) {
                    Memory::Vector<NodePtr> arguments{value.get()};
                    value = Ptr(Node::makeFunction("diff", arguments));
                }
                continue;
            }
            break;
        }
        return value;
    }

    Ptr parseAtom() {
        switch (m_token.kind) {
        case TokenKind::Number: {
            Ptr value;
            if (m_token.isInteger) {
                value = Ptr(Node::makeInteger(m_token.integerValue));
            } else {
                value = makeNumber(m_token.floatNumerator, m_token.floatDenominator);
            }
            advance();
            return value;
        }
        case TokenKind::LeftParen: {
            advance();
            Ptr inner = parseAssignment();
            if (!inner) return Ptr();
            if (m_token.kind != TokenKind::RightParen) {
                fail("缺少右括号");
                return Ptr();
            }
            advance();
            return inner;
        }
        case TokenKind::LeftBracket: {
            // 方括号列表：[a, b, c] → list(a, b, c)
            advance();
            Memory::Vector<Ptr> elements;
            if (m_token.kind != TokenKind::RightBracket) {
                for (;;) {
                    Ptr element = parseAssignment();
                    if (!element) return Ptr();
                    elements.push_back(element);
                    if (m_token.kind == TokenKind::Operator && m_token.text == ",") {
                        advance();
                        continue;
                    }
                    break;
                }
            }
            if (m_token.kind != TokenKind::RightBracket) {
                fail("列表缺少右方括号 ]");
                return Ptr();
            }
            advance();
            Memory::Vector<NodePtr> raw;
            raw.reserve(elements.size());
            for (const Ptr &element : elements) raw.push_back(element.get());
            return Ptr(Node::makeFunction("list", raw));
        }
        case TokenKind::Bar: {
            advance();
            Ptr inner = parseAssignment();
            if (!inner) return Ptr();
            if (m_token.kind != TokenKind::Bar) {
                fail("缺少配对的竖线 |");
                return Ptr();
            }
            advance();
            Memory::Vector<NodePtr> arguments{inner.get()};
            return Ptr(Node::makeFunction("abs", arguments));
        }
        case TokenKind::Identifier: return parseIdentifier();
        default:
            fail("此处需要一个表达式");
            return Ptr();
        }
    }

    Ptr parseIdentifier() {
        const std::string name = m_token.text;
        advance();

        // 常数
        if (name == "pi" || name == "PI" || name == "Pi") return Ptr(Node::makeConstant(ConstantId::Pi));
        if (name == "e" || name == "E") return Ptr(Node::makeConstant(ConstantId::E));
        if (name == "oo" || name == "inf" || name == "infinity") return Ptr(Node::makeConstant(ConstantId::Infinity));
        // `gamma` 现在是 Γ 函数，欧拉常数改由 EulerGamma / eulergamma 表示
        if (name == "EulerGamma" || name == "eulergamma") return Ptr(Node::makeConstant(ConstantId::EulerGamma));
        // 复数单位暂不作为一等公民：只有大写 `I` 保留为符号 `I`，
        // 小写 `i` 是普通变量（求和下标等场景大量使用，不能被偷换）。
        if (name == "I") return Ptr(Node::makeSymbol("I"));

        // `define(f(x), ...)`：第一个参数是"函数签名"，需要按调用解析。
        // 这个开关只在这个括号内打开，不影响后续参数与其它位置。
        // 必须放在通用函数调用分支之前，否则会被 builtin 分支按普通参数解析。
        if (m_options.allowSignature && m_token.kind == TokenKind::LeftParen &&
            (name == "define" || name == "let")) {
            advance();
            ++m_signatureDepth;
            Ptr signature = parseAssignment();
            --m_signatureDepth;
            if (!signature) return Ptr();
            Memory::Vector<Ptr> arguments{signature};
            while (m_token.kind == TokenKind::Operator && m_token.text == ",") {
                advance();
                Ptr argument = parseAssignment();
                if (!argument) return Ptr();
                arguments.push_back(argument);
            }
            if (m_token.kind != TokenKind::RightParen) {
                fail("define 缺少右括号");
                return Ptr();
            }
            advance();
            Memory::Vector<NodePtr> raw;
            raw.reserve(arguments.size());
            for (const Ptr &argument : arguments) raw.push_back(argument.get());
            return Ptr(Node::makeFunction("define", raw));
        }

        // 函数调用：内置名字，或用户已登记的定义名（`f(2)` 不再读作 f 乘 2）
        if (m_token.kind == TokenKind::LeftParen && (isFunctionName(name) || isUserFunction(name))) {
            advance();
            Memory::Vector<Ptr> arguments;
            if (m_token.kind != TokenKind::RightParen) {
                for (;;) {
                    Ptr argument = parseAssignment();
                    if (!argument) return Ptr();
                    arguments.push_back(argument);
                    if (m_token.kind == TokenKind::Operator && m_token.text == ",") {
                        advance();
                        continue;
                    }
                    break;
                }
            }
            if (m_token.kind != TokenKind::RightParen) {
                fail("函数调用缺少右括号");
                return Ptr();
            }
            advance();

            // 语法糖：sqrt → ^(1/2)，cbrt → ^(1/3)，ln → log
            if (name == "sqrt" && arguments.size() == 1) {
                return makePow(arguments.front(), makeNumber(BigInt(1), BigInt(2)));
            }
            if (name == "cbrt" && arguments.size() == 1) {
                return makePow(arguments.front(), makeNumber(BigInt(1), BigInt(3)));
            }
            const char *canonical = name == "ln" ? "log" : name.c_str();
            Memory::Vector<NodePtr> raw;
            raw.reserve(arguments.size());
            for (const Ptr &argument : arguments) raw.push_back(argument.get());
            return Ptr(Node::makeFunction(canonical, raw));
        }

        // 未知函数调用：仅初值条件、define 签名这类场景需要（`y(0)` 是"函数在 0 点取值"，
        // `f(x)` 在 define 的第一个参数里是"函数签名"，而不是 y 乘 0 / f 乘 x）。
        // 默认关闭，避免 `x(1)` 这类写法被误读成函数。
        if ((m_allowUnknownCalls || m_options.allowAnyCall ||
             (m_options.allowSignature && m_signatureDepth > 0)) &&
            m_token.kind == TokenKind::LeftParen) {
            advance();
            Memory::Vector<Ptr> arguments;
            if (m_token.kind != TokenKind::RightParen) {
                for (;;) {
                    Ptr argument = parseAssignment();
                    if (!argument) return Ptr();
                    arguments.push_back(argument);
                    if (m_token.kind == TokenKind::Operator && m_token.text == ",") {
                        advance();
                        continue;
                    }
                    break;
                }
            }
            if (m_token.kind != TokenKind::RightParen) {
                fail("函数调用缺少右括号");
                return Ptr();
            }
            advance();
            Memory::Vector<NodePtr> raw;
            raw.reserve(arguments.size());
            for (const Ptr &argument : arguments) raw.push_back(argument.get());
            return Ptr(Node::makeFunction(name.c_str(), raw));
        }

        // 未声明的名字按符号处理（`sqrt` 单独出现时也是符号）
        return Ptr(Node::makeSymbol(name));
    }

    /** @brief 已知函数名表：决定 `name(...)` 是调用还是乘法。 */
    static bool isFunctionName(const std::string &name) {
        static const char *kNames[] = {
            "sin", "cos", "tan", "cot", "sec", "csc", "asin", "acos", "atan", "acot",
            "sinh", "cosh", "tanh", "asinh", "acosh", "atanh",
            "exp", "log", "ln", "sqrt", "cbrt", "abs", "sign", "floor", "ceil", "round",
            "factorial", "gamma", "max", "min", "gcd", "lcm",
            // 命令名：让 `expand((x+1)^3)` 走函数调用而不是当成 expand*(x+1)^3
            "simplify", "simp", "expand", "factor", "together", "trigsimp", "trigs",
            "subs", "substitute", "replace", "eval",
            "define", "let", "undef", "undefine", "forget", "defs", "definitions",
            "latex", "diff", "differentiate", "derivative", "integrate", "int", "integral",
            "limit", "lim", "taylor", "series", "solve", "groebner", "gb", "resultant", "res",
            "critical", "extremum", "criticalpoints",
            "list", "equal", "range",
            // 级数与特殊函数
            "sum", "summation", "converge", "converges", "bernoulli", "guess",
            "specialsimplify", "specialsimp", "spsimp",
            "gamma", "factorial", "factorial2", "doublefactorial", "beta", "binomial",
            "harmonic", "erf", "erfc", "erfi", "Ei", "li",
            "besselj", "bessely", "besseli", "besselk",
            // 线性代数
            "det", "determinant", "inv", "inverse", "invert", "transpose", "trans",
            "trace", "tr", "matmul", "dot", "rank", "rref", "charpoly",
            "eigenvals", "eigenvalues", "identity", "eye",
            // 常微分方程
            "dsolve", "ode", "Derivative", "odeclassify",
        };
        for (const char *candidate : kNames) {
            if (name == candidate) return true;
        }
        return false;
    }

    /** @brief 名字是否是用户已登记的定义函数（由解析选项提供判据）。 */
    bool isUserFunction(const std::string &name) const {
        if (!m_options.isUserFunction) return false;
        return m_options.isUserFunction(name);
    }

    const std::string &m_text;
    bool m_allowUnknownCalls = false;
    /** @brief 用户定义相关的解析选项。 */
    ParseOptions m_options;
    /** @brief 正在解析 define 的签名参数（>0 时裸符号+括号按调用解析）。 */
    int m_signatureDepth = 0;
    std::size_t m_position = 0;
    Token m_token;
    std::string m_error;
    std::size_t m_errorPosition = 0;
};

} // namespace

bool isReservedConstantName(const std::string &name) {
    return name == "pi" || name == "PI" || name == "Pi" || name == "e" || name == "E" || name == "oo" ||
           name == "inf" || name == "infinity" || name == "EulerGamma" || name == "eulergamma";
}

namespace {

/**
 * @brief 判断输入是否以 `dsolve(` / `ode(` 开头。
 *
 * 这两个命令的参数里可能带初值条件 `y(0) = 1`，需要把裸符号后跟括号读成函数调用；
 * 其它场景保持"隐式乘法"的直觉读法。
 */
bool startsWithOdeCommand(const std::string &text) {
    std::size_t index = 0;
    while (index < text.size() && std::isspace(static_cast<unsigned char>(text[index]))) ++index;
    for (const char *name : {"dsolve", "ode"}) {
        const std::size_t length = std::strlen(name);
        if (text.compare(index, length, name) != 0) continue;
        std::size_t after = index + length;
        while (after < text.size() && std::isspace(static_cast<unsigned char>(text[after]))) ++after;
        if (after < text.size() && text[after] == '(') return true;
    }
    return false;
}

} // namespace

ParseResult parseExpression(const std::string &text) {
    Parser parser(text, startsWithOdeCommand(text));
    return parser.run();
}

ParseResult parseExpression(const std::string &text, const ParseOptions &options) {
    Parser parser(text, options, startsWithOdeCommand(text));
    return parser.run();
}

ParseResult parseInitialConditionText(const std::string &text) {
    Parser parser(text, /*allowUnknownCalls=*/true);
    return parser.run();
}

} // namespace Cas::Detail

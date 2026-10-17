#include "cas/engine.h"

#include "cas/detail/calculus.h"
#include "cas/detail/multivariate.h"
#include "cas/detail/parser.h"
#include "cas/detail/polynomial.h"
#include "cas/detail/printer.h"
#include "cas/detail/simplify.h"
#include "cas/detail/solver.h"
#include "cas/detail/trigonometry.h"

#include <chrono>
#include <cstdio>
#include <cstring>

namespace Cas {
namespace {

using namespace Cas::Detail;

/** @brief 计时器：作用域结束把毫秒数写回 Result。 */
class Timer {
public:
    explicit Timer(double *target) : m_target(target) {}
    ~Timer() {
        const auto end = std::chrono::steady_clock::now();
        *m_target = std::chrono::duration<double, std::milli>(end - m_start).count();
    }

private:
    double *m_target;
    std::chrono::steady_clock::time_point m_start = std::chrono::steady_clock::now();
};

/** @brief 去掉首尾空白。 */
std::string trim(const std::string &text) {
    std::size_t begin = 0, end = text.size();
    while (begin < end && std::isspace(static_cast<unsigned char>(text[begin]))) ++begin;
    while (end > begin && std::isspace(static_cast<unsigned char>(text[end - 1]))) --end;
    return text.substr(begin, end - begin);
}

/** @brief 取函数名（`f(...)` 形式），非函数返回空串。 */
std::string functionName(const NodePtr &node) {
    if (!node || node->kind() != NodeKind::Function || !node->name()) return std::string();
    return node->name();
}

/** @brief 收集带名参数的键值（`diff(f, x)` 里的 x、`limit(f, x, 0)` 里的 0）。 */
struct Arguments {
    Memory::Vector<Ptr> values;
    /** @brief 取第 index 个实参，越界返回空。 */
    Ptr at(std::size_t index) const { return index < values.size() ? values[index] : Ptr(); }
};

} // namespace

/* ------------------------------------------------------------------ 构造 */

Engine::Engine() { Allocator::initialize(); }

Engine::~Engine() = default;

/* ------------------------------------------------------------ 结果封装 */

namespace {

/** @brief 由表达式构造 Result。 */
Result makeResult(const Ptr &expression, bool unicode) {
    Result result;
    if (!expression) {
        result.message = "内部错误：表达式为空";
        return result;
    }
    PrintOptions options;
    options.useUnicode = unicode;
    result.text = toText(expression.get(), options);
    result.latex = toLatex(expression.get(), options);
    result.ok = true;
    return result;
}

/** @brief 构造失败 Result。 */
Result makeError(const std::string &message, std::size_t position = 0) {
    Result result;
    result.message = message;
    result.position = position;
    return result;
}

/** @brief 由符号编号构造节点。 */
Ptr idx(std::uint32_t id) { return Ptr(Node::makeSymbolById(id)); }

/** @brief 把表达式节点渲染成可再次解析的文本（用于把参数当字符串回填）。 */
std::string render(const Ptr &expression) { return expression ? toText(expression.get()) : std::string(); }

/** @brief 取参数里的变量名，缺省为 "x"。
 *
 * 只有当参数确实是符号（`x`、`y`）时才当作变量名；其他形态（数值、表达式）
 * 一律回落为默认变量，避免把 `limit(f, x, 0)` 的第三个参数 `0` 误当成变量名。
 */
std::string variableName(const Ptr &argument) {
    if (!argument) return "x";
    if (argument->kind() == NodeKind::Symbol && argument->name()) return argument->name();
    return "x";
}

} // namespace

/* ------------------------------------------------------------ 基本变换 */

Result Engine::simplify(const std::string &expression) {
    Result result;
    Timer timer(&result.elapsedMs);
    const ParseResult parsed = parseExpression(expression);
    if (!parsed.ok) return makeError(parsed.message, parsed.position);
    result = makeResult(Detail::simplify(parsed.expression), m_unicodeOutput);
    if (result.ok) m_lastResult = result.text;
    return result;
}

Result Engine::expand(const std::string &expression) {
    Result result;
    Timer timer(&result.elapsedMs);
    const ParseResult parsed = parseExpression(expression);
    if (!parsed.ok) return makeError(parsed.message, parsed.position);
    result = makeResult(Detail::expand(parsed.expression), m_unicodeOutput);
    if (result.ok) m_lastResult = result.text;
    return result;
}

Result Engine::factor(const std::string &expression) {
    Result result;
    Timer timer(&result.elapsedMs);
    const ParseResult parsed = parseExpression(expression);
    if (!parsed.ok) return makeError(parsed.message, parsed.position);
    result = makeResult(Detail::factor(parsed.expression), m_unicodeOutput);
    if (result.ok) m_lastResult = result.text;
    return result;
}

Result Engine::together(const std::string &expression) {
    Result result;
    Timer timer(&result.elapsedMs);
    const ParseResult parsed = parseExpression(expression);
    if (!parsed.ok) return makeError(parsed.message, parsed.position);
    result = makeResult(Detail::together(parsed.expression), m_unicodeOutput);
    if (result.ok) m_lastResult = result.text;
    return result;
}

Result Engine::trigSimplify(const std::string &expression) {
    Result result;
    Timer timer(&result.elapsedMs);
    const ParseResult parsed = parseExpression(expression);
    if (!parsed.ok) return makeError(parsed.message, parsed.position);
    result = makeResult(Detail::trigSimplify(parsed.expression), m_unicodeOutput);
    if (result.ok) m_lastResult = result.text;
    return result;
}

/* ------------------------------------------------------------ 微积分 */

namespace {

/** @brief 由变量名取符号编号。 */
std::uint32_t symbolOf(const std::string &name) { return internSymbol(trim(name)); }

} // namespace

Result Engine::differentiate(const std::string &expression, const std::string &variable) {
    Result result;
    Timer timer(&result.elapsedMs);
    const ParseResult parsed = parseExpression(expression);
    if (!parsed.ok) return makeError(parsed.message, parsed.position);
    const std::uint32_t id = symbolOf(variable);
    result = makeResult(Detail::differentiate(parsed.expression, id), m_unicodeOutput);
    if (result.ok) m_lastResult = result.text;
    return result;
}

Result Engine::integrate(const std::string &expression, const std::string &variable) {
    Result result;
    Timer timer(&result.elapsedMs);
    const ParseResult parsed = parseExpression(expression);
    if (!parsed.ok) return makeError(parsed.message, parsed.position);
    const std::uint32_t id = symbolOf(variable);
    Ptr antiderivative = Detail::integrate(parsed.expression, id);
    if (!antiderivative) return makeError("无法求出原函数（超出当前支持的积分规则集）");
    result = makeResult(antiderivative, m_unicodeOutput);
    if (result.ok) m_lastResult = result.text;
    return result;
}

Result Engine::integrateDefinite(const std::string &expression, const std::string &lower, const std::string &upper,
                                 const std::string &variable) {
    Result result;
    Timer timer(&result.elapsedMs);
    const ParseResult body = parseExpression(expression);
    if (!body.ok) return makeError(body.message, body.position);
    const ParseResult lowerValue = parseExpression(lower);
    if (!lowerValue.ok) return makeError("下限无法解析：" + lowerValue.message, lowerValue.position);
    const ParseResult upperValue = parseExpression(upper);
    if (!upperValue.ok) return makeError("上限无法解析：" + upperValue.message, upperValue.position);

    const std::uint32_t id = symbolOf(variable);
    Ptr value = Detail::definiteIntegral(body.expression, id, lowerValue.expression, upperValue.expression);
    if (!value) return makeError("定积分无法求值（无法求出原函数或上下限代入失败）");
    // 上下限代入常留下 cos(pi) 之类的特殊角，最后再过一遍三角化简
    Ptr cleaned = Detail::trigSimplify(value);
    if (cleaned) value = cleaned;
    result = makeResult(value, m_unicodeOutput);
    if (result.ok) m_lastResult = result.text;
    return result;
}

Result Engine::limit(const std::string &expression, const std::string &point, const std::string &variable) {
    Result result;
    Timer timer(&result.elapsedMs);
    const ParseResult body = parseExpression(expression);
    if (!body.ok) return makeError(body.message, body.position);
    const ParseResult pointValue = parseExpression(point);
    if (!pointValue.ok) return makeError("趋近点无法解析：" + pointValue.message, pointValue.position);

    const std::uint32_t id = symbolOf(variable);
    Ptr value = Detail::limit(body.expression, id, pointValue.expression);
    if (!value) return makeError("极限无法确定（超出当前支持的极限策略）");
    result = makeResult(value, m_unicodeOutput);
    if (result.ok) m_lastResult = result.text;
    return result;
}

Result Engine::taylor(const std::string &expression, const std::string &point, int order,
                      const std::string &variable) {
    Result result;
    Timer timer(&result.elapsedMs);
    if (order < 0) return makeError("展开阶数必须非负");
    const ParseResult body = parseExpression(expression);
    if (!body.ok) return makeError(body.message, body.position);
    const ParseResult pointValue = parseExpression(point);
    if (!pointValue.ok) return makeError("展开点无法解析：" + pointValue.message, pointValue.position);

    const std::uint32_t id = symbolOf(variable);
    Ptr value = Detail::taylorExpand(body.expression, id, pointValue.expression, order);
    if (!value) return makeError("Taylor 展开失败");
    result = makeResult(value, m_unicodeOutput);
    if (result.ok) m_lastResult = result.text;
    return result;
}

/* -------------------------------------------------------------- 求解 */

namespace {

/** @brief 把等式节点拆成左右两侧；非等式左=expr、右=0。 */
bool splitEquation(const Ptr &expression, Ptr *left, Ptr *right) {
    if (!expression) return false;
    if (functionName(expression.get()) == "equal" && expression->size() == 2) {
        *left = Ptr::borrow(expression->operand(0));
        *right = Ptr::borrow(expression->operand(1));
        return true;
    }
    *left = expression;
    *right = Ptr(Node::makeInteger(BigInt(0)));
    return true;
}

/** @brief 把解列表格式化为 `x = a, x = b` 或 `{a, b}`。 */
std::string formatSolutions(const Memory::Vector<Solution> &solutions, std::uint32_t variable) {
    if (solutions.empty()) return "无解（或在当前支持的范围内无法求出解析解）";
    const char *name = symbolName(variable);
    std::string text = name ? name : "x";
    text += " = ";
    for (std::size_t i = 0; i < solutions.size(); ++i) {
        if (i) text += ", " + std::string(name ? name : "x") + " = ";
        text += toText(solutions[i].value.get());
    }
    return text;
}

} // namespace

Result Engine::solve(const std::string &equation, const std::string &variable) {
    Result result;
    Timer timer(&result.elapsedMs);
    const ParseResult parsed = parseExpression(equation);
    if (!parsed.ok) return makeError(parsed.message, parsed.position);

    const std::uint32_t id = symbolOf(variable);
    Ptr left, right;
    if (!splitEquation(parsed.expression, &left, &right)) return makeError("方程形式无法识别");

    Memory::Vector<Solution> solutions;
    if (!Detail::solve(left.get(), right.get(), id, &solutions)) {
        // 返回 false 有两种情形：恒等式，或者无法识别
        if (Detail::simplify(makeAdd(left.share(), makeNegate(right)))->isZero()) {
            result.ok = true;
            result.text = "恒等式（任意值都满足）";
            result.latex = result.text;
            return result;
        }
        return makeError("无法求解该方程（超出多项式/简单超越方程的求解范围）");
    }

    result.ok = true;
    result.text = formatSolutions(solutions, id);
    result.latex = result.text;
    if (!solutions.empty()) {
        // 单个解时同步给出纯表达式，方便继续参与计算
        m_lastResult = toText(solutions.front().value.get());
    }
    return result;
}

Result Engine::solveSystem(const std::vector<std::string> &equations, const std::vector<std::string> &variables) {
    Result result;
    Timer timer(&result.elapsedMs);
    if (equations.empty() || variables.empty()) return makeError("方程组或变量表为空");

    Memory::Vector<std::uint32_t> ids;
    ids.reserve(variables.size());
    for (const std::string &name : variables) ids.push_back(symbolOf(name));

    Memory::Vector<Memory::Vector<Ptr>> rows;
    for (const std::string &line : equations) {
        const ParseResult parsed = parseExpression(line);
        if (!parsed.ok) return makeError("方程解析失败：" + parsed.message, parsed.position);
        Ptr left, right;
        if (!splitEquation(parsed.expression, &left, &right)) return makeError("方程形式无法识别");

        Ptr difference = Detail::simplify(makeAdd(left.share(), makeNegate(right)));
        Ptr expanded = Detail::expand(difference);

        // 逐个变量取系数；每次都先把该变量之外的部分剥掉，避免符号互相干扰
        Memory::Vector<Ptr> row;
        row.reserve(ids.size() + 1);
        Ptr remainder = expanded; // 始终维护"尚未取系数的部分"
        for (std::size_t i = 0; i < ids.size(); ++i) {
            UnivariatePoly poly;
            if (!toPolynomial(remainder.get(), ids[i], &poly)) {
                return makeError("当前只支持线性方程组（含非线性项）");
            }
            if (poly.degree() > 1) return makeError("当前只支持线性方程组（次数超过 1）");
            Ptr coefficient = poly.degree() == 1 ? poly.coefficient(1) : Ptr(Node::makeInteger(BigInt(0)));
            row.push_back(coefficient);
            // 剥掉这一项，剩下的只含后续变量与常数
            remainder = Detail::simplify(
                makeAdd(remainder.share(), makeNegate(makeMul(coefficient.share(), idx(ids[i])))));
        }
        // 剩余部分即常数项 c，方程形如 Σ aᵢxᵢ + c = 0 → 右端项取 -c
        row.push_back(makeNegate(remainder));
        rows.push_back(row);
    }

    Memory::Vector<Ptr> solution;
    if (!solveLinearSystem(rows, ids, &solution)) return makeError("方程组无唯一解（奇异或不相容）");

    std::string text;
    for (std::size_t i = 0; i < ids.size(); ++i) {
        if (i) text += ", ";
        const char *name = symbolName(ids[i]);
        text += std::string(name ? name : "x") + " = " + toText(solution[i].get());
    }
    result.ok = true;
    result.text = text;
    result.latex = text;
    return result;
}

Result Engine::solveInequality(const std::string &inequality, const std::string &variable) {
    Result result;
    Timer timer(&result.elapsedMs);

    // 手工切出关系符（解析器把 <、> 当非法字符，这里在文本层处理）
    static const char *kRelations[] = {"<=", ">=", "<", ">", "!="};
    std::string leftText, rightText, relation;
    for (const char *candidate : kRelations) {
        const std::size_t at = inequality.find(candidate);
        if (at == std::string::npos) continue;
        relation = candidate;
        leftText = inequality.substr(0, at);
        rightText = inequality.substr(at + std::strlen(candidate));
        break;
    }
    if (relation.empty()) return makeError("未找到关系符（支持 <、<=、>、>=、!=）");

    const ParseResult left = parseExpression(leftText);
    if (!left.ok) return makeError("左边无法解析：" + left.message, left.position);
    const ParseResult right = parseExpression(rightText);
    if (!right.ok) return makeError("右边无法解析：" + right.message, right.position);

    const std::uint32_t id = symbolOf(variable);
    std::string description;
    if (!Detail::solveInequality(left.expression.get(), right.expression.get(), id, relation.c_str(), &description)) {
        return makeError("当前只支持线性与二次不等式");
    }
    result.ok = true;
    result.text = description;
    result.latex = description;
    return result;
}

Result Engine::gcd(const std::string &left, const std::string &right, const std::string &variable) {
    Result result;
    Timer timer(&result.elapsedMs);
    const ParseResult a = parseExpression(left);
    if (!a.ok) return makeError(a.message, a.position);
    const ParseResult b = parseExpression(right);
    if (!b.ok) return makeError(b.message, b.position);
    result = makeResult(polyGcd(a.expression.get(), b.expression.get(), symbolOf(variable)), m_unicodeOutput);
    return result;
}

Result Engine::resultant(const std::string &left, const std::string &right, const std::string &variable) {
    Result result;
    Timer timer(&result.elapsedMs);
    const ParseResult a = parseExpression(left);
    if (!a.ok) return makeError(a.message, a.position);
    const ParseResult b = parseExpression(right);
    if (!b.ok) return makeError(b.message, b.position);
    result = makeResult(polyResultant(a.expression.get(), b.expression.get(), symbolOf(variable)), m_unicodeOutput);
    return result;
}

Result Engine::groebner(const std::vector<std::string> &polynomials, const std::vector<std::string> &variables) {
    Result result;
    Timer timer(&result.elapsedMs);
    if (polynomials.empty() || variables.empty()) return makeError("生成元或变量表为空");

    Memory::Vector<std::uint32_t> ids;
    for (const std::string &name : variables) ids.push_back(symbolOf(name));

    Memory::Vector<Ptr> generators;
    for (const std::string &line : polynomials) {
        const ParseResult parsed = parseExpression(line);
        if (!parsed.ok) return makeError("多项式解析失败：" + parsed.message, parsed.position);
        generators.push_back(Detail::expand(parsed.expression));
    }

    const Memory::Vector<Ptr> basis = groebnerBasis(generators, ids);
    std::string text = "{";
    for (std::size_t i = 0; i < basis.size(); ++i) {
        if (i) text += ", ";
        text += toText(basis[i].get());
    }
    text += "}";
    result.ok = true;
    result.text = text;
    result.latex = text;
    return result;
}

Result Engine::toLatexResult(const std::string &expression) {
    Result result;
    Timer timer(&result.elapsedMs);
    const ParseResult parsed = parseExpression(expression);
    if (!parsed.ok) return makeError(parsed.message, parsed.position);
    PrintOptions options;
    options.useUnicode = false;
    result.latex = toLatex(parsed.expression.get(), options);
    result.text = toText(parsed.expression.get());
    result.ok = true;
    return result;
}

/* --------------------------------------------- 节点直传（保结构）版本 */

Result Engine::derivativeOfNode(const Detail::Ptr &expression, const std::string &variable) {
    Result result;
    Timer timer(&result.elapsedMs);
    if (!expression) return makeError("diff 需要一个表达式参数");
    result = makeResult(Detail::differentiate(expression, symbolOf(variable)), m_unicodeOutput);
    if (result.ok) m_lastResult = result.text;
    return result;
}

Result Engine::limitOfNode(const Detail::Ptr &expression, const Detail::Ptr &point, const std::string &variable) {
    Result result;
    Timer timer(&result.elapsedMs);
    if (!expression || !point) return makeError("limit 需要表达式与趋近点");
    Ptr value = Detail::limit(expression, symbolOf(variable), point);
    if (!value) return makeError("极限无法确定（超出当前支持的极限策略）");
    result = makeResult(value, m_unicodeOutput);
    if (result.ok) m_lastResult = result.text;
    return result;
}

Result Engine::taylorOfNode(const Detail::Ptr &expression, const Detail::Ptr &point, int order,
                            const std::string &variable) {
    Result result;
    Timer timer(&result.elapsedMs);
    if (!expression || !point) return makeError("taylor 需要表达式与展开点");
    if (order < 0) return makeError("展开阶数必须非负");
    Ptr value = Detail::taylorExpand(expression, symbolOf(variable), point, order);
    if (!value) return makeError("Taylor 展开失败");
    result = makeResult(value, m_unicodeOutput);
    if (result.ok) m_lastResult = result.text;
    return result;
}

/* -------------------------------------------------------------- 命令分发 */

Result Engine::evaluate(const std::string &input) {
    const std::string text = trim(input);
    if (text.empty()) return makeError("输入为空");

    // 先按纯表达式解析；顶层是 `name(...)` 且 name 是命令时直接分发
    const ParseResult parsed = parseExpression(text);
    if (parsed.ok) {
        // 顶层命令无条件分发：命令名本身就是语法，不需要额外识别
        if (parsed.expression->kind() == NodeKind::Function && parsed.expression->size() >= 1) {
            Result dispatched = dispatch(functionName(parsed.expression.get()), parsed.expression.get(), text);
            // dispatch 只在"不认识这个函数名"时返回空 message 且 ok=false
            if (dispatched.ok || !dispatched.message.empty()) return dispatched;
        }
        Result simplified;
        Timer timer(&simplified.elapsedMs);
        // 先做三角/特殊角折叠，再做通用化简；三角化简失败时回退到纯代数化简
        Ptr value = Detail::trigSimplify(parsed.expression);
        if (!value) value = Detail::simplify(parsed.expression);
        simplified = makeResult(value, m_unicodeOutput);
        if (simplified.ok) m_lastResult = simplified.text;
        return simplified;
    }

    // 解析失败但看起来像不等式，交给不等式求解
    if (text.find('<') != std::string::npos || text.find('>') != std::string::npos) {
        return solveInequality(text);
    }
    return makeError(parsed.message, parsed.position);
}

std::string Engine::helpText() {
    return
        "命令（变量默认 x，可省略）：\n"
        "  simplify(expr)              化简\n"
        "  expand(expr)                展开\n"
        "  factor(expr)                因式分解\n"
        "  together(expr)              通分\n"
        "  trigsimp(expr)              三角化简\n"
        "  diff(expr [, var])          求导\n"
        "  integrate(expr [, var])     不定积分\n"
        "  integrate(expr, a, b [, var]) 定积分\n"
        "  limit(expr, point [, var])  极限（point 可为 oo）\n"
        "  taylor(expr, point, n)      Taylor 展开\n"
        "  solve(equation [, var])     解方程，equation 形如 x^2-1=0\n"
        "  solve([eq1, eq2], [x, y])   解线性方程组\n"
        "  groebner([p1, p2], [x, y])  Groebner 基\n"
        "  gcd(p, q [, var])           多项式 GCD\n"
        "  resultant(p, q [, var])     Resultant\n"
        "  latex(expr)                 输出 LaTeX\n"
        "\n"
        "直接输入表达式则做化简；支持 2x、3(x+1)、(x+1)(x-1)、x^2、sin(x) 等写法。\n"
        "REPL 中：help 显示帮助，quit 退出，unicode on/off 切换输出字符集。";
}

/* -------------------------------------------------------------- 命令分发 */

Result Engine::dispatch(const std::string &name, const Detail::Node *call, const std::string &source) {
    if (!call) return makeError("");

    const std::size_t count = call->size();
    auto argument = [call](std::size_t index) -> Ptr {
        return index < call->size() ? Ptr::borrow(call->operand(index)) : Ptr();
    };
    auto text = [&argument](std::size_t index) { return render(argument(index)); };

    if (name == "simplify" || name == "simp") return simplify(text(0));
    if (name == "expand") return expand(text(0));
    if (name == "factor") return factor(text(0));
    if (name == "together") return together(text(0));
    if (name == "trigsimp" || name == "trigSimplify" || name == "trigs") return trigSimplify(text(0));
    if (name == "latex") return toLatexResult(text(0));

    if (name == "diff" || name == "differentiate" || name == "derivative") {
        if (count == 0) return makeError("diff 需要至少一个参数");
        return derivativeOfNode(argument(0), count > 1 ? variableName(argument(1)) : "x");
    }
    if (name == "integrate" || name == "int" || name == "integral") {
        if (count == 0) return makeError("integrate 需要至少一个参数");
        // 两种定积分写法都要支持：
        //   integrate(f, a, b)      变量默认 x
        //   integrate(f, x, a, b)   变量在第二个参数
        if (count >= 4) {
            return integrateDefinite(text(0), text(2), text(3), variableName(argument(1)));
        }
        if (count == 3) {
            // 第二个参数若是纯符号，则它是变量、缺上限；否则按 a、b 处理
            if (argument(1) && argument(1)->kind() == NodeKind::Symbol) {
                return integrate(text(0), variableName(argument(1)));
            }
            return integrateDefinite(text(0), text(1), text(2), "x");
        }
        return integrate(text(0), count > 1 ? variableName(argument(1)) : "x");
    }
    if (name == "limit" || name == "lim") {
        // 支持三种写法：
        //   limit(f, x, 0)   表达式、变量、趋近点（数学课本写法）
        //   limit(f, 0, x)   表达式、趋近点、变量
        //   limit(f, x=0)    表达式、带变量的方程
        if (count == 2) {
            if (Ptr relation = argument(1)) {
                if (functionName(relation.get()) == "equal" && relation->size() == 2) {
                    return limitOfNode(argument(0), Ptr::borrow(relation->operand(1)),
                                       variableName(Ptr::borrow(relation->operand(0))));
                }
            }
            return makeError("limit 需要表达式与趋近点；可写 limit(f, x, 0) 或 limit(f, x=0)");
        }
        if (Ptr second = argument(1)) {
            // 第二个参数是符号 → 说明用了“变量在中间”的写法
            if (second->kind() == NodeKind::Symbol) {
                return limitOfNode(argument(0), argument(2), variableName(second));
            }
            return limitOfNode(argument(0), second, variableName(argument(2)));
        }
        return makeError("limit 参数不完整");
    }
    if (name == "taylor" || name == "series") {
        if (count < 3) return makeError("taylor 需要表达式、展开点与阶数");
        long long order = 0;
        if (!argument(2)->isSmallInteger(&order)) return makeError("Taylor 展开阶数必须是整数");
        return taylor(text(0), text(1), static_cast<int>(order), count > 3 ? variableName(argument(3)) : "x");
    }
    if (name == "solve") {
        if (count == 0) return makeError("solve 需要至少一个参数");
        // solve([e1, e2], [x, y]) —— 方程组
        if (argument(0)->kind() == NodeKind::Function && std::strcmp(argument(0)->name(), "list") == 0) {
            std::vector<std::string> equations;
            for (NodePtr element : argument(0)->operands()) equations.push_back(toText(element));
            std::vector<std::string> variables;
            if (Ptr list = argument(1)) {
                for (NodePtr element : list->operands()) variables.push_back(toText(element));
            }
            return solveSystem(equations, variables);
        }
        return solve(text(0), count > 1 ? variableName(argument(1)) : "x");
    }
    if (name == "groebner" || name == "gb") {
        if (count < 2) return makeError("groebner 需要多项式列表与变量列表");
        std::vector<std::string> polynomials, variables;
        for (NodePtr element : argument(0)->operands()) polynomials.push_back(toText(element));
        for (NodePtr element : argument(1)->operands()) variables.push_back(toText(element));
        return groebner(polynomials, variables);
    }
    if (name == "gcd") {
        if (count < 2) return makeError("gcd 需要两个多项式");
        return gcd(text(0), text(1), count > 2 ? variableName(argument(2)) : "x");
    }
    if (name == "resultant" || name == "res") {
        if (count < 2) return makeError("resultant 需要两个多项式");
        return resultant(text(0), text(1), count > 2 ? variableName(argument(2)) : "x");
    }

    // 不是命令：返回空 ok 空 message，让 evaluate 走默认化简
    (void)source;
    Result unknown;
    unknown.ok = false;
    unknown.message.clear();
    return unknown;
}

/* -------------------------------------------------------------- 内存 */

MemoryStats Engine::memoryStats() const { return Allocator::statistics(); }

void Engine::releaseMemory() { Allocator::initialize(); }

} // namespace Cas

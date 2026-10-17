#include "cas/engine.h"

#include "cas/detail/calculus.h"
#include "cas/detail/matrix.h"
#include "cas/detail/multivariate.h"
#include "cas/detail/ode.h"
#include "cas/detail/parser.h"
#include "cas/detail/polynomial.h"
#include "cas/detail/printer.h"
#include "cas/detail/series.h"
#include "cas/detail/simplify.h"
#include "cas/detail/special.h"
#include "cas/detail/solver.h"
#include "cas/detail/trigonometry.h"

#include <chrono>
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

/* -------------------------------------------------------------- 矩阵 */

/** @brief 是否为 list(...) 节点。 */
bool isListNode(const Ptr &node) {
    return node && node->kind() == NodeKind::Function && node->name() && std::strcmp(node->name(), "list") == 0;
}

/** @brief 矩阵结果在文本渲染之上附加一个结构化标记，供嵌套调用识别。 */
const char kMatrixMark = '\x1F';

bool parseMatrixText(const std::string &text, Detail::SymbolicMatrix *out);

/**
 * @brief 把 `[[1,2],[3,4]]` 或 `list(list(1,2),list(3,4))` 解析成矩阵。
 *
 * 也接受扁平列表 `[1,2,3,4]`（此时需要期望的行数才能定形状，这里一律当 1 行）。
 */
bool toMatrix(const Ptr &node, Detail::SymbolicMatrix *out) {
    if (!node) return false;
    if (!isListNode(node)) {
        // 上一级命令的渲染结果（带矩阵标记）也能被内层识别
        if (node->kind() == NodeKind::Symbol && node->name() && *node->name() == kMatrixMark) {
            return parseMatrixText(node->name(), out);
        }
        // 未求值的命令调用（如 det(inv(...))）不是矩阵：必须先把内层算出来
        if (node->kind() == NodeKind::Function && node->name() &&
            std::strcmp(node->name(), "list") != 0) {
            return false;
        }
        // 单个标量视为 1×1
        *out = Detail::SymbolicMatrix(1, 1);
        out->set(0, 0, node);
        return true;
    }
    if (node->size() == 0) return false;

    // 行形态：每个元素都是 list
    const bool nested = isListNode(Ptr::borrow(node->operand(0)));
    if (nested) {
        const std::size_t rows = node->size();
        std::size_t columns = 0;
        for (NodePtr child : node->operands()) {
            Ptr row(Ptr::borrow(child).share());
            if (!isListNode(row)) return false;
            if (columns == 0) columns = row->size();
            else if (row->size() != columns) return false; // 各行长度必须一致
        }
        if (columns == 0) return false;
        Memory::Vector<NodePtr> raw;
        raw.reserve(rows * columns);
        for (NodePtr child : node->operands()) {
            Ptr row(Ptr::borrow(child).share());
            for (NodePtr cell : row->operands()) raw.push_back(cell);
        }
        Memory::Vector<Ptr> entries;
        entries.reserve(raw.size());
        for (NodePtr cell : raw) entries.push_back(Ptr::borrow(cell).share());
        *out = Detail::SymbolicMatrix::fromRows(rows, columns, entries);
        return !out->isEmpty();
    }

    // 扁平列表 → 单行矩阵
    Memory::Vector<Ptr> entries;
    entries.reserve(node->size());
    for (NodePtr cell : node->operands()) entries.push_back(Ptr::borrow(cell).share());
    *out = Detail::SymbolicMatrix::fromRows(1, entries.size(), entries);
    return !out->isEmpty();
}

/** @brief 矩阵渲染成 `[[a, b], [c, d]]`。 */
std::string renderMatrix(const Detail::SymbolicMatrix &matrix, bool unicode) {
    std::string text = "[";
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        if (row > 0) text += ", ";
        text += "[";
        for (std::size_t column = 0; column < matrix.columns(); ++column) {
            if (column > 0) text += ", ";
            const Ptr &entry = matrix.at(row, column);
            text += entry ? toText(entry.get(), PrintOptions{unicode}) : std::string("0");
        }
        text += "]";
    }
    text += "]";
    return text;
}

/** @brief 矩阵的 LaTeX 形式，用 \begin{pmatrix}。 */
std::string renderMatrixLatex(const Detail::SymbolicMatrix &matrix) {
    std::string text = "\\begin{pmatrix}";
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        if (row > 0) text += " \\\\ ";
        for (std::size_t column = 0; column < matrix.columns(); ++column) {
            if (column > 0) text += " & ";
            const Ptr &entry = matrix.at(row, column);
            text += entry ? toLatex(entry.get(), PrintOptions{}) : std::string("0");
        }
    }
    text += "\\end{pmatrix}";
    return text;
}

/** @brief 用矩阵构造一个成功的 Result。 */
Result makeMatrixResult(const Detail::SymbolicMatrix &matrix, bool unicode) {
    Result result;
    result.ok = true;
    result.text = kMatrixMark + renderMatrix(matrix, unicode);
    result.latex = renderMatrixLatex(matrix);
    return result;
}

/** @brief 判断节点是否为数值零（用于逆矩阵的奇异检测）。 */
bool isNumericZero(const Ptr &node) {
    return node && node->isNumber() && node->numerator().isZero();
}

/**
 * @brief 去掉文本里所有的矩阵标记。
 *
 * 标记 `\x1F` 不是合法字符，解析器见到它会报"非法字符"。把矩阵结果喂给外层命令前
 * 需要先剥掉它——标记只用于 `evaluate` 层判断结果形态，不参与表达式语法。
 */
std::string stripMatrixMark(const std::string &text) {
    std::string clean;
    clean.reserve(text.size());
    for (char ch : text) {
        if (ch != kMatrixMark) clean.push_back(ch);
    }
    return clean;
}

/** @brief 从渲染文本还原矩阵；可以带矩阵标记，也可以是形如 `[[...]]` 的纯文本。 */
bool parseMatrixText(const std::string &text, Detail::SymbolicMatrix *out) {
    if (text.empty()) return false;
    const std::size_t begin = (text[0] == kMatrixMark) ? 1 : 0;
    if (text.size() <= begin) return false;
    const ParseResult parsed = parseExpression(text.substr(begin));
    if (!parsed.ok) return false;
    return toMatrix(parsed.expression, out);
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

Result Engine::criticalPoints(const std::string &expression, const std::string &variable) {
    Result result;
    Timer timer(&result.elapsedMs);
    const ParseResult body = parseExpression(expression);
    if (!body.ok) return makeError(body.message, body.position);

    const std::uint32_t id = symbolOf(variable);
    Ptr derivative = Detail::differentiate(body.expression, id);
    if (!derivative || !Detail::simplify(derivative)) return makeError("无法求导，无法定位驻点");

    Memory::Vector<Detail::CriticalPoint> points;
    if (!Detail::criticalPoints(body.expression, id, &points)) {
        return makeError("在 f'(x) = 0 上没有找到实驻点（或当前求解器无法解出）");
    }

    std::string text;
    char buffer[512];
    for (std::size_t index = 0; index < points.size(); ++index) {
        const Detail::CriticalPoint &point = points[index];
        if (!point.location) continue;
        PrintOptions options;
        options.useUnicode = m_unicodeOutput;
        const std::string location = toText(point.location.get(), options);
        if (point.kind == Detail::CriticalKind::Unknown) {
            std::snprintf(buffer, sizeof(buffer), "%s = %s（%s）", variable.c_str(), location.c_str(),
                          Detail::criticalKindName(point.kind));
        } else {
            const std::string value =
                point.value ? toText(point.value.get(), options) : std::string("?");
            std::snprintf(buffer, sizeof(buffer), "%s = %s（%s，f = %s）", variable.c_str(),
                          location.c_str(), Detail::criticalKindName(point.kind), value.c_str());
        }
        if (!text.empty()) text += "; ";
        text += buffer;
    }
    if (text.empty()) return makeError("未得到可用驻点");

    result.ok = true;
    result.text = text;
    // 与收敛/分类等「文字型结果」一致：LaTeX 用 \text{...} 包住，避免空字符串。
    result.latex = "\\text{" + text + "}";
    m_lastResult = result.text;
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

Result Engine::dsolve(const std::string &equation, const std::string &unknown,
                      const std::string &variable, const std::vector<std::string> &conditions) {
    Result result;
    Timer timer(&result.elapsedMs);

    const ParseResult parsed = parseExpression(equation);
    if (!parsed.ok) return makeError(parsed.message, parsed.position);

    Ptr left;
    Ptr right;
    if (!splitEquation(parsed.expression, &left, &right)) return makeError("方程形式无法识别");

    // 整理成 F = 0
    Ptr residual = Detail::simplify(makeAdd(left.share(), makeNegate(right)));
    if (!residual) return makeError("方程整理失败");

    const std::uint32_t unknownId = symbolOf(unknown);
    const std::uint32_t variableId = symbolOf(variable);

    // 初值条件也要单独解析
    Memory::Vector<Ptr> parsedConditions;
    for (const std::string &condition : conditions) {
        const ParseResult conditionParsed = parseInitialConditionText(condition);
        if (!conditionParsed.ok) return makeError("初值条件解析失败：" + conditionParsed.message,
                                                  conditionParsed.position);
        parsedConditions.push_back(conditionParsed.expression);
    }

    std::string method;
    std::string message;
    Ptr solution = Detail::solveOdeWithInitialValues(residual, unknownId, variableId, parsedConditions,
                                                     &method, &message);
    if (!solution) return makeError(message.empty() ? "无法求解该微分方程" : message);

    const char *unknownName = symbolName(unknownId);
    const std::string name = unknownName ? unknownName : unknown;

    result.ok = true;
    // 显式解（线性/二阶常系数）返回的是 `y = …` 的右端
    if (method.find("线性") != std::string::npos || method.find("二阶") != std::string::npos) {
        result.text = name + " = " + toText(solution.get());
        result.latex = name + " = " + toLatex(solution.get());
        m_lastResult = toText(solution.get());
        return result;
    }

    // 其余方法给出隐式解：写成 `… = 0`，其中已包含积分常数 C1
    result.text = "隐式解（" + method + "）：" + toText(solution.get()) + " = 0";
    result.latex = result.text;
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

/* ------------------------------------------------------------ 级数 */

namespace {

/**
 * @brief 把 `i=1..n` 解析出的 range(i, 1, n) 拆成求和上下限。
 * @return 成功时填充 limits；失败返回 false 并写入 message。
 */
bool parseRange(const Detail::Ptr &range, Detail::SumLimits *limits, std::string *message) {
    if (!range || range->kind() != NodeKind::Function || !range->name()) {
        *message = "需要区间写法 i=1..n（上限可写 oo 表示无穷），例如 sum(1/i^2, i=1..oo)";
        return false;
    }
    const std::string name = range->name();
    if (name == "range" && range->size() == 3) {
        const Ptr variable = Ptr::borrow(range->operand(0));
        if (variable->kind() != NodeKind::Symbol) {
            *message = "区间写法应为 变量=下界..上界";
            return false;
        }
        limits->variable = variable->symbolId();
        limits->lower = Detail::simplify(Ptr::borrow(range->operand(1)));
        limits->upper = Detail::simplify(Ptr::borrow(range->operand(2)));
        return true;
    }
    if (name == "equal" && range->size() == 2) {
        *message = "上限缺失：请写成 i=1..n";
        return false;
    }
    *message = "需要区间写法 i=1..n（上限可写 oo 表示无穷）";
    return false;
}

} // namespace

Result Engine::sumCommand(const Detail::Ptr &term, const Detail::Ptr &range, const Detail::Ptr &variable) {
    Result result;
    Timer timer(&result.elapsedMs);
    if (!term) return makeError("sum 需要被求和的通项");

    Detail::SumLimits limits;
    std::string message;
    if (!parseRange(range, &limits, &message)) return makeError(message);
    if (variable && variable->kind() == NodeKind::Symbol) limits.variable = variable->symbolId();

    Ptr value = Detail::summation(term, limits);
    if (!value) {
        return makeError("无法求出闭式；该和式超出当前支持的规则集（多项式、几何级数、p 级数与有限展开）");
    }
    result = makeResult(value, m_unicodeOutput);
    if (result.ok) m_lastResult = result.text;
    return result;
}

Result Engine::convergeCommand(const Detail::Ptr &term, const Detail::Ptr &range) {
    Result result;
    Timer timer(&result.elapsedMs);
    if (!term) return makeError("converge 需要级数的通项");

    Detail::SumLimits limits;
    std::string message;
    if (!parseRange(range, &limits, &message)) return makeError(message);
    if (!limits.upper || limits.upper->kind() != NodeKind::Constant ||
        limits.upper->constant() != ConstantId::Infinity) {
        return makeError("converge 只对无穷级数有意义，上限请写 oo");
    }

    Detail::Convergence verdict = Detail::Convergence::Unknown;
    Ptr value = Ptr();
    const bool decided = Detail::seriesConvergence(term, limits, &verdict, &value);

    std::string text;
    if (!decided || verdict == Detail::Convergence::Unknown) {
        text = "无法判定敛散性";
    } else if (verdict == Detail::Convergence::Divergent) {
        text = "发散";
    } else {
        text = "收敛";
        if (value) text += "，和为 " + render(value);
    }

    result.text = text;
    result.latex = "\\text{" + text + "}";
    result.ok = true;
    m_lastResult = text;
    return result;
}

/* -------------------------------------------------------------- 命令分发 */

/** @brief 可作为子表达式出现的矩阵命令（结果需要以矩阵文本形式回填）。 */
bool isMatrixCommand(const std::string &name) {
    return name == "det" || name == "determinant" || name == "inv" || name == "inverse" ||
           name == "invert" || name == "transpose" || name == "trans" || name == "trace" ||
           name == "tr" || name == "matmul" || name == "dot" || name == "identity" || name == "eye" ||
           name == "rref" || name == "charpoly" || name == "eigenvals" || name == "eigenvalues" ||
           name == "rank";
}

/**
 * @brief 可作为子表达式出现的标量命令（返回普通表达式文本）。
 *
 * `sum(i, i=1..n)`、`bernoulli(6)` 这类命令既能单独使用，也会嵌在
 * `simplify(...)`、`sum(...) + sum(...)` 里；必须允许内层求值，
 * 否则整条表达式只会原样回显。
 */
bool isScalarCommand(const std::string &name) {
    return name == "sum" || name == "summation" || name == "converge" || name == "converges" ||
           name == "bernoulli" || name == "guess" || name == "specialsimplify" ||
           name == "specialsimp" || name == "spsimp";
}

/** @brief 可作为子表达式出现、需要内层求值的命令（矩阵 + 标量）。 */
bool isResolvableCommand(const std::string &name) {
    return isMatrixCommand(name) || isScalarCommand(name);
}

/** @brief `range` 只是区间写法 `i=1..n` 的内部表示，不是真正的求值命令。 */
bool isRangeMarker(const std::string &name) { return name == "range" || name == "equal"; }

/**
 * @brief 表达式树里是否还留着未求值的命令调用。
 *
 * `resolveCommands` 会在某个内层命令求值失败时停下来；此时改写后的文本里仍会
 * 留下 `命令名(...)` 形态的片段，重新解析再分发只会得到同样的失败。这里用于在
 * 重新求值前判断"是否真的已经做到底了"，避免无意义的重解析。
 */
bool hasPendingCommand(const Ptr &node) {
    if (!node) return false;
    if (node->kind() == NodeKind::Function && node->name()) {
        const std::string name = node->name();
        // `range/equal` 是区间写法的内部节点，不是"未求值的命令"
        if (!isRangeMarker(name) && isResolvableCommand(name)) return true;
    }
    for (NodePtr child : node->operands()) {
        if (hasPendingCommand(Ptr::borrow(child))) return true;
    }
    return false;
}

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
        // 内层也可能有矩阵命令：det(inv(A))、matmul(A,B) 等
        std::string rewritten;
        const bool hit = resolveCommands(parsed.expression, &rewritten);
        if (!hit && hasPendingCommand(parsed.expression)) {
            // 内层命令确实存在，但求值失败（典型：inv 遇到奇异矩阵）。
            // 直接化简只会原样回显命令，用户看不出哪里错了，这里给出明确诊断。
            return makeError("内层矩阵命令无法求值；请检查参数（例如 inv 的矩阵是否可逆）");
        }
        if (hit && !rewritten.empty() && rewritten != text) {
            Result result;
            Timer timer(&result.elapsedMs);

            // 整个表达式被内层命令替换成了矩阵：保持矩阵渲染，不再走化简
            if (rewritten[0] == kMatrixMark) {
                Detail::SymbolicMatrix matrix;
                if (parseMatrixText(rewritten, &matrix)) {
                    result = makeMatrixResult(matrix, m_unicodeOutput);
                    m_lastResult = result.text;
                    return result;
                }
            }

            // 标量/混合结果：重新解析求值。若改写文本里还有命令残留，说明是求值失败，直接报错。
            const std::string plain = stripMatrixMark(rewritten);
            const ParseResult again = parseExpression(plain);
            if (again.ok && hasPendingCommand(again.expression)) {
                return makeError("内层矩阵命令无法求值：" + plain);
            }
            result = evaluate(plain);
            if (result.ok) m_lastResult = result.text;
            return result;
        }
        Result simplified;
        Timer timer(&simplified.elapsedMs);
        // 先做三角/特殊角折叠，再做通用化简；三角化简失败时回退到纯代数化简
        Ptr value = Detail::trigSimplify(parsed.expression);
        if (!value) value = Detail::simplify(parsed.expression);
        // 最后折叠特殊函数（Γ(n)、BesselJ(0,0) 之类）
        value = Detail::simplifySpecial(value);
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
        "  critical(expr [, var])      驻点与极值判定\n"
        "  solve(equation [, var])     解方程，equation 形如 x^2-1=0\n"
        "  solve([eq1, eq2], [x, y])   解线性方程组\n"
        "  groebner([p1, p2], [x, y])  Groebner 基\n"
        "  gcd(p, q [, var])           多项式 GCD\n"
        "  resultant(p, q [, var])     Resultant\n"
        "  latex(expr)                 输出 LaTeX\n"
        "\n"
        "线性代数（矩阵写作 [[1,2],[3,4]]，结果可与绝对值命令嵌套）：\n"
        "  det(A)                      行列式\n"
        "  inv(A)                      逆矩阵（奇异时报错）\n"
        "  transpose(A)                转置\n"
        "  trace(A)                    迹（仅方阵）\n"
        "  matmul(A, B)                矩阵乘法\n"
        "  rank(A) / rref(A)           秩 / 行最简形\n"
        "  charpoly(A [, var])         特征多项式\n"
        "  eigenvals(A)                特征值\n"
        "  identity(n)                 n 阶单位阵\n"
        "\n"
        "特殊函数（自动折叠已知值，未知值保持符号形式）：\n"
        "  gamma(z)                    Γ 函数（整数、半整数、极点）\n"
        "  factorial(n) / n!           阶乘\n"
        "  factorial2(n) / n!!         双阶乘\n"
        "  beta(a, b)                  B 函数\n"
        "  binomial(n, k)              组合数\n"
        "  harmonic(n [, r])           谐波数 H_n / 广义 H_n^(r)\n"
        "  erf(z) / erfc(z) / erfi(z)  误差函数\n"
        "  Ei(z) / li(z)               指数积分、对数积分\n"
        "  besselj(n, z) 等            Bessel 函数 J/Y/I/K\n"
        "  bernoulli(n)                伯努利数 B_n\n"
        "  specialsimplify(expr)       折叠表达式中的特殊函数\n"
        "\n"
        "级数求和（区间写作 i=1..n，上限 oo 表示无穷）：\n"
        "  sum(term, i=1..n)           幂和、p 级数、几何、算术-几何、裂项\n"
        "  converge(term, i=1..oo)     判别敛散性并给出收敛时的和\n"
        "  guess([a1, a2, ...])        由前几项猜通项公式\n"
        "\n"
        "常微分方程（未知函数默认 y，自变量默认 x；撇号 y' / y'' 表示导数）：\n"
        "  dsolve(eq)                  求解 ODE，eq 形如 y' = y、y'' + y = 0\n"
        "  dsolve(eq, y)               指定未知函数\n"
        "  dsolve(eq, y, t)            指定自变量\n"
        "  dsolve(eq, y(0)=1)          带初值条件，定出积分常数\n"
        "  一阶：可分离变量、一阶线性、齐次型（u = y/x）、伯努利方程\n"
        "  二阶：常系数线性（齐次三情形 + 待定系数法非齐次，含共振时乘 x）\n"
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
    if (name == "specialsimplify" || name == "specialsimp" || name == "spsimp") {
        const ParseResult parsed = parseExpression(text(0));
        if (!parsed.ok) return makeError(parsed.message, parsed.position);
        Ptr value = Detail::simplifySpecial(Detail::simplify(parsed.expression));
        return makeResult(value, m_unicodeOutput);
    }
    if (name == "sum" || name == "summation") return sumCommand(argument(0), argument(1), argument(2));
    if (name == "converge" || name == "converges") return convergeCommand(argument(0), argument(1));
    if (name == "bernoulli") {
        long long index = 0;
        if (!argument(0) || !argument(0)->isSmallInteger(&index) || index < 0) {
            return makeError("bernoulli 需要一个非负整数下标");
        }
        Ptr value = Detail::bernoulliNumber(index);
        if (!value) return makeError("伯努利数下标过大");
        return makeResult(value, m_unicodeOutput);
    }
    if (name == "guess") {
        const Ptr list = argument(0);
        if (!isListNode(list)) return makeError("guess 需要一个序列，如 guess([1,4,9,16,25])");
        Memory::Vector<Ptr> terms;
        for (NodePtr cell : list->operands()) {
            const ParseResult parsed = parseExpression(render(Ptr::borrow(cell)));
            terms.push_back(parsed.ok ? Detail::simplify(parsed.expression) : Ptr::borrow(cell));
        }
        const std::string indexName = count > 1 ? variableName(argument(1)) : "n";
        Ptr guess = Detail::guessSequence(terms, internSymbol(indexName));
        if (!guess) return makeError("无法从给定项识别通项");
        return makeResult(guess, m_unicodeOutput);
    }
    if (name == "expand") return expand(text(0));
    if (name == "factor") return factor(text(0));
    if (name == "together") return together(text(0));
    if (name == "trigsimp" || name == "trigSimplify" || name == "trigs") return trigSimplify(text(0));
    if (name == "latex") return toLatexResult(text(0));

    // 参数里套着未求值的命令（det(inv(A))）：交给 evaluate 的内层求值流程处理
    if (isMatrixCommand(name)) {
        for (std::size_t index = 0; index < count; ++index) {
            const Ptr child = argument(index);
            if (child && child->kind() == NodeKind::Function && child->name() &&
                isMatrixCommand(child->name())) {
                Result deferred;
                deferred.ok = false;
                deferred.message.clear();
                return deferred;
            }
        }
    }

    if (name == "det" || name == "determinant") {
        if (count == 0) return makeError("det 需要一个矩阵参数，如 det([[1,2],[3,4]])");
        Detail::SymbolicMatrix matrix;
        if (!toMatrix(argument(0), &matrix)) return makeError("无法解析矩阵；请写成 [[1,2],[3,4]] 形式");
        if (!matrix.isSquare()) return makeError("行列式只对非空方阵有定义");
        Result result;
        Timer timer(&result.elapsedMs);
        result = makeResult(Detail::determinant(matrix), m_unicodeOutput);
        if (result.ok) m_lastResult = result.text;
        return result;
    }
    if (name == "inv" || name == "inverse" || name == "invert") {
        if (count == 0) return makeError("inv 需要一个方阵参数，如 inv([[1,2],[3,4]])");
        Detail::SymbolicMatrix matrix;
        if (!toMatrix(argument(0), &matrix)) return makeError("无法解析矩阵；请写成 [[1,2],[3,4]] 形式");
        if (!matrix.isSquare()) return makeError("逆矩阵只对方阵有定义");
        Detail::SymbolicMatrix inverted;
        const Ptr determinantValue = Detail::inverse(matrix, &inverted);
        if (!determinantValue) return makeError("无法解析矩阵；请写成 [[1,2],[3,4]] 形式");
        if (isNumericZero(determinantValue)) return makeError("矩阵不可逆（行列式为 0）");
        Result result = makeMatrixResult(inverted, m_unicodeOutput);
        if (result.ok) m_lastResult = result.text;
        return result;
    }
    if (name == "transpose" || name == "trans") {
        if (count == 0) return makeError("transpose 需要一个矩阵参数");
        Detail::SymbolicMatrix matrix;
        if (!toMatrix(argument(0), &matrix)) return makeError("无法解析矩阵；请写成 [[1,2],[3,4]] 形式");
        Result result = makeMatrixResult(Detail::transpose(matrix), m_unicodeOutput);
        if (result.ok) m_lastResult = result.text;
        return result;
    }
    if (name == "trace" || name == "tr") {
        if (count == 0) return makeError("trace 需要一个方阵参数");
        Detail::SymbolicMatrix matrix;
        if (!toMatrix(argument(0), &matrix)) return makeError("无法解析矩阵；请写成 [[1,2],[3,4]] 形式");
        if (!matrix.isSquare()) return makeError("迹只对方阵有定义");
        Result result = makeResult(Detail::trace(matrix), m_unicodeOutput);
        if (result.ok) m_lastResult = result.text;
        return result;
    }
    if (name == "matmul" || name == "dot") {
        if (count < 2) return makeError("matmul 需要两个矩阵参数");
        Detail::SymbolicMatrix left, right;
        if (!toMatrix(argument(0), &left) || !toMatrix(argument(1), &right)) {
            return makeError("无法解析矩阵；请写成 [[1,2],[3,4]] 形式");
        }
        Detail::SymbolicMatrix product = Detail::multiplyMatrices(left, right);
        if (product.isEmpty()) return makeError("维度不匹配：左矩阵列数必须等于右矩阵行数");
        Result result = makeMatrixResult(product, m_unicodeOutput);
        if (result.ok) m_lastResult = result.text;
        return result;
    }
    if (name == "identity" || name == "eye") {
        if (count == 0) return makeError("identity 需要一个阶数，如 identity(3)");
        const Ptr &size = argument(0);
        if (!size || !size->isInteger() || size->numerator().isNegative()) {
            return makeError("identity 的阶数必须是非负整数");
        }
        if (!size->numerator().fitsLongLong()) return makeError("identity 的阶数过大");
        const auto order = static_cast<std::size_t>(size->numerator().toLongLong());
        if (order == 0) return makeError("identity 的阶数必须大于 0");
        Result result = makeMatrixResult(Detail::identityMatrix(order), m_unicodeOutput);
        if (result.ok) m_lastResult = result.text;
        return result;
    }
    if (name == "rank") {
        if (count == 0) return makeError("rank 需要一个矩阵参数");
        Detail::SymbolicMatrix matrix;
        if (!toMatrix(argument(0), &matrix)) return makeError("无法解析矩阵；请写成 [[1,2],[3,4]] 形式");
        std::size_t value = 0;
        if (!Detail::rank(matrix, &value)) return makeError("无法计算秩");
        Result result = makeResult(makeNumber(BigInt(static_cast<long long>(value)), BigInt(1)),
                                   m_unicodeOutput);
        if (result.ok) m_lastResult = result.text;
        return result;
    }
    if (name == "rref") {
        if (count == 0) return makeError("rref 需要一个矩阵参数");
        Detail::SymbolicMatrix matrix;
        if (!toMatrix(argument(0), &matrix)) return makeError("无法解析矩阵；请写成 [[1,2],[3,4]] 形式");
        Detail::SymbolicMatrix reduced;
        if (!Detail::reducedRowEchelon(matrix, &reduced, nullptr)) return makeError("无法化简");
        Result result = makeMatrixResult(reduced, m_unicodeOutput);
        if (result.ok) m_lastResult = result.text;
        return result;
    }
    if (name == "charpoly") {
        if (count == 0) return makeError("charpoly 需要一个方阵参数");
        Detail::SymbolicMatrix matrix;
        if (!toMatrix(argument(0), &matrix)) return makeError("无法解析矩阵；请写成 [[1,2],[3,4]] 形式");
        if (!matrix.isSquare()) return makeError("特征多项式只对方阵有定义");
        const std::uint32_t variable = internSymbol(count > 1 ? variableName(argument(1)) : "x");
        Result result;
        Timer timer(&result.elapsedMs);
        result = makeResult(Detail::characteristicPolynomial(matrix, variable), m_unicodeOutput);
        if (result.ok) m_lastResult = result.text;
        return result;
    }
    if (name == "eigenvals" || name == "eigenvalues") {
        if (count == 0) return makeError("eigenvals 需要一个方阵参数");
        Detail::SymbolicMatrix matrix;
        if (!toMatrix(argument(0), &matrix)) return makeError("无法解析矩阵；请写成 [[1,2],[3,4]] 形式");
        if (!matrix.isSquare()) return makeError("特征值只对方阵有定义");
        Memory::Vector<Ptr> values;
        if (!Detail::eigenvalues(matrix, &values)) return makeError("特征值无法求出（需要高次方程求根）");
        std::string joined;
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (i > 0) joined += ", ";
            joined += render(values[i]);
        }
        Result result;
        result.ok = true;
        result.text = "{" + joined + "}";
        result.latex = "\\{" + joined + "\\}";
        m_lastResult = result.text;
        return result;
    }

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
    if (name == "critical" || name == "extremum" || name == "criticalpoints") {
        if (count == 0) return makeError("critical 需要表达式，如 critical(x^2-1)");
        return criticalPoints(text(0), count > 1 ? variableName(argument(1)) : "x");
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
    if (name == "dsolve" || name == "ode") {
        if (count == 0) return makeError("dsolve 需要一个微分方程，如 dsolve(y' = x*y)");

        // 参数划分：第 1 个是方程；之后的参数依次是「未知函数」「自变量」和若干初值条件。
        // 初值条件形如 `y(0) = 1`（解析后是 Function("y", 0) 与右端的等式），
        // 未知函数/自变量则是裸符号，靠 elementKind 区分。
        std::string unknownName = "y";
        std::string variableArg = "x";
        std::vector<std::string> conditions;
        bool unknownTaken = false;
        bool variableTaken = false;
        for (std::size_t index = 1; index < count; ++index) {
            const Ptr child = argument(index);
            const std::string rendered = text(index);
            if (rendered.find('=') != std::string::npos) {
                conditions.push_back(rendered);
                continue;
            }
            // 裸符号：先当未知函数，再当自变量
            if (child && child->kind() == NodeKind::Symbol) {
                if (!unknownTaken) {
                    unknownName = variableName(child);
                    unknownTaken = true;
                    continue;
                }
                if (!variableTaken) {
                    variableArg = variableName(child);
                    variableTaken = true;
                    continue;
                }
            }
            conditions.push_back(rendered);
        }
        return dsolve(text(0), unknownName, variableArg, conditions);
    }
    if (name == "odeclassify") {
        if (count == 0) return makeError("odeclassify 需要一个微分方程");
        Result result;
        Timer timer(&result.elapsedMs);
        const ParseResult parsed = parseExpression(text(0));
        if (!parsed.ok) return makeError(parsed.message, parsed.position);
        Ptr left;
        Ptr right;
        if (!splitEquation(parsed.expression, &left, &right)) return makeError("方程形式无法识别");
        Ptr residual = Detail::simplify(makeAdd(left.share(), makeNegate(right)));
        const std::uint32_t unknownId = symbolOf(count > 1 ? variableName(argument(1)) : "y");
        const int order = Detail::derivativeOrder(residual, unknownId);
        std::string description;
        std::string message;
        Ptr solution = Detail::solveOde(residual, unknownId,
                                         symbolOf(count > 2 ? variableName(argument(2)) : "x"),
                                         &description, &message);
        result.ok = true;
        if (solution) {
            result.text = "阶数 " + std::to_string(order) + "；可识别类型：" + description;
        } else {
            result.text = "阶数 " + std::to_string(order) + "；无法识别：" + message;
        }
        result.latex = result.text;
        return result;
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

/* ---------------------------------------------------------- 内层命令 */

bool Engine::resolveCommands(const Ptr &node, std::string *out) {
    if (!node || !out) return false;

    // 非函数节点：递归处理子项，任一子项被改写了就按运算符拼回去。
    // 例如 `sum(i, i=1..10) + sum(i^2, i=1..10)` 的根是 Add，
    // 必须钻进去才能把两个求和命令都求值。
    if (node->kind() != NodeKind::Function) {
        const std::size_t count = node->size();
        if (count == 0) return false;
        Memory::Vector<std::string> rewritten;
        bool changed = false;
        for (std::size_t i = 0; i < count; ++i) {
            std::string piece;
            if (resolveCommands(Ptr::borrow(node->operand(i)), &piece)) {
                rewritten.push_back(piece);
                changed = true;
            } else {
                rewritten.push_back(render(Ptr::borrow(node->operand(i))));
            }
        }
        if (!changed) return false;
        switch (node->kind()) {
        case NodeKind::Add:
        case NodeKind::Mul: {
            const char *op = node->kind() == NodeKind::Add ? " + " : "*";
            std::string built;
            for (std::size_t i = 0; i < rewritten.size(); ++i) {
                if (i > 0) built += op;
                built += rewritten[i];
            }
            *out = built;
            return true;
        }
        case NodeKind::Pow: {
            if (rewritten.size() != 2) return false;
            *out = "(" + rewritten[0] + ")^(" + rewritten[1] + ")";
            return true;
        }
        default:
            return false;
        }
    }

    const std::string name = node->name() ? node->name() : std::string();

    // `sum`/`converge` 的区间实参 `i=1..n` 会被解析成内部的 `range(i,1,n)`，
    // 它不是普通表达式，render 之后既不能回写成 `i=1..n`、重新 dispatch 也认不出。
    // 所以这两个命令直接整棵调用渲染后分发，不做实参级重建。
    if (name == "sum" || name == "summation" || name == "converge" || name == "converges") {
        const std::string callText = render(node);
        const ParseResult parsedCall = parseExpression(stripMatrixMark(callText));
        if (!parsedCall.ok || parsedCall.expression->kind() != NodeKind::Function) return false;
        const Result outcome = dispatch(name, parsedCall.expression.get(), callText);
        if (!outcome.ok) return false;
        *out = outcome.text;
        return true;
    }

    // 先把子节点逐个求值成文本，顺带收集
    Memory::Vector<std::string> parts;
    bool childChanged = false;
    for (std::size_t i = 0; i < node->size(); ++i) {
        const Ptr child = Ptr::borrow(node->operand(i));
        std::string rendered;
        if (resolveCommands(child, &rendered)) {
            parts.push_back(rendered);
            childChanged = true;
        } else {
            parts.push_back(render(child));
        }
    }

    if (!isResolvableCommand(name)) {
        if (!childChanged) return false;
        // 普通函数：按"函数名(实参...)"重建，实参用改写后的文本
        std::string rebuilt = name + std::string("(");
        for (std::size_t i = 0; i < parts.size(); ++i) {
            if (i > 0) rebuilt += ", ";
            rebuilt += parts[i];
        }
        rebuilt += ")";
        *out = rebuilt;
        return true;
    }

    // 该节点本身就是命令：用"文本实参"重建一次调用，交给 dispatch
    std::string callText = name + std::string("(");
    for (std::size_t i = 0; i < parts.size(); ++i) {
        if (i > 0) callText += ", ";
        callText += parts[i];
    }
    callText += ")";

    // 子节点已经在上面递归处理过了，这里不再递归整体文本（否则会无限套娃）
    // 直接分发到该命令本身（不能再走 evaluate，否则会回到这里形成递归）
    const ParseResult parsedCall = parseExpression(stripMatrixMark(callText));
    if (!parsedCall.ok || parsedCall.expression->kind() != NodeKind::Function) return false;
    const Result outcome = dispatch(name, parsedCall.expression.get(), callText);
    if (!outcome.ok) return false;

    // 矩阵结果的渲染文本自带标记，原样回填：标记是"这段文本代表矩阵"的凭据，
    // 上层的 toMatrix 与最终的 evaluate 都靠它区分矩阵与普通列表。
    *out = outcome.text;
    return true;
}

/* -------------------------------------------------------------- 内存 */

MemoryStats Engine::memoryStats() const { return Allocator::statistics(); }

void Engine::releaseMemory() { Allocator::initialize(); }

} // namespace Cas

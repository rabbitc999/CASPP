#include "cas/detail/ode.h"

#include "cas/detail/matrix.h"
#include "cas/detail/calculus.h"
#include "cas/detail/polynomial.h"
#include "cas/detail/printer.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace Cas::Detail {

namespace {

/* ------------------------------------------------------------ 基础工具 */

Ptr integer(long long value) { return makeNumber(BigInt(value), BigInt(1)); }

/** @brief 构造符号节点。 */
Ptr symbolNode(std::uint32_t symbol) { return Ptr(Node::makeSymbolById(symbol)); }

/** @brief 该节点是否就是指定符号。 */
bool isSymbol(const Ptr &node, std::uint32_t symbol) {
    return node && node->kind() == NodeKind::Symbol && node->symbolId() == symbol;
}

/** @brief 表达式是否含某个符号（含作为导数实参的情形）。 */
bool mentions(const Ptr &expression, std::uint32_t symbol) {
    if (!expression) return false;
    if (expression->kind() == NodeKind::Symbol) return expression->symbolId() == symbol;
    for (NodePtr child : expression->operands()) {
        if (mentions(Ptr::borrow(child), symbol)) return true;
    }
    return false;
}

/**
 * @brief 按子树结构替换：把 expression 里所有与 pattern 相等的子树换成 replacement。
 *
 * 与 calculus.cpp 里的同名内部实现一致，这里独立一份是因为那个在匿名命名空间内。
 */
Ptr replaceSubtreeLocal(const Ptr &expression, const Ptr &pattern, const Ptr &replacement) {
    if (!expression || !pattern || !replacement) return Ptr();
    if (expression->equals(pattern.get())) return replacement;
    if (expression->kind() == NodeKind::Symbol || expression->kind() == NodeKind::Integer ||
        expression->kind() == NodeKind::Rational || expression->kind() == NodeKind::Constant) {
        return expression;
    }
    Memory::Vector<Ptr> holders;
    Memory::Vector<NodePtr> children;
    bool changed = false;
    for (NodePtr child : expression->operands()) {
        Ptr original = Ptr::borrow(child);
        Ptr next = replaceSubtreeLocal(original, pattern, replacement);
        if (!next) return Ptr();
        if (!next->equals(original.get())) changed = true;
        holders.push_back(next);
        children.push_back(next.get());
    }
    if (!changed) return expression;
    switch (expression->kind()) {
    case NodeKind::Add:
        return makeAddMany(children);
    case NodeKind::Mul:
        return makeMulMany(children);
    case NodeKind::Pow:
        if (children.size() != 2) return Ptr();
        return makePow(holders[0], holders[1]);
    case NodeKind::Function: {
        const char *name = expression->name();
        if (!name) return Ptr();
        return Ptr(Node::makeFunction(name, children));
    }
    default:
        return expression;
    }
}

/**
 * @brief 剥掉表达式的数值系数，返回剩下的「形状」部分。
 *
 * `-e^x` 在外形上是 `Mul(-1, E^x)`，`3*sin(x)` 是 `Mul(3, sin(x))`；判形状前
 * 需要把这些数值因子去掉。若整个表达式就是数值，返回空。
 */
Ptr numericCoefficientStripped(const Ptr &expression) {
    if (!expression) return Ptr();
    if (expression->isNumber()) return Ptr();
    if (expression->kind() != NodeKind::Mul) return expression;
    Memory::Vector<NodePtr> kept;
    Memory::Vector<Ptr> keepers;
    for (NodePtr child : expression->operands()) {
        Ptr factor = Ptr::borrow(child);
        if (factor->isNumber()) continue;
        keepers.push_back(factor);
        kept.push_back(factor.get());
    }
    if (kept.empty()) return Ptr();
    if (kept.size() == expression->size()) return expression;
    return makeMulMany(kept);
}

/** @brief 把表达式摊平成一串加法项，`holders` 负责保活到调用方使用完毕。 */void collectAddTerms(const Ptr &expression, Memory::Vector<NodePtr> *pieces, Memory::Vector<Ptr> *holders) {
    if (!expression) return;
    if (expression->kind() == NodeKind::Add) {
        for (NodePtr child : expression->operands()) {
            collectAddTerms(Ptr::borrow(child), pieces, holders);
        }
        return;
    }
    holders->push_back(expression);
    pieces->push_back(expression.get());
}

/** @brief 表达式子树里是否出现与 `atom` 结构相等的节点。 */
bool mentionsNode(const Ptr &expression, const Ptr &atom) {
    if (!expression || !atom) return false;
    if (expression->equals(atom.get())) return true;
    for (NodePtr child : expression->operands()) {
        if (mentionsNode(Ptr::borrow(child), atom)) return true;
    }
    return false;
}

/**
 * @brief 把表达式按 `atom` 收集成线性组合：成功时满足 `expr = coefficient·atom + rest`，
 *        其中 `rest` 不含 `atom`。
 *
 * 只处理 atom 以"整块因子"形式出现的情形（`atom`、`c·atom`、`a+b·atom`）。
 * 这是从 `y' = f(x)·y + g(x)` 抽出 `f`、`g` 的关键。
 *
 * @param coefficient 输出 atom 的系数。
 * @param rest        输出与 atom 无关的剩余部分（可为 0）。
 */
/** @brief 子树里是否有 Add 节点（用于判断是否需要先展开）。 */
bool containsAdd(const Ptr &expression) {
    if (!expression) return false;
    if (expression->kind() == NodeKind::Add) return true;
    for (NodePtr child : expression->operands()) {
        if (containsAdd(Ptr::borrow(child))) return true;
    }
    return false;
}

bool splitLinear(const Ptr &expression, const Ptr &atom, Ptr *coefficient, Ptr *rest) {
    if (!expression || !atom || !coefficient || !rest) return false;

    // 化简常把 `-y + x` 写成 `-(y - x)`，这种嵌套结构里 atom 被包在 Add 内层，
    // 直接按 Mul 拆会认不出。先展开一次再拆，形态就统一了。
    if (expression->kind() != NodeKind::Add && containsAdd(expression)) {
        Ptr expanded = expand(expression, 20000);
        if (expanded && !expanded->equals(expression.get())) {
            return splitLinear(expanded, atom, coefficient, rest);
        }
    }

    if (expression->equals(atom.get())) {
        *coefficient = integer(1);
        *rest = integer(0);
        return true;
    }

    // 整棵表达式完全不含 atom（例如 `a·y'' + c·y = 0` 里的 `c·y` 对 y' 而言）：
    // 视为「系数为 0、整项落入余项」，否则缺项方程会被判成不支持。
    if (!mentionsNode(expression, atom)) {
        *coefficient = integer(0);
        *rest = expression;
        return true;
    }

    if (expression->kind() == NodeKind::Mul) {
        Memory::Vector<NodePtr> others;
        Memory::Vector<Ptr> holders;
        bool found = false;
        for (NodePtr child : expression->operands()) {
            Ptr factor = Ptr::borrow(child);
            if (!found && factor->equals(atom.get())) {
                found = true;
                continue;
            }
            holders.push_back(factor);
            others.push_back(factor.get());
        }
        if (!found) return false;
        *coefficient = others.empty() ? integer(1) : makeMulMany(others);
        *rest = integer(0);
        return true;
    }

    if (expression->kind() == NodeKind::Add) {
        Memory::Vector<NodePtr> coefficientParts;
        Memory::Vector<NodePtr> restParts;
        Memory::Vector<Ptr> holders;
        for (NodePtr child : expression->operands()) {
            Ptr term = Ptr::borrow(child);
            // 完全不含 atom 的项：系数为 0，整项进 rest
            if (!mentionsNode(term, atom)) {
                holders.push_back(term);
                restParts.push_back(term.get());
                continue;
            }
            Ptr part;
            Ptr remainder;
            if (!splitLinear(term, atom, &part, &remainder)) return false;
            holders.push_back(part);
            holders.push_back(remainder);
            if (!(part->isNumber() && part->numerator().isZero())) coefficientParts.push_back(part.get());
            if (!(remainder->isNumber() && remainder->numerator().isZero())) restParts.push_back(remainder.get());
        }
        *coefficient = coefficientParts.empty() ? integer(0) : makeAddMany(coefficientParts);
        *rest = restParts.empty() ? integer(0) : makeAddMany(restParts);
        return true;
    }

    return false;
}

/* ------------------------------------------------------------ 导数节点 */

/** @brief 求导阶数的最大值，只看 `Derivative(y, x)` 这种节点。 */
int highestOrderOf(const Ptr &expression, std::uint32_t unknown) {
    int highest = 0;
    if (!expression) return highest;
    std::uint32_t variable = 0;
    int level = 1;
    if (Ptr target = asDerivative(expression, &variable, &level)) {
        if (isSymbol(target, unknown) && level > highest) highest = level;
    }
    for (NodePtr child : expression->operands()) {
        const int childOrder = highestOrderOf(Ptr::borrow(child), unknown);
        if (childOrder > highest) highest = childOrder;
    }
    return highest;
}

/** @brief 把 `y' = rhs` 的右端解出来；失败返回空。 */
Ptr isolateFirstDerivative(const Ptr &equation, std::uint32_t unknown, std::uint32_t variable) {
    const Ptr derivative = makeDerivative(symbolNode(unknown), variable, 1);
    if (!derivative) return Ptr();

    Ptr coefficient;
    Ptr rest;
    if (!splitLinear(equation, derivative, &coefficient, &rest)) return Ptr();
    // rest 里出现 y 是正常的（y' = y 就是这种情形），但它不能再含 y' 或更高阶导数，
    // 否则方程对 y' 不是线性的。
    if (mentionsNode(rest, derivative)) return Ptr();
    if (highestOrderOf(rest, unknown) > 0) return Ptr();
    if (coefficient->isNumber() && coefficient->numerator().isZero()) return Ptr();

    // coefficient·y' + rest = 0 ⇒ y' = -rest / coefficient
    return simplify(makeDivide(makeNegate(rest), coefficient));
}

/* ------------------------------------------------------------ 各方法 */

/**
 * @brief 可分离变量：`y' = f(x)·g(y)`。
 *
 * 把右端拆成「含 y 的因子」与「不含 y 的因子」，若含 y 的部分里不再出现 x，
 * 则可分离：`∫dy/g(y) = ∫f(x)dx`。返回隐式解左端减右端。
 */
Ptr trySeparable(const Ptr &rhs, std::uint32_t unknown, std::uint32_t variable, std::string *description) {
    Memory::Vector<Ptr> yFactors;
    Memory::Vector<Ptr> xFactors;
    Memory::Vector<Ptr> holders;

    auto collect = [&](const Ptr &expression) {
        if (expression->kind() == NodeKind::Mul) {
            for (NodePtr child : expression->operands()) {
                Ptr factor = Ptr::borrow(child);
                holders.push_back(factor);
                if (mentions(factor, unknown)) {
                    yFactors.push_back(factor);
                } else {
                    xFactors.push_back(factor);
                }
            }
            return;
        }
        holders.push_back(expression);
        if (mentions(expression, unknown)) {
            yFactors.push_back(expression);
        } else {
            xFactors.push_back(expression);
        }
    };
    collect(rhs);
    // 右端完全不含 y：y' = f(x)，属于可分离的特例（g(y) ≡ 1）
    Memory::Vector<NodePtr> yNodes;
    Memory::Vector<NodePtr> xNodes;
    for (const Ptr &factor : yFactors) yNodes.push_back(factor.get());
    for (const Ptr &factor : xFactors) xNodes.push_back(factor.get());
    Ptr yPart = yNodes.empty() ? integer(1) : makeMulMany(yNodes);
    Ptr xPart = xNodes.empty() ? integer(1) : makeMulMany(xNodes);
    if (!yPart || !xPart) return Ptr();
    if (mentions(yPart, variable)) return Ptr();

    Ptr left = integrate(simplify(makeDivide(integer(1), yPart)), unknown);
    if (!left) return Ptr();
    Ptr right = integrate(xPart, variable);
    if (!right) return Ptr();

    if (description) *description = "可分离变量";
    // 隐式解：∫dy/g(y) − ∫f(x)dx − C1 = 0
    const Ptr constant = symbolNode(internSymbol("C1"));
    return simplify(makeAdd(makeAdd(left, makeNegate(right)), makeNegate(constant)));
}

/**
 * @brief 一阶线性：`y' = P(x)·y + Q(x)`，解 `y = e^{∫P}(∫Q·e^{-∫P} + C1)`。
 *
 * 注意右端里可能还带 y'（`y' + y = x` 整理后是 `y' = -y + x`，这里已消掉），
 * 但更稳妥的做法是把 y' 先替换成独立符号再做线性拆分，避免 y' 被误当作 y。
 */
Ptr tryLinearFirstOrder(const Ptr &rhs, std::uint32_t unknown, std::uint32_t variable,
                        std::string *description) {
    (void)variable;
    // 右端若仍含 y' 的更高阶导数，就不可能是一阶线性
    if (highestOrderOf(rhs, unknown) > 0) return Ptr();

    Ptr p;
    Ptr q;
    if (!splitLinear(rhs, symbolNode(unknown), &p, &q)) return Ptr();
    // p、q 必须与 y 无关（q 里允许含 x，p 里允许含 x）
    if (mentions(p, unknown) || mentions(q, unknown)) return Ptr();

    Ptr pIntegral = integrate(p, variable);
    if (!pIntegral) return Ptr();
    const Ptr exponential = Ptr(Node::makeConstant(ConstantId::E));
    Ptr growth = makePow(exponential, pIntegral);
    if (!growth) return Ptr();
    Ptr decay = makePow(exponential, makeNegate(pIntegral));
    if (!decay) return Ptr();

    Ptr inner = integrate(simplify(makeMul(q, decay)), variable);
    if (!inner) return Ptr();

    const Ptr constant = symbolNode(internSymbol("C1"));
    Ptr solution = simplify(makeMul(growth, makeAdd(inner, constant)));
    if (!solution) return Ptr();

    if (description) *description = "一阶线性（积分因子）";
    // 返回 `solution - y`，由调用方拼成 y = ...
    return makeAdd(solution, makeNegate(symbolNode(unknown)));
}

/**
 * @brief 把表达式按给定的 x、y 取值做数值求值（用于齐次性判定）。
 * @return 成功返回 true 并写入 `value`；出现无法数值化的节点时返回 false。
 */
bool numericValueAt(const Ptr &expression, std::uint32_t variable, std::uint32_t unknown,
                    double xValue, double yValue, double *value) {
    if (!expression) return false;
    switch (expression->kind()) {
    case NodeKind::Integer:
    case NodeKind::Rational:
        *value = expression->numerator().toDouble() / expression->denominator().toDouble();
        return true;
    case NodeKind::Symbol:
        if (expression->symbolId() == variable) {
            *value = xValue;
            return true;
        }
        if (expression->symbolId() == unknown) {
            *value = yValue;
            return true;
        }
        return false;
    case NodeKind::Constant:
        if (expression->constant() == ConstantId::E) {
            *value = 2.71828182845904523536;
            return true;
        }
        if (expression->constant() == ConstantId::Pi) {
            *value = 3.14159265358979323846;
            return true;
        }
        return false;
    case NodeKind::Add: {
        double total = 0.0;
        for (NodePtr child : expression->operands()) {
            double piece = 0.0;
            if (!numericValueAt(Ptr::borrow(child), variable, unknown, xValue, yValue, &piece)) return false;
            total += piece;
        }
        *value = total;
        return true;
    }
    case NodeKind::Mul: {
        double product = 1.0;
        for (NodePtr child : expression->operands()) {
            double piece = 0.0;
            if (!numericValueAt(Ptr::borrow(child), variable, unknown, xValue, yValue, &piece)) return false;
            product *= piece;
        }
        *value = product;
        return true;
    }
    case NodeKind::Pow: {
        double base = 0.0;
        double exponent = 0.0;
        if (!numericValueAt(Ptr::borrow(expression->operand(0)), variable, unknown, xValue, yValue, &base))
            return false;
        if (!numericValueAt(Ptr::borrow(expression->operand(1)), variable, unknown, xValue, yValue, &exponent))
            return false;
        *value = std::pow(base, exponent);
        return true;
    }
    case NodeKind::Function: {
        if (!expression->name() || expression->size() != 1) return false;
        double argument = 0.0;
        if (!numericValueAt(Ptr::borrow(expression->operand(0)), variable, unknown, xValue, yValue, &argument))
            return false;
        const std::string name = expression->name();
        if (name == "exp") *value = std::exp(argument);
        else if (name == "log" || name == "ln") *value = std::log(argument);
        else if (name == "sin") *value = std::sin(argument);
        else if (name == "cos") *value = std::cos(argument);
        else if (name == "tan") *value = std::tan(argument);
        else if (name == "sqrt") *value = std::sqrt(argument);
        else return false;
        return true;
    }
    default:
        return false;
    }
}

/**
 * @brief 判断右端是否是 `y/x` 的 0 次齐次函数。
 *
 * 用两组成比例的取值比较：`F(2,3)` 与 `F(4,6)` 相等即认为齐次。
 * 比符号化替换更稳——符号替换后还要约分公共因子，化简层未必做得到。
 */
bool looksHomogeneous(const Ptr &rhs, std::uint32_t unknown, std::uint32_t variable) {
    if (!mentions(rhs, unknown) || !mentions(rhs, variable)) return false;
    double first = 0.0;
    double second = 0.0;
    if (!numericValueAt(rhs, variable, unknown, 2.0, 3.0, &first)) return false;
    if (!numericValueAt(rhs, variable, unknown, 4.0, 6.0, &second)) return false;
    if (!std::isfinite(first) || !std::isfinite(second)) return false;
    const double scale = std::max(1.0, std::abs(first));
    return std::abs(first - second) <= 1e-9 * scale;
}

/**
 * @brief 齐次型：`y' = F(y/x)`，令 `u = y/x`（即 `y = x·u`），
 *        得 `x·u' = F(u) − u`，分离变量后为 `∫du/(F(u)−u) = ln|x| + C1`。
 */
Ptr tryHomogeneous(const Ptr &rhs, std::uint32_t unknown, std::uint32_t variable, std::string *description) {
    const std::uint32_t ratioSymbol = internSymbol("u");
    Ptr xSymbol = symbolNode(variable);
    Ptr ySymbol = symbolNode(unknown);
    Ptr uSymbol = symbolNode(ratioSymbol);

    // 判定 + 构造 F(u)：先用数值比例检验 0 次齐次，再取 x = 1、y = u 得到 F(u)。
    // 这样避开了"符号替换后还要约公共因子"的麻烦（(x²+y²)/(xy) 就是这种情形）。
    if (!looksHomogeneous(rhs, unknown, variable)) return Ptr();
    const Ptr candidate = rhs;
    Ptr substituted = simplify(substitute(candidate, variable, integer(1)));
    if (!substituted) return Ptr();
    substituted = simplify(substitute(substituted, unknown, uSymbol));
    if (!substituted) return Ptr();
    if (mentions(substituted, variable)) return Ptr();
    if (!mentions(substituted, ratioSymbol)) return Ptr();

    // F(u) − u 常是 (u²+1)/u − u 这种未通分的形式，先通分成 1/u 再积分
    Ptr difference = makeAdd(substituted, makeNegate(uSymbol));
    Ptr combined = together(difference);
    Ptr denominator = simplify(combined ? combined : difference);
    if (!denominator) return Ptr();
    if (denominator->isNumber() && denominator->numerator().isZero()) return Ptr();

    Ptr left = integrate(simplify(makeDivide(integer(1), denominator)), ratioSymbol);
    if (!left) return Ptr();
    Ptr right = integrate(makeDivide(integer(1), xSymbol), variable);
    if (!right) return Ptr();
    // ∫du/(F(u)−u) = ∫dx/x + C，常数项统一记作 C1（与其它解法的记号一致）
    Ptr constant = symbolNode(internSymbol("C1"));
    if (description) *description = "齐次型（换元 u = y/x）";
    // 把 u 回代成 y/x，否则答案里会留一个用户不认识的 u
    Ptr implicit = simplify(makeAdd(makeAdd(left, makeNegate(right)), constant));
    if (!implicit) return Ptr();
    return substitute(implicit, ratioSymbol, makeDivide(ySymbol, xSymbol));
}

/** @brief 多项式系数（整数次幂）的指数：`y^p` 返回 p，否则返回 false。 */
bool powerOfUnknown(const Ptr &expression, std::uint32_t unknown, long long *power, Ptr *coefficient) {
    if (expression->kind() == NodeKind::Pow && isSymbol(Ptr::borrow(expression->operand(0)), unknown) &&
        Ptr::borrow(expression->operand(1))->isSmallInteger(power)) {
        *coefficient = integer(1);
        return true;
    }
    if (expression->kind() == NodeKind::Mul) {
        Memory::Vector<NodePtr> others;
        Memory::Vector<Ptr> holders;
        bool found = false;
        for (NodePtr child : expression->operands()) {
            Ptr factor = Ptr::borrow(child);
            if (!found && factor->kind() == NodeKind::Pow &&
                isSymbol(Ptr::borrow(factor->operand(0)), unknown) &&
                Ptr::borrow(factor->operand(1))->isSmallInteger(power)) {
                found = true;
                continue;
            }
            holders.push_back(factor);
            others.push_back(factor.get());
        }
        if (!found) return false;
        *coefficient = others.empty() ? integer(1) : makeMulMany(others);
        return true;
    }
    return false;
}

/**
 * @brief 伯努利方程：`y' = P(x)·y + Q(x)·y^n`。
 *
 * 换元 `v = y^(1-n)` 化为一阶线性：`v' = (1-n)P·v + (1-n)Q`。
 */
Ptr tryBernoulli(const Ptr &incomingRhs, std::uint32_t unknown, std::uint32_t variable,
                 std::string *description) {
    // 化简有时把 `y - x·y²` 写成 `-(x·y² - y)`，外层是 Mul 而不是 Add；
    // 先展开一次，形态才统一。
    Ptr rhs = incomingRhs;
    if (rhs->kind() != NodeKind::Add && containsAdd(rhs)) {
        Ptr expanded = expand(rhs, 20000);
        if (expanded) rhs = expanded;
    }
    if (rhs->kind() != NodeKind::Add) return Ptr();

    long long degree = 0;
    Ptr q = Ptr();
    Ptr p = Ptr();
    bool haveQ = false;
    bool haveP = false;

    for (NodePtr child : rhs->operands()) {
        Ptr term = Ptr::borrow(child);
        long long termPower = 0;
        Ptr coefficient;
        if (powerOfUnknown(term, unknown, &termPower, &coefficient)) {
            if (haveQ || termPower < 2) return Ptr();
            degree = termPower;
            q = coefficient;
            haveQ = true;
            continue;
        }
        Ptr linear;
        Ptr rest;
        if (!splitLinear(term, symbolNode(unknown), &linear, &rest)) return Ptr();
        if (mentions(linear, unknown) || mentions(rest, unknown)) return Ptr();
        // 完全不含 y 的独立项（例如 `y' = x² + y²` 里的 `x²`）说明右端不是
        // `P·y + Q·yⁿ` 的形状，直接拒绝，不能当成「P = 0」混过去。
        if (!(rest->isNumber() && rest->numerator().isZero())) return Ptr();
        if (haveP) return Ptr();
        p = linear;
        haveP = true;
    }
    if (!haveQ || !haveP) return Ptr();

    const long long factorValue = 1 - degree;
    const Ptr scale = integer(factorValue);
    Ptr pPrime = simplify(makeMul(scale, p));
    Ptr qPrime = simplify(makeMul(scale, q));

    Ptr pIntegral = integrate(pPrime, variable);
    if (!pIntegral) return Ptr();
    const Ptr exponential = Ptr(Node::makeConstant(ConstantId::E));
    Ptr growth = makePow(exponential, pIntegral);
    Ptr decay = makePow(exponential, makeNegate(pIntegral));
    if (!growth || !decay) return Ptr();
    Ptr inner = integrate(simplify(makeMul(qPrime, decay)), variable);
    if (!inner) return Ptr();
    const Ptr constant = symbolNode(internSymbol("C1"));
    Ptr vSolution = simplify(makeMul(growth, makeAdd(inner, constant)));
    if (!vSolution) return Ptr();

    // v = y^(1-n) ⇒ y = v^(1/(1-n))
    Ptr recovered = makePow(vSolution, makeNumber(BigInt(1), BigInt(factorValue)));
    if (!recovered) return Ptr();
    if (description) *description = "伯努利方程（换元 v = y^(1-n)）";
    return makeAdd(recovered, makeNegate(symbolNode(unknown)));
}

} // namespace

/* ---------------------------------------------------------- 导数节点 */

Ptr makeDerivative(const Ptr &function, std::uint32_t variable, int order) {
    if (!function) return Ptr();
    Ptr variableNode = symbolNode(variable);
    Ptr orderNode = integer(order);
    // 阶数总是显式写出（`Derivative(y, x, 1)`），与解析器对 `y'` 的产物保持一致，
    // 否则 equals 比较会因参数个数不同而失败。
    Memory::Vector<NodePtr> arguments{function.get(), variableNode.get(), orderNode.get()};
    Memory::Vector<Ptr> holders{variableNode, orderNode};
    return Ptr(Node::makeFunction("Derivative", arguments));
}

Ptr asDerivative(const Ptr &node, std::uint32_t *variable, int *order) {
    if (!node || node->kind() != NodeKind::Function) return Ptr();
    if (!node->name() || std::string(node->name()) != "Derivative") return Ptr();
    if (node->size() < 2) return Ptr();
    Ptr variableNode = Ptr::borrow(node->operand(1));
    if (variableNode->kind() != NodeKind::Symbol) return Ptr();
    if (variable) *variable = variableNode->symbolId();
    int level = 1;
    if (node->size() >= 3) {
        long long value = 1;
        if (Ptr::borrow(node->operand(2))->isSmallInteger(&value)) level = static_cast<int>(value);
    }
    if (order) *order = level;
    return Ptr::borrow(node->operand(0));
}

bool containsUnknown(const Ptr &expression, std::uint32_t unknown) {
    if (!expression) return false;
    if (expression->kind() == NodeKind::Symbol) return expression->symbolId() == unknown;
    std::uint32_t variable = 0;
    if (Ptr target = asDerivative(expression, &variable, nullptr)) {
        if (isSymbol(target, unknown)) return true;
    }
    for (NodePtr child : expression->operands()) {
        if (containsUnknown(Ptr::borrow(child), unknown)) return true;
    }
    return false;
}

int derivativeOrder(const Ptr &equation, std::uint32_t unknown) {
    if (!equation) return -1;
    if (!containsUnknown(equation, unknown)) return -1;
    return highestOrderOf(equation, unknown);
}

/* ---------------------------------------------------- 二阶常系数 */

/**
 * @brief 把 `a·y'' + b·y' + c·y = f(x)` 的系数抽出来。
 *
 * 只处理系数为常数（不含 x）的情形。
 * @return 成功返回 true，并写入 a、b、c 与右端 f。
 */
bool splitConstantCoefficients(const Ptr &equation, std::uint32_t unknown, std::uint32_t variable,
                               Ptr *a, Ptr *b, Ptr *c, Ptr *forcing) {
    const Ptr second = makeDerivative(symbolNode(unknown), variable, 2);
    const Ptr first = makeDerivative(symbolNode(unknown), variable, 1);
    if (!second || !first) return false;

    Ptr coefficientA;
    Ptr restAfterA;
    Ptr expandedEquation = expand(equation, 20000);
    const Ptr source = expandedEquation ? expandedEquation : equation;
    if (!splitLinear(source, second, &coefficientA, &restAfterA)) return false;

    Ptr coefficientB;
    Ptr restAfterB;
    if (!splitLinear(restAfterA, first, &coefficientB, &restAfterB)) return false;

    Ptr coefficientC;
    Ptr restAfterC;
    if (!splitLinear(restAfterB, symbolNode(unknown), &coefficientC, &restAfterC)) return false;

    // 残留部分必须与 y 无关（否则不是常系数线性方程）
    if (containsUnknown(restAfterC, unknown)) return false;
    if (mentions(coefficientA, variable) || mentions(coefficientB, variable) || mentions(coefficientC, variable))
        return false;

    *a = simplify(coefficientA);
    *b = simplify(coefficientB);
    *c = simplify(coefficientC);
    // 方程是以 `a·y'' + b·y' + c·y + rest = 0` 的形式给出的，
    // 所以右端 f(x) = -rest，符号必须翻一次，否则特解会整体取反。
    Ptr negatedRest = makeNegate(restAfterC);
    if (!negatedRest) return false;
    *forcing = simplify(negatedRest);
    return *a && *b && *c && *forcing;
}

/** @brief 解二次特征方程 `a·r² + b·r + c = 0`，返回判别式与根。 */
struct CharacteristicRoots {
    Ptr discriminant;   ///< Δ = b² − 4ac
    Ptr root1;          ///< (−b + √Δ)/(2a)
    Ptr root2;          ///< (−b − √Δ)/(2a)
    bool complexPair;   ///< Δ < 0 时为共轭复根（此时 root1/root2 的实部与虚部分别给出）
    Ptr realPart;       ///< 复根情形：实部 −b/(2a)
    Ptr imaginaryPart;  ///< 复根情形：虚部 √(−Δ)/(2a)
};

CharacteristicRoots solveCharacteristic(const Ptr &a, const Ptr &b, const Ptr &c) {
    CharacteristicRoots roots;
    Ptr discriminant = simplify(makeAdd(makePow(b, integer(2)), makeNegate(makeMul(integer(4), makeMul(a, c)))));
    roots.discriminant = discriminant;
    Ptr twoA = makeMul(integer(2), a);
    Ptr negB = makeNegate(b);

    const bool negative = discriminant->isNumber() && discriminant->numerator().isNegative();
    roots.complexPair = negative;
    if (negative) {
        roots.realPart = simplify(makeDivide(negB, twoA));
        Ptr magnitude = makePow(makeNegate(discriminant), makeNumber(BigInt(1), BigInt(2)));
        roots.imaginaryPart = simplify(makeDivide(magnitude, twoA));
        roots.root1 = roots.realPart;
        roots.root2 = roots.imaginaryPart;
        return roots;
    }
    Ptr root = makePow(discriminant, makeNumber(BigInt(1), BigInt(2)));
    roots.root1 = simplify(makeDivide(makeAdd(negB, root), twoA));
    roots.root2 = simplify(makeDivide(makeAdd(negB, makeNegate(root)), twoA));
    return roots;
}

/** @brief 齐次解：按特征根的三种情形拼出 `C1·… + C2·…`。 */
Ptr homogeneousSolution(const Ptr &a, const Ptr &b, const Ptr &c, std::uint32_t variable,
                        std::string *caseDescription) {
    const CharacteristicRoots roots = solveCharacteristic(a, b, c);
    Ptr x = symbolNode(variable);
    Ptr c1 = symbolNode(internSymbol("C1"));
    Ptr c2 = symbolNode(internSymbol("C2"));

    if (roots.complexPair) {
        // y = e^{αx}(C1 cos βx + C2 sin βx)
        Ptr exponential = makePow(Ptr(Node::makeConstant(ConstantId::E)),
                                  makeMul(roots.realPart, x));
        Ptr cosine = Ptr(Node::makeFunction("cos", Memory::Vector<NodePtr>{
                                                  makeMul(roots.imaginaryPart, x).get()}));
        Ptr sine = Ptr(Node::makeFunction("sin", Memory::Vector<NodePtr>{
                                                makeMul(roots.imaginaryPart, x).get()}));
        if (caseDescription) *caseDescription = "二阶常系数齐次（共轭复根）";
        return simplify(makeMul(exponential, makeAdd(makeMul(c1, cosine), makeMul(c2, sine))));
    }

    Ptr first = makePow(Ptr(Node::makeConstant(ConstantId::E)), makeMul(roots.root1, x));
    if (roots.root1->equals(roots.root2.get())) {
        // 重根：y = (C1 + C2·x)e^{rx}
        if (caseDescription) *caseDescription = "二阶常系数齐次（重根）";
        return simplify(makeMul(makeAdd(c1, makeMul(c2, x)), first));
    }
    Ptr second = makePow(Ptr(Node::makeConstant(ConstantId::E)), makeMul(roots.root2, x));
    if (caseDescription) *caseDescription = "二阶常系数齐次（两个相异实根）";
    return simplify(makeAdd(makeMul(c1, first), makeMul(c2, second)));
}

/**
 * @brief 非齐次特解（待定系数法）。
 *
 * 支持的自由项形状：
 * - 多项式 `P(x)`
 * - `A·e^{kx}`
 * - `A·cos(kx)`、`A·sin(kx)`
 * - 多项式 × 指数（`P(x)e^{kx}`）
 *
 * 用"试探函数代入、比较系数"实现：把试探式写成符号组合，代入
 * `a·y'' + b·y' + c·y` 后与 forcing 逐项比对，解出待定常数。
 *
 * @return 找不到合适试探式时返回空。
 */
Ptr particularSolution(const Ptr &a, const Ptr &b, const Ptr &c, const Ptr &forcing,
                       std::uint32_t variable) {
    if (!forcing) return Ptr();
    Ptr simplifiedForcing = simplify(forcing);
    if (!simplifiedForcing) return Ptr();
    // 右端为 0：无需求特解
    if (simplifiedForcing->isZero()) return integer(0);
    if (mentions(simplifiedForcing, internSymbol("y"))) return Ptr();

    Ptr x = symbolNode(variable);

    // 试探式的"形状"由 forcing 决定：取 forcing 里的函数部分与多项式部分
    // 这里按常见三类分别处理，系数用未定符号 A、B 表示，代入后解方程。

    auto undeterminedCoefficients = [&](const Memory::Vector<Ptr> &trialTerms,
                                        const Memory::Vector<std::uint32_t> &names) -> Ptr {
        // trial = Σ A_i · trialTerms[i]
        Memory::Vector<Ptr> coefficientSymbols;
        for (std::uint32_t name : names) coefficientSymbols.push_back(symbolNode(name));
        Memory::Vector<NodePtr> trialPieces;
        Memory::Vector<Ptr> trialHolders;  // 保活：trialPieces 里只是裸指针
        for (std::size_t i = 0; i < trialTerms.size(); ++i) {
            Ptr product = makeMul(coefficientSymbols[i], trialTerms[i]);
            if (!product) return Ptr();
            trialHolders.push_back(product);
            trialPieces.push_back(product.get());
        }
        Ptr trial = trialPieces.empty() ? Ptr() : makeAddMany(trialPieces);
        if (!trial) return Ptr();

        // plug = a·trial'' + b·trial' + c·trial − forcing
        // 这里刻意只做 expand 不做 simplify：simplify 会把 A_j 连同其基函数一起
        // 合并成更紧的形式（如 `A0 + A1*x` 压成单个节点序列），按基函数逐项
        // 比对时反而不好拆。expand 保留「若干项相加」的形态，正好逐项对齐。
        Ptr firstDerivative = differentiate(trial, variable);
        Ptr secondDerivative = differentiate(firstDerivative, variable);
        if (!firstDerivative || !secondDerivative) return Ptr();
        Ptr raw = makeAdd(makeAdd(makeMul(a, secondDerivative), makeMul(b, firstDerivative)),
                          makeAdd(makeMul(c, trial), makeNegate(simplifiedForcing)));
        if (!raw) return Ptr();
        Ptr plug = expand(raw, 20000);
        if (!plug) plug = raw;

        // 解线性方程组。
        // plug 恒为 0，而各 trialTerm 彼此线性无关，所以「每个 trialTerm 的系数
        // 都为 0」与「把 A_j 视作自变量、plug 作为线性式恒为 0」等价。
        //
        // 做法：把 plug 摊平成加法项，逐项「除以某个 trialTerm」；除尽且结果不含 x
        // 的那部分就是该 trialTerm 的系数（是各 A_j 的线性组合）。于是每个 trialTerm
        // 给出一条方程，交给符号矩阵求解器解 A_j。
        const std::size_t unknownCount = coefficientSymbols.size();

        Memory::Vector<NodePtr> pieces;
        Memory::Vector<Ptr> holders;
        collectAddTerms(plug, &pieces, &holders);

        // 比对方法：把 plug 摊平成一串加法项，对每个基函数 basis_r，
        // 收集「除以 basis_r 后不含 x」的项——这些就是 basis_r 的系数（各 A_j 的线性组合）。
        // 各基函数线性无关 ⇒ 每个系数都必须为 0，于是每个 basis_r 给出一条方程。
        //
        // 基函数表同时纳入两类：
        // 1) 试探式自身的 trialTerm（乘过 x 的共振项也在此列）；
        // 2) plug 里出现的「不含任何待定系数的项」——它们代表试探式覆盖不到的
        //    基函数，必须令其系数为 0 才不会漏掉方程。
        Memory::Vector<Ptr> basisList;
        auto containsBasis = [&](const Ptr &candidate) {
            for (const Ptr &existing : basisList) {
                if (existing->equals(candidate.get())) return true;
            }
            return false;
        };
        for (const Ptr &basis : trialTerms) {
            if (!containsBasis(basis)) basisList.push_back(basis);
        }
        for (std::size_t k = 0; k < pieces.size(); ++k) {
            Ptr term = Ptr::borrow(pieces[k]);
            bool hasUndetermined = false;
            for (const Ptr &symbol : coefficientSymbols) {
                if (mentions(term, symbol->symbolId())) {
                    hasUndetermined = true;
                    break;
                }
            }
            if (hasUndetermined) continue;
            // 含 symbol 之外的自由符号（理论上不该有）的项跳过
            if (!containsBasis(term)) basisList.push_back(term);
        }
        // 试探式里通常含常数基（多项式情形必有）；若没有，plug 残留的常数项会被漏掉
        bool hasConstantBasis = false;
        for (const Ptr &basis : basisList) {
            if (basis->isNumber() && basis->isOne()) hasConstantBasis = true;
        }
        if (!hasConstantBasis) basisList.push_back(integer(1));

        const std::size_t rowCount = basisList.size();
        Memory::Vector<Ptr> entries;
        Memory::Vector<Ptr> rhsEntries;
        for (std::size_t r = 0; r < rowCount; ++r) {
            const Ptr &basis = basisList[r];
            Memory::Vector<Ptr> row(unknownCount, integer(0));
            Ptr constantColumn = integer(0);
            for (std::size_t k = 0; k < pieces.size(); ++k) {
                Ptr term = Ptr::borrow(pieces[k]);
                Ptr quotient = simplify(makeDivide(term, basis));
                if (!quotient) return Ptr();
                // 商里仍含 x ⇒ 这一项不属于 basis_r
                if (mentions(quotient, variable)) continue;
                bool touched = false;
                for (std::size_t j = 0; j < unknownCount; ++j) {
                    const std::uint32_t id = coefficientSymbols[j]->symbolId();
                    if (!mentions(quotient, id)) continue;
                    Ptr linear;
                    Ptr remainder;
                    if (!splitLinear(quotient, coefficientSymbols[j], &linear, &remainder)) return Ptr();
                    if (!remainder->isZero()) return Ptr();
                    if (mentions(linear, id)) return Ptr();
                    row[j] = simplify(makeAdd(row[j], linear)).share();
                    touched = true;
                }
                if (!touched) constantColumn = simplify(makeAdd(constantColumn, quotient)).share();
            }
            for (std::size_t j = 0; j < unknownCount; ++j) {
                if (!row[j]) return Ptr();
                entries.push_back(row[j]);
            }
            rhsEntries.push_back(simplify(makeNegate(constantColumn)));
        }

        SymbolicMatrix coefficients = SymbolicMatrix::fromRows(rowCount, unknownCount, entries);
        SymbolicMatrix rightHandSide = SymbolicMatrix::fromRows(rowCount, 1, rhsEntries);
        Memory::Vector<Ptr> solutionValues;
        bool consistent = false;
        if (!solveLinearSystemMatrix(coefficients, rightHandSide, &solutionValues, &consistent)) {
            // 退化（试探式与齐次解共振）时本函数放弃
            return Ptr();
        }
        if (solutionValues.size() != unknownCount) return Ptr();
        Memory::Vector<Ptr> solved;
        for (const Ptr &value : solutionValues) {
            if (!value || mentions(value, variable)) return Ptr();
            solved.push_back(simplify(value));
        }
        if (solved.empty()) return Ptr();
        Memory::Vector<NodePtr> resultPieces;
        Memory::Vector<Ptr> resultHolders;
        for (std::size_t i = 0; i < trialTerms.size(); ++i) {
            Ptr product = makeMul(solved[i], trialTerms[i]);
            if (!product) return Ptr();
            resultHolders.push_back(product);
            resultPieces.push_back(product.get());
        }
        return simplify(makeAddMany(resultPieces));
    };

    // 情形 1：纯多项式
    UnivariatePoly polynomial;
    if (toPolynomial(simplifiedForcing.get(), variable, &polynomial) && !polynomial.isZero()) {
        Memory::Vector<Ptr> terms;
        Memory::Vector<std::uint32_t> names;
        for (std::size_t k = 0; k < polynomial.coefficients.size(); ++k) {
            terms.push_back(makePow(x, integer(static_cast<long long>(k))));
            names.push_back(internSymbol("A" + std::to_string(k)));
        }
        return undeterminedCoefficients(terms, names);
    }

    // 情形 2：A·e^{kx}（k 可为数值）
    //
    // 先剥掉可能存在的数值系数（`-e^x` 的外形是 Mul(-1, Pow(E,x))），
    // 再判断核心是否为 e^{kx}。
    if (Ptr core = numericCoefficientStripped(simplifiedForcing)) {
        const bool isExponential =
            core->kind() == NodeKind::Pow &&
            Ptr::borrow(core->operand(0))->kind() == NodeKind::Constant &&
            Ptr::borrow(core->operand(0))->constant() == ConstantId::E;
        if (isExponential) {
            Ptr exponent = Ptr::borrow(core->operand(1));
            Ptr kCoefficient;
            Ptr kRest;
            if (splitLinear(exponent, x, &kCoefficient, &kRest)) {
                if (kRest->isNumber() && kRest->numerator().isZero() && !mentions(kCoefficient, variable)) {
                    // 试探式 A·e^{kx}；若 e^{kx} 本身是齐次解（k 是特征根），
                    // 需要依次乘 x（重根时还要再乘一次）。
                    Ptr exponential = core;
                    Ptr trialTerm = exponential;
                    for (int power = 0; power < 2; ++power) {
                        Ptr first = differentiate(trialTerm, variable);
                        Ptr second = differentiate(first, variable);
                        if (!first || !second) return Ptr();
                        Ptr applied = simplify(makeAdd(
                            makeAdd(makeMul(a, second), makeMul(b, first)), makeMul(c, trialTerm)));
                        if (!applied) return Ptr();
                        if (!applied->isZero()) {
                            return undeterminedCoefficients(
                                Memory::Vector<Ptr>{trialTerm},
                                Memory::Vector<std::uint32_t>{internSymbol("A")});
                        }
                        trialTerm = makeMul(x, trialTerm);
                        if (!trialTerm) return Ptr();
                    }
                    return Ptr();  // 三重根不是二阶方程能出现的情形
                }
            }
        }
    }

    // 情形 3：cos(kx) 或 sin(kx)（可带数值系数）
    if (Ptr core = numericCoefficientStripped(simplifiedForcing)) {
        if (core->kind() == NodeKind::Function && core->size() == 1) {
            const char *name = core->name();
            const bool trigonometric = name && (std::string(name) == "cos" || std::string(name) == "sin");
            if (trigonometric) {
                Ptr argument = Ptr::borrow(core->operand(0));
                Ptr kCoefficient;
                Ptr kRest;
                if (splitLinear(argument, x, &kCoefficient, &kRest)) {
                    if (kRest->isNumber() && kRest->numerator().isZero() && !mentions(kCoefficient, variable)) {
                        Memory::Vector<NodePtr> cosineArgs{argument.get()};
                        Memory::Vector<NodePtr> sineArgs{argument.get()};
                        Ptr cosine = Ptr(Node::makeFunction("cos", cosineArgs));
                        Ptr sine = Ptr(Node::makeFunction("sin", sineArgs));

                        // 共振判定：cos(kx)、sin(kx) 若满足齐次方程，试探式整体乘 x
                        Ptr scale = integer(1);
                        Ptr probe = cosine;
                        Ptr first = differentiate(probe, variable);
                        Ptr second = differentiate(first, variable);
                        if (!first || !second) return Ptr();
                        Ptr applied = simplify(makeAdd(
                            makeAdd(makeMul(a, second), makeMul(b, first)), makeMul(c, probe)));
                        if (!applied) return Ptr();
                        if (applied->isZero()) {
                            scale = x;
                            // 重根时 cos/sin 会各共振一次，二阶方程最多再乘一次 x
                            Ptr scaledCosine = makeMul(x, cosine);
                            Ptr scaledFirst = differentiate(scaledCosine, variable);
                            Ptr scaledSecond = differentiate(scaledFirst, variable);
                            if (!scaledFirst || !scaledSecond) return Ptr();
                            Ptr twiceApplied = simplify(makeAdd(
                                makeAdd(makeMul(a, scaledSecond), makeMul(b, scaledFirst)),
                                makeMul(c, scaledCosine)));
                            if (twiceApplied && twiceApplied->isZero()) {
                                scale = makePow(x, integer(2));
                            }
                        }
                        if (!scale) return Ptr();
                        Ptr trialCosine = scale->isOne() ? cosine : makeMul(scale, cosine);
                        Ptr trialSine = scale->isOne() ? sine : makeMul(scale, sine);
                        if (!trialCosine || !trialSine) return Ptr();
                        return undeterminedCoefficients(
                            Memory::Vector<Ptr>{trialCosine, trialSine},
                            Memory::Vector<std::uint32_t>{internSymbol("A"), internSymbol("B")});
                    }
                }
            }
        }
    }

    return Ptr();
}

/* ---------------------------------------------------------- 对外入口 */

Ptr solveOde(const Ptr &equation, std::uint32_t unknown, std::uint32_t variable,
             std::string *description, std::string *message) {
    if (!equation) {
        if (message) *message = "方程为空的";
        return Ptr();
    }
    const int order = derivativeOrder(equation, unknown);
    if (order <= 0) {
        if (message) *message = "方程里没有出现 y 或其导数";
        return Ptr();
    }
    if (order > 2) {
        if (message) *message = "只支持一阶与二阶常微分方程";
        return Ptr();
    }
    if (order == 2) {
        Ptr a;
        Ptr b;
        Ptr c;
        Ptr forcing;
        if (!splitConstantCoefficients(equation, unknown, variable, &a, &b, &c, &forcing)) {
            if (message) *message = "二阶方程目前只支持常系数线性形式 a·y'' + b·y' + c·y = f(x)";
            return Ptr();
        }
        std::string homogeneousCase;
        Ptr homogeneous = homogeneousSolution(a, b, c, variable, &homogeneousCase);
        if (!homogeneous) {
            if (message) *message = "特征方程求解失败";
            return Ptr();
        }
        Ptr particular = particularSolution(a, b, c, forcing, variable);
        if (!particular) {
            if (message) {
                *message = "齐次部分已解出（" + homogeneousCase + "），但非齐次项的形状超出待定系数法范围：" +
                           toText(forcing.get());
            }
            // 只给出齐次解仍然有价值
            if (description) *description = homogeneousCase + "；非齐次特解未求出";
            Ptr negatedUnknown = makeNegate(symbolNode(unknown));
            return makeAdd(homogeneous, negatedUnknown);
        }
        if (description) *description = homogeneousCase;
        Ptr negatedUnknown = makeNegate(symbolNode(unknown));
        Ptr combined = makeAdd(homogeneous, particular);
        return makeAdd(combined, negatedUnknown);
    }

    // 化简常把 `y - x·y²` 写成 `-(x·y² - y)`，外层是 Mul 而不是 Add；
    // 各识别器都按"若干项相加"的结构工作，这里统一展开一次。
    Ptr rhs = isolateFirstDerivative(equation, unknown, variable);
    if (!rhs) {
        if (message) *message = "无法把方程整理成 y' = f(x, y) 的形式（可能不是对 y' 线性的）";
        return Ptr();
    }
    if (rhs->kind() != NodeKind::Add && containsAdd(rhs)) {
        Ptr expanded = expand(rhs, 20000);
        if (expanded) rhs = expanded;
    }

    // 顺序：线性最通用且解式最干净；伯努利要求右端是"y 的一次项 + y^n 项"的特定结构，
    // 与线性互斥；齐次与可分离放最后。
    if (Ptr solution = tryLinearFirstOrder(rhs, unknown, variable, description)) return solution;
    if (Ptr solution = tryBernoulli(rhs, unknown, variable, description)) return solution;
    if (Ptr solution = tryHomogeneous(rhs, unknown, variable, description)) return solution;
    if (Ptr solution = trySeparable(rhs, unknown, variable, description)) return solution;

    if (message) *message = "该方程超出当前支持的一阶类型（可分离/线性/齐次/伯努利）";
    return Ptr();
}

/* ---------------------------------------------------------- 初值问题 */

/**
 * @brief 解析一条初值条件，取出导数阶数、取值点与右端值。
 *
 * 接受 `y(c) = v` 与 `y'(c) = v`（解析器把 `y(0)` 建成函数节点、`y'(0)` 建成
 * 对 `Derivative` 的调用，这里两种形状都认）。
 *
 * @return 解析成功返回 true。
 */
bool parseInitialCondition(const Ptr &condition, std::uint32_t unknown, std::uint32_t *order,
                           Ptr *point, Ptr *value) {
    if (!condition || !order || !point || !value) return false;

    // 条件总是 `lhs = rhs` 的形状
    Ptr left;
    Ptr right;
    if (condition->kind() == NodeKind::Function && condition->name() &&
        std::string(condition->name()) == "equal") {
        if (condition->size() != 2) return false;
        left = Ptr::borrow(condition->operand(0));
        right = Ptr::borrow(condition->operand(1));
    } else {
        return false;
    }
    if (!left || !right) return false;

    const char *unknownName = symbolName(unknown);
    if (!unknownName) return false;

    // 形如 `y(0)`：函数名就是未知函数
    if (left->kind() == NodeKind::Function && left->name() && left->size() == 1 &&
        std::string(left->name()) == std::string(unknownName)) {
        *order = 0;
        *point = Ptr::borrow(left->operand(0));
        *value = right;
        return static_cast<bool>(*point);
    }

    // 形如 `Derivative(y(0), x, n)`：解析器对 `y'(0)` 的产物
    if (left->kind() == NodeKind::Function && left->name() && left->size() == 3 &&
        std::string(left->name()) == "Derivative") {
        Ptr inner = Ptr::borrow(left->operand(0));
        Ptr orderNode = Ptr::borrow(left->operand(2));
        long long parsedOrder = 0;
        if (!orderNode || !orderNode->isSmallInteger(&parsedOrder) || parsedOrder < 0) return false;
        if (inner->kind() == NodeKind::Function && inner->name() && inner->size() == 1 &&
            std::string(inner->name()) == std::string(unknownName)) {
            *order = static_cast<std::uint32_t>(parsedOrder);
            *point = Ptr::borrow(inner->operand(0));
            *value = right;
            return static_cast<bool>(*point);
        }
        return false;
    }
    return false;
}

Ptr solveOdeWithInitialValues(const Ptr &equation, std::uint32_t unknown, std::uint32_t variable,
                              const Memory::Vector<Ptr> &conditions, std::string *description,
                              std::string *message) {
    if (!equation) {
        if (message) *message = "方程为空的";
        return Ptr();
    }

    const int order = derivativeOrder(equation, unknown);
    if (order <= 0) {
        if (message) *message = "方程里没有出现 y 或其导数";
        return Ptr();
    }
    if (order > 2) {
        if (message) *message = "只支持一阶与二阶常微分方程";
        return Ptr();
    }

    // 无初值时直接给通解，但要统一成「y = …」的右端：显式解补上 `+y`，
    // 隐式解原样返回（调用方会按 `= 0` 打印）。
    if (conditions.empty()) {
        std::string generalMethod;
        Ptr general = solveOde(equation, unknown, variable, &generalMethod, message);
        if (!general) return Ptr();
        if (description) *description = generalMethod;
        const bool asExplicit = generalMethod.find("线性") != std::string::npos ||
                                generalMethod.find("二阶") != std::string::npos;
        if (!asExplicit) return general;
        return simplify(makeAdd(general, symbolNode(unknown)));
    }

    std::string method;
    Ptr general = solveOde(equation, unknown, variable, &method, message);
    if (!general) return Ptr();

    // 通解有两种形态：
    // - 显式（线性 / 二阶常系数）：返回的是 `… - y`，补回 y 才是 `y = …`
    // - 隐式（可分离 / 齐次型 / 伯努利）：本身就是 `F(x, y) = 0`
    const bool explicitForm = method.find("线性") != std::string::npos || method.find("二阶") != std::string::npos;
    if (!explicitForm) {
        if (message) *message = "隐式解（" + method + "）里 y 与 x 缠在一起，无法直接套用初值条件";
        return Ptr();
    }

    Ptr unknownSymbol = symbolNode(unknown);
    // solveOde 的显式解形态是 `… - y`，因此这里加回 `+y` 得到 `y = …`
    Ptr rightHandSide = simplify(makeAdd(general, unknownSymbol));
    if (!rightHandSide) return Ptr();

    // 显式解的右端不应再含未知函数；若含，说明求解器内部把某个条件和方程混在一起了
    if (containsUnknown(rightHandSide, unknown)) {
        if (message) *message = "求得的显式解仍含 y，无法代入初值条件";
        return Ptr();
    }

    // 条件是逐条代入的，因此二阶方程要先给 `y(x0)`、再给 `y'(x0)`；
    // 每条条件直接解出一个还没定的常数。
    Memory::Vector<Ptr> constants;
    constants.push_back(symbolNode(internSymbol("C1")));
    if (order >= 2) constants.push_back(symbolNode(internSymbol("C2")));

    std::size_t conditionCount = 0;
    for (const Ptr &condition : conditions) {
        std::uint32_t conditionOrder = 0;
        Ptr point;
        Ptr value;
        if (!parseInitialCondition(condition, unknown, &conditionOrder, &point, &value)) {
            if (message) *message = "无法识别初值条件：" + toText(condition.get());
            return Ptr();
        }
        if (static_cast<int>(conditionOrder) > order) {
            if (message) *message = "初值条件的阶数超过了方程的阶数";
            return Ptr();
        }

        // 取 y 的 conditionOrder 阶导数，再在取值点求值
        Ptr target = rightHandSide;
        for (std::uint32_t k = 0; k < conditionOrder; ++k) {
            target = differentiate(target, variable);
            if (!target) return Ptr();
        }
        Ptr substituted = substitute(target, variable, point);
        if (!substituted) return Ptr();
        Ptr residual = simplify(makeAdd(substituted, makeNegate(value)));
        if (!residual) return Ptr();

        // residual 是各常数的线性式：取第一个还没确定的常数的系数与其余部分
        bool solvedOne = false;
        for (const Ptr &constant : constants) {
            if (!mentions(rightHandSide, constant->symbolId())) continue;  // 已解出
            Ptr linear;
            Ptr rest;
            if (!splitLinear(residual, constant, &linear, &rest)) continue;
            if (linear->isNumber() && linear->numerator().isZero()) continue;

            // 其余项（把该常数置 0 后的残差）作为方程右端
            Ptr constantPart = substitute(residual, constant->symbolId(), integer(0));
            if (!constantPart) return Ptr();
            constantPart = simplify(constantPart);
            if (!constantPart) return Ptr();

            Ptr solvedValue = simplify(makeDivide(makeNegate(constantPart), linear));
            if (!solvedValue || mentions(solvedValue, variable)) return Ptr();

            rightHandSide = simplify(substitute(rightHandSide, constant->symbolId(), solvedValue));
            if (!rightHandSide) return Ptr();
            solvedOne = true;
            break;
        }
        if (!solvedOne) {
            if (message) *message = "初值条件无法确定积分常数（条件可能与方程不相容）";
            return Ptr();
        }
        ++conditionCount;
    }

    if (description) *description = method + "（已代入初值）";
    (void)conditionCount;
    return rightHandSide;
}



} // namespace Cas::Detail

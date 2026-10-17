// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#include "cas/detail/solver.h"
#include "cas/detail/calculus.h"
#include "cas/detail/multivariate.h"
#include "cas/detail/polynomial.h"
#include "cas/detail/printer.h"

#include <algorithm>
#include <cstdio>
#include <cstring>

namespace Cas::Detail {
namespace {

/** @brief 快捷构造整数。 */
Ptr integer(long long value) { return Ptr(Node::makeInteger(BigInt(value))); }

/**
 * @brief 代回验证：把解代入方程，检查两边是否数值相等。
 *
 * 只在两边都化为数值时才判真；含未定符号的情形视为"未验证"，不算失败。
 */
bool verifySolution(const NodePtr &left, const NodePtr &right, std::uint32_t variable, const Ptr &value) {
    Ptr leftValue = substitute(Ptr::borrow(left), variable, value);
    Ptr rightValue = substitute(Ptr::borrow(right), variable, value);
    if (!leftValue || !rightValue) return false;
    if (!leftValue->isNumber() || !rightValue->isNumber()) return true;
    return leftValue->numerator() * rightValue->denominator() ==
           rightValue->numerator() * leftValue->denominator();
}

/** @brief 二次方程的精确解：x = (-b ± √Δ) / (2a)。 */
Memory::Vector<Solution> solveQuadratic(const UnivariatePoly &polynomial) {
    Memory::Vector<Solution> solutions;
    const Ptr a = polynomial.coefficient(2);
    const Ptr b = polynomial.coefficient(1);
    const Ptr c = polynomial.coefficient(0);
    Ptr discriminant = simplify(makeAdd(makeMul(b, b), makeNegate(makeMul(integer(4), makeMul(a, c)))));
    if (!discriminant) return solutions;

    // 判别式为完全平方的有理数：给出两个有理根
    if (discriminant->isNumber()) {
        const BigInt numerator = discriminant->numerator();
        const BigInt denominator = discriminant->denominator();
        if (numerator.isNegative()) return solutions; // 复根：不返回
        if (numerator.isPerfectSquare() && denominator.isPerfectSquare()) {
            const BigInt root = numerator.isqrt() / denominator.isqrt();
            if (root * root == numerator / denominator) {
                Ptr twoA = makeMul(integer(2), a);
                Solution first;
                first.value = simplify(makeDivide(makeAdd(makeNegate(b), makeNumber(root, BigInt(1))), twoA));
                first.verified = true;
                solutions.push_back(first);
                if (!root.isZero()) {
                    Solution second;
                    second.value =
                        simplify(makeDivide(makeAdd(makeNegate(b), makeNumber(root.negated(), BigInt(1))), twoA));
                    second.verified = true;
                    solutions.push_back(second);
                }
                return solutions;
            }
        }
    }

    // 无理根：保留 √Δ 形式
    Ptr root = makePow(discriminant, makeNumber(BigInt(1), BigInt(2)));
    Ptr twoA = makeMul(integer(2), a);
    Solution first;
    first.value = simplify(makeDivide(makeAdd(makeNegate(b), root), twoA));
    solutions.push_back(first);
    Solution second;
    second.value = simplify(makeDivide(makeAdd(makeNegate(b), makeNegate(root)), twoA));
    solutions.push_back(second);
    return solutions;
}

/** @brief 多项式求根：低次走公式，高次走有理根 + 降次。 */
Memory::Vector<Solution> solvePolynomial(const UnivariatePoly &polynomial, std::uint32_t variable) {
    Memory::Vector<Solution> solutions;
    if (polynomial.isZero()) return solutions; // 恒等式：任何值都是解，交由调用方处理

    switch (polynomial.degree()) {
    case 0:
        return solutions; // 非零常数 = 0：无解
    case 1: {
        Ptr a = polynomial.coefficient(1), b = polynomial.coefficient(0);
        Solution solution;
        solution.value = simplify(makeDivide(makeNegate(b), a));
        solution.verified = true;
        solutions.push_back(solution);
        return solutions;
    }
    case 2:
        return solveQuadratic(polynomial);
    default: break;
    }

    // 高次：尝试有理根，逐个剥离
    UnivariatePoly working = polynomial;
    for (int guard = 0; guard < 64 && working.degree() >= 1; ++guard) {
        // 常数项为 0 → x = 0
        if (working.coefficient(0)->isZero()) {
            Solution zero;
            zero.value = integer(0);
            zero.verified = true;
            solutions.push_back(zero);
            Memory::Vector<Ptr> shifted(working.coefficients.begin() + 1, working.coefficients.end());
            working.coefficients = shifted;
            working.trim();
            continue;
        }

        // 试小整数根
        long long foundRoot = 0;
        bool found = false;
        for (long long candidate = 1; candidate <= 512 && !found; ++candidate) {
            for (int sign : {1, -1}) {
                const long long value = candidate * sign;
                Ptr evaluated = substitute(fromPolynomial(working), variable, integer(value));
                if (evaluated && evaluated->isZero()) {
                    foundRoot = value;
                    found = true;
                    break;
                }
            }
        }
        if (!found) break;

        Solution solution;
        solution.value = integer(foundRoot);
        solution.verified = true;
        solutions.push_back(solution);

        // 综合除法降次
        UnivariatePoly divisor, quotient, remainder;
        divisor.variable = variable;
        divisor.coefficients.push_back(integer(-foundRoot));
        divisor.coefficients.push_back(integer(1));
        if (!polyDivide(working, divisor, &quotient, &remainder)) break;
        if (!remainder.isZero()) break;
        working = quotient;
    }

    // 剩余二次部分
    if (working.degree() == 2) {
        const auto tail = solveQuadratic(working);
        for (const Solution &solution : tail) solutions.push_back(solution);
    } else if (working.degree() == 1) {
        Ptr a = working.coefficient(1), b = working.coefficient(0);
        Solution solution;
        solution.value = simplify(makeDivide(makeNegate(b), a));
        solution.verified = true;
        solutions.push_back(solution);
    }
    return solutions;
}

/** @brief 去重并按大小排序（能比较数值时）。 */
void normalizeSolutions(Memory::Vector<Solution> *solutions) {
    Memory::Vector<Solution> unique;
    for (const Solution &solution : *solutions) {
        if (!solution.value) continue;
        bool duplicate = false;
        for (const Solution &existing : unique) {
            if (existing.value->equals(solution.value.get())) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) unique.push_back(solution);
    }
    std::stable_sort(unique.begin(), unique.end(), [](const Solution &left, const Solution &right) {
        if (left.value->isNumber() && right.value->isNumber()) {
            const BigInt l = left.value->numerator() * right.value->denominator();
            const BigInt r = right.value->numerator() * left.value->denominator();
            return l.compare(r) < 0;
        }
        return false;
    });
    *solutions = unique;
}

/* ------------------------------------------------ 超越方程的简单情形 */

/**
 * @brief 处理 e^(kx+b) = c、ln(kx+b) = c、sin(...) = c 这类可反解结构。
 */
bool solveTranscendental(const NodePtr &left, const NodePtr &right, std::uint32_t variable,
                         Memory::Vector<Solution> *solutions) {
    // 归一为 f(x) = 0 形式：left - right
    Ptr difference = simplify(makeAdd(Ptr::borrow(left), makeNegate(Ptr::borrow(right))));
    if (!difference) return false;

    // 形式一：加法里有一项是函数调用，其余为常数 → 移项后反解
    Memory::Vector<NodePtr> terms;
    if (difference->kind() == NodeKind::Add) {
        for (NodePtr child : difference->operands()) terms.push_back(child);
    } else {
        terms.push_back(difference.get());
    }

    NodePtr target = nullptr;
    Ptr rest = integer(0);
    for (NodePtr term : terms) {
        // 判断该项是否"含变量的函数调用"
        bool hasFunction = false;
        if (term->kind() == NodeKind::Function) {
            hasFunction = dependsOn(term, variable);
        } else if (term->kind() == NodeKind::Mul) {
            for (NodePtr factor : term->operands()) {
                if (factor->kind() == NodeKind::Function && dependsOn(factor, variable)) hasFunction = true;
            }
        }
        if (hasFunction && !target) {
            target = term;
        } else {
            rest = makeAdd(rest, Ptr::borrow(term));
        }
    }
    if (!target) return false;

    // 取出函数部分
    NodePtr functionNode = target;
    Ptr coefficient = integer(1);
    if (target->kind() == NodeKind::Mul) {
        for (NodePtr factor : target->operands()) {
            if (factor->kind() == NodeKind::Function && dependsOn(factor, variable)) {
                functionNode = factor;
            } else {
                coefficient = makeMul(coefficient, Ptr::borrow(factor));
            }
        }
    }

    // 目标：function(inner) = -rest / coefficient
    Ptr targetValue = simplify(makeDivide(makeNegate(rest), coefficient));
    if (!targetValue) return false;

    const char *name = functionNode->name();
    Ptr inner = Ptr::borrow(functionNode->operand(0));

    // 反函数映射
    Ptr innerValue;
    if (std::strcmp(name, "exp") == 0 ||
        (functionNode->kind() == NodeKind::Pow &&
         functionNode->operand(0)->kind() == NodeKind::Constant &&
         functionNode->operand(0)->constant() == ConstantId::E)) {
        // e^u = c → u = ln c
        if (targetValue->isNumber() && !targetValue->numerator().isNegative()) {
            Memory::Vector<NodePtr> arguments;
            arguments.push_back(targetValue.get());
            innerValue = Ptr(Node::makeFunction("log", arguments));
        }
    } else if (std::strcmp(name, "log") == 0) {
        // ln u = c → u = e^c
        innerValue = makePow(Ptr(Node::makeConstant(ConstantId::E)), targetValue);
    } else if (std::strcmp(name, "sin") == 0 || std::strcmp(name, "cos") == 0) {
        // 三角方程：主值 + 周期性说明，这里只给主值分支
        const char *inverse = std::strcmp(name, "sin") == 0 ? "asin" : "acos";
        Memory::Vector<NodePtr> arguments;
        arguments.push_back(targetValue.get());
        innerValue = Ptr(Node::makeFunction(inverse, arguments));
    }

    if (!innerValue) return false;

    // 最后解线性方程 inner = innerValue
    Memory::Vector<Solution> innerSolutions;
    if (!solve(inner.get(), innerValue.get(), variable, &innerSolutions)) return false;
    for (const Solution &solution : innerSolutions) solutions->push_back(solution);
    return true;
}

} // namespace

bool solve(const NodePtr &left, const NodePtr &right, std::uint32_t variable,
           Memory::Vector<Solution> *solutions) {
    solutions->clear();
    if (!left || !right) return false;

    // 归一为 f(x) = 0
    Ptr difference = simplify(makeAdd(Ptr::borrow(left), makeNegate(Ptr::borrow(right))));
    if (!difference) return false;

    // 恒等式：0 = 0
    if (difference->isZero()) return false;

    // 多项式路径
    Ptr expanded = expand(difference);
    UnivariatePoly polynomial;
    if (toPolynomial(expanded.get(), variable, &polynomial) && !polynomial.isZero()) {
        *solutions = solvePolynomial(polynomial, variable);
        normalizeSolutions(solutions);
        for (Solution &solution : *solutions) {
            if (!solution.verified) solution.verified = verifySolution(left, right, variable, solution.value);
        }
        return true;
    }

    // 换元路径：x^2 视作新变量（双二次等）
    {
        // 构造临时变量名并替换 x^2 → t
        Ptr t(Node::makeSymbol("_t"));
        const std::uint32_t tId = t->symbolId();
        Ptr replaced = substitute(expanded, variable, makePow(Ptr(Node::makeSymbolById(variable)), integer(2)));
        UnivariatePoly inT;
        if (toPolynomial(replaced.get(), tId, &inT) && inT.degree() >= 1 && inT.degree() <= 2) {
            Memory::Vector<Solution> tSolutions;
            const Ptr zero = integer(0);
            if (solve(replaced.get(), zero.get(), tId, &tSolutions) && !tSolutions.empty()) {
                bool allResolved = true;
                for (const Solution &tSolution : tSolutions) {
                    // x^2 = t → x = ±√t
                    if (tSolution.value->isNumber() && tSolution.value->isNegative()) continue; // 复数根跳过
                    Ptr root = makePow(tSolution.value, makeNumber(BigInt(1), BigInt(2)));
                    Memory::Vector<Solution> xSolutions;
                    const Ptr xSquared = makePow(Ptr(Node::makeSymbolById(variable)), integer(2));
                    solve(xSquared.get(), root.get(), variable, &xSolutions);
                    for (const Solution &xSolution : xSolutions) solutions->push_back(xSolution);
                }
                if (allResolved && !solutions->empty()) {
                    normalizeSolutions(solutions);
                    return true;
                }
            }
        }
    }

    // 超越方程路径
    if (solveTranscendental(left, right, variable, solutions)) {
        normalizeSolutions(solutions);
        return true;
    }
    return false;
}

bool solveInequality(const NodePtr &left, const NodePtr &right, const std::uint32_t variable,
                     const char *relation, std::string *description) {
    Ptr difference = simplify(makeAdd(Ptr::borrow(left), makeNegate(Ptr::borrow(right))));
    if (!difference) return false;

    UnivariatePoly polynomial;
    if (!toPolynomial(expand(difference).get(), variable, &polynomial)) return false;

    Ptr name(Node::makeSymbolById(variable));
    const std::string symbol = toText(name.get());

    if (polynomial.degree() == 1) {
        const Ptr a = polynomial.coefficient(1), b = polynomial.coefficient(0);
        if (!a->isNumber()) return false;
        // a x + b (relation) 0  →  x (relation') -b/a（a < 0 时方向翻转）
        Ptr threshold = simplify(makeDivide(makeNegate(b), a));
        bool flip = a->isNegative();
        std::string op = relation;
        if (flip) {
            if (op == "<") op = ">";
            else if (op == ">") op = "<";
            else if (op == "<=") op = ">=";
            else if (op == ">=") op = "<=";
        }
        *description = symbol + " " + op + " " + toText(threshold.get());
        return true;
    }

    if (polynomial.degree() == 2) {
        Memory::Vector<Solution> roots = solveQuadratic(polynomial);
        if (roots.size() != 2) return false;
        const Ptr a = polynomial.coefficient(2);
        const bool opening = !a->isNegative();
        const std::string r1 = toText(roots[0].value.get());
        const std::string r2 = toText(roots[1].value.get());
        const bool less = std::strcmp(relation, "<") == 0 || std::strcmp(relation, "<=") == 0;
        const bool inclusive = std::strcmp(relation, "<=") == 0 || std::strcmp(relation, ">=") == 0;
        const char *inner = inclusive ? "\\le" : "<";
        const char *outer = inclusive ? "\\ge" : ">";
        if (opening) {
            *description = less ? r1 + " " + inner + " " + symbol + " " + inner + " " + r2
                                : symbol + " " + outer + " " + r2 + "\\ \\text{or}\\ " + symbol + " " + inner + " " + r1;
        } else {
            *description = less ? symbol + " " + outer + " " + r2 + "\\ \\text{or}\\ " + symbol + " " + inner + " " + r1
                                : r1 + " " + inner + " " + symbol + " " + inner + " " + r2;
        }
        return true;
    }
    return false;
}

bool solveLinearSystem(const Memory::Vector<Memory::Vector<Ptr>> &equations,
                       const Memory::Vector<std::uint32_t> &variables, Memory::Vector<Ptr> *solution) {
    solution->clear();
    const std::size_t rows = equations.size();
    const std::size_t columns = variables.size();
    if (rows == 0 || columns == 0) return false;

    // 输入约定：每行是增广矩阵的一行，前 columns 个是系数，第 columns 个是右端项。
    // 系数必须是数值：含符号（如 `x*y = 6` 里 x 的系数是 y）说明方程组其实是非线性的，
    // 按符号做消元会给出 `x = 6/y` 这种"把未知数当参数"的假解，必须拒绝。
    Memory::Vector<Memory::Vector<Ptr>> matrix = equations;
    for (auto &row : matrix) {
        if (row.size() != columns + 1) return false;
        row.resize(columns + 1);
        for (Ptr &entry : row) {
            if (!entry) entry = integer(0);
            if (!entry->isNumber()) return false;
        }
    }

    // 高斯-约当消元
    std::size_t pivotRow = 0;
    Memory::Vector<std::size_t> pivotColumns;
    for (std::size_t column = 0; column < columns && pivotRow < rows; ++column) {
        // 选绝对值最大的主元（数值情形更稳；符号情形取首个非零）
        std::size_t selected = rows;
        for (std::size_t row = pivotRow; row < rows; ++row) {
            if (!matrix[row][column]->isZero()) {
                selected = row;
                break;
            }
        }
        if (selected == rows) continue;

        std::swap(matrix[pivotRow], matrix[selected]);
        const Ptr pivot = matrix[pivotRow][column];
        for (std::size_t c = column; c <= columns; ++c) {
            matrix[pivotRow][c] = simplify(makeDivide(matrix[pivotRow][c], pivot));
        }
        // 归一化后回代：把所有其它行的该列消成 0
        for (std::size_t row = 0; row < rows; ++row) {
            if (row == pivotRow) continue;
            const Ptr factor = matrix[row][column];
            if (factor->isZero()) continue;
            for (std::size_t c = column; c <= columns; ++c) {
                matrix[row][c] = simplify(makeAdd(matrix[row][c], makeNegate(makeMul(factor, matrix[pivotRow][c]))));
            }
        }
        pivotColumns.push_back(column);
        ++pivotRow;
    }

    // 矛盾行：0 = 非零
    for (std::size_t row = pivotRow; row < rows; ++row) {
        bool allZero = true;
        for (std::size_t c = 0; c < columns; ++c) {
            if (!matrix[row][c]->isZero()) {
                allZero = false;
                break;
            }
        }
        if (allZero && !matrix[row][columns]->isZero()) return false;
    }
    if (pivotColumns.size() < columns) return false; // 欠定

    solution->resize(columns);
    for (std::size_t i = 0; i < pivotColumns.size(); ++i) (*solution)[pivotColumns[i]] = matrix[i][columns];
    // 未覆盖的变量补 0（理论上不会发生，因为上面的欠定检查）
    for (std::size_t i = 0; i < columns; ++i) {
        if (!(*solution)[i]) (*solution)[i] = integer(0);
    }
    return true;
}

/* ------------------------------------------------ 多项式方程组（非线性） */

namespace {

/** @brief 分支上限：超过就停止展开，如实报告「解过多」。 */
constexpr std::size_t kMaxBranches = 64;

/** @brief 已确定的赋值（变量编号 → 值）。 */
using Assignments = Memory::Vector<std::pair<std::uint32_t, Ptr>>;

/** @brief 三角消元的递归状态。 */
struct TriangularState {
    Memory::Vector<Ptr> equations;         ///< 当前仍生效的 f = 0 列表
    Memory::Vector<std::uint32_t> pending; ///< 尚未确定的变量
    Memory::Vector<std::uint32_t> allVariables; ///< 全部变量（Groebner 变量表必须覆盖）
    Assignments fixed;                     ///< 已确定的赋值
    Memory::Vector<std::uint32_t> free;    ///< 判定为自由变量（消不出方程）的变量
    Memory::Vector<SystemSolution> *out = nullptr; ///< 输出
    bool overflow = false;                 ///< 分支数超限
};

/** @brief 多项式是否只含变量表中列出的符号，且系数为有理数。 */
bool onlyInVariables(const Ptr &polynomial, const Memory::Vector<std::uint32_t> &variables) {
    const Memory::Vector<std::uint32_t> free = freeSymbols(polynomial.get());
    for (std::uint32_t symbol : free) {
        bool allowed = false;
        for (std::uint32_t candidate : variables) {
            if (symbol == candidate) {
                allowed = true;
                break;
            }
        }
        if (!allowed) return false;
    }
    // Groebner 基只在有理数域上保证终止，符号系数直接拒绝
    MultivariatePoly check;
    Memory::Vector<std::uint32_t> ordered = variables;
    return toMultivariate(polynomial.get(), ordered, &check);
}

/**
 * @brief 过滤方程列表：丢掉恒等式，遇到非零常数则报告矛盾。
 * @return 出现矛盾（非零常数 = 0）返回 true
 */
bool reduceEquations(Memory::Vector<Ptr> *equations) {
    Memory::Vector<Ptr> kept;
    kept.reserve(equations->size());
    for (const Ptr &equation : *equations) {
        if (!equation) continue;
        Ptr simplified = simplify(equation);
        if (!simplified) continue;
        if (simplified->isZero()) continue;      // 0 = 0：恒等式，丢弃
        if (simplified->isNumber()) return true; // 非零常数 = 0：矛盾
        kept.push_back(simplified);
    }
    *equations = kept;
    return false;
}

/** @brief 解单变量多项式，返回去重后的所有实根。 */
Memory::Vector<Ptr> rootsOf(const Ptr &polynomial, std::uint32_t variable) {
    Memory::Vector<Ptr> roots;
    UnivariatePoly poly;
    if (!toPolynomial(polynomial.get(), variable, &poly)) return roots;
    if (poly.isZero()) return roots;

    // 一次方程允许符号系数（`y = 1 - x` 的通解形式），高次不行：
    // `x^2 + y^2 - 1 = 0` 关于 y 的根是 ±√(1-x²)，那是分支函数而不是单个解，
    // 贸然返回会把参数解搞成伪解。高次时只接受纯数值系数。
    if (poly.degree() > 1) {
        for (const Ptr &coefficient : poly.coefficients) {
            if (!coefficient) return roots;
            if (!coefficient->isNumber()) return roots;
        }
    }

    for (const Solution &solution : solvePolynomial(poly, variable)) {
        if (!solution.value) continue;
        const std::string key = toText(solution.value.get());
        bool duplicate = false;
        for (const Ptr &existing : roots) {
            if (toText(existing.get()) == key) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) roots.push_back(solution.value);
    }
    return roots;
}

/** @brief 从待解变量表中移除一个变量（不论它排在哪个位置）。 */
void dropVariable(Memory::Vector<std::uint32_t> *pending, std::uint32_t variable) {
    for (std::size_t i = 0; i < pending->size(); ++i) {
        if ((*pending)[i] == variable) {
            pending->erase(pending->begin() + static_cast<std::ptrdiff_t>(i));
            return;
        }
    }
}

void solveRecursive(TriangularState &state);

/**
 * @brief 解出 `pending.back()` 这一个变量，对每个根代回后继续递归。
 *
 * 求该变量的方程有两个来源，按代价从低到高尝试：
 *  1. 方程组里本来就有只含该变量的方程 → 直接求根；
 *  2. 否则算一次 Groebner 基。该变量排在变量表**末尾**（最主元），
 *     所以基里必然出现一个只含它的多项式——这就是消元的结果。
 *
 * 每得到一个根就代回其它方程：代回后出现非零常数说明这个根不满足别的方程，
 * 直接剪枝。所以「非根」的中间值不会扩散成分支。
 */
void solveRecursive(TriangularState &state) {
    if (state.overflow) return;

    if (state.pending.empty()) {
        // 变量都解完了。剩下的方程有两种可能：
        //  1. equations 为空       —— 干净地解完整；
        //  2. equations 非空       —— 这些方程在已求出的根上退化成 0（冗余，例如
        //                            x+y=1 与 2x+2y=2），不代表矛盾，因为矛盾早在
        //                            代回时被识别成非零常数并剪枝了。
        // 两种都要产出一个分支，否则整个方程组会被误报成「无解」。
        SystemSolution solution;
        solution.complete = state.free.empty();
        solution.assignments = state.fixed;
        if (state.out->size() < kMaxBranches) {
            state.out->push_back(solution);
        } else {
            state.overflow = true;
        }
        return;
    }

    if (state.equations.empty()) {
        // 方程用完了但变量还有剩：参数曲面 / 欠定，如实报告
        SystemSolution solution;
        solution.complete = false;
        solution.assignments = state.fixed;
        if (state.out->size() < kMaxBranches) {
            state.out->push_back(solution);
        } else {
            state.overflow = true;
        }
        return;
    }

    // 取 **首位** 变量来解。
    // `pending` 的首位对应 `variables` 中最早的那个变量，也就是 Lex 消元中最先被
    // 消净、最后才被解出的那个——此时其它变量都已代入成常数，它的系数必然不含
    // 其它符号，单变量求根才能成功。
    // 反过来取末位（`back()`）会先解最主元的变量，此刻其它变量还没定，
    // 方程形如 `y = 1 - x`，系数含 x，被 `toPolynomial` 判为非常数而求根失败，
    // 于是 `x+y=1` 这种再简单不过的方程反而报「无解」。
    const std::uint32_t variable = state.pending.front();
    const Memory::Vector<std::uint32_t> savedPending = state.pending;
    const Memory::Vector<Ptr> savedEquations = state.equations;

    // ---- 1. 方程组里是否已有只含该变量的方程
    Memory::Vector<Ptr> candidates;
    for (const Ptr &equation : savedEquations) {
        if (onlyInVariables(equation, {variable})) {
            candidates.push_back(equation);
            break;
        }
    }

    // ---- 2. 没有就用 Groebner 基消元
    if (candidates.empty()) {
        // 变量序的语义（实测 `compareMonomial` + Lex）：变量表**末位**最主元，
        // 消元时被消净的是**首位**那个变量。
        // 例：Lex `[x, y]` 对 {x^2+y^2-25, xy-12} 给出
        //     {y^3/12 + x - 25y/12, y^4 - 25y^2 + 144}  ← 只剩 y，x 被消净。
        // 本轮的变量要成为"剩下的那个"，所以必须排在**末位**；
        // 排首位会被消掉，基里就找不到只含它的多项式。
        //
        // 变量表必须覆盖**所有**变量（state.allVariables），不能只用 pending：
        // 方程里可能还带着已跳过、尚未代回的变量（如欠定时留下的 x），
        // 变量表缺了它，toMultivariate 会因"未声明符号"把生成元全部丢掉，
        // 基就成了空的，分支被误判成矛盾。
        Memory::Vector<std::uint32_t> ordered;
        ordered.reserve(state.allVariables.size());
        for (std::uint32_t other : state.allVariables) {
            if (other != variable) ordered.push_back(other);
        }
        ordered.push_back(variable);

        // Lex 更直接给出消元结果；次数不高时计算量也可接受。
        Memory::Vector<Ptr> basis = groebnerBasis(savedEquations, ordered, MonomialOrder::Lex);
        if (basis.empty()) {
            // Lex 失败（极端情况）时退回 Grevlex 再试一次
            basis = groebnerBasis(savedEquations, ordered, MonomialOrder::Grevlex);
        }
        if (basis.empty()) return; // 生成理想的方程组矛盾，本分支无解
        for (const Ptr &element : basis) {
            if (onlyInVariables(element, {variable})) {
                candidates.push_back(element);
                break;
            }
        }
    }

    // ---- 2.5 消元不成（方程数少于未知量，基里没有只含该变量的元素）：
    // 直接从原方程里解出该变量，留下含参表达式。
    // 这一步是欠定方程组能给出参数解的关键：
    // `x+y=1` 关于 x 求解得 `x = 1-y`，y 自然成为自由变量；
    // 若在这里放弃，整组方程就会被误判成「无解」。
    if (candidates.empty()) {
        for (const Ptr &equation : savedEquations) {
            UnivariatePoly poly;
            if (!toPolynomial(equation.get(), variable, &poly)) continue;
            if (poly.degree() != 1) continue;
            // 首项系数为 0 时这个方程其实不含该变量，跳过
            const Ptr leading = poly.coefficient(1);
            if (!leading || leading->isZero()) continue;
            candidates.push_back(equation);
            break;
        }
    }

    // ---- 3. 求根
    Memory::Vector<Ptr> roots;
    for (const Ptr &candidate : candidates) {
        roots = rootsOf(candidate, variable);
        if (!roots.empty()) break;
    }

    if (roots.empty()) {
        if (candidates.empty()) {
            // 消不出单变量方程：该变量是自由变量，先解更靠前的变量。
            // 必须登记到 free 里，否则分支会被当成"完整解"输出，
            // 格式化阶段就只能打印 `x = ?`。
            dropVariable(&state.pending, variable);
            state.free.push_back(variable);
            solveRecursive(state);
            state.free.pop_back();
            state.pending = savedPending;
        }
        // 有候选多项式但无实根（如 y^2+1=0 只有复根）：本分支无解
        return;
    }

    // ---- 4. 逐根代回，递归解剩下的变量
    for (const Ptr &root : roots) {
        if (state.overflow) return;

        Memory::Vector<Ptr> next;
        bool contradiction = false;
        for (const Ptr &equation : savedEquations) {
            // 先展开再化简：代回符号根（如 x = 1-y）会得到 `y-(y-1)-1` 这类
            // 结构上还没合并的式子，只 simplify 不够，会残留一个恒等式，
            // 后续把它当真实约束去求根，分支就再也走不下去。
            Ptr substituted = simplify(expand(substitute(equation, variable, root)));
            if (!substituted) {
                contradiction = true;
                break;
            }
            if (substituted->isZero()) continue;
            if (substituted->isNumber()) {
                contradiction = true; // 这个根不满足别的方程
                break;
            }
            next.push_back(substituted);
        }
        if (contradiction) continue;

        state.equations = next;
        state.pending = savedPending;
        dropVariable(&state.pending, variable);
        state.fixed.push_back({variable, root});
        solveRecursive(state);
        state.fixed.pop_back();
    }

    // 恢复状态，供上层（或下一个根）使用
    state.equations = savedEquations;
    state.pending = savedPending;
}

} // namespace

bool solveSystemNonlinear(const Memory::Vector<Ptr> &polynomials, const Memory::Vector<std::uint32_t> &variables,
                          Memory::Vector<SystemSolution> *solutions) {
    solutions->clear();
    if (polynomials.empty() || variables.empty()) return false;

    // 规范化输入：展开、确认是有理系数多项式
    Memory::Vector<Ptr> equations;
    equations.reserve(polynomials.size());
    for (const Ptr &polynomial : polynomials) {
        if (!polynomial) continue;
        Ptr expanded = expand(polynomial);
        if (!expanded) return false;
        if (!onlyInVariables(expanded, variables)) return false;
        equations.push_back(expanded);
    }

    if (reduceEquations(&equations)) return true; // 矛盾：无解，正常返回空表
    if (equations.empty()) {
        // 全是恒等式：任意取值都是解
        SystemSolution solution;
        solution.complete = false;
        solutions->push_back(solution);
        return true;
    }

    TriangularState state;
    state.equations = equations;
    state.pending = variables;
    state.allVariables = variables;
    state.out = solutions;
    solveRecursive(state);

    if (state.overflow) {
        solutions->clear();
        return false;
    }
    return true;
}

} // namespace Cas::Detail

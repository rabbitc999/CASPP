// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#include "cas/detail/solver.h"
#include "cas/detail/calculus.h"
#include "cas/detail/polynomial.h"
#include "cas/detail/printer.h"

#include <algorithm>
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
    Memory::Vector<Memory::Vector<Ptr>> matrix = equations;
    for (auto &row : matrix) {
        if (row.size() != columns + 1) return false;
        row.resize(columns + 1);
        for (Ptr &entry : row) {
            if (!entry) entry = integer(0);
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

} // namespace Cas::Detail

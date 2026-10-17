/**
 * @file matrix.cpp
 * @brief 符号矩阵实现。
 */

#include "cas/detail/matrix.h"
#include "cas/detail/polynomial.h"
#include "cas/detail/simplify.h"
#include "cas/detail/solver.h"

#include <cstring>

namespace Cas::Detail {
namespace {

/* ------------------------------------------------------------ 小工具 */

Ptr zeroValue() { return makeNumber(BigInt(0), BigInt(1)); }
Ptr oneValue() { return makeNumber(BigInt(1), BigInt(1)); }

/** @brief 判断节点是否为数值零。 */
bool isZeroValue(const Ptr &node) {
    if (!node || !node->isNumber()) return false;
    return node->numerator().isZero();
}

/** @brief 判断节点是否为数值一。 */
bool isOneValue(const Ptr &node) {
    if (!node || !node->isNumber()) return false;
    return node->denominator().isOne() && node->numerator().isOne();
}

/** @brief 交换两行。 */
void swapRows(SymbolicMatrix *matrix, std::size_t a, std::size_t b) {
    if (a == b) return;
    for (std::size_t j = 0; j < matrix->columns(); ++j) {
        Ptr left = matrix->at(a, j);
        Ptr right = matrix->at(b, j);
        matrix->set(a, j, right);
        matrix->set(b, j, left);
    }
}

/**
 * @brief 分数消元（Gauss-Jordan），就地化为 RREF。
 *
 * 用有理数除法，元素都是符号时也能工作（结果为有理式）。
 */
bool fractionReducedEchelon(SymbolicMatrix *matrix, Memory::Vector<std::size_t> *pivots) {
    const std::size_t rowCount = matrix->rows();
    const std::size_t columnCount = matrix->columns();
    std::size_t pivotRow = 0;

    for (std::size_t column = 0; column < columnCount && pivotRow < rowCount; ++column) {
        // 找首个非零元作为主元
        std::size_t chosen = rowCount;
        for (std::size_t row = pivotRow; row < rowCount; ++row) {
            Ptr value = simplify(matrix->at(row, column));
            if (!value || isZeroValue(value)) continue;
            chosen = row;
            break;
        }
        if (chosen == rowCount) continue; // 该列无主元
        swapRows(matrix, pivotRow, chosen);

        // 归一化主元行
        Ptr pivot = simplify(matrix->at(pivotRow, column));
        if (!pivot) return false;
        if (!isOneValue(pivot)) {
            for (std::size_t j = column; j < columnCount; ++j) {
                Ptr cell = simplify(matrix->at(pivotRow, j));
                Ptr scaled = simplify(makeDivide(cell, pivot));
                matrix->set(pivotRow, j, scaled ? scaled : cell);
            }
        }

        // 消去其余各行在该列的元素
        for (std::size_t row = 0; row < rowCount; ++row) {
            if (row == pivotRow) continue;
            Ptr factor = simplify(matrix->at(row, column));
            if (!factor || isZeroValue(factor)) continue;
            for (std::size_t j = column; j < columnCount; ++j) {
                Ptr cell = simplify(matrix->at(row, j));
                Ptr pivotCell = simplify(matrix->at(pivotRow, j));
                Ptr product = simplify(makeMul(factor, pivotCell));
                Ptr updated = simplify(makeAdd(cell, makeNegate(product)));
                matrix->set(row, j, updated ? updated : cell);
            }
        }

        if (pivots) pivots->push_back(column);
        ++pivotRow;
    }
    return true;
}

} // namespace

/* -------------------------------------------------------- 类成员实现 */

SymbolicMatrix::SymbolicMatrix(std::size_t rows, std::size_t columns)
    : m_rows(rows), m_columns(columns) {
    m_entries.reserve(rows * columns);
    for (std::size_t i = 0; i < rows * columns; ++i) m_entries.push_back(zeroValue());
}

SymbolicMatrix SymbolicMatrix::fromRows(std::size_t rows, std::size_t columns,
                                        const Memory::Vector<Ptr> &entries) {
    SymbolicMatrix matrix(rows, columns);
    if (entries.size() != rows * columns) return SymbolicMatrix();
    for (std::size_t i = 0; i < entries.size(); ++i) {
        matrix.m_entries[i] = entries[i] ? simplify(entries[i]) : zeroValue();
    }
    return matrix;
}

const Ptr &SymbolicMatrix::at(std::size_t row, std::size_t column) const {
    static const Ptr empty;
    if (row >= m_rows || column >= m_columns) return empty;
    return m_entries[row * m_columns + column];
}

void SymbolicMatrix::set(std::size_t row, std::size_t column, const Ptr &value) {
    if (row >= m_rows || column >= m_columns) return;
    m_entries[row * m_columns + column] = value;
}

bool SymbolicMatrix::isZero() const {
    for (const Ptr &entry : m_entries) {
        Ptr simplified = simplify(entry);
        if (!simplified || !isZeroValue(simplified)) return false;
    }
    return true;
}

bool SymbolicMatrix::isIdentity() const {
    if (!isSquare()) return false;
    for (std::size_t row = 0; row < m_rows; ++row) {
        for (std::size_t column = 0; column < m_columns; ++column) {
            Ptr value = simplify(at(row, column));
            if (row == column) {
                if (!value || !isOneValue(value)) return false;
            } else {
                if (!value || !isZeroValue(value)) return false;
            }
        }
    }
    return true;
}

bool SymbolicMatrix::equals(const SymbolicMatrix &other) const {
    if (m_rows != other.m_rows || m_columns != other.m_columns) return false;
    for (std::size_t i = 0; i < m_entries.size(); ++i) {
        Ptr left = simplify(m_entries[i]);
        Ptr right = simplify(other.m_entries[i]);
        if (!left || !right) return false;
        if (!left->equals(right.get())) return false;
    }
    return true;
}

/* ------------------------------------------------------------ 基本运算 */

SymbolicMatrix transpose(const SymbolicMatrix &matrix) {
    if (matrix.isEmpty()) return SymbolicMatrix();
    SymbolicMatrix result(matrix.columns(), matrix.rows());
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        for (std::size_t column = 0; column < matrix.columns(); ++column) {
            result.set(column, row, matrix.at(row, column));
        }
    }
    return result;
}

SymbolicMatrix addMatrices(const SymbolicMatrix &left, const SymbolicMatrix &right) {
    if (left.rows() != right.rows() || left.columns() != right.columns()) return SymbolicMatrix();
    SymbolicMatrix result(left.rows(), left.columns());
    for (std::size_t row = 0; row < left.rows(); ++row) {
        for (std::size_t column = 0; column < left.columns(); ++column) {
            Ptr sum = simplify(makeAdd(left.at(row, column), right.at(row, column)));
            result.set(row, column, sum ? sum : zeroValue());
        }
    }
    return result;
}

SymbolicMatrix subtractMatrices(const SymbolicMatrix &left, const SymbolicMatrix &right) {
    if (left.rows() != right.rows() || left.columns() != right.columns()) return SymbolicMatrix();
    SymbolicMatrix result(left.rows(), left.columns());
    for (std::size_t row = 0; row < left.rows(); ++row) {
        for (std::size_t column = 0; column < left.columns(); ++column) {
            Ptr difference = simplify(makeAdd(left.at(row, column),
                                              makeNegate(right.at(row, column))));
            result.set(row, column, difference ? difference : zeroValue());
        }
    }
    return result;
}

SymbolicMatrix multiplyMatrices(const SymbolicMatrix &left, const SymbolicMatrix &right) {
    if (left.columns() != right.rows()) return SymbolicMatrix();
    SymbolicMatrix result(left.rows(), right.columns());
    for (std::size_t row = 0; row < left.rows(); ++row) {
        for (std::size_t column = 0; column < right.columns(); ++column) {
            Ptr accumulator = zeroValue();
            for (std::size_t k = 0; k < left.columns(); ++k) {
                Ptr product = makeMul(left.at(row, k), right.at(k, column));
                accumulator = simplify(makeAdd(accumulator, product));
            }
            result.set(row, column, accumulator ? accumulator : zeroValue());
        }
    }
    return result;
}

SymbolicMatrix scaleMatrix(const Ptr &scalar, const SymbolicMatrix &matrix) {
    if (!scalar) return SymbolicMatrix();
    SymbolicMatrix result(matrix.rows(), matrix.columns());
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        for (std::size_t column = 0; column < matrix.columns(); ++column) {
            Ptr scaled = simplify(makeMul(scalar, matrix.at(row, column)));
            result.set(row, column, scaled ? scaled : zeroValue());
        }
    }
    return result;
}

SymbolicMatrix identityMatrix(std::size_t size) {
    SymbolicMatrix result(size, size);
    for (std::size_t i = 0; i < size; ++i) result.set(i, i, oneValue());
    return result;
}

/* -------------------------------------------------------------- 行列式 */

/**
 * @brief Bareiss 分数消元求行列式。
 *
 * 每一步用上一轮的主元做精确整除：
 *   a[i][j] = (pivot * a[i][j] - a[i][k] * a[k][j]) / previous
 * 中间结果始终在同一个域内，不会像 Gauss 消元那样分数膨胀。
 * 参考 sympy `poly.matrices.dense.ddm_irref_den` 的主循环。
 */
Ptr determinant(const SymbolicMatrix &matrix) {
    if (!matrix.isSquare()) return Ptr();
    const std::size_t size = matrix.rows();
    if (size == 0) return oneValue();

    SymbolicMatrix work = matrix;
    Ptr previous = oneValue();
    int sign = 1;

    for (std::size_t k = 0; k + 1 < size; ++k) {
        // 选主元：优先取非零元，必要时交换行
        std::size_t pivotRow = size;
        for (std::size_t row = k; row < size; ++row) {
            Ptr candidate = simplify(work.at(row, k));
            if (candidate && !isZeroValue(candidate)) {
                pivotRow = row;
                break;
            }
        }
        if (pivotRow == size) return zeroValue(); // 整列为零
        if (pivotRow != k) {
            swapRows(&work, k, pivotRow);
            sign = -sign;
        }

        Ptr pivot = simplify(work.at(k, k));
        if (!pivot) return Ptr();

        for (std::size_t i = k + 1; i < size; ++i) {
            for (std::size_t j = k + 1; j < size; ++j) {
                Ptr leftTerm = simplify(makeMul(pivot, work.at(i, j)));
                Ptr rightTerm = simplify(makeMul(work.at(i, k), work.at(k, j)));
                Ptr difference = simplify(makeAdd(leftTerm, makeNegate(rightTerm)));
                Ptr updated = difference;
                if (!isOneValue(previous)) {
                    updated = simplify(makeDivide(difference, previous));
                }
                if (!updated) return Ptr();
                work.set(i, j, updated);
            }
            work.set(i, k, zeroValue());
        }
        previous = pivot;
    }

    Ptr result = simplify(work.at(size - 1, size - 1));
    if (!result) return Ptr();
    return sign < 0 ? simplify(makeNegate(result)) : result;
}

Ptr trace(const SymbolicMatrix &matrix) {
    if (!matrix.isSquare()) return Ptr();
    Ptr accumulator = zeroValue();
    for (std::size_t i = 0; i < matrix.rows(); ++i) {
        accumulator = simplify(makeAdd(accumulator, matrix.at(i, i)));
    }
    return accumulator;
}

bool reducedRowEchelon(const SymbolicMatrix &matrix, SymbolicMatrix *out,
                       Memory::Vector<std::size_t> *pivots) {
    if (!out || matrix.isEmpty()) return false;
    SymbolicMatrix work = matrix;
    Memory::Vector<std::size_t> columns;
    if (!fractionReducedEchelon(&work, &columns)) return false;
    *out = work;
    if (pivots) *pivots = columns;
    return true;
}

bool rank(const SymbolicMatrix &matrix, std::size_t *out) {
    if (!out || matrix.isEmpty()) return false;
    SymbolicMatrix reduced;
    Memory::Vector<std::size_t> pivots;
    if (!reducedRowEchelon(matrix, &reduced, &pivots)) return false;
    *out = pivots.size();
    return true;
}

/* ---------------------------------------------------------------- 逆 */

Ptr inverse(const SymbolicMatrix &matrix, SymbolicMatrix *out) {
    if (!out || !matrix.isSquare()) return Ptr();
    const std::size_t size = matrix.rows();

    // 增广 [A | I]，消元后左半就是 A^-1
    SymbolicMatrix augmented(size, size * 2);
    for (std::size_t row = 0; row < size; ++row) {
        for (std::size_t column = 0; column < size; ++column) {
            augmented.set(row, column, matrix.at(row, column));
        }
        augmented.set(row, size + row, oneValue());
    }

    Memory::Vector<std::size_t> pivots;
    if (!fractionReducedEchelon(&augmented, &pivots)) return Ptr();
    // 只看左半部分的主元：右半（单位阵）几乎总能提供主元，
    // 若统统计数会把奇异矩阵误判为可逆。
    std::size_t leftPivots = 0;
    for (std::size_t index = 0; index < pivots.size(); ++index) {
        if (pivots[index] < size) ++leftPivots;
    }
    if (leftPivots != size) return zeroValue(); // 奇异矩阵

    SymbolicMatrix result(size, size);
    for (std::size_t row = 0; row < size; ++row) {
        for (std::size_t column = 0; column < size; ++column) {
            result.set(row, column, augmented.at(row, size + column));
        }
    }
    *out = result;
    return determinant(matrix); // 返回行列式，便于调用方区分"奇异"与"失败"
}

/* ------------------------------------------------------- 特征多项式 */

/**
 * @brief Faddeev–LeVerrier 递推求特征多项式。
 *
 *   M_0 = 0, c_0 = 1
 *   M_k = A·M_{k-1} + c_{k-1}·I
 *   c_k = -tr(A·M_k) / k
 * 则 p(λ) = λ^n + c_1 λ^{n-1} + … + c_n。
 */
Ptr characteristicPolynomial(const SymbolicMatrix &matrix, std::uint32_t variable) {
    if (!matrix.isSquare()) return Ptr();
    const std::size_t size = matrix.rows();
    if (size == 0) return oneValue();

    Memory::Vector<Ptr> coefficients; // c_0 … c_n
    coefficients.push_back(oneValue());

    SymbolicMatrix accumulated(size, size); // M_0 = 0
    for (std::size_t k = 1; k <= size; ++k) {
        // M_k = A·M_{k-1} + c_{k-1}·I
        SymbolicMatrix product = multiplyMatrices(matrix, accumulated);
        for (std::size_t i = 0; i < size; ++i) {
            Ptr diagonal = simplify(makeAdd(product.at(i, i), coefficients[k - 1]));
            product.set(i, i, diagonal ? diagonal : zeroValue());
        }
        accumulated = product;

        // c_k = -tr(A·M_k) / k
        SymbolicMatrix productWithA = multiplyMatrices(matrix, accumulated);
        Ptr total = simplify(trace(productWithA));
        Ptr divided = simplify(makeDivide(simplify(makeNegate(total)),
                                          makeNumber(BigInt(static_cast<long long>(k)), BigInt(1))));
        if (!divided) return Ptr();
        coefficients.push_back(divided);
    }

    // 组装 λ^n + c_1 λ^(n-1) + … + c_n
    Ptr polynomial = zeroValue();
    Ptr variableNode = Ptr(Node::makeSymbolById(variable));
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        const std::size_t power = size - i;
        Ptr term = coefficients[i];
        if (power > 0) term = makeMul(term, makePow(variableNode, makeNumber(BigInt(static_cast<long long>(power)), BigInt(1))));
        polynomial = simplify(makeAdd(polynomial, term));
    }
    return polynomial;
}

/* ----------------------------------------------------------- 特征值 */

bool eigenvalues(const SymbolicMatrix &matrix, Memory::Vector<Ptr> *out) {
    if (!out || !matrix.isSquare()) return false;
    out->clear();

    // 特征多项式用独立变量 λ，避免与矩阵元素里的符号冲突
    const std::uint32_t lambda = internSymbol("_lambda");
    Ptr characteristic = characteristicPolynomial(matrix, lambda);
    if (!characteristic) return false;

    Memory::Vector<Solution> solutions;
    Ptr zero = zeroValue(); // 必须具名持有，否则 get() 得到悬垂指针
    if (!solve(characteristic.get(), zero.get(), lambda, &solutions)) return false;
    for (const Solution &solution : solutions) {
        if (solution.value) out->push_back(solution.value);
    }
    return !out->empty();
}

/* -------------------------------------------------------- 解线性方程 */

bool solveLinearSystemMatrix(const SymbolicMatrix &coefficients, const SymbolicMatrix &rightHandSide,
                             Memory::Vector<Ptr> *out, bool *consistent) {
    if (!out) return false;
    if (consistent) *consistent = true;
    if (coefficients.isEmpty() || rightHandSide.rows() != coefficients.rows()) return false;
    if (rightHandSide.columns() != 1 && rightHandSide.rows() != 1) return false;

    const std::size_t rowCount = coefficients.rows();
    const std::size_t columnCount = coefficients.columns();

    // 增广矩阵
    SymbolicMatrix augmented(rowCount, columnCount + 1);
    for (std::size_t row = 0; row < rowCount; ++row) {
        for (std::size_t column = 0; column < columnCount; ++column) {
            augmented.set(row, column, coefficients.at(row, column));
        }
        const Ptr &rhs = rightHandSide.columns() == 1 ? rightHandSide.at(row, 0) : rightHandSide.at(0, row);
        augmented.set(row, columnCount, rhs);
    }

    SymbolicMatrix reduced;
    Memory::Vector<std::size_t> pivots;
    if (!reducedRowEchelon(augmented, &reduced, &pivots)) return false;

    // 矛盾行：系数全零但右端非零
    for (std::size_t row = 0; row < rowCount; ++row) {
        bool allZero = true;
        for (std::size_t column = 0; column < columnCount; ++column) {
            Ptr value = simplify(reduced.at(row, column));
            if (value && !isZeroValue(value)) {
                allZero = false;
                break;
            }
        }
        if (allZero) {
            Ptr rhs = simplify(reduced.at(row, columnCount));
            if (rhs && !isZeroValue(rhs)) {
                if (consistent) *consistent = false;
                return false;
            }
        }
    }

    // 主元列数少于未知数个数 → 无穷多解，自由变量取 0
    out->clear();
    out->reserve(columnCount);
    for (std::size_t column = 0; column < columnCount; ++column) out->push_back(zeroValue());

    for (std::size_t index = 0; index < pivots.size(); ++index) {
        const std::size_t pivotColumn = pivots[index];
        if (pivotColumn >= columnCount) continue; // 增广列不属于未知数
        Ptr value = simplify(reduced.at(index, columnCount));
        (*out)[pivotColumn] = value ? value : zeroValue();
    }
    return true;
}

} // namespace Cas::Detail

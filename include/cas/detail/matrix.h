#ifndef CAS_DETAIL_MATRIX_H
#define CAS_DETAIL_MATRIX_H

/**
 * @file matrix.h
 * @brief 符号矩阵与线性代数。
 *
 * 元素是通用的符号表达式（`Ptr`），因此矩阵可以同时承载数值与含参矩阵。
 * 行列式走 Bareiss 分数消元：每一步都做精确整除，中间结果不会像普通
 * Gauss 消元那样产生分数膨胀——这一点参考了 sympy 的 `ddm_irref_den`。
 *
 * 约定：所有返回 `Ptr` 的函数在失败时返回空 `Ptr`；所有接收 `Ptr` 的函数在
 * 参数为空时立即返回空/失败，不做断言。
 */

#include "cas/detail/node.h"

#include <cstddef>
#include <cstdint>
#include <string>

namespace Cas::Detail {

/** @brief 二维符号矩阵（行主序稠密存储）。 */
class SymbolicMatrix {
public:
    SymbolicMatrix() = default;
    /** @brief 构造 rows × columns 的全零矩阵。 */
    SymbolicMatrix(std::size_t rows, std::size_t columns);

    /** @brief 从行主序元素列表构造；元素个数必须等于 rows*columns。 */
    static SymbolicMatrix fromRows(std::size_t rows, std::size_t columns,
                                   const Memory::Vector<Ptr> &entries);

    std::size_t rows() const { return m_rows; }
    std::size_t columns() const { return m_columns; }
    bool isEmpty() const { return m_rows == 0 || m_columns == 0; }
    bool isSquare() const { return m_rows == m_columns && m_rows > 0; }

    /** @brief 取元素；越界返回空 Ptr。 */
    const Ptr &at(std::size_t row, std::size_t column) const;
    /** @brief 写元素；越界忽略。 */
    void set(std::size_t row, std::size_t column, const Ptr &value);

    /** @brief 是否全部元素都是零。 */
    bool isZero() const;
    /** @brief 是否为单位阵。 */
    bool isIdentity() const;
    /** @brief 结构相等（元素逐项 `equals`）。 */
    bool equals(const SymbolicMatrix &other) const;

private:
    std::size_t m_rows = 0;
    std::size_t m_columns = 0;
    Memory::Vector<Ptr> m_entries;
};

/* ------------------------------------------------------------ 基本运算 */

SymbolicMatrix transpose(const SymbolicMatrix &matrix);
SymbolicMatrix addMatrices(const SymbolicMatrix &left, const SymbolicMatrix &right);
SymbolicMatrix subtractMatrices(const SymbolicMatrix &left, const SymbolicMatrix &right);
SymbolicMatrix multiplyMatrices(const SymbolicMatrix &left, const SymbolicMatrix &right);
/** @brief 标量乘：每个元素乘以 scalar。 */
SymbolicMatrix scaleMatrix(const Ptr &scalar, const SymbolicMatrix &matrix);
/** @brief 单位阵。 */
SymbolicMatrix identityMatrix(std::size_t size);

/* ------------------------------------------------------------ 不变量 */

/**
 * @brief 行列式（Bareiss 分数消元）。
 *
 * 非方阵、或元素无法化成数/符号表达式时返回空。
 */
Ptr determinant(const SymbolicMatrix &matrix);

/** @brief 迹：主对角线元素之和。 */
Ptr trace(const SymbolicMatrix &matrix);

/**
 * @brief 化简后的行阶梯形（RREF）与主元列号。
 *
 * 内部用分数消元（允许有理数中间结果），返回矩阵的每个元素都已 `simplify`。
 */
bool reducedRowEchelon(const SymbolicMatrix &matrix, SymbolicMatrix *out,
                       Memory::Vector<std::size_t> *pivots);

/** @brief 秩：主元个数。 */
bool rank(const SymbolicMatrix &matrix, std::size_t *out);

/**
 * @brief 逆矩阵。
 *
 * 用伴随矩阵除以行列式，行列式为 0 时返回空。
 */
Ptr inverse(const SymbolicMatrix &matrix, SymbolicMatrix *out);

/**
 * @brief 特征多项式 det(λI − A)。
 *
 * 通过 Faddeev–LeVerrier 递推构造，返回以 `variable` 为变量的表达式。
 */
Ptr characteristicPolynomial(const SymbolicMatrix &matrix, std::uint32_t variable);

/**
 * @brief 特征值（实根与有理根）。
 *
 * 对特征多项式做有理根分解，逐因子求根；无法求出的因子会被跳过。
 */
bool eigenvalues(const SymbolicMatrix &matrix, Memory::Vector<Ptr> *out);

/**
 * @brief 解 A·x = b。
 *
 * 返回唯一解向量；无解或无穷多解时返回 false，并通过 `consistent` 区分。
 */
bool solveLinearSystemMatrix(const SymbolicMatrix &coefficients, const SymbolicMatrix &rightHandSide,
                             Memory::Vector<Ptr> *out, bool *consistent);

} // namespace Cas::Detail

#endif // CAS_DETAIL_MATRIX_H

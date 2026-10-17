#ifndef CAS_DETAIL_MULTIVARIATE_H
#define CAS_DETAIL_MULTIVARIATE_H

/**
 * @file multivariate.h
 * @brief 多元多项式与 Groebner 基。
 *
 * 表示法：`Monomial` 是"指数向量 + 有理系数"，`MultivariatePoly` 是单项式集合。
 * 系数限定为**有理数**——这是 Groebner 基算法能保证终止的前提（域上才能做除法），
 * 符号系数会退化成多项式环，需要额外的消元序，超出本工程范围。
 *
 * Groebner 基采用 Buchberger 算法 + 乘积判据（Buchberger 第一判据）跳过冗余 S-多项式。
 * 变量序默认按符号编号升序对应 x1 < x2 < ...（分次逆字典序 grevlex，与 Singular 默认一致）。
 */

#include "cas/detail/simplify.h"

namespace Cas::Detail {

/** @brief 单项式：指数向量（与 `variables` 对齐）+ 有理系数。 */
struct Monomial {
    Memory::Vector<long long> exponents; ///< 各变量的指数
    BigInt numerator;                    ///< 系数分子
    BigInt denominator;                  ///< 系数分母（恒正）
};

/** @brief 多元多项式：变量的有序表 + 单项式列表。 */
struct MultivariatePoly {
    Memory::Vector<std::uint32_t> variables; ///< 变量符号编号，按字典序
    Memory::Vector<Monomial> terms;          ///< 单项式

    /** @brief 是否为零多项式。 */
    bool isZero() const { return terms.empty(); }
};

/** @brief 把表达式转成多元多项式（系数必须是有理数）。 */
bool toMultivariate(const NodePtr &expression, const Memory::Vector<std::uint32_t> &variables,
                    MultivariatePoly *out);

/** @brief 多元多项式转回表达式。 */
Ptr fromMultivariate(const MultivariatePoly &polynomial);

/** @brief 多元多项式展开（等价于 `expand` 的多变量版本）。 */
Ptr expandMultivariate(const Ptr &expression, std::size_t termLimit = 100000);

/* --------------------------------------------------------- Groebner 基 */

/** @brief 单项式序：分次逆字典（grevlex），与多数代数系统的默认序一致。 */
enum class MonomialOrder {
    Lex,     ///< 纯字典序，便于消元
    Grevlex, ///< 分次逆字典序，计算最快
};

/**
 * @brief 计算 Groebner 基。
 * @param generators 生成元表达式
 * @param variables  变量表（顺序即消元优先级）
 * @param order      单项式序
 * @return 约化 Groebner 基（首一化，按序降序）
 */
Memory::Vector<Ptr> groebnerBasis(const Memory::Vector<Ptr> &generators,
                                  const Memory::Vector<std::uint32_t> &variables,
                                  MonomialOrder order = MonomialOrder::Grevlex);

/** @brief 多元多项式 GCD（仅在两元、且为互素单变量情形时保证精确）。 */
Ptr multivariateGcd(const Ptr &left, const Ptr &right, const Memory::Vector<std::uint32_t> &variables);

} // namespace Cas::Detail

#endif // CAS_DETAIL_MULTIVARIATE_H

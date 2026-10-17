// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#ifndef CAS_DETAIL_BIGINT_H
#define CAS_DETAIL_BIGINT_H

/**
 * @file bigInt.h
 * @brief 任意精度整数。
 *
 * 表示：绝对值以 2^32 为基的小端 limb 数组，符号单独存；|v| <= 2^61 的整数走内联存储，
 * 不触碰内存池——绝大多数符号表达式里的整数都在这个范围。
 *
 * 乘法：小块用教科书算法，中块用 Karatsuba，大块走三段分解 + NTT。
 * 除法：Knuth 算法 D（16 位半字长除），商向下取整、余数符号跟随除数，与 Python 语义一致。
 */

#include "config.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace Cas::Detail {

class BigInt {
public:
    BigInt() = default;
    BigInt(long long value);
    BigInt(int value) : BigInt(static_cast<long long>(value)) {}
    explicit BigInt(const char *text, int base = 10);
    explicit BigInt(const std::string &text, int base = 10);

    bool isZero() const;
    bool isNegative() const { return m_negative; }
    bool isOne() const;
    bool isOdd() const;
    int sign() const;

    std::string toString() const;
    std::string toStringHex() const;
    double toDouble() const;
    std::size_t bitLength() const;
    /** @brief 幅度占用的 32 位 limb 个数。 */
    std::size_t limbCount() const;
    /** @brief 取第 index 个 limb（小端，越界返回 0）。 */
    std::uint32_t limbAt(std::size_t index) const;

    BigInt abs() const;
    BigInt negated() const;
    BigInt operator-() const { return negated(); }

    BigInt &operator+=(const BigInt &other);
    BigInt &operator-=(const BigInt &other);
    BigInt &operator*=(const BigInt &other);
    BigInt &operator/=(const BigInt &other);
    BigInt &operator%=(const BigInt &other);

    BigInt operator+(const BigInt &other) const;
    BigInt operator-(const BigInt &other) const;
    BigInt operator*(const BigInt &other) const;
    BigInt operator/(const BigInt &other) const;
    BigInt operator%(const BigInt &other) const;

    BigInt operator<<(int bits) const;
    BigInt operator>>(int bits) const;

    int compare(const BigInt &other) const;
    bool operator==(const BigInt &other) const;
    bool operator!=(const BigInt &other) const { return !(*this == other); }
    bool operator<(const BigInt &other) const { return compare(other) < 0; }
    bool operator<=(const BigInt &other) const { return compare(other) <= 0; }
    bool operator>(const BigInt &other) const { return compare(other) > 0; }
    bool operator>=(const BigInt &other) const { return compare(other) >= 0; }

    /** @brief 是否可用 long long 精确取回。 */
    bool fitsLongLong() const { return m_small; }
    /** @brief 取低 64 位（按补码语义解释）；大数时截断，不抛异常。 */
    long long toLongLong() const;

    /** @brief 向下整除与取余：q = floor(a/b)，r = a - q*b（r 与 b 同号）。 */
    static void divMod(const BigInt &a, const BigInt &b, BigInt *quotient, BigInt *remainder);
    /** @brief 非负模。 */
    BigInt mod(const BigInt &modulus) const;
    /** @brief 模逆；不互素返回 0。 */
    BigInt modInverse(const BigInt &modulus) const;
    static BigInt power(const BigInt &base, unsigned long long exponent);
    static BigInt gcd(const BigInt &a, const BigInt &b);
    BigInt isqrt() const;
    bool isPerfectSquare() const;
    BigInt nthRoot(unsigned n) const;
    bool isPerfectPower(unsigned n) const;

private:
    using Limbs = std::vector<std::uint32_t>;

    void normalize();
    void setFromLongLong(long long value);
    Limbs toLimbs() const;
    static void limbsToBig(const Limbs &limbs, bool negative, BigInt *out);
    static int compareAbs(const BigInt &a, const BigInt &b);
    static void addAbs(const Limbs &a, const Limbs &b, Limbs *out);
    static void subAbs(const Limbs &a, const Limbs &b, Limbs *out); // 要求 |a| >= |b|
    static Limbs mulAbs(const Limbs &a, const Limbs &b);
    static void divModAbs(const Limbs &a, const Limbs &b, Limbs *quotient, Limbs *remainder);

    bool m_negative = false;
    bool m_small = true;
    long long m_smallValue = 0;
    Limbs m_limbs;
};

} // namespace Cas::Detail

#endif // CAS_DETAIL_BIGINT_H

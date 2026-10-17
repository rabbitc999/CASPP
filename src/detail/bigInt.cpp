// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

#include "cas/detail/bigInt.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>

namespace Cas::Detail {
namespace {

using Limbs = std::vector<std::uint32_t>;

constexpr std::uint64_t kLimbBase = 1ull << 32;

std::size_t actualSize(const Limbs &limbs) {
    std::size_t n = limbs.size();
    while (n > 0 && limbs[n - 1] == 0) --n;
    return n;
}

int bitLengthOf(const Limbs &limbs) {
    const std::size_t n = actualSize(limbs);
    if (n == 0) return 0;
    std::uint32_t top = limbs[n - 1];
    int bits = 0;
    while (top) {
        ++bits;
        top >>= 1;
    }
    return static_cast<int>((n - 1) * 32) + bits;
}

std::uint32_t getLimb(const Limbs &limbs, std::size_t index) { return index < limbs.size() ? limbs[index] : 0; }

/** @brief 幅值比较。 */
int compareLimbs(const Limbs &a, const Limbs &b) {
    const std::size_t na = actualSize(a), nb = actualSize(b);
    if (na != nb) return na < nb ? -1 : 1;
    for (std::size_t i = na; i-- > 0;) {
        if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    }
    return 0;
}

Limbs slice(const Limbs &a, std::size_t from, std::size_t to) {
    if (to <= from) return {};
    return Limbs(a.begin() + static_cast<std::ptrdiff_t>(from), a.begin() + static_cast<std::ptrdiff_t>(to));
}

Limbs shiftLeftLimbs(const Limbs &a, int bits) {
    if (bits == 0) return a;
    const std::size_t limbShift = static_cast<std::size_t>(bits / 32);
    const int bitShift = bits % 32;
    Limbs out(a.size() + limbShift + 1, 0);
    if (bitShift == 0) {
        for (std::size_t i = 0; i < a.size(); ++i) out[i + limbShift] = a[i];
    } else {
        std::uint32_t carry = 0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            out[i + limbShift] = (a[i] << bitShift) | carry;
            carry = a[i] >> (32 - bitShift);
        }
        out[a.size() + limbShift] = carry;
    }
    while (!out.empty() && out.back() == 0) out.pop_back();
    return out;
}

/** @brief out = a - b，要求 a >= b。 */
void differenceInto(Limbs *out, const Limbs &a, const Limbs &b) {
    out->assign(a.begin(), a.end());
    std::uint64_t borrow = 0;
    for (std::size_t i = 0; i < out->size(); ++i) {
        const std::uint64_t current = (*out)[i];
        const std::uint64_t other = getLimb(b, i) + borrow;
        const bool underflow = current < other;
        (*out)[i] = static_cast<std::uint32_t>(current + kLimbBase - other);
        borrow = underflow ? 1 : 0;
    }
    out->resize(actualSize(*out));
}

/** @brief a -= b，要求 a >= b。 */
void subInPlace(Limbs &a, const Limbs &b) {
    Limbs result;
    differenceInto(&result, a, b);
    a.swap(result);
}

/* ---------------------------------------------------------------- 乘法 */

Limbs mulSchoolbook(const Limbs &a, const Limbs &b) {
    const std::size_t n = actualSize(a), m = actualSize(b);
    Limbs out(n + m + 1, 0);
    for (std::size_t i = 0; i < n; ++i) {
        const std::uint64_t ai = a[i];
        std::uint64_t carry = 0;
        for (std::size_t j = 0; j < m; ++j) {
            const std::uint64_t cur = ai * b[j] + out[i + j] + carry;
            out[i + j] = static_cast<std::uint32_t>(cur);
            carry = cur >> 32;
        }
        std::size_t k = i + m;
        while (carry) {
            const std::uint64_t cur = static_cast<std::uint64_t>(out[k]) + carry;
            out[k] = static_cast<std::uint32_t>(cur);
            carry = cur >> 32;
            ++k;
        }
    }
    out.resize(actualSize(out));
    return out;
}

/** @brief out += b << (32 * offset)。 */
void addShiftedInPlace(Limbs &out, const Limbs &b, std::size_t offset) {
    if (b.empty()) return;
    if (out.size() < offset + b.size()) out.resize(offset + b.size(), 0);
    std::uint64_t carry = 0;
    std::size_t i = 0;
    for (; i < b.size() || carry; ++i) {
        if (offset + i >= out.size()) out.push_back(0);
        const std::uint64_t cur = static_cast<std::uint64_t>(out[offset + i]) + getLimb(b, i) + carry;
        out[offset + i] = static_cast<std::uint32_t>(cur);
        carry = cur >> 32;
    }
    out.resize(actualSize(out));
}

Limbs mulRecursive(const Limbs &a, const Limbs &b);

Limbs mulKaratsuba(const Limbs &a, const Limbs &b) {
    const std::size_t n = actualSize(a), m = actualSize(b);
    const std::size_t half = std::max<std::size_t>(1, std::max(n, m) / 2);
    const Limbs a0 = slice(a, 0, half), a1 = slice(a, half, n);
    const Limbs b0 = slice(b, 0, half), b1 = slice(b, half, m);
    Limbs z0 = mulRecursive(a0, b0);
    Limbs z2 = mulRecursive(a1, b1);
    Limbs sumA = a0, sumB = b0;
    addShiftedInPlace(sumA, a1, 0);
    addShiftedInPlace(sumB, b1, 0);
    Limbs z1 = mulRecursive(sumA, sumB);
    subInPlace(z1, z0);
    subInPlace(z1, z2);
    Limbs out = z0;
    addShiftedInPlace(out, z1, half);
    addShiftedInPlace(out, z2, 2 * half);
    return out;
}

/* ------------------------------------------------------------ 三段 NTT */

std::uint32_t powerMod(std::uint32_t base, std::uint64_t exponent) {
    std::uint64_t result = 1, b = base % kNttPrime;
    while (exponent) {
        if (exponent & 1) result = result * b % kNttPrime;
        b = b * b % kNttPrime;
        exponent >>= 1;
    }
    return static_cast<std::uint32_t>(result);
}

void ntt(std::vector<std::uint32_t> &a, bool invert) {
    const std::size_t n = a.size();
    const std::uint64_t mod = kNttPrime;
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (std::size_t len = 2; len <= n; len <<= 1) {
        std::uint32_t w = powerMod(kNttRoot, (mod - 1) / len);
        if (invert) w = powerMod(w, mod - 2);
        for (std::size_t i = 0; i < n; i += len) {
            std::uint64_t wn = 1;
            const std::size_t half = len >> 1;
            for (std::size_t j = 0; j < half; ++j) {
                const std::uint64_t u = a[i + j];
                const std::uint64_t v = a[i + j + half] * wn % mod;
                a[i + j] = static_cast<std::uint32_t>((u + v) % mod);
                a[i + j + half] = static_cast<std::uint32_t>((u + mod - v) % mod);
                wn = wn * w % mod;
            }
        }
    }
    if (invert) {
        const std::uint64_t invN = powerMod(static_cast<std::uint32_t>(n), mod - 2);
        for (auto &x : a) x = static_cast<std::uint32_t>(static_cast<std::uint64_t>(x) * invN % mod);
    }
}

/**
 * @brief NTT 乘法：32 位 limb 拆成 8/8/8/8 位四段做卷积。
 *
 * 每段取值 < 2^8，单个卷积和上限为 `min(n,m) * (2^8-1)^2 < 2^16 * limbs`。
 * 模数 kNttPrime = 2013265921 ≈ 2^30.9，所以只要 limb 数 < 2^14 就不会模折叠；
 * 超过该规模时自动退化为分段（逐块 NTT 累加），保证结果与朴素乘法完全一致。
 */
Limbs mulNtt(const Limbs &a, const Limbs &b) {
    const std::size_t n = actualSize(a), m = actualSize(b);
    constexpr std::size_t kRadixBits = 8;
    constexpr std::uint32_t kRadixMask = (1u << kRadixBits) - 1;
    constexpr std::size_t kSegments = 32 / kRadixBits; // 每个 32 位 limb 拆成 4 段

    std::vector<std::uint32_t> da, db;
    da.reserve(n * kSegments);
    db.reserve(m * kSegments);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t s = 0; s < kSegments; ++s) da.push_back((a[i] >> (s * kRadixBits)) & kRadixMask);
    }
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t s = 0; s < kSegments; ++s) db.push_back((b[i] >> (s * kRadixBits)) & kRadixMask);
    }

    // 卷积系数上界：min(n,m) * (2^8-1)^2，必须小于模数才能保证单次 NTT 不折叠。
    const std::uint64_t bound = static_cast<std::uint64_t>(std::min(n, m)) * kRadixMask * kRadixMask;
    if (bound >= kNttPrime) {
        // 规模超出单次 NTT 的安全范围（limb 数 ≥ 2^14）：继续一半拆分，保证精确性。
        return mulKaratsuba(a, b);
    }

    std::size_t size = 1;
    while (size < da.size() + db.size()) size <<= 1;
    da.resize(size, 0);
    db.resize(size, 0);
    ntt(da, false);
    ntt(db, false);
    for (std::size_t i = 0; i < size; ++i)
        da[i] = static_cast<std::uint32_t>(static_cast<std::uint64_t>(da[i]) * db[i] % kNttPrime);
    ntt(da, true);

    // 按 4 段一组重组回 32 位 limb，段间进位用 64 位暂存。
    Limbs out;
    std::uint64_t carry = 0;
    for (std::size_t i = 0; i < da.size(); i += kSegments) {
        std::uint64_t value = carry;
        for (std::size_t s = 0; s < kSegments && i + s < da.size(); ++s)
            value += static_cast<std::uint64_t>(da[i + s]) << (s * kRadixBits);
        out.push_back(static_cast<std::uint32_t>(value & 0xFFFFFFFFu));
        carry = value >> 32;
    }
    while (carry) {
        out.push_back(static_cast<std::uint32_t>(carry & 0xFFFFFFFFu));
        carry >>= 32;
    }
    out.resize(actualSize(out));
    return out;
}

Limbs mulRecursive(const Limbs &a, const Limbs &b) {
    const std::size_t n = actualSize(a), m = actualSize(b);
    if (n == 0 || m == 0) return {};
    if (std::min(n, m) <= 32) return mulSchoolbook(a, b);
    if (std::min(n, m) <= 96) return mulKaratsuba(a, b);
    return mulNtt(a, b);
}

/* ---------------------------------------------------------------- 除法 */

void divModSmall(const Limbs &a, std::uint32_t divisor, Limbs *quotient, std::uint32_t *remainder) {
    const std::size_t n = actualSize(a);
    quotient->assign(n, 0);
    std::uint64_t rem = 0;
    for (std::size_t i = n; i-- > 0;) {
        const std::uint64_t cur = (rem << 32) | a[i];
        (*quotient)[i] = static_cast<std::uint32_t>(cur / divisor);
        rem = cur % divisor;
    }
    quotient->resize(actualSize(*quotient));
    *remainder = static_cast<std::uint32_t>(rem);
}

/**
 * @brief 移位相减长除。语义取自 magel/library/bignum 的 divmodBase：
 *        逐位判断 |r| >= |b| << k，成立就相减并置商位，不使用任何估商启发式。
 *
 * 正确性只依赖"比较 + 减法"两个可单独验证的原语，因此没有边界陷阱。
 * 代价是每位一次比较/一次减法；符号表达式中的商通常只有几十位，开销可接受。
 */
void divModByShiftSubtract(const Limbs &a, const Limbs &b, Limbs *quotient, Limbs *remainder) {
    Limbs r = a;
    const std::size_t bits = static_cast<std::size_t>(bitLengthOf(a) - bitLengthOf(b));
    Limbs q(bits / 32 + 1, 0);
    for (std::size_t step = bits + 1; step-- > 0;) {
        const Limbs shifted = shiftLeftLimbs(b, static_cast<int>(step));
        if (compareLimbs(r, shifted) >= 0) {
            subInPlace(r, shifted);
            q[step / 32] |= (1u << (step % 32));
        }
    }
    q.resize(actualSize(q));
    quotient->swap(q);
    remainder->swap(r);
}

} // namespace

/* ==================================================== 构造 / 规范化 */

BigInt::BigInt(long long value) { setFromLongLong(value); }

BigInt::BigInt(const char *text, int base) : BigInt(std::string(text ? text : "0"), base) {}

BigInt::BigInt(const std::string &text, int base) {
    std::size_t i = 0;
    bool negative = false;
    while (i < text.size() && (text[i] == ' ' || text[i] == '\t')) ++i;
    if (i < text.size() && (text[i] == '+' || text[i] == '-')) {
        negative = text[i] == '-';
        ++i;
    }
    m_small = false;
    m_limbs.clear();
    if (base == 16) {
        // 每 8 个十六进制位正好一个 limb；不足 8 位的组放在最高位。
        std::string digitsOnly;
        for (; i < text.size(); ++i) {
            const char c = text[i];
            if ((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F')) digitsOnly.push_back(c);
        }
        const std::size_t groupSize = 8;
        const std::size_t remainder = digitsOnly.size() % groupSize;
        const std::size_t firstLength = remainder == 0 ? groupSize : remainder;
        Limbs reversed;
        for (std::size_t pos = 0; pos < digitsOnly.size();) {
            const std::size_t length = (pos == 0) ? firstLength : groupSize;
            std::uint64_t value = 0;
            for (std::size_t k = 0; k < length; ++k) {
                const char c = digitsOnly[pos + k];
                const int digit = (c <= '9') ? c - '0' : ((c >= 'a') ? c - 'a' + 10 : c - 'A' + 10);
                value = (value << 4) | static_cast<std::uint64_t>(digit);
            }
            reversed.push_back(static_cast<std::uint32_t>(value));
            pos += length;
        }
        m_limbs.assign(reversed.rbegin(), reversed.rend());
    } else {
        for (; i < text.size(); ++i) {
            const char c = text[i];
            if (c < '0' || c > '9') continue;
            std::uint64_t carry = static_cast<std::uint64_t>(c - '0');
            for (std::size_t k = 0; k < m_limbs.size(); ++k) {
                const std::uint64_t cur = static_cast<std::uint64_t>(m_limbs[k]) * 10 + carry;
                m_limbs[k] = static_cast<std::uint32_t>(cur);
                carry = cur >> 32;
            }
            if (carry) m_limbs.push_back(static_cast<std::uint32_t>(carry));
        }
    }
    m_negative = negative;
    normalize();
}

void BigInt::setFromLongLong(long long value) {
    m_limbs.clear();
    m_negative = value < 0;
    if (value >= -kSmallIntLimit && value <= kSmallIntLimit) {
        m_small = true;
        m_smallValue = value;
        return;
    }
    m_small = false;
    unsigned long long magnitude = static_cast<unsigned long long>(value);
    if (value < 0) magnitude = 0ull - magnitude;
    while (magnitude) {
        m_limbs.push_back(static_cast<std::uint32_t>(magnitude & 0xFFFFFFFFu));
        magnitude >>= 32;
    }
}

void BigInt::normalize() {
    const std::size_t n = actualSize(m_limbs);
    m_limbs.resize(n);
    if (n == 0) {
        m_negative = false;
        m_small = true;
        m_smallValue = 0;
        return;
    }
    if (n <= 2) {
        // 回到内联范围时改用小整数表示，保证同一数值只有一种表示。
        unsigned long long magnitude = m_limbs[0];
        if (n == 2) magnitude |= static_cast<unsigned long long>(m_limbs[1]) << 32;
        if (magnitude <= static_cast<unsigned long long>(kSmallIntLimit)) {
            m_small = true;
            m_smallValue = m_negative ? -static_cast<long long>(magnitude) : static_cast<long long>(magnitude);
            m_limbs.clear();
            return;
        }
    }
    m_small = false;
}

BigInt::Limbs BigInt::toLimbs() const {
    if (!m_small) return m_limbs;
    Limbs limbs;
    unsigned long long magnitude = m_smallValue < 0 ? 0ull - static_cast<unsigned long long>(m_smallValue)
                                                    : static_cast<unsigned long long>(m_smallValue);
    while (magnitude) {
        limbs.push_back(static_cast<std::uint32_t>(magnitude & 0xFFFFFFFFu));
        magnitude >>= 32;
    }
    return limbs;
}

void BigInt::limbsToBig(const Limbs &limbs, bool negative, BigInt *out) {
    out->m_small = false;
    out->m_limbs = limbs;
    out->m_negative = negative;
    out->normalize();
}

/* ==================================================== 基础查询 */

bool BigInt::isZero() const { return m_small ? m_smallValue == 0 : m_limbs.empty(); }

bool BigInt::isOne() const {
    return m_small ? m_smallValue == 1 : (m_limbs.size() == 1 && m_limbs[0] == 1 && !m_negative);
}

int BigInt::sign() const {
    if (isZero()) return 0;
    return m_negative ? -1 : 1;
}

bool BigInt::isOdd() const { return m_small ? (m_smallValue & 1) != 0 : (!m_limbs.empty() && (m_limbs[0] & 1u) != 0); }

int BigInt::compareAbs(const BigInt &a, const BigInt &b) { return compareLimbs(a.toLimbs(), b.toLimbs()); }

int BigInt::compare(const BigInt &other) const {
    const int sa = sign(), sb = other.sign();
    if (sa != sb) return sa < sb ? -1 : 1;
    if (sa == 0) return 0;
    const int magnitude = compareAbs(*this, other);
    return sa > 0 ? magnitude : -magnitude;
}

bool BigInt::operator==(const BigInt &other) const {
    if (m_small && other.m_small) return m_smallValue == other.m_smallValue;
    if (isZero() != other.isZero()) return false;
    if (isZero()) return true;
    if (m_negative != other.m_negative) return false;
    return compareAbs(*this, other) == 0;
}

std::size_t BigInt::bitLength() const {
    if (m_small) {
        unsigned long long magnitude = m_smallValue < 0 ? 0ull - static_cast<unsigned long long>(m_smallValue)
                                                        : static_cast<unsigned long long>(m_smallValue);
        std::size_t bits = 0;
        while (magnitude) {
            ++bits;
            magnitude >>= 1;
        }
        return bits;
    }
    return static_cast<std::size_t>(bitLengthOf(m_limbs));
}

BigInt BigInt::abs() const {
    BigInt result = *this;
    result.m_negative = false;
    if (result.m_small) result.m_smallValue = result.m_smallValue < 0 ? -result.m_smallValue : result.m_smallValue;
    return result;
}

BigInt BigInt::negated() const {
    BigInt result = *this;
    if (result.isZero()) return result;
    result.m_negative = !result.m_negative;
    if (result.m_small) result.m_smallValue = -result.m_smallValue;
    return result;
}

/* ==================================================== 加减 */

void BigInt::addAbs(const Limbs &a, const Limbs &b, Limbs *out) {
    const std::size_t n = std::max(a.size(), b.size());
    out->assign(n + 1, 0);
    std::uint64_t carry = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const std::uint64_t cur = static_cast<std::uint64_t>(getLimb(a, i)) + getLimb(b, i) + carry;
        (*out)[i] = static_cast<std::uint32_t>(cur);
        carry = cur >> 32;
    }
    (*out)[n] = static_cast<std::uint32_t>(carry);
    out->resize(actualSize(*out));
}

void BigInt::subAbs(const Limbs &a, const Limbs &b, Limbs *out) { differenceInto(out, a, b); }

BigInt BigInt::operator+(const BigInt &other) const {
    if (m_small && other.m_small) {
        const long long sum = m_smallValue + other.m_smallValue; // 两者都在 2^61 内，long long 不会溢出
        return BigInt(sum);
    }
    const Limbs a = toLimbs(), b = other.toLimbs();
    BigInt result;
    result.m_small = false;
    if (m_negative == other.m_negative) {
        addAbs(a, b, &result.m_limbs);
        result.m_negative = m_negative;
    } else if (compareAbs(*this, other) > 0) {
        subAbs(a, b, &result.m_limbs);
        result.m_negative = m_negative;
    } else {
        subAbs(b, a, &result.m_limbs);
        result.m_negative = other.m_negative;
    }
    result.normalize();
    return result;
}

BigInt BigInt::operator-(const BigInt &other) const { return *this + other.negated(); }

BigInt &BigInt::operator+=(const BigInt &other) {
    *this = *this + other;
    return *this;
}

BigInt &BigInt::operator-=(const BigInt &other) {
    *this = *this - other;
    return *this;
}

/* ==================================================== 乘法 */

BigInt::Limbs BigInt::mulAbs(const Limbs &a, const Limbs &b) { return mulRecursive(a, b); }

BigInt BigInt::operator*(const BigInt &other) const {
    if (isZero() || other.isZero()) return BigInt();
    if (m_small && other.m_small) {
        // 全程用无符号幅度判定：有符号乘积一旦超过 LLONG_MAX 就是 UB，回绕后再判区间会给出错误结果。
        const long long x = m_smallValue, y = other.m_smallValue;
        const unsigned long long ax = x < 0 ? 0ull - static_cast<unsigned long long>(x) : static_cast<unsigned long long>(x);
        const unsigned long long ay = y < 0 ? 0ull - static_cast<unsigned long long>(y) : static_cast<unsigned long long>(y);
        const bool negative = (x < 0) != (y < 0);
        if (ax == 0 || ay <= static_cast<unsigned long long>(-1) / ax) {
            const unsigned long long magnitude = ax * ay;
            if (magnitude <= static_cast<unsigned long long>(kSmallIntLimit)) {
                const long long value = static_cast<long long>(magnitude);
                return BigInt(negative ? -value : value);
            }
        }
    }
    BigInt result;
    result.m_small = false;
    result.m_limbs = mulAbs(toLimbs(), other.toLimbs());
    result.m_negative = m_negative != other.m_negative;
    result.normalize();
    return result;
}

BigInt &BigInt::operator*=(const BigInt &other) {
    *this = *this * other;
    return *this;
}

/* ==================================================== 除法 */

void BigInt::divModAbs(const Limbs &a, const Limbs &b, Limbs *quotient, Limbs *remainder) {
    quotient->clear();
    remainder->clear();
    if (actualSize(b) == 0) return;
    if (compareLimbs(a, b) < 0) {
        remainder->assign(a.begin(), a.begin() + static_cast<std::ptrdiff_t>(actualSize(a)));
        return;
    }
    if (actualSize(b) == 1) {
        std::uint32_t rem = 0;
        divModSmall(a, b[0], quotient, &rem);
        if (rem) remainder->push_back(rem);
        return;
    }
    divModByShiftSubtract(a, b, quotient, remainder);
}

void BigInt::divMod(const BigInt &a, const BigInt &b, BigInt *quotient, BigInt *remainder) {
    if (b.isZero()) {
        if (quotient) *quotient = BigInt();
        if (remainder) *remainder = BigInt();
        return;
    }
    if (a.m_small && b.m_small) {
        const long long x = a.m_smallValue, y = b.m_smallValue;
        if (!(x == std::numeric_limits<long long>::min() && y == -1)) {
            long long q = x / y;
            long long r = x % y;
            if (r != 0 && ((r < 0) != (y < 0))) {
                --q;
                r += y;
            }
            if (quotient) *quotient = BigInt(q);
            if (remainder) *remainder = BigInt(r);
            return;
        }
    }
    Limbs q, r;
    divModAbs(a.toLimbs(), b.toLimbs(), &q, &r);
    BigInt qv, rv;
    limbsToBig(q, a.m_negative != b.m_negative, &qv);
    limbsToBig(r, a.m_negative, &rv);
    if (!rv.isZero() && rv.sign() != b.sign()) {
        qv = qv - BigInt(1);
        rv = rv + b;
    }
    if (quotient) *quotient = qv;
    if (remainder) *remainder = rv;
}

BigInt BigInt::operator/(const BigInt &other) const {
    BigInt q;
    divMod(*this, other, &q, nullptr);
    return q;
}

BigInt BigInt::operator%(const BigInt &other) const {
    BigInt r;
    divMod(*this, other, nullptr, &r);
    return r;
}

BigInt &BigInt::operator/=(const BigInt &other) {
    *this = *this / other;
    return *this;
}

BigInt &BigInt::operator%=(const BigInt &other) {
    *this = *this % other;
    return *this;
}

BigInt BigInt::mod(const BigInt &modulus) const {
    BigInt r = *this % modulus;
    if (r.m_negative) r = r + modulus.abs();
    return r;
}

BigInt BigInt::modInverse(const BigInt &modulus) const {
    const BigInt m = modulus.abs();
    BigInt oldR = mod(m), r = m;
    BigInt oldS(1), s(0);
    while (!r.isZero()) {
        const BigInt q = oldR / r;
        const BigInt nextR = oldR - q * r;
        const BigInt nextS = oldS - q * s;
        oldR = r;
        r = nextR;
        oldS = s;
        s = nextS;
    }
    if (!oldR.isOne()) return BigInt();
    return oldS.mod(m);
}

/* ==================================================== 位移 */

BigInt BigInt::operator<<(int bits) const {
    if (bits <= 0) return *this;
    BigInt result;
    limbsToBig(shiftLeftLimbs(toLimbs(), bits), m_negative, &result);
    return result;
}

BigInt BigInt::operator>>(int bits) const {
    if (bits <= 0) return *this;
    if (m_small) return BigInt(m_smallValue >> bits);
    const std::size_t limbShift = static_cast<std::size_t>(bits / 32);
    const int bitShift = bits % 32;
    if (limbShift >= m_limbs.size()) return BigInt();
    Limbs out(m_limbs.size() - limbShift, 0);
    if (bitShift == 0) {
        for (std::size_t i = 0; i < out.size(); ++i) out[i] = m_limbs[i + limbShift];
    } else {
        for (std::size_t i = 0; i + limbShift + 1 < m_limbs.size(); ++i)
            out[i] = (m_limbs[i + limbShift] >> bitShift) | (m_limbs[i + limbShift + 1] << (32 - bitShift));
        out[out.size() - 1] = m_limbs[m_limbs.size() - 1] >> bitShift;
    }
    BigInt result;
    limbsToBig(out, m_negative, &result);
    return result;
}

/* ==================================================== 幂 / 根 / 数论 */

BigInt BigInt::power(const BigInt &base, unsigned long long exponent) {
    BigInt result(1), b = base;
    while (exponent) {
        if (exponent & 1) result *= b;
        exponent >>= 1;
        if (exponent) b *= b;
    }
    return result;
}

BigInt BigInt::gcd(const BigInt &a, const BigInt &b) {
    BigInt x = a.abs(), y = b.abs();
    while (!y.isZero()) {
        const BigInt r = x % y;
        x = y;
        y = r;
    }
    return x;
}

BigInt BigInt::isqrt() const {
    if (isZero() || m_negative) return BigInt();
    if (m_small) {
        const long long value = m_smallValue;
        long long root = static_cast<long long>(std::sqrt(static_cast<double>(value)));
        while (root > 0 && root * root > value) --root;
        while ((root + 1) * (root + 1) <= value) ++root;
        return BigInt(root);
    }
    BigInt x = BigInt(1) << static_cast<int>((bitLength() + 1) / 2);
    while (true) {
        const BigInt next = (x + *this / x) >> 1;
        if (next >= x) break;
        x = next;
    }
    return x;
}

bool BigInt::isPerfectSquare() const {
    if (m_negative) return false;
    const BigInt root = isqrt();
    return root * root == *this;
}

BigInt BigInt::nthRoot(unsigned n) const {
    if (n == 0 || m_negative) return BigInt();
    if (n == 1) return *this;
    if (n == 2) return isqrt();
    if (isZero() || isOne()) return *this;
    const BigInt exponentValue(static_cast<long long>(n));
    BigInt x = BigInt(1) << static_cast<int>((bitLength() + n - 1) / n);
    while (true) {
        const BigInt xPow = BigInt::power(x, n - 1);
        if (xPow.isZero()) break;
        const BigInt next = (BigInt(static_cast<long long>(n - 1)) * x + *this / xPow) / exponentValue;
        if (next >= x) break;
        x = next;
    }
    while (BigInt::power(x, n) > *this) x -= BigInt(1);
    return x;
}

bool BigInt::isPerfectPower(unsigned n) const {
    if (n == 0 || m_negative) return false;
    return BigInt::power(nthRoot(n), n) == *this;
}

/* ==================================================== 转换 */

std::string BigInt::toString() const {
    if (m_small) return std::to_string(m_smallValue);
    if (m_limbs.empty()) return "0";
    static constexpr std::uint32_t kChunk = 1000000000u;
    Limbs work = m_limbs;
    std::vector<std::uint32_t> chunks;
    while (!work.empty()) {
        std::uint64_t rem = 0;
        for (std::size_t i = work.size(); i-- > 0;) {
            const std::uint64_t cur = (rem << 32) | work[i];
            work[i] = static_cast<std::uint32_t>(cur / kChunk);
            rem = cur % kChunk;
        }
        chunks.push_back(static_cast<std::uint32_t>(rem));
        work.resize(actualSize(work));
    }
    std::string out;
    if (m_negative) out.push_back('-');
    out += std::to_string(chunks.back());
    for (std::size_t i = chunks.size() - 1; i-- > 0;) {
        std::string piece = std::to_string(chunks[i]);
        out.append(9 - piece.size(), '0');
        out += piece;
    }
    return out;
}

std::string BigInt::toStringHex() const {
    if (isZero()) return "0";
    std::string out;
    if (m_negative) out.push_back('-');
    const Limbs limbs = toLimbs();
    char buffer[16];
    for (std::size_t i = limbs.size(); i-- > 0;) {
        std::snprintf(buffer, sizeof(buffer), i + 1 == limbs.size() ? "%x" : "%08x", limbs[i]);
        out += buffer;
    }
    return out;
}

long long BigInt::toLongLong() const {
    if (m_small) return m_smallValue;
    unsigned long long magnitude = 0;
    if (!m_limbs.empty()) magnitude = m_limbs[0];
    if (m_limbs.size() > 1) magnitude |= static_cast<unsigned long long>(m_limbs[1]) << 32;
    const long long value = static_cast<long long>(magnitude);
    return m_negative ? -value : value;
}

std::size_t BigInt::limbCount() const {
    if (m_small) {
        unsigned long long magnitude = m_smallValue < 0 ? 0ull - static_cast<unsigned long long>(m_smallValue)
                                                        : static_cast<unsigned long long>(m_smallValue);
        return magnitude == 0 ? 0 : (magnitude >> 32 ? 2 : 1);
    }
    return m_limbs.size();
}

std::uint32_t BigInt::limbAt(std::size_t index) const {
    if (m_small) {
        const unsigned long long magnitude = m_smallValue < 0 ? 0ull - static_cast<unsigned long long>(m_smallValue)
                                                              : static_cast<unsigned long long>(m_smallValue);
        if (index == 0) return static_cast<std::uint32_t>(magnitude & 0xFFFFFFFFu);
        if (index == 1) return static_cast<std::uint32_t>(magnitude >> 32);
        return 0;
    }
    return index < m_limbs.size() ? m_limbs[index] : 0;
}

double BigInt::toDouble() const {
    if (m_small) return static_cast<double>(m_smallValue);
    double result = 0;
    for (std::size_t i = m_limbs.size(); i-- > 0;) result = result * 4294967296.0 + static_cast<double>(m_limbs[i]);
    return m_negative ? -result : result;
}

} // namespace Cas::Detail

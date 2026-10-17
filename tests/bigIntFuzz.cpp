// SPDX-License-Identifier: MIT
// Copyright (c) 2026 小兔糖糖OvO

// 大数随机对拍：cstdlib 只做本地一致性校验（q*b+r == a、r < |b|），再覆盖 NTT 乘法路径。
#include "cas/detail/bigInt.h"

#include <cstdio>
#include <random>
#include <string>

using Cas::Detail::BigInt;

int main() {
    std::mt19937_64 rng(20261017);
    int badDivision = 0, badMultiply = 0, trials = 0;

    for (int trial = 0; trial < 2000; ++trial) {
        const int digitsA = 2 + static_cast<int>(rng() % 60);
        const int digitsB = 1 + static_cast<int>(rng() % 30);
        auto randomNumber = [&](int digits) {
            std::string text;
            text.push_back(static_cast<char>('1' + rng() % 9));
            for (int i = 1; i < digits; ++i) text.push_back(static_cast<char>('0' + rng() % 10));
            return BigInt(text);
        };
        const BigInt a = randomNumber(digitsA);
        const BigInt b = randomNumber(digitsB);
        const BigInt q = a / b, r = a % b;
        ++trials;
        if (q * b + r != a || r >= b || r.sign() < 0) {
            if (badDivision < 4) {
                std::printf("BAD div a=%s b=%s q=%s r=%s\n", a.toString().c_str(), b.toString().c_str(),
                            q.toString().c_str(), r.toString().c_str());
            }
            ++badDivision;
        }

        // 乘法与除法互逆：用 (a+b)(a-b) 展开式校验
        if ((a + b) * (a - b) != a * a - b * b) {
            if (badMultiply < 3) {
                std::printf("BAD mul a=%s b=%s\n", a.toString().c_str(), b.toString().c_str());
                std::printf("     lhs=%s\n", ((a + b) * (a - b)).toString().c_str());
                std::printf("     rhs=%s\n", (a * a - b * b).toString().c_str());
            }
            ++badMultiply;
        }
    }

    // 强制覆盖 NTT 路径（>96 limb）
    const BigInt big1 = BigInt::power(BigInt(3), 7000);
    const BigInt big2 = BigInt::power(BigInt(7), 5000);
    const BigInt product = big1 * big2;
    if (product / big2 != big1 || product % big2 != BigInt(0)) {
        std::printf("BAD big: product/big2 == big1 ? %d, product%%big2 == 0 ? %d\n",
                    product / big2 == big1, product % big2 == BigInt(0));
        ++badMultiply;
    }
    if (BigInt(product.toString()) != product) {
        std::printf("BAD big: decimal roundtrip\n");
        ++badMultiply;
    }
    if (BigInt(product.toStringHex(), 16) != product) {
        std::printf("BAD big: hex roundtrip\n");
        ++badMultiply;
    }

    std::printf("除法对拍 %d 例，失败 %d；乘法校验失败 %d\n", trials, badDivision, badMultiply);
    return (badDivision == 0 && badMultiply == 0) ? 0 : 1;
}

// 数值层与分配器的独立自检：在接入表达式树之前先跑通这一层。
#include "cas/allocator.h"
#include "cas/detail/bigInt.h"
#include "cas/memory.h"

#include <cstdio>
#include <cstring>
#include <string>

using Cas::Detail::BigInt;

namespace {

int failures = 0;

void check(bool condition, const char *what) {
    if (!condition) {
        ++failures;
        std::printf("  [FAIL] %s\n", what);
    }
}

std::string big(long long v) { return BigInt(v).toString(); }

void testSmall() {
    check(big(0) == "0", "zero");
    check((BigInt(7) + BigInt(-3)).toString() == "4", "7 + -3");
    check((BigInt(3) - BigInt(10)).toString() == "-7", "3 - 10");
    check((BigInt(-6) * BigInt(7)).toString() == "-42", "(-6)*7");
    check((BigInt(-7) / BigInt(2)).toString() == "-4", "floor(-7/2) == -4");
    check((BigInt(-7) % BigInt(2)).toString() == "1", "-7 % 2 == 1");
    check((BigInt(7) % BigInt(-2)).toString() == "-1", "7 % -2 == -1");
}

void testBig() {
    // 20! = 2432902008176640000，超过 2^61，走大数路径
    BigInt f(1);
    for (int i = 2; i <= 30; ++i) f *= BigInt(static_cast<long long>(i));
    check(f.toString() == "265252859812191058636308480000000", "30! 十进制");

    const std::string digits = f.toString();
    check(BigInt(digits).toString() == digits, "十进制往返");
    check(BigInt(f.toStringHex(), 16) == f, "十六进制往返");

    // (a+b)(a-b) == a^2 - b^2，覆盖 NTT 与 Karatsuba 路径
    BigInt a("123456789012345678901234567890123456789012345678901234567890");
    BigInt b("987654321098765432109876543210987654321098765432109876543210");
    check((a + b) * (a - b) == a * a - b * b, "(a+b)(a-b) = a^2-b^2");

    // NTT 路径：各 120 个 limb，乘法结果再除以原数应还原
    BigInt huge = BigInt::power(BigInt(7), 400);
    BigInt other = BigInt::power(BigInt(11), 300);
    check(huge * other / other == huge, "(7^400 * 11^300) / 11^300");

    BigInt q, r;
    BigInt::divMod(huge * other + BigInt(12345), huge, &q, &r);
    check(q == other && r == BigInt(12345), "大数带余除法");

    check((BigInt(-100) / BigInt(7)).toString() == "-15", "floor(-100/7)");
    check((BigInt(-100) % BigInt(7)).toString() == "5", "-100 mod 7");
}

void testRoots() {
    check(BigInt(144).isqrt().toString() == "12", "isqrt(144)");
    check(BigInt(145).isqrt().toString() == "12", "isqrt(145)");
    check(BigInt(1000000000000000000LL).isqrt().toString() == "1000000000", "isqrt(1e18)");
    check(BigInt::power(BigInt(123456789), 2).isPerfectSquare(), "大数完全平方");
    check(!BigInt::power(BigInt(123456789), 2).negated().isPerfectSquare(), "负数非完全平方");
    check(BigInt(1000).nthRoot(3).toString() == "10", "1000 的立方根");
    check(BigInt::power(BigInt(97), 5).isPerfectPower(5), "97^5 是完全五次幂");
    check(BigInt::gcd(BigInt::power(BigInt(6), 20), BigInt::power(BigInt(6), 12)).toString() ==
              BigInt::power(BigInt(6), 12).toString(),
          "gcd(6^20, 6^12)");
    check(BigInt(17).modInverse(BigInt(3120)).toString() == "2753", "17 在 3120 下的逆元");
}

void testPool() {
    using namespace Cas::Detail;
    check(Cas::Allocator::initialize(), "池初始化");
    const auto before = Cas::Allocator::statistics();
    {
        Cas::Memory::Vector<int> values;
        for (int i = 0; i < 1000; ++i) values.push_back(i);
        Cas::Memory::String text(1000, 'x');
        Cas::Memory::Map<int, int> table;
        for (int i = 0; i < 200; ++i) table[i] = i * i;
        check(values[999] == 999 && text.size() == 1000 && table[7] == 49, "池内容器可用");
    }
    const auto after = Cas::Allocator::statistics();
    check(after.liveBlocks == before.liveBlocks, "容器析构后块数归零");
    check(after.liveBytes == before.liveBytes, "容器析构后字节归零");
    check(after.totalBlocks > before.totalBlocks, "累计分配次数增长");
    check(after.peakBytes > 0, "记录峰值");
    check(Cas::Allocator::quiescent(), "池已静默");
}

void testPoolSemantics() {
    using Cas::Detail::MemoryPool;
    alignas(16) unsigned char buffer[4096];
    MemoryPool pool(buffer, sizeof(buffer));
    check(pool.ready(), "外部缓冲区初始化");

    void *a = pool.allocate(100, 16);
    void *b = pool.allocate(200, 16);
    void *c = pool.allocate(300, 16);
    check(a && b && c, "三块分配成功");
    check(reinterpret_cast<std::uintptr_t>(a) % 16 == 0, "16 字节对齐");
    std::memset(a, 0xAB, 100);
    std::memset(b, 0xCD, 200);
    check(pool.statistics().liveBlocks == 3, "活动块计数");

    check(pool.release(b), "释放中间块");
    check(!pool.release(b), "重复释放被拒绝");
    int stackValue = 0;
    check(!pool.release(&stackValue), "池外地址被拒绝");

    // 释放 a 后应与空前缀合并，再能容纳接近原前缀的大块
    check(pool.release(a), "释放首块");
    void *big1 = pool.allocate(280, 8);
    check(big1 != nullptr, "合并后可分配");
    check(pool.release(big1) && pool.release(c), "清理");

    void *grown = pool.allocate(64, 8);
    void *grown2 = pool.reallocate(grown, 900);
    check(grown2 != nullptr, "扩容成功");
    check(!pool.reallocate(grown, 900), "扩容后旧地址失效");
    check(pool.release(grown2), "释放扩容块");
    check(pool.statistics().liveBlocks == 0, "全部释放");
    check(pool.statistics().liveBytes == 0, "字节归零");
    check(!pool.allocate(100000, 8), "容量不足返回空指针");
}

} // namespace

int main() {
    std::printf("allocate/数值层自检\n");
    std::fflush(stdout);
    testSmall(); std::printf("small ok\n"); std::fflush(stdout);
    testBig(); std::printf("big ok\n"); std::fflush(stdout);
    testRoots(); std::printf("roots ok\n"); std::fflush(stdout);
    testPool(); std::printf("pool ok\n"); std::fflush(stdout);
    testPoolSemantics(); std::printf("poolsem ok\n"); std::fflush(stdout);
    if (failures == 0) {
        std::printf("全部通过\n");
        return 0;
    }
    std::printf("%d 项失败\n", failures);
    return 1;
}

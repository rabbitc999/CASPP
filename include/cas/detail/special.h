#ifndef CAS_DETAIL_SPECIAL_H
#define CAS_DETAIL_SPECIAL_H

/**
 * @file special.h
 * @brief 特殊函数：Gamma / Beta / 误差函数 / 指数积分 / Bessel / 组合数 / 谐波数。
 *
 * 设计原则与 sympy 的 `functions.special` 一致：**每个函数在自己的 `eval` 里做
 * 尽可能多的精确折叠，做不了就原样留下**。区别在于本工程不做自动求值猜测——
 * 只有在数学上无争议的恒等式才折叠（`Gamma(n) = (n-1)!`、`Gamma(1/2) = √π`、
 * `erf(0) = 0`），其余一律保留为函数节点，由用户显式调用 `specialsimplify` 或
 * 求值命令处理。
 *
 * 所有函数都是"符号优先"的：参数是数值时才走数值路径（整数阶乘、半整数 Gamma），
 * 是符号时走恒等式重写（`Gamma(z+1) = z·Gamma(z)`、`BesselJ(-n, z) = (-1)^n J(n, z)`）。
 */

#include "cas/detail/node.h"
#include "cas/detail/simplify.h"

namespace Cas::Detail {

/* ------------------------------------------------------------ Gamma 族 */

/**
 * @brief Gamma 函数 Γ(z)。
 *
 * 折叠规则（按优先级）：
 *  - `Γ(n)`，n 正整数 → `(n-1)!`
 *  - `Γ(1/2 + n)`，n 非负整数 → `(2n)!/(4^n n!) · √π`
 *  - `Γ(1/2 - n)`，n 非负整数 → `(-4)^n n!/(2n)! · √π`
 *  - `Γ(0)`、`Γ(负整数)` → 发散，返回未求值的 `Γ(z)` 并置 `out->diverges`
 */
Ptr gamma(const Ptr &argument);

/** @brief 阶乘 n!，n 为非负整数时精确求值；其余情形转成 `Γ(n+1)`。 */
Ptr factorial(const Ptr &argument);

/** @brief 双阶乘 n!!：`n(n-2)(n-4)…`。 */
Ptr doubleFactorial(const Ptr &argument);

/**
 * @brief Beta 函数 B(a, b)。
 *
 * 整数/半整数组合用 `B(a,b) = Γ(a)Γ(b)/Γ(a+b)` 精确求值；
 * 其余情形保留并尝试 `B(a,b) = B(b,a)` 的规范排序。
 */
Ptr beta(const Ptr &first, const Ptr &second);

/** @brief 二项式系数 C(n, k)，n 为非负整数且 k 在 [0, n] 内时精确求值。 */
Ptr binomial(const Ptr &n, const Ptr &k);

/**
 * @brief 谐波数 H(n) 与广义谐波数 H(n, r) = Σ_{k=1..n} 1/k^r。
 *
 * n 为非负整数时展开成有理数（上限 10000，防爆炸）；n 为符号时用递推恒等式。
 */
Ptr harmonic(const Ptr &n, const Ptr &order);

/* ------------------------------------------------------- 误差函数族 */

/** @brief 误差函数 erf(z) = 2/√π ∫₀^z e^(-t²) dt。 */
Ptr erf(const Ptr &argument);
/** @brief 互补误差函数 erfc(z) = 1 - erf(z)。 */
Ptr erfc(const Ptr &argument);
/** @brief 虚误差函数 erfi(z) = -i·erf(i·z)，实数域上用级数定义。 */
Ptr erfi(const Ptr &argument);

/* ---------------------------------------------------- 指数积分与对数积分 */

/** @brief 指数积分 Ei(x) = -∫_{-x}^∞ e^(-t)/t dt。 */
Ptr expIntegralEi(const Ptr &argument);
/** @brief 对数积分 li(x) = Ei(ln x)。 */
Ptr logarithmicIntegral(const Ptr &argument);

/* ------------------------------------------------------------- Bessel */

/**
 * @brief 第一类 Bessel 函数 J_ν(z)。
 *
 * 折叠：`J_ν(0)` = 1（ν=0）/ 0（ν 正整数）；整数阶用 `J_{-n} = (-1)^n J_n` 规范化；
 * 半整数阶保留（不做初等化）。
 */
Ptr besselJ(const Ptr &order, const Ptr &argument);
/** @brief 第二类 Bessel 函数 Y_ν(z)。 */
Ptr besselY(const Ptr &order, const Ptr &argument);
/** @brief 修正 Bessel 函数 I_ν(z)。 */
Ptr besselI(const Ptr &order, const Ptr &argument);
/** @brief 修正 Bessel 函数 K_ν(z)。 */
Ptr besselK(const Ptr &order, const Ptr &argument);

/* ------------------------------------------------------------- 统一入口 */

/**
 * @brief 特殊函数的总折叠入口。
 *
 * 输入一个 Function 节点（或任意表达式），自底向上把可折叠的特殊函数替换成
 * 精确值，再做一次通用化简。这是 `specialsimplify` 命令的实现核心，也被
 * `simplify` 内部的数值折叠阶段调用。
 */
Ptr simplifySpecial(const Ptr &expression);

/**
 * @brief 尝试对单个函数调用做折叠。
 * @return 折叠结果；无法折叠时返回 `Ptr()`（空指针），调用方应保留原节点。
 */
Ptr evaluateSpecialCall(const std::string &name, const Memory::Vector<Ptr> &arguments);

} // namespace Cas::Detail

#endif // CAS_DETAIL_SPECIAL_H

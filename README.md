# caspp —— C++17 代数符号计算器

一个从零实现的计算机代数系统（CAS）核心。符号表达式树、大整数运算、多项式算法、
微积分、方程求解与 Groebner 基，全部基于自研内存池分配。

参考了 `gitee.com/rabbitc999/cas` 的内存分配器设计与 sympy 的算法组织方式，
命名风格沿用 Qt 式驼峰（`Cas::Engine`、`calculate()`、`memoryStats()`）。

## 构建与运行

```bash
mingw32-make            # 构建 build/caspp 与全部测试
mingw32-make test       # 运行全部测试（等价于 make check）
mingw32-make run        # 启动 REPL
mingw32-make bench      # 跑一组较重的表达式
mingw32-make clean      # 清理 build/
```

需要 g++ 支持 C++17。产物全部落在 `build/`，源码树保持干净。

### 命令行

```bash
build/caspp -c "diff(x^2*sin(x), x)"      # 单条求值
build/caspp -q -c "factor(x^2-1)"         # 只输出结果，不打印提示
build/caspp                               # 进入 REPL
```

参数：`-c/--command`（可重复）、`-q/--quiet`、`-u/--unicode`、`-v/--version`、`-h/--help`。

### REPL

```
> simplify(2x + 3x)
5*x
  [LaTeX] 5 x
  [0.00 ms]
> help
> mem
> unicode on
> exit
```

`%` 引用上一次结果。

## 支持的运算

| 类别 | 命令 | 示例 |
| --- | --- | --- |
| 化简 | `simplify` / 直接输入 | `simplify(2x/2)` → `x` |
| 展开 | `expand` | `expand((x+1)^3)` → `x^3 + 3*x^2 + 3*x + 1` |
| 因式分解 | `factor` | `factor(2x^2+7x+3)` → `(x + 3)*(2*x + 1)` |
| 通分 | `together` | `together(1/x + 1/(x+1))` |
| 三角化简 | `trigsimp` | `trigsimp(sin(x)^2+cos(x)^2)` → `1` |
| 求导 | `diff` | `diff(sin(x^2), x)` → `2*x*cos(x^2)` |
| 不定积分 | `integrate(f, x)` | `integrate(x^3*exp(x), x)` |
| 定积分 | `integrate(f, x, a, b)` 或 `integrate(f, a, b)` | `integrate(1/(x^2+1), x, 0, 1)` → `1/4*pi` |
| 极限 | `limit(f, x, p)` / `limit(f, p, x)` / `limit(f, x=p)` | `limit(sin(x)/x, x, 0)` → `1` |
| Taylor 展开 | `taylor(f, 0, n)` | `taylor(cos(x), 0, 6)` |
| 解方程 | `solve(eq, x)` | `solve(x^2-5x+6=0, x)` → `x = 2, x = 3` |
| 解方程组 | `solve([eq1, eq2], [x, y])` | 高斯-约当消元 |
| 不等式 | 直接输入 `>` `<` | `2x+1<3` → `x < 1` |
| 多项式 GCD | `gcd(p, q)` | |
| 结式 | `resultant(p, q, x)` | |
| Groebner 基 | `groebner([f, g], [x, y])` | `{y^2 - 1/2, x - y}` |
| LaTeX | `latex(expr)` | `\frac{1}{x^{1}}` |
| 内存统计 | `mem` | 池容量/活动字节/峰值 |

### 积分能力

- 基本初等函数（幂、指数、对数、三角、双曲、反三角）
- 线性换元（自动识别内层导数为常数）
- 分部积分（LIATE 式的 `uRank` 优先级，递归深度上限 6）
- 有理函数：分母一次 → 多项式除法 + 对数；分母二次 → 按判别式分三种情形，得到对数组合或 `atan`
- 数值上更复杂的被积函数会被明确判定为「超出当前支持的积分规则集」，不会返回错误答案

### 极限策略

依次尝试：直接代入 → 有理式约分（多项式 GCD）→ 洛必达 → 分子分母各自 Taylor 展开到 8 阶取首个非零项之商。
`±∞` 方向单独处理分子分母次数比较。

## 架构

```
include/cas/
  allocator.h        公开的分配器接口与 MemoryStats
  engine.h           门面：Result / Engine（全部对外 API）
  memory.h           PoolAllocator<T> 与容器别名
  detail/
    config.h         池大小、小整数上限、NTT 参数
    memoryPool.h     首次适配链表池
    bigInt.h         任意精度整数
    node.h           表达式树节点与引用计数 Ptr
    simplify.h       makeAdd/makeMul/makePow/makeNegate/simplify
    polynomial.h     一元多项式与 expand/factor/gcd/resultant/together
    calculus.h       diff/integrate/limit/taylor/substitute
    solver.h         solve/solveInequality/solveLinearSystem
    trigonometry.h   trigsimp 及相关重写规则
    multivariate.h   多元多项式与 Groebner 基
    parser.h         递归下降解析器
    printer.h        文本与 LaTeX 输出
```

### 内存模型

不使用全局 GC，改用**引用计数 + 不可变节点 + RAII `Ptr` 守卫**：

- 节点创建时 `m_refs = 1`，由返回它的 `Ptr` 持有
- `Ptr` 拷贝即 `retain`，析构即 `release`，计数归零时归还池
- `Ptr::borrow(NodePtr)` 用于临时观察，不改变计数
- 所有节点分配都走自研内存池（首次适配 + 块拆分 + 相邻空闲块双向合并），
  池未就绪时自动退回 `::operator new`

这样避免了参考项目里 mark-sweep GC 的悬垂风险，同时保留了池化带来的分配性能。

### 大整数

基 `2^32` 小端 limb，符号独立存储，`±2^61` 以内走内联表示。

- 乘法分段：≤32 limb 教科书算法 → ≤96 limb Karatsuba → 更大规模走三段分解 + NTT（模数 `2013265921`，原根 31，每 limb 拆 4 段 × 8 位）
- 除法用移位相减长除（参考 magel 的 `divmodBase`）
- 精确开方、k 次整数根（二分）、GCD

## 测试

五套回归测试，`make test` 全绿：

| 测试 | 内容 |
| --- | --- |
| `bigIntFuzz` | 除法对拍 2000 例 + 乘法交叉校验 |
| `numericTest` | 数值运算与内存池语义 |
| `expressionTest` | 表达式树、化简、打印 |
| `polynomialTest` | 展开/因式分解/通分/GCD |
| `calculusTest` | 118 项：解析、化简、展开、分解、求导、积分、极限、Taylor、求解、不等式、三角、Groebner、LaTeX、内存 |

## 已知边界

- 积分不对全部有理函数做完整部分分式；分母三次及以上、或分母不可约高次时返回「超出支持的规则集」
- 求解只覆盖一次、二次方程与线性方程组，高次方程依赖有理根定理
- 三角化简是规则驱动的有限轮不动点迭代，不做任意恒等式的通搜
- 输出保持精确有理数，不自动数值化

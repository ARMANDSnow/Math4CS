---
course: Math4CS
chapter: NA2
topic: 线性系统与 LU 分解 Linear Systems and LU Decomposition
teacher: 应凯
date: 2026-06-19
tags: [Math4CS, 数值算法, 线性方程组, 高斯消元, LU分解, 三角矩阵]
---

# NA2 线性系统与 LU 分解（Linear Systems and LU Decomposition）

> 来源：应凯老师 `2-LinearSystemAndLU.pptx`（教材 Solomon《Numerical Algorithms》第 2–3 章）· 考纲第 2 章
> 一句话定位：怎么又快又稳地解线性方程组 $Ax=b$。核心招数是把矩阵拆成"两个三角矩阵相乘"（LU 分解），因为三角矩阵最好解。
> **考纲地位：LU 分解会考计算。** 这是个手算大题候选，务必练熟整套流程。

**本章导航**
1. 这章解决什么问题
2. 从零铺垫：线性方程组与三角矩阵
3. 可解性：什么时候有唯一解
4. 高斯消元 + 三种"行操作矩阵"
5. LU 分解：把消元过程"记录"下来
6. 用 LU 解方程：前代 + 回代
7. 手算例题（HW4.2 原题，全流程）
8. 考点雷达
9. 易错点 & 陷阱
10. 本章速查卡
11. 练习题（自测）

---

## 1. 这章解决什么问题

无数实际问题最后都归结为解一个线性方程组：

$$Ax=b$$

其中 $A$ 是已知矩阵、$b$ 是已知向量、$x$ 是要求的未知量。问题是：$A$ 可能很大，直接求逆 $x=A^{-1}b$ 又慢又不稳。

这章的核心思想：**有些矩阵天生好解**（三角矩阵、对角矩阵），那就想办法把一般的 $A$ **拆成好解的矩阵相乘**，再分步解。LU 分解就是把 $A$ 拆成"下三角 $L$ × 上三角 $U$"。

---

## 2. 从零铺垫：线性方程组与三角矩阵

### 2.1 线性方程组的矩阵写法

三个方程三个未知数，比如

$$\begin{cases}2x+4y=2\\3x+5y=4\end{cases}\quad\Longleftrightarrow\quad \underbrace{\begin{bmatrix}2&4\\3&5\end{bmatrix}}_{A}\underbrace{\begin{bmatrix}x\\y\end{bmatrix}}_{x}=\underbrace{\begin{bmatrix}2\\4\end{bmatrix}}_{b}$$

矩阵 $A$ 的每一行 = 一个方程的系数。

### 2.2 为什么"三角矩阵"最好解

**上三角矩阵**（对角线下方全是 0）：

$$\begin{bmatrix}2&4\\0&-1\end{bmatrix}\begin{bmatrix}x\\y\end{bmatrix}=\begin{bmatrix}2\\1\end{bmatrix}$$

最后一行只剩一个未知数：$-1\cdot y=1\Rightarrow y=-1$。把它代回上一行，又只剩一个未知数 $x$。这种"从最后一行往上一个个代回去"叫**回代（back-substitution）**。

**下三角矩阵**同理，从第一行往下解，叫**前代（forward-substitution）**。

直觉：三角矩阵之所以好解，是因为它让未知数"一个一个暴露出来"，每步只需解一个一元方程。**整个数值线代的策略，就是想方设法把矩阵变成三角形 / 对角形 / 正交形。**

---

## 3. 可解性：什么时候有唯一解

线性系统 $Ax=b$（$A\in\mathbb{R}^{m\times n}$）有三种情形：

| 形状 | 名称 | 解的情况 |
|---|---|---|
| $m=n$ 方阵且满秩 | 适定 | **唯一解** |
| $m>n$ 高瘦矩阵（方程多于未知数） | 超定 Over-determined | 一般**无解**，退而求最小二乘（见 [[NA3-最小二乘与参数回归]]） |
| $m<n$ 矮胖矩阵（方程少于未知数） | 欠定 Under-determined | **无穷多解**，需附加条件（如最小范数，见 [[NA3-最小二乘与参数回归]]） |

PPT 原话：没有矮胖系统能有唯一解；高瘦系统会遇到无解的右端项。**本章先聚焦最干净的情形：$A$ 是 $n\times n$ 可逆方阵，有唯一解。**

---

## 4. 高斯消元 + 三种"行操作矩阵"

高斯消元就是初中解方程组的"加减消元法"，只不过我们要把每一步操作写成**矩阵乘法**，这样才能"记录"下来做成 LU。三种基本行操作各对应一种矩阵：

### 4.1 置换矩阵 Permutation（换行）

交换两行。例：交换第 1、2 行的置换矩阵

$$P=\begin{bmatrix}0&1\\1&0\end{bmatrix}$$

用处：当主元（pivot，即对角线上用来消元的元素）是 0 时，换一行上来。逆就是它自己（再换回去）。

### 4.2 缩放矩阵 Scaling（某行乘常数）

把第 $k$ 行整体乘 $a_k$，对应对角矩阵 $S=\mathrm{diag}(a_1,\dots,a_n)$。逆 = 各对角元取倒数。

### 4.3 消元矩阵 Elimination（核心）⭐

用第 $i$ 行的倍数去消掉第 $j$ 行某个位置。这是 LU 的关键。

例：要用第 1 行消掉第 2 行的首元（把 $(2,1)$ 位置变 0），乘子（multiplier）$m_{21}=\dfrac{a_{21}}{a_{11}}$，消元矩阵：

$$E=\begin{bmatrix}1&0\\-m_{21}&1\end{bmatrix}$$

**它的逆只需把乘子变号**（因为"减回去"就是"加回来"）：

$$E^{-1}=\begin{bmatrix}1&0\\+m_{21}&1\end{bmatrix}$$

> 这个"逆=乘子变号"的性质，是 LU 里 $L$ 能轻松写出来的根本原因。

### 4.4 消元的两半：前向消元 + 回代

- **前向消元（Forward Elimination）**：从上往下，逐列把主元下方的元素消成 0，最后把 $A$ 变成上三角 $U$。
- **回代（Back-Substitution）**：从下往上，逐个解出未知数。

**计算量**：每选一个主元要对其下方各行做消元，总共约 $\dfrac{n^3}{3}$ 次运算，即 $O(n^3)$。这是高斯消元的代价。

---

## 5. LU 分解：把消元过程"记录"下来

### 5.1 核心思想

前向消元把 $A$ 变成上三角 $U$，过程就是一串消元矩阵左乘：

$$E_k\cdots E_2 E_1\,A=U$$

把这些 $E$ 全挪到右边（取逆）：

$$A=\underbrace{E_1^{-1}E_2^{-1}\cdots E_k^{-1}}_{L}\;U$$

**奇妙之处**：这一堆消元矩阵的逆乘起来，恰好是一个**下三角矩阵 $L$**，而且它的对角线下方元素**就是各步的消元乘子 $m_{ij}$**（直接抄下来，不用算！）。于是：

$$\boxed{A=LU}$$

- $L$（Lower）：下三角，对角线全为 1，下方填消元乘子 $m_{ij}$。
- $U$（Upper）：上三角，就是前向消元的最终结果。

### 5.2 为什么要 LU（动机）

如果只解一次 $Ax=b$，直接高斯消元就行，何必分解？因为实际中**常常同一个 $A$、要解很多不同的 $b$**（$Ax=b_1, Ax=b_2,\dots$）。

- 高斯消元每换一个 $b$ 都要重做 $O(n^3)$。
- LU **只分解一次** $O(n^3)$，之后每个 $b$ 只需两次三角求解，各 $O(n^2)$，便宜得多。

> 一句话：**LU = 把 $O(n^3)$ 的力气一次性花掉，换来后续每个 $b$ 都只花 $O(n^2)$。**

### 5.3 紧凑存储（了解即可）

$L$ 的对角线恒为 1 不用存，$U$ 的下三角是 0 不用存——于是 $L$、$U$ 可以**挤在同一个矩阵里**存，省一半内存。考试一般不细考，知道有这回事即可。

---

## 6. 用 LU 解方程：前代 + 回代

有了 $A=LU$，解 $Ax=b$ 拆成两步三角求解：

$$Ax=b\;\Rightarrow\;LUx=b\;\Rightarrow\;\begin{cases}\textbf{① 前代：}\;Ly=b\quad(\text{解出中间量 }y=Ux)\\[2pt]\textbf{② 回代：}\;Ux=y\quad(\text{解出最终 }x)\end{cases}$$

两步都是三角系统，各 $O(n^2)$。流程：先用下三角 $L$ 从上往下解出 $y$，再用上三角 $U$ 从下往上解出 $x$。

---

## 7. 手算例题（HW4.2 原题，全流程）

> **题**：用高斯消元解下列方程组，写出每一步对应的消元矩阵，并把左边矩阵分解为 $A=LU$。
> $$\begin{bmatrix}2&4\\3&5\end{bmatrix}\begin{bmatrix}x\\y\end{bmatrix}=\begin{bmatrix}2\\4\end{bmatrix}$$

**第 1 步：前向消元，消去 $(2,1)$ 位置的 3**

乘子 $m_{21}=\dfrac{a_{21}}{a_{11}}=\dfrac{3}{2}$。做"第 2 行 $-\dfrac32\times$第 1 行"。对应消元矩阵：

$$E_1=\begin{bmatrix}1&0\\-\tfrac32&1\end{bmatrix}$$

作用到 $A$：

$$E_1A=\begin{bmatrix}1&0\\-\tfrac32&1\end{bmatrix}\begin{bmatrix}2&4\\3&5\end{bmatrix}=\begin{bmatrix}2&4\\\;3-\tfrac32\cdot2&\;5-\tfrac32\cdot4\end{bmatrix}=\begin{bmatrix}2&4\\0&-1\end{bmatrix}=U$$

**第 2 步：写出 $L=E_1^{-1}$（乘子变号）**

$$L=E_1^{-1}=\begin{bmatrix}1&0\\\tfrac32&1\end{bmatrix},\qquad U=\begin{bmatrix}2&4\\0&-1\end{bmatrix}$$

**验证 $LU=A$**：

$$\begin{bmatrix}1&0\\\tfrac32&1\end{bmatrix}\begin{bmatrix}2&4\\0&-1\end{bmatrix}=\begin{bmatrix}2&4\\3&\;6-1\end{bmatrix}=\begin{bmatrix}2&4\\3&5\end{bmatrix}=A\;\checkmark$$

**第 3 步：前代 $Ly=b$**（解中间量 $y$）

$$\begin{bmatrix}1&0\\\tfrac32&1\end{bmatrix}\begin{bmatrix}y_1\\y_2\end{bmatrix}=\begin{bmatrix}2\\4\end{bmatrix}\;\Rightarrow\;\begin{cases}y_1=2\\[2pt]\tfrac32\cdot2+y_2=4\Rightarrow y_2=1\end{cases}\;\Rightarrow\;y=\begin{bmatrix}2\\1\end{bmatrix}$$

**第 4 步：回代 $Ux=y$**（解最终 $x$）

$$\begin{bmatrix}2&4\\0&-1\end{bmatrix}\begin{bmatrix}x\\y\end{bmatrix}=\begin{bmatrix}2\\1\end{bmatrix}\;\Rightarrow\;\begin{cases}-1\cdot y=1\Rightarrow y=-1\\[2pt]2x+4(-1)=2\Rightarrow x=3\end{cases}$$

**最终解**：

$$\boxed{x=3,\quad y=-1}$$

验算：$\begin{bmatrix}2&4\\3&5\end{bmatrix}\begin{bmatrix}3\\-1\end{bmatrix}=\begin{bmatrix}6-4\\9-5\end{bmatrix}=\begin{bmatrix}2\\4\end{bmatrix}\;\checkmark$

---

## 8. 考点雷达

- 考纲明确：**LU 分解会考计算**——这是手算大题候选。要熟练"前向消元写 $U$ → 乘子变号写 $L$ → 前代解 $y$ → 回代解 $x$"全流程。
- 必练：给一个 $2\times2$ 或 $3\times3$ 的 $A$，要求写每步消元矩阵 + 给出 $L,U$ +（可能）解方程。
- 概念可能被问：为什么用 LU（多个右端项复用分解）；高斯消元的复杂度 $O(n^3)$；三角矩阵为什么好解。
- 横向联系：对称正定矩阵有更省的 **Cholesky 分解** $A=LL^\top$（见 [[NA3-最小二乘与参数回归]]）；正交化分解见 [[NA5-列空间与QR]]。

---

## 9. 易错点 & 陷阱

1. **$L$、$U$ 谁是谁**：$L$ 下三角（Lower，对角线全 1，放乘子），$U$ 上三角（Upper，消元结果）。别写反。
2. **乘子符号**：消元矩阵 $E$ 里乘子是 **$-m$**（做减法），而 $L=E^{-1}$ 里乘子是 **$+m$**（变号）。$L$ 里直接填正的乘子 $m_{ij}=a_{ij}/a_{jj}$。
3. **$L$ 的对角线**：恒为 1，别填别的。
4. **顺序别乱**：先前代 $Ly=b$ 求 $y$，再回代 $Ux=y$ 求 $x$。把 $L$、$U$ 用反或顺序颠倒会全错。
5. **主元为 0**：若某步主元是 0，不能直接消元，要先用置换矩阵换行（这叫选主元 pivoting），否则除以 0。
6. **多步消元的 $L$**：$3\times3$ 时有多个乘子，$L$ 的下三角各位置 $(i,j)$ 填的是"用第 $j$ 行消第 $i$ 行"的乘子 $m_{ij}$，对号入座别填错位置。

---

## 10. 本章速查卡

**目标**：解 $Ax=b$（$A$ 为 $n\times n$ 可逆方阵）。

**LU 分解**：$A=LU$，$L$ 下三角（对角线 1，下方填消元乘子 $m_{ij}=a_{ij}/a_{jj}$），$U$ 上三角（前向消元结果）。

**消元矩阵**：消第 $i$ 行用 $E$（含 $-m$），其逆含 $+m$；这些逆乘起来就是 $L$。

**解方程两步**：
1. 前代 $Ly=b$ → 求中间量 $y$（从上往下）
2. 回代 $Ux=y$ → 求最终 $x$（从下往上）

**复杂度**：分解 $O(n^3)$；之后每个右端项 $b$ 各 $O(n^2)$。这是用 LU（而非每次重做消元）的理由。

**可解性**：方阵满秩→唯一解；超定（高瘦）→最小二乘；欠定（矮胖）→最小范数解。

**手算口诀**：前向消元写 $U$ → 乘子变号写 $L$ → 前代解 $y$ → 回代解 $x$ → 代回验算。

---

## 11. 练习题（自测）

> 仿期末 / 作业风格，全新设定与数字，与上文例题、原作业均不重复。建议先独立做完，再展开核对答案。

### 练习 1 · LU 分解（$3\times3$ 手算 $L,U$）（仿真题 Q7，⭐⭐⭐，约 8 分）

对下列矩阵做 LU 分解。要求写出前向消元各步的乘子，给出下三角 $L$（对角线全 1）与上三角 $U$，并验证 $LU=A$。

$$A=\begin{pmatrix}2&-1&1\\4&1&4\\-2&10&9\end{pmatrix}$$

<details>
<summary>📖 参考答案</summary>

**第 1 步：用第 1 行消去第 1 列下方两元。** 主元 $a_{11}=2$。

乘子 $m_{21}=\dfrac{a_{21}}{a_{11}}=\dfrac{4}{2}=2$，$m_{31}=\dfrac{a_{31}}{a_{11}}=\dfrac{-2}{2}=-1$。

- 第 2 行 $-\,2\times$第 1 行：$[\,4-4,\;1-(-2),\;4-2\,]=[\,0,\;3,\;2\,]$
- 第 3 行 $-(-1)\times$第 1 行（即第 3 行 $+$第 1 行）：$[\,-2+2,\;10+(-1),\;9+1\,]=[\,0,\;9,\;10\,]$

得中间矩阵

$$\begin{pmatrix}2&-1&1\\0&3&2\\0&9&10\end{pmatrix}$$

**第 2 步：用新第 2 行消去第 2 列下方。** 主元为 $3$。

乘子 $m_{32}=\dfrac{9}{3}=3$。第 3 行 $-\,3\times$第 2 行：$[\,0,\;9-9,\;10-6\,]=[\,0,\;0,\;4\,]$。得

$$U=\begin{pmatrix}2&-1&1\\0&3&2\\0&0&4\end{pmatrix}$$

**第 3 步：把乘子对号入座填进 $L$**（$L$ 对角线全 1，位置 $(i,j)$ 填 $m_{ij}$）：

$$L=\begin{pmatrix}1&0&0\\2&1&0\\-1&3&1\end{pmatrix}$$

**验证 $LU=A$**：

$$\begin{pmatrix}1&0&0\\2&1&0\\-1&3&1\end{pmatrix}\begin{pmatrix}2&-1&1\\0&3&2\\0&0&4\end{pmatrix}=\begin{pmatrix}2&-1&1\\4&-2+3&2+2\\-2&1+9&-1+6+4\end{pmatrix}=\begin{pmatrix}2&-1&1\\4&1&4\\-2&10&9\end{pmatrix}=A\;\checkmark$$

</details>

### 练习 2 · 用 LU 解方程组（前代 + 回代）（仿 HW4.2，⭐⭐，约 7 分）

已知某矩阵 $A$ 的 LU 分解为

$$L=\begin{pmatrix}1&0&0\\3&1&0\\2&-1&1\end{pmatrix},\qquad U=\begin{pmatrix}1&2&-1\\0&-3&5\\0&0&2\end{pmatrix}$$

利用该分解（不要重新消元）求解 $Ax=b$，其中 $b=\begin{pmatrix}1\\1\\2\end{pmatrix}$。请分别写出前代得到的中间量 $y$ 与回代得到的最终解 $x$。

<details>
<summary>📖 参考答案</summary>

$Ax=b\Rightarrow LUx=b$。令 $y=Ux$，分两步三角求解。

**第 1 步：前代 $Ly=b$（从上往下，$L$ 对角线全 1）。**

$$\begin{pmatrix}1&0&0\\3&1&0\\2&-1&1\end{pmatrix}\begin{pmatrix}y_1\\y_2\\y_3\end{pmatrix}=\begin{pmatrix}1\\1\\2\end{pmatrix}$$

- $y_1=1$
- $3y_1+y_2=1\Rightarrow y_2=1-3=-2$
- $2y_1-y_2+y_3=2\Rightarrow 2-(-2)+y_3=2\Rightarrow y_3=2-4=-2$

$$y=\begin{pmatrix}1\\-2\\-2\end{pmatrix}$$

**第 2 步：回代 $Ux=y$（从下往上）。**

$$\begin{pmatrix}1&2&-1\\0&-3&5\\0&0&2\end{pmatrix}\begin{pmatrix}x_1\\x_2\\x_3\end{pmatrix}=\begin{pmatrix}1\\-2\\-2\end{pmatrix}$$

- $2x_3=-2\Rightarrow x_3=-1$
- $-3x_2+5x_3=-2\Rightarrow -3x_2-5=-2\Rightarrow -3x_2=3\Rightarrow x_2=-1$
- $x_1+2x_2-x_3=1\Rightarrow x_1-2+1=1\Rightarrow x_1=2$

$$\boxed{x=\begin{pmatrix}2\\-1\\-1\end{pmatrix}}$$

**验算**：可先算出 $A=LU=\begin{pmatrix}1&2&-1\\3&3&2\\2&7&-5\end{pmatrix}$，则

$$Ax=\begin{pmatrix}1&2&-1\\3&3&2\\2&7&-5\end{pmatrix}\begin{pmatrix}2\\-1\\-1\end{pmatrix}=\begin{pmatrix}2-2+1\\6-3-2\\4-7+5\end{pmatrix}=\begin{pmatrix}1\\1\\2\end{pmatrix}=b\;\checkmark$$

</details>

### 练习 3 · 高斯消元与选主元（pivoting）（仿真题 Q7 概念延伸，⭐⭐，约 6 分）

对下列矩阵做前向消元：

$$A=\begin{pmatrix}2&1&-1\\4&2&3\\-2&3&1\end{pmatrix}$$

(1) 消去第 1 列下方两元后，观察第 2 个主元（$(2,2)$ 位置）出现什么问题？
(2) 该怎么处理才能继续消元？请写出对应的置换矩阵 $P$，并完成消元得到上三角 $U$。
(3) 这个 $A$ 是否可逆（即 $Ax=b$ 是否有唯一解）？请用 $\det A$ 说明，并指出"主元出现 0"是否就等于"矩阵不可逆"。

<details>
<summary>📖 参考答案</summary>

**(1) 第 1 列消元。** 主元 $a_{11}=2$，乘子 $m_{21}=\dfrac{4}{2}=2$，$m_{31}=\dfrac{-2}{2}=-1$。

- 第 2 行 $-\,2\times$第 1 行：$[\,0,\;2-2,\;3+2\,]=[\,0,\;0,\;5\,]$
- 第 3 行 $+$第 1 行：$[\,0,\;3+1,\;1-1\,]=[\,0,\;4,\;0\,]$

$$\begin{pmatrix}2&1&-1\\0&\mathbf{0}&5\\0&4&0\end{pmatrix}$$

**问题**：第 2 个主元 $(2,2)=0$，无法用它作除数继续消元（会除以 0）。

**(2) 处理：选主元（pivoting）——把下方非零行换上来。** 第 3 行的 $(3,2)=4\neq0$，故交换第 2、3 行，对应置换矩阵

$$P=\begin{pmatrix}1&0&0\\0&0&1\\0&1&0\end{pmatrix}$$

换行后矩阵已是上三角（$(3,2)$ 位置恰为 0，$m_{32}=0$ 无需再消）：

$$U=\begin{pmatrix}2&1&-1\\0&4&0\\0&0&5\end{pmatrix}$$

（即得到带选主元的分解 $PA=LU$，其中 $L=\begin{pmatrix}1&0&0\\-1&1&0\\2&0&1\end{pmatrix}$，可自行验证 $LU=PA$。）

**(3) 可逆性。** 沿第 1 列展开：

$$\det A=2\,(2\cdot1-3\cdot3)-1\,(4\cdot1-3\cdot(-2))+(-1)\,(4\cdot3-2\cdot(-2))$$
$$=2(2-9)-1(4+6)-1(12+4)=2(-7)-10-16=-14-10-16=-40\neq0$$

所以 $A$ 可逆，$Ax=b$ 有唯一解。

**结论**：消元中"某个主元变成 0"**并不**意味着矩阵不可逆——本题 $\det A=-40\neq0$ 却出现了零主元。它只说明当前行顺序不适合直接消元，换行（选主元）后即可正常进行。真正不可逆的标志是消元到底后 $U$ 的对角线出现 0（即主元乘积 $=0$）。本题 $U$ 对角线为 $2,4,5$，乘积 $40$，与 $|\det A|$ 一致。

</details>

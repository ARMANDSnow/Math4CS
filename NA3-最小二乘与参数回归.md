---
course: Math4CS
chapter: NA3
topic: 设计与分析线性系统 / 最小二乘与参数回归 Least Squares & Regression
teacher: A · 应凯
date: 2026-06-19
tags: [Math4CS, 数值算法, 最小二乘, 回归, 正规方程, 正则化, Tikhonov, Cholesky]
---

# NA3 最小二乘与参数回归（Least Squares & Regression）

> 来源：A 老师 `3-DesignAnalyLinearSys.pptx`（教材 Solomon 第 3–4 章）· 考纲第 3 章
> 一句话定位：给一堆（带噪声的）数据，找最符合它们的模型参数。当数据多于参数（超定）无法精确拟合时，就退而求"误差平方和最小"——这就是最小二乘。
> **考纲地位：大题（参数估计 + 最小二乘）。** 旧真题 Q7(1) 图像去模糊 6 分直接考。Cholesky 分解"不考计算"但概念/关系会考。

**本章导航**
1. 这章解决什么问题
2. 从零铺垫：回归 = 用基函数拟合数据
3. 恰定情形：n 个实验 n 个参数（多项式回归 / Vandermonde）
4. 超定情形：最小二乘 min‖Ax−b‖²
5. 正规方程 Normal Equations（核心推导）
6. Gram 矩阵 AᵀA 的性质
7. 欠定情形：最小范数解（HW4.4）
8. 正则化 Tikhonov（旧真题 Q7(1)）
9. Cholesky 分解（概念）
10. 手算例题汇总
11. 考点雷达
12. 易错点 & 陷阱
13. 本章速查卡

---

## 1. 这章解决什么问题

科学实验/机器学习里最常见的任务：**我有一组数据点，想找一个函数去拟合它们。**

- 如果数据没噪声、且实验次数刚好等于参数个数 → 解一个普通线性方程组（NA2 的 LU 就行）。
- 但现实里数据**有噪声**，而且我们往往**做很多次实验**（数据点远多于参数）。这时想让函数**精确穿过每个点**既不可能也不明智（会过拟合噪声）。
- 折中办法：让函数**整体上离所有点最近**——具体说，让"残差平方和"最小。这就是**最小二乘（Least Squares）**。

---

## 2. 从零铺垫：回归 = 用基函数拟合数据

### 2.1 回归的基本设定

我们假设数据背后有个"黑箱"函数 $f$，输入 $x$ 输出 $y$。我们**猜**它的形式，比如线性 $f(x)=ax$、或多项式 $f(x)=a_0+a_1x+a_2x^2$。每次实验得到一对 $(x_i,y_i)$，我们要**反推参数**（$a$、$a_0,a_1,a_2\dots$）。

### 2.2 关键思想：把模型写成"基函数的线性组合"

PPT 反复强调的核心（"write f as a linear combination of basis functions"）：

$$f(x)=c_1\phi_1(x)+c_2\phi_2(x)+\cdots+c_n\phi_n(x)$$

- $\phi_j(x)$：**基函数**（你选的、已知的函数，如 $1,x,x^2,\dots$）。
- $c_j$：**待求参数**（未知数）。

妙处：虽然 $f$ 对 $x$ 可能是非线性的（如含 $x^2$），但它**对参数 $c_j$ 是线性的**！所以代入数据后会得到一个关于参数的**线性方程组**——又回到了我们会解的问题。

---

## 3. 恰定情形：n 个实验 n 个参数

做 $n$ 次实验、有 $n$ 个参数时，得到 $n$ 个方程 $n$ 个未知数，凑成方阵系统 $Ac=y$。

**多项式回归**的例子：用 $f(x)=a_0+a_1x+\cdots+a_{n-1}x^{n-1}$ 拟合 $n$ 个点，系数矩阵是著名的 **Vandermonde 矩阵**：

$$\underbrace{\begin{bmatrix}1&x_1&x_1^2&\cdots&x_1^{n-1}\\1&x_2&x_2^2&\cdots&x_2^{n-1}\\\vdots&&&&\vdots\\1&x_n&x_n^2&\cdots&x_n^{n-1}\end{bmatrix}}_{\text{Vandermonde}}\begin{bmatrix}a_0\\a_1\\\vdots\\a_{n-1}\end{bmatrix}=\begin{bmatrix}y_1\\y_2\\\vdots\\y_n\end{bmatrix}$$

每行把一个数据点 $x_i$ 代进各基函数。解这个方阵系统就得到系数。

> 问题（PPT 提出）：为什么非得做恰好 $n$ 次实验？如果 $y$ 有测量误差怎么办？精确穿过每个噪声点 = **过拟合**。这就引出超定 + 最小二乘。

---

## 4. 超定情形：最小二乘 min‖Ax−b‖²

当实验次数 $m$ **大于**参数个数 $n$（$m>n$，高瘦矩阵 $A\in\mathbb{R}^{m\times n}$），方程比未知数多，一般**无解**（没有哪组参数能让函数同时精确穿过所有点）。

退一步：找一组参数 $x$，让**预测值 $Ax$ 和观测值 $b$ 的差距最小**。用欧氏距离的平方衡量：

$$\boxed{\min_{x}\ \lVert Ax-b\rVert_2^2}$$

这就是**超定最小二乘问题**。$\lVert Ax-b\rVert_2^2=\sum_i(\text{第 }i\text{ 个预测}-\text{第 }i\text{ 个观测})^2$ 就是残差平方和。

直觉：找不到穿过所有点的线，那就找一条"总体上离所有点最近"的线。

---

## 5. 正规方程 Normal Equations（核心推导）

怎么解 $\min_x\lVert Ax-b\rVert_2^2$？把目标展开求导即可。

**第 1 步：展开二次型**

$$\lVert Ax-b\rVert_2^2=(Ax-b)^\top(Ax-b)=x^\top A^\top A x-2b^\top A x+b^\top b$$

（中间用到 $x^\top A^\top b=b^\top Ax$，因为它是标量、等于自身转置，所以两项合并成 $-2b^\top Ax$。）

**第 2 步：对 $x$ 求梯度并令其为 0**

$$\nabla_x=2A^\top A x-2A^\top b=0$$

**第 3 步：得到正规方程**

$$\boxed{A^\top A\,x=A^\top b}$$

若 $A$ 列满秩，$A^\top A$ 可逆，解为：

$$x=(A^\top A)^{-1}A^\top b$$

> 几何直觉：最小二乘解让残差向量 $Ax-b$ **垂直于** $A$ 的列空间（投影）。"正规（normal）"就是"正交/垂直"的意思——残差和列空间正交，正是 $A^\top(Ax-b)=0$ 这个条件。这条几何直觉在 [[NA5-列空间与QR]] 会用 QR 重新解释一遍。

---

## 6. Gram 矩阵 AᵀA 的性质

正规方程的核心矩阵 $A^\top A$ 叫 **Gram 矩阵**，两条关键性质：

1. **对称**：$(A^\top A)^\top=A^\top A$。
2. **半正定（PSD）**：对任意 $x$，$x^\top(A^\top A)x=(Ax)^\top(Ax)=\lVert Ax\rVert_2^2\ge0$。
   若 $A$ 列满秩，则 $Ax\ne0$（$x\ne0$ 时），所以 $A^\top A$ 还是**正定（PD）**、可逆。

> 这两条很重要：对称正定（SPD）矩阵能用更省的 **Cholesky 分解**（§9）来解正规方程，而不必走一般 LU。

---

## 7. 欠定情形：最小范数解（HW4.4）

反过来，若参数比方程多（$m<n$，矮胖矩阵），$Ax=b$ 有**无穷多解**。怎么挑唯一一个？常用办法：选**范数最小**的那个解（"最简洁"的解）。

**最小范数问题**：

$$\min_x\ \lVert x\rVert_2\quad\text{s.t. }Ax=b$$

**结论（HW4.4 要证）**：当 $A$ 行满秩时，最小范数解是

$$\boxed{x_s=A^\top(AA^\top)^{-1}b}$$

**证明**（HW4.4 完整作答）：

先验证 $x_s$ 可行：$Ax_s=AA^\top(AA^\top)^{-1}b=b$ ✓。

任取另一可行解 $x$（满足 $Ax=b$），把它写成 $x=(x-x_s)+x_s$，展开范数平方：

$$\lVert x\rVert_2^2=\lVert(x-x_s)+x_s\rVert_2^2=\lVert x-x_s\rVert_2^2+\lVert x_s\rVert_2^2+2(x-x_s)^\top x_s$$

关键是证**交叉项为 0**（即 $x-x_s$ 与 $x_s$ 正交）：

$$(x-x_s)^\top x_s=(x-x_s)^\top A^\top(AA^\top)^{-1}b=\big[A(x-x_s)\big]^\top(AA^\top)^{-1}b=\big[\underbrace{Ax}_{b}-\underbrace{Ax_s}_{b}\big]^\top(\cdots)=0$$

于是

$$\lVert x\rVert_2^2=\lVert x-x_s\rVert_2^2+\lVert x_s\rVert_2^2\ge\lVert x_s\rVert_2^2$$

等号当且仅当 $x=x_s$。所以 $x_s$ 是**唯一**的最小范数解。$\blacksquare$

> 记忆：超定（高瘦）→ 正规方程 $x=(A^\top A)^{-1}A^\top b$；欠定（矮胖）→ 最小范数 $x_s=A^\top(AA^\top)^{-1}b$。**两个公式互为"镜像"**：一个是 $A^\top A$，一个是 $AA^\top$。

---

## 8. 正则化 Tikhonov（旧真题 Q7(1)）

### 8.1 为什么要正则化

两种麻烦会让最小二乘"翻车"：

- **过拟合**：精确拟合噪声数据，模型在新数据上很差。
- **病态**：$A^\top A$ 接近奇异（条件数极大），解对噪声极度敏感。PPT 那个例子（输入 $1,1$ 和 $1,1.00001$ 几乎一样）就会让解爆炸。

解决办法：在目标里加一个"**惩罚解太大**"的项，逼解平滑/变小。这就是 **Tikhonov 正则化**。

### 8.2 Tikhonov 正则化的形式与解

$$\boxed{\min_x\ \lVert Ax-b\rVert_2^2+\lambda\lVert x\rVert_2^2}\qquad(\lambda>0\text{ 为正则化参数})$$

求梯度令其为 0：$2A^\top Ax-2A^\top b+2\lambda x=0$，得**正则化正规方程**：

$$\boxed{(A^\top A+\lambda I)\,x=A^\top b}$$

加的 $\lambda I$ 让矩阵一定可逆、条件数变好。$\lambda$ 越大，解越小越平滑（但偏差越大）；$\lambda$ 越小，越接近原始最小二乘。

> 常见正则化谱系（cheat sheet 可写）：Tikhonov / Ridge（$\lambda\lVert x\rVert_2^2$，让解平滑）、Lasso（$\beta\lVert x\rVert_1$，让解稀疏）、Elastic Net（两者都加）。本课重点是 Tikhonov。SVD 视角的正则化见 [[NA7-奇异值分解SVD]]。

---

## 9. Cholesky 分解（概念）

> 考纲：**Cholesky 不考计算**，但概念和"与 QR 的关系"会考（旧真题 Q7(2)）。所以这节重在理解，计算例题放在 §10 供参考。

### 9.1 是什么

当矩阵 $A$ **对称正定（SPD）** 时（正规方程的 $A^\top A$ 正是 SPD），它可以分解成

$$\boxed{A=LL^\top}$$

其中 $L$ 是下三角矩阵。这是 LU 在对称正定情形下的"特化版"：因为对称，只需存一个 $L$（而不是 $L$ 和 $U$ 两个），**省一半内存、省一半计算**。

### 9.2 计算公式（逐元素）

$$l_{jj}=\sqrt{a_{jj}-\sum_{k<j}l_{jk}^2},\qquad l_{ij}=\frac{1}{l_{jj}}\Big(a_{ij}-\sum_{k<j}l_{ik}l_{jk}\Big)\ (i>j)$$

直觉：对角元开方，非对角元"减掉已算部分再除以对角元"。

### 9.3 用处

解 SPD 系统 $Ax=b$（如正规方程）：分解 $A=LL^\top$ 后，前代 $Ly=b$ + 回代 $L^\top x=y$，和 LU 一样两步三角求解。它是数值上最稳、最省的解法之一。

---

## 10. 手算例题汇总

### 例题 1（旧真题 Q7(1)，图像去模糊，6 分）

> **题**：拍到一张模糊灰度图 $\vec x_0\in\mathbb{R}^p$（$p$=像素数），模糊过程是线性变换、记为矩阵 $G$，想恢复清晰图 $\vec x\in\mathbb{R}^p$。
> (i) 用最小二乘表达去模糊问题；(ii) 为提高稳定性加入 Tikhonov 正则化项。

**解答**

(i) 模糊模型是"清晰图经过 $G$ 得到模糊图"，即 $G\vec x\approx\vec x_0$。让模型输出和观测的模糊图差距最小：

$$\boxed{\min_{\vec x}\ \lVert G\vec x-\vec x_0\rVert_2^2}$$

(ii) 图像去模糊往往病态（$G$ 条件数大），直接解会放大噪声。加 Tikhonov 项：

$$\boxed{\min_{\vec x}\ \lVert G\vec x-\vec x_0\rVert_2^2+\lambda\lVert \vec x\rVert_2^2}$$

（对应解 $(G^\top G+\lambda I)\vec x=G^\top\vec x_0$。题目只要求写表达式即可。）

### 例题 2（HW4.4，最小范数解）

见 §7 完整证明：$x_s=A^\top(AA^\top)^{-1}b$，用"$x=(x-x_s)+x_s$ 展开 + 交叉项正交"证明它范数最小。

### 例题 3（HW4.3，Cholesky 分解，理解用——考纲说不考计算）

> **题**：对对称正定矩阵 $A$ 做 Cholesky 分解 $A=LL^\top$。
> $$A=\begin{bmatrix}4&-2&4\\-2&5&-4\\4&-4&14\end{bmatrix}$$

逐元素套 §9.2 公式：

- $l_{11}=\sqrt4=2$
- $l_{21}=a_{21}/l_{11}=-2/2=-1$，$\;l_{31}=a_{31}/l_{11}=4/2=2$
- $l_{22}=\sqrt{a_{22}-l_{21}^2}=\sqrt{5-1}=2$
- $l_{32}=(a_{32}-l_{31}l_{21})/l_{22}=(-4-(2)(-1))/2=-1$
- $l_{33}=\sqrt{a_{33}-l_{31}^2-l_{32}^2}=\sqrt{14-4-1}=3$

$$L=\begin{bmatrix}2&0&0\\-1&2&0\\2&-1&3\end{bmatrix}$$

验证 $LL^\top=A$（自行乘开即可，对角线 $4,5,14$、非对角 $-2,4,-4$ 都对得上）。

---

## 11. 考点雷达

- 考纲：**大题 = 参数估计 + 最小二乘**。要会：① 把拟合问题写成 $\min\lVert Ax-b\rVert^2$；② 推导正规方程 $A^\top Ax=A^\top b$；③ 加 Tikhonov 项 $(A^\top A+\lambda I)x=A^\top b$。
- 旧真题 **Q7(1)（6 分）= 图像去模糊**，正是"写最小二乘 + 加 Tikhonov"，几乎可套模板。
- **Cholesky 不考计算**，但要会说它是什么（SPD 的 $A=LL^\top$、省一半）、以及与 QR 的关系（在 [[NA5-列空间与QR]] 证）。
- 最小范数解 $x_s=A^\top(AA^\top)^{-1}b$（HW4.4）属于经典证明题，"展开 + 交叉项正交"的套路要会。
- 横向联系：正规方程的几何（残差⊥列空间）+ QR 解法见 [[NA5-列空间与QR]]；SVD/伪逆解最小二乘见 [[NA7-奇异值分解SVD]]；正定矩阵判定见 [[Ch4-优化]] 的 PSD 部分。

---

## 12. 易错点 & 陷阱

1. **正规方程方向**：是 $A^\top A x=A^\top b$，别漏了左乘的 $A^\top$；解是 $(A^\top A)^{-1}A^\top b$ 不是 $A^{-1}b$（$A$ 不是方阵，没有逆）。
2. **超定 vs 欠定公式混用**：高瘦超定用 $A^\top A$；矮胖欠定用 $AA^\top$。记反就全错。
3. **Tikhonov 加在哪**：惩罚项是 $+\lambda\lVert x\rVert_2^2$（惩罚解的大小），对应矩阵里 $+\lambda I$。别加到 $b$ 上。
4. **展开二次型漏因子**：$\lVert Ax-b\rVert^2=x^\top A^\top Ax-2b^\top Ax+b^\top b$，中间项系数是 **−2**（两项合并），别写成 −1。
5. **最小范数证明的关键步**：必须证 $(x-x_s)\perp x_s$（交叉项为 0），靠的是 $A(x-x_s)=b-b=0$。这一步是得分点，别跳。
6. **Cholesky 要求**：只有**对称正定**矩阵才有 Cholesky；对角元开方时根号内必须 $>0$，否则不是 SPD。

---

## 13. 本章速查卡

**回归**：模型 = 基函数线性组合 $f=\sum_j c_j\phi_j$，对参数线性 → 线性方程组。多项式回归 → Vandermonde 矩阵。

**超定最小二乘**：$\min_x\lVert Ax-b\rVert_2^2$
- 展开：$x^\top A^\top Ax-2b^\top Ax+b^\top b$
- **正规方程**：$A^\top Ax=A^\top b\;\Rightarrow\;x=(A^\top A)^{-1}A^\top b$
- 几何：残差 $Ax-b\perp$ 列空间

**欠定最小范数**：$\min\lVert x\rVert$ s.t. $Ax=b\;\Rightarrow\;x_s=A^\top(AA^\top)^{-1}b$

**Tikhonov 正则化**：$\min\lVert Ax-b\rVert_2^2+\lambda\lVert x\rVert_2^2\;\Rightarrow\;(A^\top A+\lambda I)x=A^\top b$

**Gram 矩阵 $A^\top A$**：对称、半正定（列满秩则正定可逆）。

**Cholesky**（SPD，不考计算但考概念）：$A=LL^\top$，$L$ 下三角，省一半内存；$l_{jj}=\sqrt{a_{jj}-\sum_{k<j}l_{jk}^2}$，$l_{ij}=(a_{ij}-\sum_{k<j}l_{ik}l_{jk})/l_{jj}$。

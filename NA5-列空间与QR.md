---
course: Math4CS
chapter: NA5
topic: 列空间与 QR 分解 Column Spaces and QR
teacher: 应凯
date: 2026-06-19
tags: [Math4CS, 数值算法, QR分解, 正交, Gram-Schmidt, Householder, 最小二乘]
---

# NA5 列空间与 QR 分解（Column Spaces and QR）

> 来源：应凯老师 `5-columnQR.pptx`（教材 Solomon 第 5 章）· 考纲第 5 章
> 一句话定位：正交矩阵是"最好解"的矩阵（转置即逆、保长度）。把一般矩阵 $A$ 正交化成 $A=QR$，就能又稳又准地解最小二乘——而且避开 [[NA3-最小二乘与参数回归]] 正规方程"条件数平方"的毛病。
> **考纲地位：Gram-Schmidt、Householder QR 都会考计算。** 必练手算。

**本章导航**
1. 这章解决什么问题
2. 从零铺垫：正交向量与正交矩阵
3. 向量投影
4. QR 分解是什么
5. Gram-Schmidt 正交化（会考计算）
6. Householder 反射（会考计算）
7. 用 QR 解最小二乘
8. Cholesky 与 QR 的关系（旧真题 Q7(2)）
9. 手算例题（HW4.6，两种方法）
10. 考点雷达
11. 易错点 & 陷阱
12. 本章速查卡
13. 练习题（自测）

---

## 1. 这章解决什么问题

[[NA3-最小二乘与参数回归]] 用正规方程 $A^\top Ax=A^\top b$ 解最小二乘，但 [[NA4-范数与条件数]] 指出 $A^\top A$ 把条件数**平方**了——数值上不稳。

本章给出更好的路子：把 $A$ 分解成 **正交矩阵 $Q$ × 上三角矩阵 $R$**（即 $A=QR$）。正交矩阵保长度、转置即逆，于是最小二乘问题化简成一个三角系统 $Rx=Q^\top b$，回代即可，**全程不碰 $A^\top A$**。

怎么把 $A$ 变正交？两种招数：**Gram-Schmidt**（几何上逐列去掉投影）和 **Householder 反射**（用镜面反射，更稳）。这两个都是考纲点名的计算考点。

---

## 2. 从零铺垫：正交向量与正交矩阵

### 2.1 正交、正交归一

- **正交（orthogonal）**：两向量垂直，内积为 0，$u^\top v=0$。
- **正交归一（orthonormal）**：一组向量两两正交、且每个长度为 1。

### 2.2 正交矩阵及其黄金性质

列向量都是正交归一的方阵 $Q$ 叫**正交矩阵**，核心性质：

$$\boxed{Q^\top Q=I\quad\Longleftrightarrow\quad Q^{-1}=Q^\top}$$

"求逆 = 转置"，几乎不花代价。还有**等距性（isometry）**：正交矩阵作用不改变长度和角度——

$$\|Qx\|_2=\|x\|_2,\qquad (Qx)^\top(Qy)=x^\top y$$

证明很简单：$\|Qx\|^2=(Qx)^\top(Qx)=x^\top Q^\top Q x=x^\top x=\|x\|^2$。

直觉：正交矩阵 = 旋转/反射，只"转动"不"拉伸"。这正是它在数值上稳的根源——**它不会放大误差**（条件数 $=1$）。

---

## 3. 向量投影

QR 的几何核心是"投影"。把向量 $a$ 投影到方向 $u$ 上：

$$\text{proj}_u(a)=\frac{a^\top u}{u^\top u}\,u$$

若 $u=q$ 已归一化（$\|q\|=1$），简化为 $\text{proj}_q(a)=(a^\top q)\,q$。

**正交分解**：$a$ 可拆成"沿 $q$ 的分量"+"垂直于 $q$ 的分量"：

$$a=\underbrace{(a^\top q)q}_{\text{平行}}+\underbrace{\big(a-(a^\top q)q\big)}_{\text{垂直}}$$

后面那个"垂直部分"就是 Gram-Schmidt 每一步要保留的东西。

---

## 4. QR 分解是什么

$$\boxed{A=QR}$$

- $Q$：列正交归一（$Q^\top Q=I$），与 $A$ 张成同一个列空间。
- $R$：上三角矩阵。

直觉：$Q$ 是把 $A$ 的列"摆正"（正交化）后的方向，$R$ 记录"怎么用这些正交方向线性组合回原来的列"。两种构造法（Gram-Schmidt / Householder）得到的 $Q,R$ 可能差正负号，但都满足 $A=QR$。

---

## 5. Gram-Schmidt 正交化（会考计算）

**思想**：逐列处理。每拿到一个新列，就**减掉它在前面所有已正交方向上的投影**，剩下的垂直部分再归一化。

**算法**：对 $A$ 的各列 $a_1,a_2,\dots,a_n$：

$$
\begin{aligned}
&q_1=\frac{a_1}{\|a_1\|}\\
&v_j=a_j-\sum_{i<j}(a_j^\top q_i)\,q_i,\qquad q_j=\frac{v_j}{\|v_j\|}\quad(j\ge2)
\end{aligned}
$$

而 $R$ 的元素就是这些投影系数：

$$R_{ij}=q_i^\top a_j\ (i<j),\qquad R_{jj}=\|v_j\|,\qquad R_{ij}=0\ (i>j)$$

**记忆**：$q$ = "去掉投影再归一化"；$R$ 上三角，对角填"归一化前的长度"，上方填"投影系数"。

> **Modified Gram-Schmidt（MGS）**：经典 GS 在浮点下会因舍入丢失正交性；MGS 把"一次性减所有投影"改成"算出一个 $q_i$ 就立刻从剩余所有列里减掉它的分量"，数值更稳。结果一样，只是减投影的顺序不同。考试知道"MGS 更稳"即可。

---

## 6. Householder 反射（会考计算）

**动机**：Gram-Schmidt 是"逐步搭正交基"；Householder 反过来——**用一系列镜面反射，把 $A$ 直接打成上三角 $R$**，反射矩阵乘起来就是 $Q$。它数值上比经典 GS 更稳。

### 6.1 Householder 反射矩阵

给定向量 $v$，反射矩阵：

$$\boxed{H=I-\frac{2\,vv^\top}{v^\top v}}$$

它把空间沿"垂直于 $v$ 的镜面"翻转。性质：**对称（$H^\top=H$）且正交（$H^\top H=I$）**，所以 $H^{-1}=H$（反射两次回到原地）。

### 6.2 怎么用它清零一列

想把某列向量 $x$ 反射成只有第一个分量非零的形式 $(\pm\|x\|,0,\dots,0)$，取

$$v=x\mp\|x\|\,e_1$$

（符号通常取与 $x_1$ 相反以避免相减相消）。这样 $Hx=(\pm\|x\|,0,\dots,0)$，一步把该列对角线下方全清零。

**整体流程**：第 1 个反射 $H_1$ 清零第 1 列对角线下方，第 2 个反射 $H_2$ 清零第 2 列……最终 $H_{n-1}\cdots H_1A=R$（上三角），于是

$$A=\underbrace{H_1H_2\cdots H_{n-1}}_{Q}\,R$$

（因为每个 $H$ 是自身的逆，挪到右边直接相乘。）

---

## 7. 用 QR 解最小二乘

这是 QR 的实战价值。要解 $\min_x\|Ax-b\|_2^2$，把 $A=QR$ 代入。对于列满秩的 reduced QR（$Q$ 列正交归一、$R$ 方阵可逆），正规方程化简为：

$$A^\top Ax=A^\top b\;\Rightarrow\;R^\top\underbrace{Q^\top Q}_{I}Rx=R^\top Q^\top b\;\Rightarrow\;\boxed{Rx=Q^\top b}$$

$R$ 是上三角，直接**回代**就能解出 $x$。

> 关键好处：全程没出现 $A^\top A$，条件数没被平方，比正规方程稳得多。这就是 PPT 反复强调"正交矩阵 preferred"的原因。

---

## 8. Cholesky 与 QR 的关系（旧真题 Q7(2)）

> 这道 6 分题把 [[NA3-最小二乘与参数回归]] 的 Cholesky 和本章的 QR 串了起来，是经典推导题。

设 $A=QR$（QR 分解），又设 $A^\top A=LL^\top$（Cholesky 分解，$L$ 下三角）。

**(i) 验证 $M\equiv A(L^\top)^{-1}$ 是正交矩阵**

$$M^\top M=\big(A(L^\top)^{-1}\big)^\top A(L^\top)^{-1}=(L^\top)^{-\top}\,A^\top A\,(L^\top)^{-1}=L^{-1}(LL^\top)(L^\top)^{-1}=I$$

（用了 $A^\top A=LL^\top$ 和 $(L^\top)^{-\top}=L^{-1}$。）所以 $M^\top M=I$，$M$ 正交。$\blacksquare$

**(ii) Cholesky 与 QR 的关系**

由 $A=QR$ 得

$$A^\top A=(QR)^\top(QR)=R^\top\underbrace{Q^\top Q}_{I}R=R^\top R$$

而 $A^\top A=LL^\top$，所以 $R^\top R=LL^\top$。$R$ 上三角 ⇒ $R^\top$ 下三角，对比得

$$\boxed{L=R^\top}\quad(\text{即 } R=L^\top,\ \text{差对角线正负号})$$

也就是说：**$A^\top A$ 的 Cholesky 下三角因子 $L$，正好是 $A$ 的 QR 分解中上三角 $R$ 的转置。** 而且 (i) 里的 $M=A(L^\top)^{-1}=A R^{-1}=Q$——$M$ 就是 QR 里的 $Q$。

---

## 9. 手算例题（HW4.6，两种方法）

> **题**：对 $A=\begin{bmatrix}1&1&1\\0&1&1\\0&1&0\end{bmatrix}$ 分别用 (1) Householder 反射 和 (2) Gram-Schmidt 做 QR 分解。
> 列向量：$a_1=(1,0,0)^\top,\ a_2=(1,1,1)^\top,\ a_3=(1,1,0)^\top$。

### 方法 1：Gram-Schmidt

**第 1 列**：$\|a_1\|=1$，$q_1=(1,0,0)^\top$。

**第 2 列**：$a_2^\top q_1=1$，

$$v_2=a_2-(a_2^\top q_1)q_1=(1,1,1)-(1,0,0)=(0,1,1),\quad \|v_2\|=\sqrt2,\quad q_2=\tfrac{1}{\sqrt2}(0,1,1)$$

**第 3 列**：$a_3^\top q_1=1$，$\;a_3^\top q_2=\tfrac{1}{\sqrt2}$，

$$v_3=a_3-(a_3^\top q_1)q_1-(a_3^\top q_2)q_2=(1,1,0)-(1,0,0)-\tfrac12(0,1,1)=(0,\tfrac12,-\tfrac12)$$

$$\|v_3\|=\tfrac{1}{\sqrt2},\quad q_3=\tfrac{1}{\sqrt2}(0,1,-1)$$

**组装**（$R_{ij}=q_i^\top a_j$）：

$$Q=\begin{bmatrix}1&0&0\\0&\tfrac{1}{\sqrt2}&\tfrac{1}{\sqrt2}\\0&\tfrac{1}{\sqrt2}&-\tfrac{1}{\sqrt2}\end{bmatrix},\qquad R=\begin{bmatrix}1&1&1\\0&\sqrt2&\tfrac{1}{\sqrt2}\\0&0&\tfrac{1}{\sqrt2}\end{bmatrix}$$

验证 $QR=A$ ✓（$\tfrac{1}{\sqrt2}\approx0.707,\ \sqrt2\approx1.414$）。

### 方法 2：Householder

**第 1 列** $a_1=(1,0,0)^\top$ 已经只有首分量非零，对角线下方本就是 0，所以 $H_1=I$（无需反射）。

**第 2 列**：取子向量 $x=(1,1)^\top$（第 2、3 行、第 2 列），要反射成 $(\sqrt2,0)^\top$。取 $v=x-\|x\|e_1=(1-\sqrt2,\,1)^\top$，算得（嵌进 3×3）

$$H_2=\begin{bmatrix}1&0&0\\0&\tfrac{1}{\sqrt2}&\tfrac{1}{\sqrt2}\\0&\tfrac{1}{\sqrt2}&-\tfrac{1}{\sqrt2}\end{bmatrix}$$

于是 $Q=H_1H_2=H_2$、$R=H_2A$，得到与 Gram-Schmidt **完全相同**的 $Q,R$（因为本题第一列特殊）。

> 提示：换一种符号约定（把 $x$ 反射到 $-\sqrt2$）会得到第 2、3 列、第 2、3 行整体变号的 $Q,R$——仍是合法 QR。考试时只要 $QR=A$、$Q$ 列正交归一、$R$ 上三角即给分。

**顺带验证 Cholesky-QR 关系**：$A^\top A=\begin{bmatrix}1&1&1\\1&3&2\\1&2&2\end{bmatrix}$，其 Cholesky $L=\begin{bmatrix}1&0&0\\1&\sqrt2&0\\1&\tfrac{1}{\sqrt2}&\tfrac{1}{\sqrt2}\end{bmatrix}=R^\top$，与 §8 结论一致。

---

## 10. 考点雷达

- 考纲：**Gram-Schmidt、Householder QR 都会考计算**——给一个 $2\times2$/$3\times3$ 矩阵要求做 QR，两种方法都可能点名。
- 必练：Gram-Schmidt 全流程（去投影→归一化→组 $R$）；Householder 反射矩阵 $H=I-2vv^\top/v^\top v$ 的构造与"清零一列"。
- 旧真题 **Q7(2)（6 分）**：验证 $M=A(L^\top)^{-1}$ 正交 + 推导 Cholesky 与 QR 关系（$L=R^\top$）。属于必拿分的推导题。
- 概念常考：正交矩阵性质（$Q^\top Q=I$、保长度）；为什么用 QR 解最小二乘（避开 $A^\top A$ 的条件数平方）。
- 横向联系：最小二乘见 [[NA3-最小二乘与参数回归]]；条件数平方见 [[NA4-范数与条件数]]；更强的正交分解 SVD 见 [[NA7-奇异值分解SVD]]。

---

## 11. 易错点 & 陷阱

1. **忘了归一化**：Gram-Schmidt 里 $q_j$ 一定要除以 $\|v_j\|$；只减投影不归一化得到的是正交但非正交归一，$R$ 会错。
2. **$R$ 元素填错**：$R_{ij}=q_i^\top a_j$（用**原列** $a_j$，不是 $v_j$），对角 $R_{jj}=\|v_j\|$。
3. **投影要减"全部"前序方向**：算 $q_3$ 时要同时减 $q_1$、$q_2$ 两个投影，别漏。
4. **Householder 符号**：$v=x\mp\|x\|e_1$ 的符号取与 $x_1$ 相反更稳；不同符号给出的 $Q,R$ 差正负号，都合法。
5. **$H$ 别忘了系数 2**：$H=I-\dfrac{2vv^\top}{v^\top v}$，漏掉 2 就不是反射（不正交）。
6. **Cholesky-QR 关系方向**：是 $L=R^\top$（Cholesky 下三角 = QR 上三角的转置），别写成 $L=R$。

---

## 12. 本章速查卡

**正交矩阵**：$Q^\top Q=I$，$Q^{-1}=Q^\top$，保长度 $\|Qx\|=\|x\|$，条件数 $=1$。

**投影**：$\text{proj}_u(a)=\dfrac{a^\top u}{u^\top u}u$；归一化时 $=(a^\top q)q$。

**QR 分解**：$A=QR$，$Q$ 列正交归一、$R$ 上三角。

**Gram-Schmidt**：$q_1=a_1/\|a_1\|$；$v_j=a_j-\sum_{i<j}(a_j^\top q_i)q_i$，$q_j=v_j/\|v_j\|$；$R_{ij}=q_i^\top a_j$，$R_{jj}=\|v_j\|$。

**Householder**：$H=I-\dfrac{2vv^\top}{v^\top v}$（对称正交）；取 $v=x\mp\|x\|e_1$ 把一列对角线下方清零；$Q=H_1H_2\cdots$，$R=H_{n-1}\cdots H_1A$。

**QR 解最小二乘**：$Rx=Q^\top b$（回代），避开 $A^\top A$ 病态。

**Cholesky-QR 关系**：$A^\top A=R^\top R=LL^\top\Rightarrow L=R^\top$；$M=A(L^\top)^{-1}=Q$ 正交。

---

## 13. 练习题（自测）

> 仿期末 / 作业风格，全新设定与数字，与上文例题、原作业均不重复。建议先独立做完，再展开核对答案。

### 练习 1 · Gram-Schmidt 正交化（仿 HW4.6，⭐⭐，约 7 分）

设
$$A=\begin{pmatrix}3&5&0\\0&0&5\\4&0&0\end{pmatrix},\qquad a_1=\begin{pmatrix}3\\0\\4\end{pmatrix},\ a_2=\begin{pmatrix}5\\0\\0\end{pmatrix},\ a_3=\begin{pmatrix}0\\5\\0\end{pmatrix}.$$
用 **Gram-Schmidt** 把 $a_1,a_2,a_3$ 正交归一，给出标准正交基 $q_1,q_2,q_3$ 与上三角矩阵 $R$，并写出 $A=QR$。

<details>
<summary>📖 参考答案</summary>

**第 1 列**：$\|a_1\|=\sqrt{3^2+4^2}=5$，
$$q_1=\frac{a_1}{\|a_1\|}=\frac{1}{5}(3,0,4)^\top.$$

**第 2 列**：投影系数 $a_2^\top q_1=\dfrac{(5)(3)+0+0}{5}=3$（即 $a_2^\top q_1=3$）。
$$v_2=a_2-(a_2^\top q_1)q_1=(5,0,0)^\top-3\cdot\tfrac15(3,0,4)^\top=\Big(5-\tfrac{9}{5},\,0,\,-\tfrac{12}{5}\Big)^\top=\Big(\tfrac{16}{5},0,-\tfrac{12}{5}\Big)^\top.$$
$\|v_2\|=\sqrt{\big(\tfrac{16}{5}\big)^2+\big(\tfrac{12}{5}\big)^2}=\sqrt{\tfrac{256+144}{25}}=\sqrt{16}=4$，
$$q_2=\frac{v_2}{4}=\frac{1}{5}(4,0,-3)^\top.$$

**第 3 列**：$a_3=(0,5,0)^\top$，而 $q_1,q_2$ 都在 $x\text{-}z$ 平面内（第二分量为 0），故
$$a_3^\top q_1=0,\qquad a_3^\top q_2=0,\qquad v_3=a_3=(0,5,0)^\top,\quad\|v_3\|=5,\quad q_3=(0,1,0)^\top.$$

**校验正交性**：$q_1^\top q_2=\tfrac{1}{25}\big[(3)(4)+(4)(-3)\big]=0$ ✓，$q_1^\top q_3=q_2^\top q_3=0$ ✓。

**组装 $R$**（$R_{ij}=q_i^\top a_j$，对角 $R_{jj}=\|v_j\|$）：
$$R_{12}=q_1^\top a_2=3,\quad R_{13}=q_1^\top a_3=0,\quad R_{23}=q_2^\top a_3=0,$$
$$Q=\begin{pmatrix}\tfrac35&\tfrac45&0\\[2pt]0&0&1\\[2pt]\tfrac45&-\tfrac35&0\end{pmatrix},\qquad R=\begin{pmatrix}5&3&0\\0&4&0\\0&0&5\end{pmatrix}.$$
逐列验证 $QR=A$：第 1 列 $5q_1=(3,0,4)^\top=a_1$；第 2 列 $3q_1+4q_2=(\tfrac95,0,\tfrac{12}{5})^\top+(\tfrac{16}{5},0,-\tfrac{12}{5})^\top=(5,0,0)^\top=a_2$；第 3 列 $5q_3=(0,5,0)^\top=a_3$。三列均吻合 ✓。

</details>

### 练习 2 · Householder 反射（仿 Q7(2ii)，⭐⭐，约 6 分）

给定向量 $x=\begin{pmatrix}1\\2\\2\end{pmatrix}$。构造 Householder 反射矩阵 $H=I-\dfrac{2vv^\top}{v^\top v}$，使 $Hx$ 只有第一个分量非零（即 $Hx=\pm\|x\|\,e_1$）。写出反射向量 $v$、矩阵 $H$ 与作用结果 $Hx$，并验证 $H$ 是正交矩阵。

<details>
<summary>📖 参考答案</summary>

**范数**：$\|x\|=\sqrt{1^2+2^2+2^2}=\sqrt9=3$。

**取符号**：$x_1=1>0$，为避免相减相消，取与 $x_1$ 相反的符号，即把 $x$ 反射到 $-\|x\|e_1=(-3,0,0)^\top$，于是
$$v=x-(-\|x\|)e_1=x+\|x\|e_1=(1+3,\,2,\,2)^\top=(4,2,2)^\top.$$
（注意 $v=2(2,1,1)^\top$，用 $v'=(2,1,1)^\top$ 得到的 $H$ 完全相同，此时 $v'^\top v'=6$。）

**系数**：$v^\top v=4^2+2^2+2^2=24$，故 $\dfrac{2}{v^\top v}=\dfrac{2}{24}=\dfrac{1}{12}$。

**外积**：
$$vv^\top=\begin{pmatrix}16&8&8\\8&4&4\\8&4&4\end{pmatrix},\qquad \frac{2vv^\top}{v^\top v}=\frac{1}{12}\begin{pmatrix}16&8&8\\8&4&4\\8&4&4\end{pmatrix}=\begin{pmatrix}\tfrac43&\tfrac23&\tfrac23\\[2pt]\tfrac23&\tfrac13&\tfrac13\\[2pt]\tfrac23&\tfrac13&\tfrac13\end{pmatrix}.$$

**反射矩阵**：
$$H=I-\frac{2vv^\top}{v^\top v}=\begin{pmatrix}1-\tfrac43&-\tfrac23&-\tfrac23\\[2pt]-\tfrac23&1-\tfrac13&-\tfrac13\\[2pt]-\tfrac23&-\tfrac13&1-\tfrac13\end{pmatrix}=\frac13\begin{pmatrix}-1&-2&-2\\-2&2&-1\\-2&-1&2\end{pmatrix}.$$

**作用结果**：
$$Hx=\frac13\begin{pmatrix}-1&-2&-2\\-2&2&-1\\-2&-1&2\end{pmatrix}\begin{pmatrix}1\\2\\2\end{pmatrix}=\frac13\begin{pmatrix}-1-4-4\\-2+4-2\\-2-2+4\end{pmatrix}=\frac13\begin{pmatrix}-9\\0\\0\end{pmatrix}=\begin{pmatrix}-3\\0\\0\end{pmatrix}=-\|x\|e_1\ \checkmark$$

**验证正交**：$H$ 对称（$H^\top=H$）。算一列范数即可，例如第 1 列 $\tfrac13(-1,-2,-2)^\top$，长度平方 $=\tfrac19(1+4+4)=1$；任两列点积为 0，如第 1、2 列 $\tfrac19\big[(-1)(-2)+(-2)(2)+(-2)(-1)\big]=\tfrac19(2-4+2)=0$。故 $H^\top H=I$，$H$ 正交。$\blacksquare$

</details>

### 练习 3 · QR 解最小二乘（仿 §7 + Q7(2)，⭐⭐⭐，约 10 分）

对超定方程组 $Ax=b$，其中
$$A=\begin{pmatrix}1&-2\\1&1\\1&1\\1&2\end{pmatrix},\qquad b=\begin{pmatrix}0\\4\\2\\4\end{pmatrix},$$
求最小二乘解 $\hat x=\arg\min_x\|Ax-b\|_2$。要求：先对 $A$ 做（reduced）Gram-Schmidt 得 $Q,R$，再用 $R\hat x=Q^\top b$ 回代求解，并给出残差 $r=b-A\hat x$ 及其范数 $\|r\|_2$。

<details>
<summary>📖 参考答案</summary>

记 $A$ 两列为 $a_1=(1,1,1,1)^\top$，$a_2=(-2,1,1,2)^\top$。

**第 1 列**：$\|a_1\|=\sqrt4=2$，$q_1=\tfrac12(1,1,1,1)^\top$。

**第 2 列**：$a_2^\top q_1=\tfrac12(-2+1+1+2)=\tfrac12\cdot2=1$，
$$v_2=a_2-(a_2^\top q_1)q_1=(-2,1,1,2)^\top-1\cdot\tfrac12(1,1,1,1)^\top=\Big(-\tfrac52,\tfrac12,\tfrac12,\tfrac32\Big)^\top.$$
$\|v_2\|=\sqrt{\tfrac{25+1+1+9}{4}}=\sqrt{\tfrac{36}{4}}=\sqrt9=3$，
$$q_2=\frac{v_2}{3}=\frac{1}{6}(-5,1,1,3)^\top.$$
（校验 $q_1^\top q_2=\tfrac{1}{12}(-5+1+1+3)=0$ ✓。）

**上三角 $R$**（$R_{11}=\|a_1\|=2,\ R_{12}=q_1^\top a_2=1,\ R_{22}=\|v_2\|=3$）：
$$R=\begin{pmatrix}2&1\\0&3\end{pmatrix}.$$

**右端 $Q^\top b$**：
$$q_1^\top b=\tfrac12(0+4+2+4)=\tfrac12\cdot10=5,\qquad q_2^\top b=\tfrac16\big[(-5)(0)+(1)(4)+(1)(2)+(3)(4)\big]=\tfrac16(0+4+2+12)=\tfrac{18}{6}=3.$$
故 $Q^\top b=(5,3)^\top$。

**回代解 $R\hat x=Q^\top b$**：
$$\begin{pmatrix}2&1\\0&3\end{pmatrix}\begin{pmatrix}\hat x_1\\\hat x_2\end{pmatrix}=\begin{pmatrix}5\\3\end{pmatrix}\Rightarrow 3\hat x_2=3\Rightarrow\hat x_2=1,\quad 2\hat x_1+1=5\Rightarrow \hat x_1=2.$$
$$\boxed{\hat x=(2,\,1)^\top}.$$

**残差**：$A\hat x=2a_1+1\cdot a_2=(2,2,2,2)^\top+(-2,1,1,2)^\top=(0,3,3,4)^\top$，
$$r=b-A\hat x=(0,4,2,4)^\top-(0,3,3,4)^\top=(0,1,-1,0)^\top,\qquad \|r\|_2=\sqrt{0+1+1+0}=\sqrt2.$$
校验最优性（残差应正交于列空间）：$r^\top a_1=0+1-1+0=0$，$r^\top a_2=0+1-1+0=0$ ✓，确为最小二乘解。

> 用正规方程 $A^\top A\hat x=A^\top b$ 复核：$A^\top A=\begin{pmatrix}4&2\\2&10\end{pmatrix}$，$A^\top b=(10,14)^\top$，解得 $\hat x=(2,1)^\top$，与 QR 法一致。

</details>

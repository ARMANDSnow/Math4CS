---
course: Math4CS
chapter: 4
topic: 优化 Optimization
date: 2026-05-09
tags: [Math4CS, 优化, 凸优化, 线性规划, 二次规划, 梯度下降, 牛顿法, 内点法]
---

# Chapter 4 优化（Optimization）—— 自学精讲版

> 来源：上海交大 *Math for CS* · Lecture 4（Xiaodong Gu）
> 整理思路：先跟着目录把全章扫一遍建立框架；然后逐节深读；最后把"应用案例"里的代码跑通。所有 PPT 出现过的关键概念、公式、例子，都在下面对应章节里展开了。
> 备注：本笔记主要用于自学，因此宁愿啰嗦也尽量不漏一步推导。

---

## 📑 目录索引

- [[#1. 什么是优化？|1 · 什么是优化？]]
- [[#2. 优化在计算机科学中无处不在|2 · 优化在计算机科学中无处不在]]
- [[#3. 本章关注的范围：Constrained · Continuous · Convex|3 · 本章关注的范围]]
- [[#4. 线性代数预备|4 · 线性代数预备]]
- [[#5. 正（半）定矩阵 PSD / PD|5 · 正（半）定矩阵]]
- [[#6. 线性规划 LP|6 · 线性规划 LP]]
- [[#7. 二次规划 QP|7 · 二次规划 QP]]
- [[#8. 二次约束二次规划 QCQP|8 · 二次约束二次规划 QCQP]]
- [[#9. 几何规划 GP|9 · 几何规划 GP]]
- [[#10. 一般优化问题的统一形式|10 · 一般优化问题统一形式]]
- [[#11. 无约束优化|11 · 无约束优化]]
- [[#12. 迭代法概览|12 · 迭代法概览]]
- [[#13. 梯度下降 Gradient Descent|13 · 梯度下降]]
- [[#14. 牛顿法 Newton's Method|14 · 牛顿法]]
- [[#15. 有约束优化与罚函数法|15 · 有约束优化与罚函数法]]
- [[#16. 指示函数与障碍思想|16 · 指示函数与障碍思想]]
- [[#17. 对数障碍函数|17 · 对数障碍函数]]
- [[#18. 中心路径|18 · 中心路径]]
- [[#19. 内点法|19 · 内点法]]
- [[#20. 应用案例 1：云计算资源分配|20 · 应用 1：云资源分配]]
- [[#21. 应用案例 2：图像分割|21 · 应用 2：图像分割]]
- [[#22. 应用案例 3：投资组合优化|22 · 应用 3：投资组合]]
- [[#23. 核心公式总表|23 · 核心公式总表]]
- [[#24. 全章学习主线|24 · 全章学习主线]]
- [[#25. 考点 / 作业最可能考的 10 件事|25 · 考点 / 作业 Top 10]]
- [[#26. 习题精解（HW3 最后两题）|26 · 习题精解（HW3 最后两题）]]

> 💡 在 Obsidian 中按 `Cmd+Opt+L` 调出右侧大纲面板，可以看到所有三级标题。

---

## 1. 什么是优化？

**Optimization，优化，本质上就是：**

> 在一堆可选方案里，找到"最好"的那个。

数学上通常写成：

$$\min_{x}\ f(x) \qquad \text{or} \qquad \max_{x}\ f(x)$$

$$\text{subject to } x \in \mathcal{X}$$

这里有 **三个核心对象**：

| 概念 | 含义 | 小白理解 |
|---|---|---|
| $x$ | 决策变量 | 你能控制的旋钮 |
| $f(x)$ | 目标函数 | 你想变大或变小的指标 |
| $\mathcal{X}$ | 可行域 / 约束集合 | 旋钮可以转到的合法范围 |

### 一个最朴素的例子：工厂造桌椅

设：

$$x_1 = \text{桌子数量}, \quad x_2 = \text{椅子数量}$$

利润（PPT 中用 \$200 / \$100，简化系数后）：

$$2x_1 + x_2$$

资源限制：

$$x_1 + 2x_2 \le 6 \quad (\text{金属})$$
$$3x_1 + x_2 \le 9 \quad (\text{木材})$$
$$x_1, x_2 \ge 0$$

那么优化问题就是：

$$\max_{x_1, x_2}\ 2x_1+x_2 \quad \text{s.t. } x_1+2x_2 \le 6,\ 3x_1+x_2 \le 9,\ x_1,x_2 \ge 0$$

> 一句话：**目标函数告诉我们"想要什么"，约束告诉我们"不能违反什么"。**

---

## 2. 优化在计算机科学中无处不在

PPT 开头一口气列了一串例子，证明优化绝不是一门"纯数学课"：

- **多点路径规划**：从多个地点中找最短路线（导航 App）
- **云-边-端资源调度**：服务器、边缘设备、终端怎么分配任务
- **操作系统调度**：CPU 先处理哪个任务
- **人工智能训练**：训练神经网络就是在最小化 loss
- **机器人轨迹优化**：怎么走最安全、最省能量
- **最小二乘回归**：找一条直线拟合数据
- **最短路径问题**：图里从 $s$ 到 $t$ 找最短路线

### 例：最小二乘回归（PPT 重点例子）

给数据：

$$(x_1, r_1), (x_2, r_2), \dots, (x_m, r_m),\quad x_i \in \mathbb{R}^n,\ r_i \in \mathbb{R}$$

希望找一个线性函数 $f(x) = w^\top x$，让预测值和真实值误差最小：

$$\min_{w}\ \lVert Xw - r\rVert_2^2$$

其中：
- $X$：输入数据矩阵（每行一个样本）
- $w$：要学习的参数向量
- $r$：真实标签
- $\lVert Xw - r\rVert_2^2$：平方误差

> ⭐ 这就是机器学习中最基本的思想：**找参数，让误差最小**。后面所有"训练"本质都是这个套路的变种。

### 例：最短路径（图优化）

给定有向图 $G=(V,E)$，每条边 $e$ 有长度 $c_e$。设：

$$x_e = \begin{cases}1 & \text{选中边 }e \\ 0 & \text{未选中}\end{cases}$$

最短路问题写成 LP（把整数松弛为 $x_e \ge 0$ 也成立）：

$$\min\ \sum_{e\in E} c_e x_e$$

$$\text{s.t.}\ \sum_{e\to v} x_e = \sum_{e\leftarrow v} x_e,\ \forall v \notin\{s,t\}$$
$$\sum_{e\leftarrow s} x_e = 1,\quad 0 \le x_e \le 1$$

> 直觉：把"流量 1"从 $s$ 注入，要求每个中间点"流入=流出"，最后一定从 $t$ 流出。让总长度最小。

---

## 3. 本章关注的范围：Constrained · Continuous · Convex

PPT 明确说，本章主要关注：

| 类型 | 含义 |
|---|---|
| **Constrained** | 有约束 |
| **Continuous** | 变量是连续值（可取小数） |
| **Convex** | 凸优化 |

### 3.1 什么是连续优化？

变量可以取小数，例如 $x = 2.37$。
对比 **整数规划（IP）**：$x \in \{0,1,2,\dots\}$，那是另一类问题（NP-hard）。

### 3.2 什么是凸优化？

**直观说，凸优化问题像一个"碗" 🥣**：
- 碗底就是最低点
- 没有乱七八糟的局部坑
- **只要找到一个局部最优点，它一定是全局最优点**

严格定义：函数 $f$ 是 **凸函数**，如果对任意 $x, y$ 与 $0 \le \theta \le 1$，都有：

$$f(\theta x + (1-\theta)y) \le \theta f(x) + (1-\theta) f(y)$$

> **小白理解**：图像上任取两个点连一条直线，**直线永远在函数曲线上方或重合**。

为什么本课只讲凸？因为：
- 凸优化通常 **有多项式时间算法**（LP 的单纯形/内点法、QP/QCQP 的内点法）
- 凸优化 **理论保证全局最优**

凸优化家族（嵌套关系）：

$$\text{LP} \subset \text{QP} \subset \text{SOCP} \subset \text{SDP} \subset \text{Cone} \subset \text{Convex}$$

整数规划 IP 不在这个谱系里——它是非凸的硬骨头。

---

## 4. 线性代数预备

优化大量使用向量、矩阵、范数、二次型，先把语言对齐。

### 4.1 向量 Vector

一个 $n$ 维向量是 $n$ 个数排成一列：

$$x = \begin{bmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{bmatrix}, \quad x^\top = (x_1, x_2, \dots, x_n)$$

`⊤` 表示 **transpose**（转置）。

**几何理解**：
- 空间中的一个点
- 或者从原点出发的一根有方向、有长度的箭头

例：$x = (3,4)^\top$ 就是平面上的点 $(3,4)$，也是从原点指向那里的箭头。

### 4.2 向量加法

对应位置相加：

$$z = x + y = (x_1+y_1, \dots, x_n+y_n)^\top$$

例：$\begin{bmatrix}1\\2\end{bmatrix} + \begin{bmatrix}3\\5\end{bmatrix} = \begin{bmatrix}4\\7\end{bmatrix}$

### 4.3 数乘

每个元素都乘那个数：

$$y = ax = (ax_1, \dots, ax_n)^\top$$

例：$3\begin{bmatrix}1\\2\end{bmatrix} = \begin{bmatrix}3\\6\end{bmatrix}$

### 4.4 内积 Inner Product

两个向量 $x, y \in \mathbb{R}^n$ 的内积：

$$x \cdot y = x^\top y = y^\top x = \sum_{i=1}^n x_i y_i$$

例：$x=(4,6,1)^\top,\ y=(1,3,1)^\top$：

$$x^\top y = 4\cdot 1 + 6\cdot 3 + 1\cdot 1 = 23$$

> 若 $x^\top y = 0$，称 $x$ 与 $y$ **正交**（orthogonal，几何上就是"垂直"）。

### 4.5 范数 Norm（向量长度）

| 名称 | 公式 | 直观 |
| --- | --- | --- |
| $\ell_2$（欧氏） | $\lVert x\rVert_2 = \sqrt{x_1^2+\dots+x_n^2}$ | 直线距离 |
| $\ell_1$ | $\lVert x\rVert_1 = \lvert x_1\rvert+\dots+\lvert x_n\rvert$ | 曼哈顿距离 |
| $\ell_p$ | $\lVert x\rVert_p = \left(\lvert x_1\rvert^p+\dots+\lvert x_n\rvert^p\right)^{1/p}$ | 一般化 |

例：$x = (3, -4)^\top$
- $\lVert x\rVert_2 = \sqrt{9+16} = 5$
- $\lVert x\rVert_1 = 3 + 4 = 7$

> **为什么后面会同时出现 $\ell_1$ 和 $\ell_2$？** 在机器学习里 $\ell_2$ 常作 loss（圆滑、好求导），$\ell_1$ 常作正则化（让解稀疏）。这是 ML 里的"L1 / L2 正则"。

### 4.6 矩阵 Matrix

$A \in \mathbb{R}^{m\times n}$ 表示 $A$ 有 $m$ 行 $n$ 列：

$$A = \begin{bmatrix} a_{11} & \cdots & a_{1n} \\ \vdots & \ddots & \vdots \\ a_{m1} & \cdots & a_{mn}\end{bmatrix}$$

$a_{ij}$：第 $i$ 行、第 $j$ 列的元素。

**矩阵也可以看成是"按列叠起来的向量组"**：

$$A = [\,a_1\ \cdots\ a_n\,]$$

或"按行叠起来"：每行是 $a^\top_j$。这个视角后面会反复出现。

### 4.7 矩阵加法

对应元素相加：$C = A+B,\ c_{ij} = a_{ij}+b_{ij}$

### 4.8 转置

行列对调。

向量：$x=\begin{bmatrix}x_1\\x_2\end{bmatrix} \Rightarrow x^\top = (x_1,\ x_2)$

矩阵：$A = \begin{bmatrix}1&2\\3&4\end{bmatrix} \Rightarrow A^\top = \begin{bmatrix}1&3\\2&4\end{bmatrix}$

### 4.9 矩阵-向量乘法（最关键）

$A \in \mathbb{R}^{m\times n},\ x \in \mathbb{R}^n$：

$$Ax = y \in \mathbb{R}^m, \quad y_i = \sum_{j=1}^{n} a_{ij} x_j$$

**两种理解视角**：
1. **行视角**：$y_i$ 是 $A$ 第 $i$ 行 $\cdot\ x$（一连串内积）
2. **列视角**：$Ax$ 是 $A$ 各列向量的**加权和**，权重就是 $x_j$
   $$Ax = x_1 a_1 + x_2 a_2 + \dots + x_n a_n$$

> 几何上，矩阵乘法对空间做：**拉伸 / 旋转 / 压缩 / 剪切**。例如 PPT 那张蒙娜丽莎被矩阵作用后变形的图，就是一次矩阵乘法。

### 4.10 二次型 $x^\top A x$

$x \in \mathbb{R}^n,\ A \in \mathbb{R}^{n\times n}$，则：

$$x^\top A x = \sum_{i,j} a_{ij} x_i x_j \in \mathbb{R}$$

**注意：结果是一个标量（数）**，不是向量。

例：$x=\begin{bmatrix}x\\y\end{bmatrix},\ A=I=\begin{bmatrix}1&0\\0&1\end{bmatrix}$

$$x^\top I x = (x\ y)\begin{bmatrix}1&0\\0&1\end{bmatrix}\begin{bmatrix}x\\y\end{bmatrix} = x^2 + y^2$$

> 二次型在 QP、最小二乘、牛顿法里反复出现。**记住它的形状：先写 $x$ 横躺，再写矩阵，再写 $x$ 竖立，最后得到一个数**。

---

## 5. 正（半）定矩阵 PSD / PD

⭐ **本章最关键的概念之一**——它是判断"我这问题是不是凸"的核心工具。

### 5.1 正半定 Positive Semi-Definite, PSD

对称矩阵 $A$ 若满足：

$$x^\top A x \ge 0,\ \forall x \in \mathbb{R}^n$$

则称 $A$ 是 **正半定**，记作 $A \succeq 0$。

### 5.2 正定 Positive Definite, PD

更强：

$$x^\top A x > 0,\ \forall x \ne 0$$

记作 $A \succ 0$。**正定 ⇒ 正半定**，反之不一定。

### 5.3 三个直观例子（配 PPT 图）

| 矩阵 $A$ | $z^\top A z$ | 类型 | 几何形状 |
|---|---|---|---|
| $\begin{bmatrix}1&0\\0&1\end{bmatrix}$ | $x^2+y^2$ | **正定** PD | 朝上的碗 🥣 |
| $\begin{bmatrix}1&0\\0&0\end{bmatrix}$ | $x^2$ | **半正定** PSD | 水槽（一个方向是平的）|
| $\begin{bmatrix}1&0\\0&-1\end{bmatrix}$ | $x^2-y^2$ | **不定** indefinite | 马鞍 🐎 |

### 5.4 为什么记 PSD？

> **PSD ⇔ 函数 $f(x)=x^\top A x$ 是凸函数 ⇔ 它的图像是个"碗" ⇔ 后面 QP/QCQP 能漂亮求解的根本原因。**

等价定义（同样有用）：
- $A$ 对称 + 所有特征值 $\lambda_i \ge 0$ ⇔ $A \succeq 0$
- $A$ 对称 + 所有特征值 $\lambda_i > 0$ ⇔ $A \succ 0$

---

## 6. 线性规划 LP

线性规划是本章第一个核心问题。

### 6.1 标准形式

**目标函数和约束都是线性的**：

$$\max_{x}\ c^\top x \quad \text{s.t. } Ax \le b,\ x \ge 0$$

或：

$$\min_{x}\ c^\top x \quad \text{s.t. } Ax \le b,\ x \ge 0$$

- $x$：变量
- $c^\top x$：线性目标
- $Ax \le b$：线性不等式约束
- $x \ge 0$：非负约束（按位 element-wise）

### 6.2 工厂生产例子（PPT 完整版）

设 $x_1$ = 桌子数量，$x_2$ = 椅子数量。简化系数后利润 $= 2x_1 + x_2$。

资源约束：

$$x_1 + 2x_2 \le 6\ (\text{金属})$$
$$3x_1 + x_2 \le 9\ (\text{木材})$$
$$x_1, x_2 \ge 0$$

矩阵形式：

$$c = \begin{bmatrix}2\\1\end{bmatrix},\ x=\begin{bmatrix}x_1\\x_2\end{bmatrix},\ A=\begin{bmatrix}1&2\\3&1\end{bmatrix},\ b=\begin{bmatrix}6\\9\end{bmatrix}$$

$$\max\ c^\top x \quad \text{s.t. } Ax \le b,\ x \ge 0$$

**最优点**位于两条资源约束线的交点：

$$\begin{cases}x_1+2x_2=6\\3x_1+x_2=9\end{cases} \Rightarrow x_1^\star=2.4,\ x_2^\star=1.8$$

最大化的目标值：$2(2.4)+1.8 = 6.6$（按真实美元 $200(2.4)+100(1.8) = 660$）。

> **小思考**：如果只造桌子（$x_2=0$），木材约束给 $x_1 \le 3$，利润 = $6$；只造椅子（$x_1=0$），金属约束给 $x_2 \le 3$，利润 = $3$。结果都不如混合的 $6.6$，说明组合方案更优。

### 6.3 一般制造问题

$n$ 种产品、$m$ 种原料：
- $x_j$：产品 $j$ 的产量
- $a_{ij}$：每单位产品 $j$ 用原料 $i$ 多少
- $b_i$：原料 $i$ 总量
- $c_j$：产品 $j$ 的单位利润

$$\max\ z = \sum_j c_j x_j$$
$$\text{s.t.}\ \sum_j a_{ij} x_j \le b_i\ (\forall i),\ x_j \ge 0$$

矩阵化即标准 LP：$\max c^\top x,\ \text{s.t. } Ax \le b,\ x \ge 0$。

### 6.4 LP 的规范形式（Canonical Form）

任何 LP 都可以转成：

$$\max\ c^\top x \quad \text{s.t. } a_i^\top x \le b_i,\ x_j \ge 0$$

**四条转换规则**：

1. **最小化 → 最大化**：$\min c^\top x \Leftrightarrow \max (-c^\top x)$
2. **$\ge$ → $\le$**：两边乘 $-1$，$a^\top x \ge b \Rightarrow -a^\top x \le -b$
3. **等式 → 两个不等式**：$a^\top x = b \Rightarrow a^\top x \le b$ 且 $-a^\top x \le -b$
4. **无符号变量 → 两个非负变量**：$x_j$ 无界 $\Rightarrow x_j = x_j^+ - x_j^-,\ x_j^+,x_j^- \ge 0$

### 6.5 LP 的几何理解（精彩！）

每个不等式 $a_i^\top x \le b_i$ 定义一个 **半空间**（half-space），方向由系数向量 $a_i$ 决定。

多个半空间相交 → **多面体**（polyhedron / 凸多边形）：
$$\{x : Ax \le b\}$$

目标函数 $c^\top x = z$ 是一族 **平行的等值线 / 等值面**。

> **优化过程**：沿着 $c$ 的方向把目标函数等值线"往外推"，直到最后一刻还能碰到可行域的那条线，就是最优值。

⭐ **重要结论**：LP 的最优解（如果存在且有界）**一定出现在多面体的某个角点（顶点）上**。

这就是 **单纯形法（Simplex Method）** 的核心思想：在角点之间跳跃，每次找一个让目标更好的角点，直到没法再改进。

> **为什么是角点？** 在多面体内部或边的中间，目标值 $c^\top x$ 一定有一个方向可以继续改进；只有在角点（多个约束同时收紧）时才"无路可走"。

---

## 7. 二次规划 QP

LP 的目标函数是线性的；QP **允许目标函数中出现二次项**。

### 7.1 QP 的形式

$$\min_x\ \tfrac12 x^\top P x + q^\top x + r \quad \text{s.t. } Ax \le b$$

**关键要求**：$P \succeq 0$（半正定）。

> **为什么？** 因为这样目标函数是凸的（一个"碗"），整个问题是 **凸 QP**，能用多项式算法（内点法）求解。如果 $P$ 不是 PSD，目标可能是马鞍，最优解可能跑到无穷远，问题就难了（非凸 QP 是 NP-hard）。

### 7.2 带约束最小二乘（PPT 例子）

数据 $(a_1,b_1),\dots,(a_m,b_m)$，希望拟合线性函数，且参数 $x_i$ 有上下界 $l_i \le x_i \le u_i$。

目标：

$$\min_x\ \lVert Ax - b\rVert_2^2$$

**展开二次型**（这一步要会做！）：

$$\lVert Ax - b\rVert_2^2 = (Ax-b)^\top(Ax-b)$$

$$= x^\top A^\top A x - x^\top A^\top b - b^\top A x + b^\top b$$

注意 $x^\top A^\top b$ 是标量，等于它的转置 $b^\top A x$，所以两项合并：

$$= x^\top A^\top A x - 2 b^\top A x + b^\top b$$

带上约束就是个标准 QP：

$$\min_x\ x^\top \underbrace{A^\top A}_{P} x - 2 b^\top A x + b^\top b$$
$$\text{s.t. } l_i \le x_i \le u_i$$

⭐ **为什么 $A^\top A$ 一定 PSD？**

$$x^\top (A^\top A) x = (Ax)^\top (Ax) = \lVert Ax\rVert_2^2 \ge 0$$

所以最小二乘永远是凸 QP，永远能解。

---

## 8. 二次约束二次规划 QCQP

**Q**uadratically **C**onstrained **Q**uadratic **P**rogramming——目标和约束都可以是二次。

### 8.1 标准形式

$$\min_x\ \tfrac12 x^\top P_0 x + q_0^\top x + r_0$$
$$\text{s.t. } \tfrac12 x^\top P_i x + q_i^\top x + r_i \le 0,\ i=1,\dots,m$$
$$Ax = b$$

若所有 $P_i \succeq 0$ → 凸 QCQP。可行域是若干 **椭球的交**（再交上仿射子空间）。

包含关系：**QP ⊂ QCQP**（QP 就是 $P_1=\dots=P_m=0$ 的特例）。

### 8.2 案例：Max Cut（最大割问题）

无向图 $G=(V,E)$，希望把顶点分成两组 $(S, V\setminus S)$，让 **跨越两边的边数最大**。

定义：

$$x_i = \begin{cases}+1 & v_i \in S\\ -1 & v_i \in V\setminus S\end{cases}$$

**关键观察**：
- 若边 $(i,j)$ 被割开（$v_i, v_j$ 不在同一边）：$x_i x_j = -1$，所以 $\dfrac{1-x_i x_j}{2} = 1$
- 若没被割开：$x_i x_j = 1$，$\dfrac{1-x_i x_j}{2} = 0$

所以 Max Cut 写成：

$$\max\ \sum_{(i,j) \in E} \frac{1-x_i x_j}{2} \quad \text{s.t. } x_i \in \{-1, +1\}$$

把 $x_i^2 = 1$ 代替 $x_i \in \{-1,+1\}$（等价的代数约束），用邻接矩阵 $G$ 写成：

$$\min\ x^\top G x \quad \text{s.t. } x_i^2 = 1$$

> ⚠️ 注意：虽然这看起来是 QCQP 的形状，但 $x_i^2 = 1$ 是 **非凸** 的等式约束（可行域只有 $2^n$ 个离散点），所以 Max Cut 实际上 **NP-hard**。课程提它是为了展示"建模技巧"，求解需要用 **SDP 松弛**（不在本章范围）。

---

## 9. 几何规划 GP

GP 是另一类特殊优化，常用于**工程设计**类问题（体积、面积、长宽比）。

### 9.1 单项式 Monomial

$$f(x) = c \cdot x_1^{a_1} x_2^{a_2} \cdots x_n^{a_n}$$

其中 $c \ge 0$，$a_i \in \mathbb{R}$（**指数可以是任意实数**，不必整数；这一点和高中"单项式"不一样）。

例：$f(x) = 3 x_1^2 x_2^{-1/2}$ 是 GP 里合法的单项式。

### 9.2 正多项式 Posynomial

**正多项式 = 若干单项式之和**。

例：$f(x) = 2 x_1 x_2 + 3 x_1^{-1} + 5 x_2^2$

### 9.3 GP 的标准形式

$$\min\ f_0(x) \quad \text{s.t. } f_i(x) \le b_i,\ h_j(x) = b_j',\ x \succ 0$$

其中 $f_i$ 是 posynomial，$h_j$ 是 monomial，$b_i,b_j' > 0$。

### 9.4 案例：手提箱设计（PPT 完整推导）

变量：高 $h$、宽 $w$、深 $d$。

- **目标**：最小化表面积（节省材料）：$\min\ 2(hw + hd + wd)$
- **体积约束**：$hwd \ge 5$，等价 $h^{-1}w^{-1}d^{-1} \le \tfrac15$
- **比例约束**：$h/w \le 2,\ h/d \le 3$
- **航空尺寸**：$h+w+d \le 7$

GP 形式：

$$\min\ 2hw + 2hd + 2wd$$
$$\text{s.t.}\ h^{-1}w^{-1}d^{-1} \le \tfrac15$$
$$\quad\ \, hw^{-1} \le 2,\ hd^{-1} \le 3$$
$$\quad\ \, h+w+d \le 7,\ h,w,d \ge 0$$

⭐ **关键技巧**：令 $\tilde h = \log h,\ \tilde w = \log w,\ \tilde d = \log d$，GP 经过对数变换变成 **凸问题**：

$$\min\ 2 e^{\tilde h+\tilde w} + 2 e^{\tilde h+\tilde d} + 2 e^{\tilde w+\tilde d}$$
$$\text{s.t.}\ e^{-\tilde h-\tilde w-\tilde d} \le \tfrac15,\ e^{\tilde h-\tilde w} \le 2,\ \dots$$

变换后的形式叫 **凸 GP**，可以用内点法求解。

---

## 10. 一般优化问题的统一形式

把前面所有种类合并成一个抽象框架：

$$\min_x\ f_0(x)$$
$$\text{s.t. } f_i(x) \le 0,\ i=1,\dots,m$$
$$\quad\ \ h_j(x) = 0,\ j=1,\dots,p$$

- $f_0$：目标函数
- $f_i$：不等式约束
- $h_j$：等式约束
- $x^\star$：最优解

**几何直观（PPT 那张图）**：
- $f_0$ 的等高线表示目标函数值
- $f_i, h_j$ 共同决定可行域 $\mathcal{X}$
- $-\nabla f_0$ 是目标函数下降最快的方向
- 最优点 $x^\star$ 通常出现在"等高线和可行域边界相切"的位置

---

## 11. 无约束优化

最简单的局面：

$$\min_x\ f(x)$$

PPT 假设 $f$ 是凸函数 + 二阶连续可微（一阶、二阶导数都存在且光滑）。

### 11.1 一维直觉

$f: \mathbb{R} \to \mathbb{R}$，最小点满足：

$$f'(x^\star) = 0$$

也就是"平坦点"（flat point）。对凸函数来说，**这就是全局最小**。

### 11.2 多维：梯度 Gradient

$f: \mathbb{R}^n \to \mathbb{R}$ 的梯度是一个向量，每一项是对相应变量的偏导：

$$\nabla_x f(x) = \begin{bmatrix} \dfrac{\partial f}{\partial x_1}(x) \\ \vdots \\ \dfrac{\partial f}{\partial x_n}(x) \end{bmatrix}$$

**直观含义**：梯度指向函数 **增长最快** 的方向；它的反方向 $-\nabla f$ 就是 **下降最快** 的方向。

⭐ **无约束凸优化的最优条件**：

$$\nabla f(x^\star) = 0$$

### 11.3 梯度计算例子（PPT 例子）

$$f(x) = 2 x_1^2 x_2 - x_1 x_3^3$$

逐项求偏导：

$$\frac{\partial f}{\partial x_1} = 4 x_1 x_2 - x_3^3$$
$$\frac{\partial f}{\partial x_2} = 2 x_1^2$$
$$\frac{\partial f}{\partial x_3} = -3 x_1 x_3^2$$

合起来：

$$\nabla f(x) = \begin{bmatrix}4 x_1 x_2 - x_3^3 \\ 2 x_1^2 \\ -3 x_1 x_3^2\end{bmatrix}$$

> 求梯度本质就是"对每个变量分别求偏导，再叠成列向量"。机械操作，多练几次就会。

---

## 12. 迭代法概览

有些问题可以直接解 $\nabla f(x^\star) = 0$ 拿到解析解（如标准最小二乘），但 **大多数复杂问题没有解析解**，要靠迭代。

迭代法产生一串点：

$$x^{(0)}, x^{(1)}, x^{(2)}, \dots$$

希望：

$$f(x^{(k)}) \to \text{OPT}$$

常见算法：
- **Gradient Descent** 梯度下降
- **Steepest Descent** 最速下降
- **Newton's Method** 牛顿法

---

## 13. 梯度下降 Gradient Descent

最基本的迭代算法。

**核心更新公式**：

$$\boxed{x_{t+1} = x_t - \eta_t\, \nabla f(x_t)}$$

- $x_t$：当前点
- $\nabla f(x_t)$：当前梯度
- $-\nabla f(x_t)$：下降方向
- $\eta_t$：**步长 / 学习率**（step size / learning rate）

> **小白比喻**：你在山上想下到谷底，每一步都摸出"哪边最陡往下"，然后朝那边迈一小步。$\eta$ 就是你迈的步子有多大。

### 13.1 步长 $\eta$ 太关键了

| 步长 | 效果 |
|---|---|
| 太小 | slow progress（爬太久） |
| 太大 | oscillations（在最低点两边来回弹） |
| 非常大 | instability（直接发散，函数值越来越大） |

> 所以神经网络训练里"调学习率"是个大学问，太小训练不动，太大爆炸。

### 13.2 梯度下降的麻烦：尺度不一致

PPT 例子：

$$f(x) = 0.01 x_1^2 + x_2^2 - 0.5 x_1 - x_2$$

$$\nabla f = (0.02 x_1 - 0.5,\ 2 x_2 - 1)$$

两个方向系数差 100 倍：$x_2$ 方向变化快，$x_1$ 方向变化慢。

结果：梯度下降在 $x_2$ 方向飞速振荡，在 $x_1$ 方向慢吞吞挪——**之字形震荡**，效率极差（PPT 里那只小青蛙就是这意思 🐸）。

> 这种"长条形碗"在工程里很常见（数据没归一化时尤其明显）。两条出路：**(a) 把数据归一化** 或 **(b) 用能自动适应曲率的二阶方法 → 牛顿法**。

---

## 14. 牛顿法 Newton's Method

牛顿法思想：**不只看一阶梯度，还看二阶曲率**。
- 梯度下降只知道"往哪边下"
- 牛顿法还知道"这个方向弯得有多厉害"，从而决定该跨多大步

### 14.1 二阶 Taylor 近似

在当前点 $x$ 附近，把 $f(x+v)$ 用二阶 Taylor 展开近似：

$$\hat f(x+v) \approx f(x) + \nabla f(x)^\top v + \tfrac12 v^\top \nabla^2 f(x)\, v$$

其中 $\nabla^2 f(x)$ 是 **Hessian 矩阵**——所有二阶偏导组成的方阵 $H_{ij} = \frac{\partial^2 f}{\partial x_i \partial x_j}$。

牛顿法选择 $v$ 让这个二阶近似最小。求 $\hat f$ 对 $v$ 的梯度并令其为 0：

$$\nabla \hat f(x+v) = \nabla f(x) + \nabla^2 f(x)\, v = 0$$

$$\Rightarrow v^\star = -[\nabla^2 f(x)]^{-1} \nabla f(x)$$

这就是 **牛顿方向**：

$$\Delta x_{\text{nt}} = -[\nabla^2 f(x)]^{-1} \nabla f(x)$$

**更新公式**：

$$\boxed{x_{t+1} = x_t - [\nabla^2 f(x_t)]^{-1} \nabla f(x_t)}$$

### 14.2 牛顿法为什么这么快？

考虑凸二次函数：

$$f(x) = \tfrac12 x^\top P x + q^\top x + r,\quad P \succ 0$$

- 梯度：$\nabla f(x) = Px + q$
- Hessian：$\nabla^2 f(x) = P$

从任意起点 $x_t$ 出发，牛顿一步：

$$x_{t+1} = x_t - P^{-1}(P x_t + q) = x_t - x_t - P^{-1} q = -P^{-1} q$$

而真正的最优解就是 $\nabla f(x^\star) = 0 \Rightarrow x^\star = -P^{-1} q$。

⭐ **结论：对任何凸二次函数，牛顿法一步收敛到最优！**

PPT 那个 $f = 0.01 x_1^2 + x_2^2 - 0.5 x_1 - x_2$ 的例子，梯度下降需要几十上百步还在震荡，牛顿法 1 步搞定。

### 14.3 代价

每步要算并求逆 Hessian。$n$ 维问题下：
- Hessian 大小 $n \times n$，需要 $O(n^2)$ 内存
- 求逆 $O(n^3)$ 计算

这在 $n = 10^9$ 的深度学习里行不通——所以 ML 实践中常用 **拟牛顿法（L-BFGS）** 或 **Adam** 这类近似方法，它们用一阶信息近似二阶。

---

## 15. 有约束优化与罚函数法

带约束的一般问题：

$$\min_x\ f_0(x)\ \text{s.t. } f_i(x) \le 0,\ Ax = b$$

PPT 提到的常见解法：

- **Penalty Function Method**（罚函数法）
  - Quadratic Penalty Method（二次罚函数）
  - **Interior Point Method**（内点法 / 对数障碍法）⭐
- **Augmented Lagrangian**（增广拉格朗日）

核心思路：**把约束塞进目标里**，转成无约束问题再用前面的算法（梯度/牛顿）解。

---

## 16. 指示函数与障碍思想

原问题：

$$\min\ f_0(x)\ \text{s.t. } f_i(x) \le 0$$

可以改写为：

$$\min\ f_0(x) + \sum_{i=1}^m I_-(f_i(x))$$

其中 **指示函数**：

$$I_-(u) = \begin{cases}0 & u \le 0\\ +\infty & u > 0\end{cases}$$

**意思**：满足约束就不罚；违反就给无穷大惩罚（强制不能违反）。

> ⚠️ 问题：$I_-$ 不可导（有跳跃），梯度法/牛顿法都用不了。

解决：**用一个光滑的函数近似它**——这就引出了对数障碍函数。

---

## 17. 对数障碍函数

PPT 给出近似：

$$\hat I_-(u) = -\mu \log(-u),\quad \text{dom} = \{u < 0\},\ \mu > 0$$

**直观理解**：
- 当 $u$ 远离 0（深入可行域内部，$u \ll 0$）：$\log(-u)$ 是个普通有限值，惩罚很小
- 当 $u \to 0^-$（迭代点逼近约束边界）：$\log(-u) \to -\infty$，所以 $-\mu \log(-u) \to +\infty$，惩罚爆炸

> 像在可行域的边界 **竖了一道无形的墙**，把迭代点关在内部。这就是"内点法"名字的来源。

### 17.1 原问题的障碍近似

原问题转化为：

$$\min\ f_0(x) + \mu\, \phi(x),\quad \phi(x) = -\sum_{i=1}^m \log(-f_i(x))$$
$$\text{s.t.}\ Ax = b\quad \text{(等式约束保留)}$$

不等式约束 **隐式** 地由 $\phi$ 强制——只要在内部，$\phi$ 有定义。

例：线性约束 $b_i - a_i^\top x \le 0$（即 $a_i^\top x \ge b_i$）的对数障碍：

$$\phi(x) = -\sum_i \log(a_i^\top x - b_i)$$

PPT 给的图（三角形可行域里的等高线像水波纹一样从中心向边界飙升）就是这个 $\phi$。

---

## 18. 中心路径

近似问题：

$$\min\ f_0(x) + \mu\, \phi(x),\ \text{s.t. } Ax = b$$

记其最优解为 $x^\star(\mu)$。

| $\mu$ | 效果 |
|---|---|
| 大 | 障碍占主导，解被推到可行域 **中心**（远离边界）—— 问题平滑、好解 |
| 小 | 障碍变弱，解逼近 **真正最优点 $x^\star$**（往往在边界角上）—— 但越来越病态、难解 |

不同 $\mu$ 对应不同的 $x^\star(\mu)$。把这些点串起来，得到一条曲线：

$$\text{central path} = \{x^\star(\mu) : \mu > 0\}$$

PPT 里那张三角形图，从中心 $x^\star(10)$ 沿着虚线一路滑到角点 $x^\star$，就是中心路径。

---

## 19. 内点法（Interior Point Method）

**核心思路**（PPT 完整算法）：

```
给定: 严格可行点 x，初始 μ，精度 ε，缩减系数 β < 1
重复:
  1. Centering 步: 用牛顿法求解
       x*(μ) = argmin_x  μ·f₀(x) + φ(x),   s.t.  Ax = b
  2. 更新:    x := x*(μ)
  3. 停止:    若 m·μ < ε，结束
  4. 收紧:    μ := β·μ      ← 让障碍变弱，解往边界靠拢
```

**直观流程**：
1. 起点放在可行域内部某处
2. 在当前 $\mu$ 下用牛顿法找一个"中心"
3. 把 $\mu$ 缩小一点（比如乘 0.1）
4. 再以上一步的解为起点继续 Newton
5. 重复 → 解沿着中心路径逼近真正的最优

> **和单纯形法的对比**：
> - 单纯形：在多面体的角点之间跳跃 → 路径走外边
> - 内点法：从可行域内部穿过去，沿中心路径接近最优 → 路径走内里
>
> 所以才叫 **interior** point method。

⭐ 这是现代 LP/QP/SDP 求解器（CVXPY、Mosek、SciPy）背后的主力算法。

---

## 20. 应用案例 1：云计算资源分配

PPT 的第一个完整案例，本质就是 LP。

### 20.1 题目

云服务商分配资源跑两类任务：
- **数据处理任务**：每个利润 \$200，需 1 单位 CPU + 3 单位内存
- **机器学习任务**：每个利润 \$100，需 2 单位 CPU + 1 单位内存
- 总资源：6,000 CPU + 9,000 内存

问：最大化利润，每类各跑多少？

### 20.2 建模

设 $x$ = 数据处理任务数，$y$ = 机器学习任务数：

$$\max\ 200 x + 100 y$$
$$\text{s.t.}\ x + 2y \le 6000\ (\text{CPU})$$
$$\quad\ \ 3x + y \le 9000\ (\text{内存})$$
$$\quad\ \ x, y \ge 0$$

### 20.3 解

两条约束交点：

$$\begin{cases}x + 2y = 6000\\ 3x + y = 9000\end{cases} \Rightarrow x = 2400,\ y = 1800$$

最大利润：$200(2400) + 100(1800) = 660{,}000$

### 20.4 Python 求解

```python
from scipy.optimize import linprog

# linprog 默认 minimize，所以系数取负号
c = [-200, -100]
A = [[1, 2],   # CPU 约束
     [3, 1]]   # 内存约束
b = [6000, 9000]
x_bounds = (0, None)

result = linprog(c, A_ub=A, b_ub=b,
                 bounds=[x_bounds, x_bounds],
                 method='highs')
print(f"x={result.x[0]:.2f}, y={result.x[1]:.2f}")
print(f"最大利润 = ${-result.fun:.2f}")
```

---

## 21. 应用案例 2：图像分割

### 21.1 题目

灰度图 $I$，每个像素 $I_{ij} \in [0,1]$（0=黑，1=白）。要把目标和背景分开：找一个 0/1 矩阵 $x$，$x_{ij}=1$ 表示该像素属于目标，$x_{ij}=0$ 表示背景。

### 21.2 建模

$$\min\ \sum_{i,j} x_{ij}(1 - 2 I_{ij}) \quad \text{s.t. } x_{ij} \in \{0,1\}$$

⭐ **为什么这个目标合理？** 看系数 $1 - 2 I_{ij}$：

| 像素亮度 | $1 - 2I$ | 优化器倾向 | 解释 |
|---|---|---|---|
| 亮 ($I \approx 1$) | $\approx -1$（负） | $x \to 1$（让 $x \cdot (-1)$ 更负）| 判为目标 |
| 暗 ($I \approx 0$) | $\approx +1$（正） | $x \to 0$（让 $x \cdot 1$ 不增加）| 判为背景 |

> 系数的符号自动告诉优化器"这个像素该归哪边"。这是建模的精彩之处——把"分类"问题翻译成数值优化。

### 21.3 松弛 Relaxation

$x_{ij} \in \{0,1\}$ 是离散约束，这是个 **整数规划**（NP-hard）。把它松弛成连续约束：

$$0 \le x_{ij} \le 1$$

得到的连续凸问题（实际上是 LP）：

$$\min\ \sum_{i,j} x_{ij}(1 - 2 I_{ij})\quad \text{s.t. } 0 \le x_{ij} \le 1$$

由于目标对每个 $x_{ij}$ 是线性的，最优解的每个 $x_{ij}$ 自动取到 0 或 1（取决于系数符号），这刚好对应分割结果。

### 21.4 Python 求解

```python
import numpy as np
import cvxpy as cp

def image_segmentation(image):
    I = image.flatten()
    n = len(I)
    x = cp.Variable(n)
    objective = cp.Minimize(cp.sum(cp.multiply(x, 1 - 2*I)))
    constraints = [x >= 0, x <= 1]
    cp.Problem(objective, constraints).solve()
    return x.value.reshape(image.shape)
```

---

## 22. 应用案例 3：投资组合优化

### 22.1 题目

银行有 A、B 两类理财产品：
- A：售价 1 万元 / 份，预期收益 2 万元
- B：售价 2 万元 / 份，预期收益 1 万元

某客户有 10 万元现金，分别买 $x_1, x_2$ 份（连续值）。客户做完调研，估算 **风险**（方差）为：

$$\text{Risk} = \tfrac12 (x_1^2 + x_2^2)\ \text{万元}$$

**综合收益 = 预期收益 − 风险**，最大化它。

### 22.2 建模

$$\max\ 2 x_1 + x_2 - \tfrac12 (x_1^2 + x_2^2)$$
$$\text{s.t. } x_1 + 2 x_2 \le 10,\ x_1, x_2 \ge 0$$

注意目标函数里 $-\tfrac12(x_1^2 + x_2^2)$ 是凹的（求最大），加上线性项依然凹 → 这是凸 QP。

### 22.3 用梯度上升 + 投影解

求梯度：

$$\nabla f = \begin{bmatrix}2 - x_1 \\ 1 - x_2\end{bmatrix}$$

更新（**因为是 max，所以朝梯度方向走**）：

$$x_1^{(k+1)} = x_1^{(k)} + \eta(2 - x_1^{(k)})$$
$$x_2^{(k+1)} = x_2^{(k)} + \eta(1 - x_2^{(k)})$$

每步更新后还要 **投影回可行域**（如果违反 $x_1 + 2 x_2 \le 10$ 或非负约束就拉回边界），迭代到收敛。

### 22.4 这个题的解析最优解（PPT 没说但值得算）

先求无约束最优：$\nabla f = 0 \Rightarrow x_1 = 2,\ x_2 = 1$

检查约束：$x_1 + 2 x_2 = 2 + 2 = 4 \le 10$ ✓ 满足

所以最优解就是 $x_1^\star = 2, x_2^\star = 1$，最大综合收益：

$$f^\star = 2(2) + 1 - \tfrac12(4 + 1) = 5 - 2.5 = 2.5\ \text{万元}$$

> ⚠️ 注意：这里约束是 $\le$ 而不是 $=$，所以 **最优解不一定花完全部 10 万**。如果题目改成"必须花完"（$x_1 + 2x_2 = 10$），最优解会换一个。

### 22.5 Python 框架

```python
import numpy as np

def f(x): return 2*x[0] + x[1] - 0.5*(x[0]**2 + x[1]**2)
def grad(x): return np.array([2 - x[0], 1 - x[1]])

def project(x):
    x = np.maximum(x, 0)                       # 非负
    if x[0] + 2*x[1] > 10:                     # 简单投影到约束边界
        # 这里只是示意，严格投影需要解 KKT
        x[0] = max(x[0], 0); x[1] = (10 - x[0]) / 2
    return x

x = np.array([0.0, 0.0])
eta = 0.1
for _ in range(1000):
    x_new = project(x + eta * grad(x))
    if np.linalg.norm(x_new - x) < 1e-6: break
    x = x_new
print(f"最优解 x = {x}, f = {f(x):.4f}")
```

---

## 23. 核心公式总表

### 优化基本形式

$$\min_x\ f(x)\ \text{s.t. } x \in \mathcal{X}$$

### 一般约束优化

$$\min_x\ f_0(x)\ \text{s.t. } f_i(x) \le 0\ (i=1..m),\ h_j(x) = 0\ (j=1..p)$$

### 内积

$$x^\top y = \sum_i x_i y_i$$

### 范数

$$\lVert x\rVert_2 = \sqrt{\textstyle\sum_i x_i^2},\quad \lVert x\rVert_1 = \sum_i \lvert x_i\rvert,\quad \lVert x\rVert_p = \left(\sum_i \lvert x_i\rvert^p\right)^{1/p}$$

### 矩阵-向量乘法

$$y = Ax,\quad y_i = \sum_{j=1}^n a_{ij} x_j$$

### PSD / PD

$$A \succeq 0 \Leftrightarrow x^\top A x \ge 0,\ \forall x$$
$$A \succ 0 \Leftrightarrow x^\top A x > 0,\ \forall x \ne 0$$

### 线性规划 LP

$$\max\ c^\top x\ \text{s.t. } Ax \le b,\ x \ge 0$$

### 二次规划 QP

$$\min\ \tfrac12 x^\top P x + q^\top x + r\ \text{s.t. } Ax \le b,\ P \succeq 0$$

### 最小二乘展开

$$\lVert Ax - b\rVert_2^2 = x^\top A^\top A x - 2 b^\top A x + b^\top b$$

### 梯度

$$\nabla f(x) = \begin{bmatrix}\partial f/\partial x_1\\ \vdots \\ \partial f/\partial x_n\end{bmatrix}$$

### 无约束凸优化最优条件

$$\nabla f(x^\star) = 0$$

### 梯度下降

$$x_{t+1} = x_t - \eta_t\, \nabla f(x_t)$$

### 牛顿法

$$x_{t+1} = x_t - [\nabla^2 f(x_t)]^{-1} \nabla f(x_t)$$

### 指示函数

$$I_-(u) = \begin{cases}0 & u \le 0 \\ +\infty & u > 0\end{cases}$$

### 对数障碍函数

$$I_-(u) \approx -\mu \log(-u),\quad u < 0$$

### 障碍法近似问题

$$\min\ f_0(x) - \mu \sum_{i=1}^m \log(-f_i(x))\ \text{s.t. } Ax = b$$

---

## 24. 全章学习主线

> **优化 = 在约束允许的范围内，让目标函数达到最好；为做到这一点，我们先用线性代数描述变量和函数，再根据目标和约束的形状把问题分成 LP/QP/QCQP/GP，最后用梯度下降、牛顿法、内点法等算法求解。**

```
线性代数 ⇒ 建模 ⇒ 问题类型 ⇒ 算法 ⇒ 应用
```

| 问题类型 | 目标 | 约束 | 关键工具 |
|---|---|---|---|
| **LP** | 线性 | 线性 | 单纯形 / 内点法 |
| **QP** | 凸二次 | 线性 | 内点法（最小二乘） |
| **QCQP** | 凸二次 | 凸二次 | 内点法（椭球交） |
| **GP** | posynomial | posynomial | 取 log → 凸 |

| 算法 | 用什么信息 | 速度 | 代价 |
|---|---|---|---|
| 梯度下降 | 一阶 $\nabla f$ | 慢，受尺度影响 | 便宜 |
| 牛顿法 | 一阶 + 二阶 Hessian | 通常一两步 | Hessian 求逆贵 |
| 内点法 | 障碍 + 牛顿子步 | 多项式 | 凸优化通用 |

---

## 25. 考点 / 作业最可能考的 10 件事

1. **会把实际问题建模成优化问题**（资源分配、生产计划、投资组合）
2. **会写 LP 的矩阵形式**：$\max c^\top x,\ \text{s.t. } Ax \le b,\ x \ge 0$
3. **理解 LP 几何意义**：半空间、多面体、**最优解在角点**
4. **会判断正定 / 正半定**：$x^\top A x \ge 0$ 或看特征值
5. **会把最小二乘展开成二次形式**：$\lVert Ax-b\rVert_2^2 = x^\top A^\top A x - 2 b^\top A x + b^\top b$
6. **会求梯度** $\nabla f(x)$（对每个变量分别求偏导）
7. **会写梯度下降更新**：$x_{t+1} = x_t - \eta_t \nabla f(x_t)$，并 **解释步长太大/太小** 的问题
8. **会写牛顿方向**：$\Delta x_{\text{nt}} = -[\nabla^2 f(x)]^{-1} \nabla f(x)$
9. **理解内点法 / 对数障碍函数** $-\mu\log(-u)$ 把约束塞进目标的思想
10. **会用 cvxpy / scipy.optimize 求解一个简单问题**（建议跑通 PPT 的三个案例）

> 这 10 个点基本覆盖了 ch4 的全部主干。把这份笔记当对照表，每个点能用自己的话讲清楚 + 推一遍公式 + 写一段代码，这章就稳了。

---

## 26. 习题精解（HW3 最后两题）

> 这两题分别考"建模 + 罚函数转换"和"牛顿法手算"，是本章最典型的两类考法。把这两题吃透，对应 [[#25. 考点 / 作业最可能考的 10 件事|考点 Top 10]] 里的第 6、7、8、9 条基本都覆盖了。

### 习题 3：相对熵最小化（4 pts）

#### 题目

两个概率分布 $x, y \in \mathbb{R}^n_{++}$ 之间的 **相对熵（KL 散度）** 定义为：

$$D(x \| y) = \sum_{k=1}^n x_k \log \frac{x_k}{y_k}$$

给定一个已知的概率分布 $y = (y_1, \dots, y_n)$，要找一个分布 $x = (x_1, \dots, x_n)$ 使得它与 $y$ 的相对熵最小，并满足约束。

#### (1) 写出优化问题（1 pt）

**关键观察**：
- $x$ 是概率分布 ⇒ 各分量之和为 1：$\sum_i x_i = 1$，即 $\mathbf{1}^\top x = 1$
- 目标函数中有 $\log x_i$ ⇒ 要求 $x_i > 0$，即 $x \in \mathbb{R}^n_{++}$（题目已给）

$$\boxed{\min_x\ \sum_{i=1}^n x_i \log \frac{x_i}{y_i} \quad \text{s.t. } \mathbf{1}^\top x = 1,\ x \in \mathbb{R}^n_{++}}$$

> 顺便提一下：这种"找最接近某个已知分布"的建模在机器学习里超常见——交叉熵 loss、变分推断、强化学习的 policy 训练，全是这个套路。

#### (2) 用罚函数转无约束（3 pts）

⭐ 题目特意说："变量在优化过程中 **不需要严格在可行域内**"——这意味着 **不能用对数障碍（内点法）**，因为内点法要求迭代点始终严格满足约束。

改用 **二次罚函数法（Quadratic Penalty Method）**：把等式约束塞进目标里，用一个二次项作为"惩罚"。

只有一个等式约束 $\mathbf{1}^\top x - 1 = 0$，加它的平方惩罚：

$$\boxed{\min_x\ \sum_{i=1}^n x_i \log \frac{x_i}{y_i} + \frac{\rho}{2}(\mathbf{1}^\top x - 1)^2}$$

其中 $\rho > 0$ 是 **罚参数**（penalty parameter）。

**直觉**：
- 若 $\mathbf{1}^\top x = 1$（满足约束）→ 罚项 = 0，不影响目标
- 若 $\mathbf{1}^\top x \ne 1$（违反约束）→ 罚项变正，目标变大，优化器被推回去
- $\rho$ 越大，约束被强制得越严；实践中常 **从小到大** 调（先让问题平滑好解，再逐步收紧）

**为什么用二次（平方），不用一次（绝对值）？**
- 平方处处光滑可导，可以直接用梯度法 / 牛顿法
- 一次（$|\mathbf{1}^\top x - 1|$）在零点不可导

**正项性 $x \in \mathbb{R}^n_{++}$ 怎么办？**
由目标函数中的 $\log x_i$ **隐式保证**：
- 当 $x_i \to 0^+$ 时，$\log(x_i/y_i) \to -\infty$，但 $x_i \log(x_i/y_i) \to 0$
- 当 $x_i \le 0$ 时，$\log x_i$ 无定义——优化器自动"不会去那里"

所以不必额外加一项罚正项性。

> **对比内点法**：内点法会写成 $\min D(x\|y) - \mu \sum_i \log x_i$，强制 $x_i > 0$ 严格保持。但题目允许 $x$ 暂时不可行，所以二次罚就够了。

---

### 习题 4：牛顿法一步收敛（3 pts）

#### 题目

无约束优化问题：

$$\min\ f(x_1, x_2) = x_1^2 + 2 x_2^2 + 2 x_1 + 3$$

给定初始点 $x_0 = (0, 1)^\top$。

#### (1) 证明 $d_0 = (-1, 0)^\top$ 是 $x_0$ 处的下降方向（1 pt）

**判据**（下降方向定义）：

$$d \text{ 是 } x \text{ 处的下降方向} \Longleftrightarrow \nabla f(x)^\top d < 0$$

> **几何意义**：$\nabla f$ 指向上升最快的方向；$d$ 和 $\nabla f$ 的夹角大于 90°（内积为负）就意味着 $d$ 至少有一个分量是沿着下坡的。

**步骤 1**：求梯度

$$\frac{\partial f}{\partial x_1} = 2 x_1 + 2, \quad \frac{\partial f}{\partial x_2} = 4 x_2$$

$$\nabla f(x) = \begin{bmatrix} 2 x_1 + 2 \\ 4 x_2 \end{bmatrix}$$

**步骤 2**：代入 $x_0 = (0, 1)$

$$\nabla f(x_0) = \begin{bmatrix} 2 \cdot 0 + 2 \\ 4 \cdot 1 \end{bmatrix} = \begin{bmatrix} 2 \\ 4 \end{bmatrix}$$

**步骤 3**：计算内积 $\nabla f(x_0)^\top d_0$

$$\nabla f(x_0)^\top d_0 = (2,\ 4) \begin{bmatrix} -1 \\ 0 \end{bmatrix} = 2 \cdot (-1) + 4 \cdot 0 = -2 < 0$$

∴ $d_0 = (-1, 0)^\top$ 是 $x_0$ 处的 **下降方向** ✅

#### (2) 用牛顿法更新 $x_0$（2 pts）

**牛顿更新公式**：

$$x_{k+1} = x_k - [\nabla^2 f(x_k)]^{-1} \nabla f(x_k)$$

**步骤 1**：求 Hessian $\nabla^2 f$

对梯度再求一次偏导：

$$\frac{\partial^2 f}{\partial x_1^2} = 2,\quad \frac{\partial^2 f}{\partial x_2^2} = 4,\quad \frac{\partial^2 f}{\partial x_1 \partial x_2} = 0$$

$$\nabla^2 f(x) = \begin{bmatrix} 2 & 0 \\ 0 & 4 \end{bmatrix}$$

> Hessian 与 $x$ 无关，是常数矩阵——因为 $f$ 是二次函数。
> 顺带：Hessian 是 PD（两个特征值 2, 4 都 > 0），所以 $f$ 是 **严格凸** 的，全局最小存在且唯一。

**步骤 2**：求 Hessian 的逆

对角矩阵的逆 = 对角元取倒数：

$$[\nabla^2 f(x_0)]^{-1} = \begin{bmatrix} 1/2 & 0 \\ 0 & 1/4 \end{bmatrix}$$

**步骤 3**：求牛顿方向 $p_0$

$$p_0 = -[\nabla^2 f(x_0)]^{-1}\, \nabla f(x_0) = -\begin{bmatrix} 1/2 & 0 \\ 0 & 1/4 \end{bmatrix} \begin{bmatrix} 2 \\ 4 \end{bmatrix} = -\begin{bmatrix} 1 \\ 1 \end{bmatrix} = \begin{bmatrix} -1 \\ -1 \end{bmatrix}$$

**步骤 4**：更新

$$x_1 = x_0 + p_0 = \begin{bmatrix} 0 \\ 1 \end{bmatrix} + \begin{bmatrix} -1 \\ -1 \end{bmatrix} = \begin{bmatrix} -1 \\ 0 \end{bmatrix}$$

$$\boxed{x_1 = (-1,\ 0)^\top}$$

#### 验证：这是不是全局最优？

令 $\nabla f(x) = 0$：

$$\begin{cases} 2 x_1 + 2 = 0 \\ 4 x_2 = 0 \end{cases} \Rightarrow x_1^\star = -1,\ x_2^\star = 0$$

正是 $x_1 = (-1, 0)^\top$！

最优值：$f(-1, 0) = (-1)^2 + 2(0)^2 + 2(-1) + 3 = 1 - 2 + 3 = 2$

⭐ **这印证了讲义 [[#14.2 牛顿法为什么这么快？|14.2]] 的结论**：对任何凸二次函数，从任意起点出发，牛顿法 **一步收敛到全局最优**。

#### 题目隐含的"小思考"

如果用 **梯度下降** 从 $x_0 = (0, 1)$ 出发要多久？

更新公式：$x_{t+1} = x_t - \eta \nabla f(x_t)$

代入：
$$x_{t+1} = \begin{bmatrix} x_1 - \eta(2 x_1 + 2) \\ x_2 - 4 \eta x_2 \end{bmatrix} = \begin{bmatrix} (1 - 2\eta) x_1 - 2 \eta \\ (1 - 4\eta) x_2 \end{bmatrix}$$

两个方向的收缩系数分别是 $1 - 2\eta$ 和 $1 - 4\eta$——尺度不一样！这就是 [[#13.2 梯度下降的麻烦：尺度不一致|梯度下降的尺度问题]] 在小例子上的体现。
- 若 $\eta = 0.5$，$x_2$ 方向系数 $1 - 2 = -1$，会震荡
- 若 $\eta = 0.25$，$x_2$ 方向一步到位，但 $x_1$ 方向系数 $1 - 0.5 = 0.5$，还要好几步

而牛顿法用 Hessian 的逆 **自动给每个方向用了正确的步长**（$\eta_1 = 1/2,\ \eta_2 = 1/4$），所以一步搞定。

---

## 📌 补充阅读 / 后续可深入方向

- [ ] 对偶理论 + KKT 条件（很多教材这里很重要，本课讲义未深入）
- [ ] 随机梯度下降 SGD 与动量法（深度学习训练）
- [ ] L-BFGS（拟牛顿法，不用算 Hessian）
- [ ] SDP 半定规划（Max-Cut 的松弛求解）
- [ ] ADMM / Proximal 方法（大规模优化 / 分布式）

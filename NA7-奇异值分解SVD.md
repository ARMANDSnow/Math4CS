---
course: Math4CS
chapter: NA7
topic: 奇异值分解 Singular Value Decomposition
teacher: A · 应凯
date: 2026-06-19
tags: [Math4CS, 数值算法, SVD, 奇异值, 伪逆, 低秩近似, EckartYoung, 正则化]
---

# NA7 奇异值分解（Singular Value Decomposition, SVD）

> 来源：A 老师 `7-svd.pptx`（教材 Solomon 第 7 章）· 考纲第 7 章
> 一句话定位：SVD 是"万能分解"——**任意**矩阵都能拆成"旋转 → 沿轴拉伸 → 旋转"。它统一了最小二乘、伪逆、低秩压缩、条件数、正则化，是整门数值课的集大成。
> **考纲地位：SVD 会考计算。** 旧真题 Q7(3) 6 分考"三矩阵的大小/性质 + 几何意义"。

**本章导航**
1. 这章解决什么问题
2. 从零铺垫：A 把单位球映成椭球
3. SVD 定义与三矩阵性质（旧真题 Q7(3)-i）
4. 几何意义：旋转-拉伸-旋转（旧真题 Q7(3)-ii）
5. 怎么算 SVD（会考计算）
6. 伪逆与解线性系统
7. Eckart-Young 低秩近似
8. 矩阵范数与 SVD
9. SVD 视角的 Tikhonov 正则化
10. 手算例题（HW5.5）
11. 考点雷达
12. 易错点 & 陷阱
13. 本章速查卡

---

## 1. 这章解决什么问题

前面几章各管一摊：LU 解方阵、QR 解最小二乘、特征值找本征方向。SVD 把它们**统一**了——它对**任意** $m\times n$ 矩阵都成立（不要求方阵、不要求对称），并一举回答：

- 矩阵的"拉伸结构"是什么？（奇异值）
- 怎么稳健解最小二乘 / 欠定系统？（伪逆）
- 怎么用更少的数据近似一个矩阵？（低秩近似、图像压缩）
- 条件数、正则化怎么从一个统一视角看？

---

## 2. 从零铺垫：A 把单位球映成椭球

理解 SVD 的关键图景（PPT 的 "Geometry of A"）：

> 取所有长度为 1 的向量（单位球面），用矩阵 $A$ 作用，得到的像是一个**椭球**。这个椭球的各**半轴长度**就是 $A$ 的**奇异值** $\sigma_1\ge\sigma_2\ge\cdots\ge0$，半轴**方向**就是左奇异向量 $u_i$。

直觉：矩阵作为变换，在某些正交方向上拉伸最猛（$\sigma_1$）、某些方向上压缩最狠（$\sigma_n$）。SVD 就是把这组"最佳拉伸方向"找出来。

**与特征值的桥**（PPT 的 Lemma）：$A^\top A$ 是对称半正定矩阵，由 [[NA6-特征值问题]] 谱定理可正交对角化，它的特征值正是 $\sigma_i^2$，特征向量正是右奇异向量 $v_i$。这就是 SVD 能算出来的根据。

---

## 3. SVD 定义与三矩阵性质（旧真题 Q7(3)-i）

$$\boxed{A=U\Sigma V^\top}\qquad A\in\mathbb{R}^{m\times n}$$

| 矩阵 | 大小 | 性质 |
|---|---|---|
| $U$ | $m\times m$ | **正交**（$U^\top U=I$），列 = **左奇异向量**，是列空间（输出空间）的正交归一基 |
| $\Sigma$ | $m\times n$ | **对角**（其余位置为 0），对角元 = **奇异值** $\sigma_1\ge\sigma_2\ge\cdots\ge\sigma_r>0=\cdots$ |
| $V$ | $n\times n$ | **正交**（$V^\top V=I$），列 = **右奇异向量**，是行空间（输入空间）的正交归一基 |

补充性质：$\sigma_i=\sqrt{\lambda_i(A^\top A)}\ge0$；非零奇异值的个数 = $\operatorname{rank}(A)=r$。

> 这正是旧真题 Q7(3)(i) 的标准答案：**写出三矩阵的大小（$m\times m$、$m\times n$、$n\times n$）和性质（$U,V$ 正交，$\Sigma$ 对角且奇异值降序非负）。**

---

## 4. 几何意义：旋转-拉伸-旋转（旧真题 Q7(3)-ii）

$A=U\Sigma V^\top$ 作用到向量 $x$ 上，从右往左读：

$$x\;\xrightarrow{\,V^\top\,}\;\text{① 旋转/反射（输入空间）}\;\xrightarrow{\,\Sigma\,}\;\text{② 沿坐标轴拉伸（各乘 }\sigma_i)\;\xrightarrow{\,U\,}\;\text{③ 旋转/反射（输出空间）}$$

| 矩阵 | 几何变换 |
|---|---|
| $V^\top$ | 正交变换：在输入空间**旋转/反射**（不变形） |
| $\Sigma$ | 沿各坐标轴**拉伸/压缩**（第 $i$ 轴乘 $\sigma_i$） |
| $U$ | 正交变换：在输出空间**旋转/反射**（不变形） |

**一句话**：任何线性变换 = 旋转 → 沿轴缩放 → 旋转。这就是 Q7(3)(ii) 的答案。

---

## 5. 怎么算 SVD（会考计算）

PPT 的 "Simple Strategy"，三步：

1. **算 $A^\top A$**（$n\times n$ 对称半正定），做特征分解。特征值 $\lambda_i$ 的平方根就是奇异值 $\sigma_i=\sqrt{\lambda_i}$（**降序**排），对应特征向量就是 $V$ 的列 $v_i$。
2. **得 $\Sigma$**：把 $\sigma_i$ 放对角线。
3. **算 $U$ 的列**：$u_i=\dfrac{1}{\sigma_i}Av_i$（对 $\sigma_i>0$）。若 $m$ 大、需要补满 $U$ 成 $m\times m$ 正交阵，再用正交补凑齐剩下的列。

> 口诀：**右边 $V$ 从 $A^\top A$ 的特征向量来，奇异值 = 特征值开方，左边 $U=Av/\sigma$。**

---

## 6. 伪逆与解线性系统

当 $A$ 不可逆（非方阵或奇异）时，用 **伪逆（Moore-Penrose pseudoinverse）** 统一处理：

$$\boxed{A^+=V\Sigma^+U^\top}$$

其中 $\Sigma^+$ 把每个非零奇异值取倒数 $1/\sigma_i$ 再转置。则 $x=A^+b$ 一举给出：

- 超定（高瘦）时 → **最小二乘解**（等价于 [[NA3-最小二乘与参数回归]] 的正规方程解，但更稳）。
- 欠定（矮胖）时 → **最小范数解**。

所以伪逆是"解一切线性系统"的统一工具，背后就是 SVD。

---

## 7. Eckart-Young 低秩近似

**外积形式**：把 SVD 拆成奇异值加权的秩 1 外积之和：

$$A=\sum_{i=1}^r\sigma_i\,u_iv_i^\top$$

**Eckart-Young 定理**：想用一个**秩 $k$** 矩阵最好地近似 $A$，最优解就是**只保留前 $k$ 个最大奇异值**那几项：

$$A_k=\sum_{i=1}^k\sigma_i\,u_iv_i^\top$$

它是所有秩 $\le k$ 矩阵里离 $A$ 最近的（2-范数和 Frobenius 范数下都最优），误差 $\|A-A_k\|_2=\sigma_{k+1}$。

> 应用：图像压缩、推荐系统、PPT 的 **Eigenfaces**（人脸用前几个奇异方向近似）。直觉：大奇异值 = 主要信息，小奇异值 = 细节/噪声，丢掉小的就压缩了。

---

## 8. 矩阵范数与 SVD

SVD 把 [[NA4-范数与条件数]] 的范数/条件数全都串起来：

$$\|A\|_2=\sigma_1\ (\text{最大奇异值}),\qquad \|A\|_F=\sqrt{\textstyle\sum_i\sigma_i^2},\qquad \kappa_2(A)=\frac{\sigma_{\max}}{\sigma_{\min}}$$

这也解释了 HW5.4：$\kappa_2(A^\top A)=\kappa_2(A)^2$，因为 $A^\top A$ 的奇异值是 $\sigma_i^2$。

---

## 9. SVD 视角的 Tikhonov 正则化

[[NA3-最小二乘与参数回归]] 的 Tikhonov 解 $(A^\top A+\lambda I)x=A^\top b$，用 SVD 展开成"**奇异值滤波**"：

$$x=\sum_i\frac{\sigma_i}{\sigma_i^2+\lambda}\,(u_i^\top b)\,v_i$$

看那个系数 $\dfrac{\sigma_i}{\sigma_i^2+\lambda}$：

- $\sigma_i$ 大（$\gg\sqrt\lambda$）：$\approx 1/\sigma_i$，和伪逆一样正常保留。
- $\sigma_i$ 小（$\ll\sqrt\lambda$）：$\approx \sigma_i/\lambda\to0$，被**压制**。

直觉：小奇异值方向最容易被噪声放大（病态来源），Tikhonov 正好把这些方向"调暗"，于是更稳。这是正则化最深刻的解释。

---

## 10. 手算例题（HW5.5）

### 例题 1（HW5.5(a)，3×3，会考计算）

> **题**：对 $A=\begin{bmatrix}0&0&1\\0&\sqrt2&0\\\sqrt3&0&0\end{bmatrix}$ 求 SVD $A=U\Sigma V^\top$。

**第 1 步：算 $A^\top A$**

$$A^\top A=\begin{bmatrix}3&0&0\\0&2&0\\0&0&1\end{bmatrix}$$

（已是对角阵，特征值一目了然。）

**第 2 步：奇异值与 $V$**

特征值 $3,2,1$，奇异值 $\sigma=\sqrt3,\sqrt2,1$（已降序）。特征向量是 $e_1,e_2,e_3$，所以

$$\Sigma=\begin{bmatrix}\sqrt3&0&0\\0&\sqrt2&0\\0&0&1\end{bmatrix},\qquad V=I=\begin{bmatrix}1&0&0\\0&1&0\\0&0&1\end{bmatrix}$$

**第 3 步：算 $U$（$u_i=Av_i/\sigma_i$）**

- $u_1=\frac{1}{\sqrt3}Ae_1=\frac{1}{\sqrt3}(0,0,\sqrt3)^\top=(0,0,1)^\top$
- $u_2=\frac{1}{\sqrt2}Ae_2=\frac{1}{\sqrt2}(0,\sqrt2,0)^\top=(0,1,0)^\top$
- $u_3=\frac{1}{1}Ae_3=(1,0,0)^\top$

$$U=\begin{bmatrix}0&0&1\\0&1&0\\1&0&0\end{bmatrix}$$

验证 $U\Sigma V^\top=U\Sigma=A$ ✓。（注：奇异向量的整体正负号可不同，只要 $U,V$ 正交、$\Sigma$ 降序、$U\Sigma V^\top=A$ 即正确。）

### 例题 2（HW5.5(b)，列向量）

> **题**：求 $A=\begin{bmatrix}-5\\3\end{bmatrix}$ 的 SVD。

$A$ 是 $2\times1$。$A^\top A=(-5)^2+3^2=34$，奇异值 $\sigma_1=\sqrt{34}$，$V=[1]$。

$$u_1=\frac{1}{\sqrt{34}}A=\frac{1}{\sqrt{34}}\begin{bmatrix}-5\\3\end{bmatrix},\quad u_2=\frac{1}{\sqrt{34}}\begin{bmatrix}3\\5\end{bmatrix}\ (\text{正交补})$$

$$U=\frac{1}{\sqrt{34}}\begin{bmatrix}-5&3\\3&5\end{bmatrix},\quad \Sigma=\begin{bmatrix}\sqrt{34}\\0\end{bmatrix},\quad V=[1]$$

验证 $U\Sigma V^\top=$ $U$ 第一列 $\times\sqrt{34}=(-5,3)^\top$ ✓。（$\sqrt{34}\approx5.831$。单列向量的 SVD 很简单：$\sigma=\|A\|$，$u_1=A/\|A\|$，$V=1$。）

### 例题 3（旧真题 Q7(3)，概念）

见 §3（三矩阵大小/性质）和 §4（几何：$V^\top$ 旋转 → $\Sigma$ 拉伸 → $U$ 旋转）。直接默写即可，6 分。

---

## 11. 考点雷达

- 考纲：**SVD 会考计算**。必练三步法（$A^\top A$ 特征分解 → $V,\sigma$ → $U=Av/\sigma$），尤其 $2\times2$/对角型/单列这种好算的。
- 旧真题 **Q7(3)（6 分）**：三矩阵大小性质 + 几何意义（旋转-拉伸-旋转），属背诵+理解题。
- 高频概念：奇异值 = $\sqrt{\lambda(A^\top A)}$；$\|A\|_2=\sigma_1$、$\kappa_2=\sigma_{\max}/\sigma_{\min}$；Eckart-Young 低秩近似；伪逆解最小二乘；Tikhonov 的奇异值滤波。
- 横向联系：基于谱定理 [[NA6-特征值问题]]；统一最小二乘 [[NA3-最小二乘与参数回归]] 与范数/条件数 [[NA4-范数与条件数]]；比 QR [[NA5-列空间与QR]] 更强的正交分解。

---

## 12. 易错点 & 陷阱

1. **奇异值是特征值开方**：$\sigma_i=\sqrt{\lambda_i(A^\top A)}$，别忘开方；且奇异值恒 $\ge0$。
2. **三矩阵大小**：$U$ 是 $m\times m$、$V$ 是 $n\times n$、$\Sigma$ 是 $m\times n$（和 $A$ 同形）。别把 $U,V$ 大小搞反。
3. **降序排列**：奇异值必须从大到小排，$V,U$ 的列也要对应同序。
4. **$U=Av/\sigma$ 别漏除 $\sigma$**：$u_i$ 要归一化，$Av_i$ 的长度正好是 $\sigma_i$。
5. **符号自由**：奇异向量可整体变号，只要 $U\Sigma V^\top=A$ 且正交性成立都对，别因为和标准答案差负号就以为错。
6. **$\|A\|_2$ vs $\|A\|_F$**：$\|A\|_2=\sigma_1$（最大奇异值），$\|A\|_F=\sqrt{\sum\sigma_i^2}$（全部平方和开方），不同。

---

## 13. 本章速查卡

**SVD**：$A=U\Sigma V^\top$，$U_{m\times m}$ 正交、$\Sigma_{m\times n}$ 对角（$\sigma_1\ge\cdots\ge0$）、$V_{n\times n}$ 正交。$\sigma_i=\sqrt{\lambda_i(A^\top A)}$，非零个数 = rank。

**几何**：$V^\top$ 旋转 → $\Sigma$ 沿轴拉伸 → $U$ 旋转（单位球 → 椭球，半轴 = 奇异值）。

**算法**：① $A^\top A$ 特征分解得 $V$、$\sigma=\sqrt\lambda$；② $\Sigma$；③ $u_i=Av_i/\sigma_i$。

**伪逆**：$A^+=V\Sigma^+U^\top$（非零 $\sigma$ 取倒数）；$x=A^+b$ = 最小二乘/最小范数解。

**Eckart-Young**：最佳秩 $k$ 近似 $A_k=\sum_{i\le k}\sigma_iu_iv_i^\top$，误差 $\|A-A_k\|_2=\sigma_{k+1}$。

**范数**：$\|A\|_2=\sigma_1$，$\|A\|_F=\sqrt{\sum\sigma_i^2}$，$\kappa_2=\sigma_{\max}/\sigma_{\min}$。

**Tikhonov（SVD 视角）**：$x=\sum_i\dfrac{\sigma_i}{\sigma_i^2+\lambda}(u_i^\top b)v_i$，小奇异值被压制 → 更稳。

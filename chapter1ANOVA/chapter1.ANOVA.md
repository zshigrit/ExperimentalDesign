# ANalysis Of VAriance (ANOVA)
## Theories

**The law of total variance for ANOVA**

```math
Var(y)=E(Var(y \mid x))+Var(E(y \mid x))
```

Note: it is not $y_{i,j}$ This is conditional on x=i; neither $\bar y..$

$$
Var(y)=\frac{\sum_i\sum_j (y_{ij}-\bar y_{..})^2}{an}
$$

$$
E(Var(y \mid x)) = \frac{\sum_i \left[\sum_j (y_{ij}-\bar y_{i.})^2/n\right]}{a} = \frac{\sum_i\sum_j (y_{ij}-\bar y_{i.})^2}{an}
$$

$$
Var(E(y \mid x)) = \frac{\sum_i (\bar y_{i.}-\bar y_{..})^2}{a}
$$

$$
\sum_i\sum_j (y_{ij}-\bar y_{..})^2 = \sum_i\sum_j (y_{ij}-\bar y_{i.})^2 + n\sum_i (\bar y_{i.}-\bar y_{..})^2
$$

$$
SS_{\text{total}}=SS_E+SS_{\text{treatment}}
$$


### Effect Model

$$
X\beta =
\begin{bmatrix}
1 & 1 & 0 & 0 \\
1 & 1 & 0 & 0 \\
1 & 0 & 1 & 0 \\
1 & 0 & 1 & 0 \\
1 & 0 & 0 & 1 \\
1 & 0 & 0 & 1
\end{bmatrix}
\begin{bmatrix}
\mu \\
\delta_1 \\
\delta_2 \\
\delta_3
\end{bmatrix}
$$

{Rank}(X) = 3 ( not full rank)


**Estimation**

$$
\hat{\beta} = (X^{T}X)^{-}X^{T}Y
$$

If $X^{T}X$ is **not invertible**, then:

- A generalized inverse is required.
- There is no unique solution.

### Projection matrix

### Estimable functions
The estimable functions are linear combinations of the rows of Xβ, e.g., ρ′Xβ.

### Test statistics


$$
F = \frac{SS_{\text{trt}}/(a-1)}{SS_E/\big(a(n-1)\big)}\sim\frac{\chi^2_{a-1}/(a-1)}{\chi^2_{a(n-1)}/\big(a(n-1)\big)}
$$
<br><br>

$$
\frac{\bar{Y}_{i\cdot} - \mu_i}{\sqrt{\sigma^2/n}}
\sim N(0,1),
\qquad
\frac{SSE}{\sigma^2}
\sim \chi^2_{a(n-1)}.
$$

Therefore,

$$
t=\frac{\displaystyle \frac{\bar{Y}_{i\cdot} - \mu_i}{\sqrt{\sigma^2/n}}}{\displaystyle \sqrt{\frac{SSE/\sigma^2}{a(n-1)}}}
=\frac{\bar{Y}_{i\cdot} - \mu_i}{\sqrt{MSE/n}}\sim t_{a(n-1)}.
$$

<br>

The **Wald** statistic is

$$
W =
(R\hat{\beta}-r)' 
\left[ R\,\widehat{\mathrm{Var}}(\hat{\beta})\,R' \right]^{-1}
(R\hat{\beta}-r).
$$

For the ANOVA model,

$$
\widehat{\mathrm{Var}}(\hat{\beta}) = \sigma^2 (X'X)^{-1}
$$

and we estimate

$$
\sigma^2
$$

by

$$
MSE = \frac{SSE}{a(n-1)}.
$$

Under the null hypothesis,

$$
W \sim \chi^2_q
$$
**The one-way ANOVA F-test is actually a Wald test for a set of linear restrictions on the treatment means.**


The ANOVA null hypothesis is

$$
H_0 : \mu_1 = \mu_2 = \cdots = \mu_a
$$

This can be written as linear restrictions

$$
R\beta = 0
$$

#### The Wald Test (Jessica's JSB https://www.dropbox.com/scl/fi/roxdm0m0cs25mxdp6ngjs/stat205_lecture2.pdf?rlkey=y2y6v3glwxt9onz53j31musw7&e=1)

$$
\beta =
\begin{pmatrix}
\beta_{1\,(p_1\times1)} \\
\beta_{2\,(p_2\times1)}
\end{pmatrix}
$$

$$
\widehat{\mathrm{Var}}(\hat{\beta}) =
\begin{pmatrix}
\widehat{\mathrm{Var}}(\hat{\beta}_1) & \widehat{\mathrm{Cov}}(\hat{\beta}_1,\hat{\beta}_2) \\
\widehat{\mathrm{Cov}}(\hat{\beta}_1,\hat{\beta}_2) & \widehat{\mathrm{Var}}(\hat{\beta}_2)
\end{pmatrix}
$$

#### Hypothesis testing

$$
H_0 : \beta_2 = 0_{p_2 \times 1}
$$

$$
H_1 : \beta_2 \ne 0_{p_2 \times 1}
$$

#### The Wald statistic

$$
W = \hat{\beta}_2^{T}\,
\widehat{\mathrm{Var}}(\hat{\beta}_2)^{-1}\,
\hat{\beta}_2
$$

Note that if \($p_2$ = 1\),

$$
W = T^2
$$

#### Distribution of the Wald statistic

Given the normality assumptions

$$
\hat{\beta}_2 \sim \mathcal{N}(\cdot,\cdot), \quad
e \sim \mathcal{N}(\cdot,\cdot)
$$

the Wald statistic follows the exact \(F\)-distribution

$$
\frac{W}{p_2} \overset{H_0}{\sim} F(p_2, n - p)
$$

#### Asymptotic result

Asymptotically, the normality assumptions can be relaxed with large sample size (Central Limit Theorem):

$$
W \overset{H_0}{\sim} \chi^2_{p_2}, \qquad n \to \infty
$$

### Unbalanced case
### Contrast
- sum of squares for a contrast


In one-way ANOVA with \(g\) groups and equal sample size \(n\), a **contrast** is defined as

$$
C = \sum_{i=1}^{g} c_i \bar{Y}_{i\cdot},
$$

where the coefficients satisfy

$$
\sum_{i=1}^{g} c_i = 0.
$$



The sum of squares associated with the contrast is

$$
F = \frac{\hat{C}^2}{\mathrm{Var}(\hat{C})}=\frac{SSc/1}{MSE}.
$$

Since

$$
\hat{C} = \sum_{i=1}^{g} c_i \bar{Y}_{i\cdot}
$$

and

$$
\mathrm{Var}(\hat{C}) = \sigma^2 \left(\frac{1}{n}\sum_{i=1}^{g} c_i^2\right),
$$

we obtain

$$
SS_C =
\frac{\left(\sum_{i=1}^{g} c_i \bar{Y}_{i\cdot}\right)^2}
{\left(\frac{1}{n}\sum_{i=1}^{g} c_i^2\right)}.
$$


#### Vector Form

Let

$$
\bar{\mathbf{Y}} =
\begin{pmatrix}
\bar{Y}_{1\cdot} \\
\bar{Y}_{2\cdot} \\
\vdots \\
\bar{Y}_{g\cdot}
\end{pmatrix},
\qquad
\mathbf{c} =
\begin{pmatrix}
c_1 \\
c_2 \\
\vdots \\
c_g
\end{pmatrix}.
$$

Then the contrast can be written as

$$
\hat{C} = \mathbf{c}^\top \bar{\mathbf{Y}}.
$$

The sum of squares becomes

$$
SS_C =
\frac{(\mathbf{c}^\top \bar{\mathbf{Y}})^2}
{\frac{1}{n}\mathbf{c}^\top \mathbf{c}}.
$$

Connection to the F Test

A single contrast corresponds to a **1-degree-of-freedom F test**:

$$
F = \frac{SS_C}{MSE}.
$$
This tests whether the linear combination of group means defined by the contrast differs from zero.
### Residual 
#### (Matrix)

$$
\mathbf{e} = \mathbf{y} - X\hat{\beta}
= \mathbf{y} - M\mathbf{y}
= (I - M)\mathbf{y}
$$

$$
E(\mathbf{e}) = (I - M)E(\mathbf{y})
= (I - M)(X\beta)
= \mathbf{0}
$$

$$
\mathrm{Var}(\mathbf{e})
= (I - M)\,\mathrm{Var}(\mathbf{y})\,(I - M)
= (I - M)(\sigma^2 I)(I - M)
$$

$$
= \sigma^2 (I - M)
= \sigma^2 \, \text{blkdiag}\!\left(I_n - \frac{1}{n}U_n\right)
$$

where

$$
U_n = \mathbf{1}_n \mathbf{1}_n^{T}
$$

$$
\mathbf{e} \sim MVN\!\left(\mathbf{0},\; \sigma^2 (I - M)\right)
$$

**Note**

$$
\lim_{n \to \infty} \mathrm{Var}(\mathbf{e}) = \sigma^2 I
$$

#### scalar
$$
\mathrm{Var}(e_{ij})
= \mathrm{Var}(Y_{ij} - \bar{Y}_{i\cdot})
= \mathrm{Var}(Y_{ij}) + \mathrm{Var}(\bar{Y}_{i\cdot}) - 2\,\mathrm{Cov}(Y_{ij}, \bar{Y}_{i\cdot})
$$

$$
= \sigma^2 + \frac{1}{n^2}\sum_{k=1}^{n} \mathrm{Var}(Y_{ik})
- 2\,\mathrm{Cov}\!\left(Y_{ij}, \frac{1}{n}\sum_{k=1}^{n} Y_{ik}\right)
$$

$$
= \sigma^2 + \frac{\sigma^2}{n}
- \frac{2}{n}\sum_{k=1}^{n} \mathrm{Cov}(Y_{ij}, Y_{ik})
$$

Note:

$$
\mathrm{Cov}(Y_{ij}, Y_{ik}) = 0 \quad \text{unless } k=j
$$

Thus

$$
= \sigma^2 + \frac{\sigma^2}{n}
- \frac{2}{n}\mathrm{Cov}(Y_{ij}, Y_{ij})
$$

$$
= \sigma^2 + \frac{\sigma^2}{n} - \frac{2}{n}\sigma^2
$$

$$
= \sigma^2\left(1 - \frac{1}{n}\right)
$$

Finally,

$$
\lim_{n \to \infty} \mathrm{Var}(e_{ij}) = \sigma^2
$$

$$
\mathrm{Cov}(e_{ij}, e_{k\ell})
= \mathrm{Cov}(Y_{ij}-\bar Y_{i\cdot},\, Y_{k\ell}-\bar Y_{k\cdot})
$$

$$
= \mathrm{Cov}(Y_{ij}, Y_{k\ell})
- \mathrm{Cov}(Y_{ij}, \bar Y_{k\cdot})
- \mathrm{Cov}(Y_{k\ell}, \bar Y_{i\cdot})
+ \mathrm{Cov}(\bar Y_{i\cdot}, \bar Y_{k\cdot})
$$

**(i) If $i \neq k$:**

$$
\mathrm{Cov}(e_{ij}, e_{k\ell}) = 0
$$

**(ii) If \(i = k\) and \(j = l\):**

$$
\mathrm{Cov}(e_{ij}, e_{ij})
= \mathrm{Var}(e_{ij})
= \sigma^2\left(1 - \frac{1}{n}\right)
$$

**(iii) If \(i = k\) and \($j \neq \ell$\):**

$$
\mathrm{Cov}(e_{ij}, e_{i\ell})
= 0
- \mathrm{Cov}\!\left(Y_{ij}, \frac{1}{n}Y_{ij}\right)
- \mathrm{Cov}\!\left(Y_{i\ell}, \frac{1}{n}Y_{i\ell}\right)
+ \frac{1}{n^2}\sum_{m=1}^{n} \mathrm{Cov}(Y_{im}, Y_{im})
$$

$$
= -\frac{1}{n}\sigma^2 - \frac{1}{n}\sigma^2 + \frac{1}{n}\sigma^2
$$

$$
= -\frac{\sigma^2}{n}
$$

**Note**

$$
\text{As } n \to \infty,\; e_{ij} \text{ and } e_{i\ell} \text{ become uncorrelated.}
$$
---
## Implementation
## Examples
## Practice
# ANalysis Of VAriance (ANOVA)
## Theories
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
F = \frac{SS_{\text{trt}}/(a-1)}{SS_E/\big(a(n-1)\big)}

\sim
\frac{\chi^2_{a-1}/(a-1)}{\chi^2_{a(n-1)}/\big(a(n-1)\big)}
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
---
## Implementation
## Examples
## Practice
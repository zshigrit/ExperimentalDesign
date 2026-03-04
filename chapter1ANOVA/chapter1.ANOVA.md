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

## Implementation
## Examples
## Practice
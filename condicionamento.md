# Condicionamento da matriz A e das matrizes de iteração

- M1 = Jacobi
- M2 = Gauss-Seidel
- M3 = Jacobi ordem 2
- M4 = Gauss-Seidel ordem 2

Condicionamento: $kappa(B) = \|B\| \cdot \|B^{-1}\|$

## A (36 x 36)

| Norma | Valor da norma | Condicionamento |
|---|---:|---:|
| Euclidiana | 3.582368e+00 | 6.452663e+00 |
| ponderada | 4.041459e+00 | 1.257914e+01 |
| Soma máxima de colunas | 5.333400e+00 | 1.365091e+01 |
| ponderada | 5.333400e+00 | 2.751369e+01 |
| Soma máxima de linhas | 5.333400e+00 | 1.499944e+01 |
| ponderada | 8.389128e+00 | 3.757540e+01 |
| Frobenius | 1.219194e+01 | 6.223314e+01 |
| ponderada | 1.340082e+01 | 8.190424e+01 |
| Máximo | 2.666700e+00 | 2.666700e+00 |
| ponderada | 2.666700e+00 | 2.666700e+00 |

> Versão ponderada: W = diag(|A[i][i]|), aplicada por W*B*W^-1.

## M1 - Jacobi (36 x 36)

| Norma | Valor da norma | Condicionamento |
|---|---:|---:|
| Euclidiana | 7.813955e-01 | inf |
| ponderada | 1.214328e+00 | inf |
| Soma máxima de colunas | 1.000000e+00 | inf |
| ponderada | 1.200100e+00 | inf |
| Soma máxima de linhas | 1.000000e+00 | inf |
| ponderada | 2.145884e+00 | inf |
| Frobenius | 1.445693e+00 | inf |
| ponderada | 2.537842e+00 | inf |
| Máximo | 1.750103e-01 | inf |
| ponderada | 4.667000e-01 | inf |

> Versão ponderada: W = diag(|A[i][i]|), aplicada por W*B*W^-1.

## M2 - Gauss-Seidel (36 x 36)

| Norma | Valor da norma | Condicionamento |
|---|---:|---:|
| Euclidiana | 1.292098e+00 | inf |
| ponderada | 1.803554e+00 | inf |
| Soma máxima de colunas | 1.790024e+00 | inf |
| ponderada | 2.385902e+00 | inf |
| Soma máxima de linhas | 2.097963e+00 | inf |
| ponderada | 3.097187e+00 | inf |
| Frobenius | 1.825918e+00 | inf |
| ponderada | 3.006851e+00 | inf |
| Máximo | 2.542313e-01 | inf |
| ponderada | 5.125323e-01 | inf |

> Versão ponderada: W = diag(|A[i][i]|), aplicada por W*B*W^-1.

## M3 - Jacobi ordem 2 (72 x 72)

| Norma | Valor da norma | Condicionamento |
|---|---:|---:|
| Euclidiana | 1.454685e+00 | inf |
| ponderada | 2.081846e+00 | inf |
| Soma máxima de colunas | 3.000000e+00 | inf |
| ponderada | 3.300150e+00 | inf |
| Soma máxima de linhas | 2.000000e+00 | inf |
| ponderada | 3.718825e+00 | inf |
| Frobenius | 7.050004e+00 | inf |
| ponderada | 7.713070e+00 | inf |
| Máximo | 1.000000e+00 | inf |
| ponderada | 1.000000e+00 | inf |

> Versão ponderada: W = diag(|A[i][i]|), aplicada por W*B*W^-1.

## M4 - Gauss-Seidel ordem 2 (72 x 72)

| Norma | Valor da norma | Condicionamento |
|---|---:|---:|
| Euclidiana | 1.833411e+00 | inf |
| ponderada | 2.694382e+00 | inf |
| Soma máxima de colunas | 3.982027e+00 | inf |
| ponderada | 5.078853e+00 | inf |
| Soma máxima de linhas | 3.443211e+00 | inf |
| ponderada | 4.993048e+00 | inf |
| Frobenius | 7.173385e+00 | inf |
| ponderada | 8.018641e+00 | inf |
| Máximo | 1.000000e+00 | inf |
| ponderada | 1.000000e+00 | inf |

> Versão ponderada: W = diag(|A[i][i]|), aplicada por W*B*W^-1.

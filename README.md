# TRFD: Trust-Region method based on Finite Differences
## Purpose

TRFD solves composite nonsmooth problems of the form

$$\min f(x)=h(F(x)), \quad$$ such that $\quad Mx \leq v,$ $\quad lb \leq x \leq ub$,

where $F:R^n \to R^m$, $h:R^m \to R$, $M$ is a matrix, $v$ is a vector, and $lb$ and $ub$ are scalars.

TRFD only relies on the function values provided for $F$. In the current implementation, the outer function $h$ can be the L1-norm, the Max-function or the Identity function (in the latter case, for smooth problems with $m=1$).

## Conditions of use

All publications using TRFD must cite:

D. Davar, G. N. Grapiglia: _TRFD: A derivative-free trust-region method based on finite differences for composite nonsmooth optimization_, arXiv preprint arXiv:2410.09165, 2024.

## How to use

To solve your own problem, you **must provide** the function Ffun that evaluates $F$.

You can **optionally** modify the file _TRFD.m_ in the section _Algorithmic parameters_, where you may modify the default values used by TRFD. Moreover, you can **optionally** modify the file _inner_solver.m_ by setting the variable _options_2_, which chooses the inner-solver used by TRFD when the p-norm is the **2-norm**.

## The command line 
 
You can run TRFD by typing:

[x_min, f_min, nf, stop, H] = TRFD (x0, Ffun, nfmax, M, v, lb, ub, instance)

**Input**

- **x0** [n x 1]: initial point
- **Ffun**: function provided by the user that computes $F$
- **nfmax** [1 x 1]: maximum number of function evaluations allowed
- **M** [... x n]: matrix for linear constraints on $x$ ([ ] if no matrix)
- **v** [... x 1]: vector for linear constraints on $x$ ([ ] if no vector)
- **lb** [1 x 1]: lower bound on $x$ ([ ] if no lower bound)
- **ub** [1 x 1]: upper bound on $x$ ([ ] if no upper bound)
- **instance** [1 x 1]: defines the function $h$, and the p-norm used by TRFD to solve the trust-region subproblems. The table below provides, for each value of **instance**, the corresponding $h$ and p-norm.

| instance | $h$ | p-norm used by TRFD|
|:----------:|:----------:|:----------:|
| **1** | L1-norm: ‖F(x)‖₁ | 1-norm |
| **2** | L1-norm: ‖F(x)‖₁ | 2-norm |
| **3** | L1-norm: ‖F(x)‖₁ | Inf-norm |
| **4** | Max-function: max(F(x)) | 1- or Inf-norm |
| **5** | Max-function: max(F(x)) | 2-norm |
| **6** | Identity function: F(x) | 1-norm |
| **7** | Identity function: F(x) | 2-norm |
| **8** | Identity function: F(x) | Inf-norm |

**Output**

- **x_min** [n x 1]: vector yielding the lowest function value found within nfmax evaluations
- **f_min** [1 x 1]: lowest function value found with nfmax evaluations
- **nf** [1 x 1]: number of function evaluations used
- **stop** [1 x 1]: an integer identifying the reason TRFD stopped:
  - stop = 1: number of function evaluations bigger than or equal to nfmax
  - stop = 2: trust-region radius less than the tolerance
  - stop = 3: approximate stationarity measure less than the tolerance
  - stop = -1: error in the execution of linprog or fmincon
- **H** [nfmax x 1]: a vector such that H(i) is the smallest value of $h(F(x))$ obtained by TRFD after i function evaluations

## Examples

As a running-example for **composite problems**, consider a basic $F$ of 2 components and 2 variables, defined by the function Ffun as

Ffun = @(x) [x(1)^(2) + (x(2)-1)^(2) + x(2) - 1; -x(1)^(2) - (x(2)-1)^(2) + x(2) + 1]

- For solving the **composite** problem min ‖F(x)‖₁ with x0 = [1; 1], nfmax = 300, without constraint and by using the 1-norm, type:

  - [x_min, f_min, nf, stop, H] = TRFD ([1; 1], Ffun, 300, [ ], [ ], [ ], [ ], 1)

- For solving the **composite** problem min max(F(x)) with x0 = [1; 1], nfmax = 300, lb = 0.1, ub = 20 and by using the 2-norm, type:

  - [x_min, f_min, nf, stop, H] = TRFD ([1; 1], Ffun, 300, [ ], [ ], 0.1, 20, 5)

---

As a running-example for **smooth problems**, consider a basic $F$ of 1 component and 2 variables, defined by the function Ffun as

Ffun = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2

- For solving the **smooth** problem min F(x) with x0 = [0.5; 0.5], nfmax = 300, M = ones(1, 2), v = 10, lb = 0.1, ub = 20 and by using the Inf-norm, type:

  - [x_min, f_min, nf, stop, H] = TRFD ([0.5; 0.5], Ffun, 300, ones(1, 2), 10, 0.1, 20, 8)

## Contact

Dânâ Davar,  Geovani Nunes Grapiglia <br>
UCLouvain, Belgium <br>
ICTEAM Institute <br>
dana.davar@uclouvain.be, geovani.grapiglia@uclouvain.be <br>

May 2025

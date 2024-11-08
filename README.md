# TRFD: Trust-Region method based on Finite Differences
## Purpose

TRFD solves composite nonsmooth problems of the form

$$\min f(x)=h(F(x)),$$ $F:R^n \to R^m$, $h:R^m \to R,$

using only the function values provided for $F$. In the current implementation, the outer function $h$ can be either the L1-norm or the Max-function.

## Conditions of use

All publications describing a work using this software must cite:

Davar D., Grapiglia G. N., "TRFD: A derivative-free trust-region method based on finite differences for composite nonsmooth optimization", arXiv preprint arXiv:2410.09165, 2024.

## How to use

To solve your own problem, you MUST SPECIFY how to evaluate $F$ in the file Ffun.m.

You can OPTIONALLY modify the file TRFD_composite.m. In the section "ALGORITHMIC PARAMETERS", you may modify the default values for the parameters to be used by TRFD.

## The command line 
 
After you specified how to evaluate $F$, you can run TRFD by typing:

[x, f_min, nf, stop, H] = TRFD_composite(x0, m, @Ffun, h, nfmax)

where x0 is the initial point, m is the number of components of $F$, @Ffun is the function that evaluates $F$, h is a number defining the outer function (h = 1 for the L1-norm, h = 2 for the Max-function), and nfmax is the maximum number of function evaluations allowed.

As a running-example, the file Ffun.m already provides a basic $F$ of 2 components and 2 variables. For example, to solve min norm(F,1) with x0 = [10;10] and nfmax = 300, you should type:

[x, f_min, nf, stop] = TRFD_composite([10;10], 2, @Ffun, 1, 300)

For solving min max(F), type:

[x, f_min, nf, stop] = TRFD_composite([10;10], 2, @Ffun, 2, 300)

## Contact

Dânâ Davar,  Geovani Nunes Grapiglia <br>
Université catholique de Louvain, Belgium <br>
Department of Mathematical Engineering <br>
dana.davar@uclouvain.be, geovani.grapiglia@uclouvain.be <br>

November 2024

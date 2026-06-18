# ec-symmetry

## Introduction

A set of Maple functions which implement the algorithm described in
[1]; namely: given a pair of real symmetric matrices with symbolic
entries and a desired arrangement of their eigenvalues, compute a
quantifier-free condition on the entries of the matrices so that their
eigenvalues are arranged in the given way.

This implementation uses the theory of symmetric polynomials. See
[ec-symmetry](https://github.com/danielprofili/ec-signature) for an
alternative algorithm which solves the same problem.

## Usage
The main function is `ECsym`.

```
# In : F, G, symbolic real symmetric matrices
# Out: Csym, Asym such that c = EC(F,G) iff c = Csym . v(Asym)
```

## Example
Let
```
  F = [ 1  -2  ]       G =  [ 1        p^2+p  p + 1 ]
      [ -2  1  ] ,          [ p^2 + p  -p     p - 1 ]
                            [ p + 1    p - 1  p     ]
```

where `p` is a parameter. Suppose we would like to find a
quantifier-free condition
on `p` so that the eigenvalue configuration of `F` and `G` is
`(1,1)`; i.e., the eigenvalues of `F` (red) and `G` (blue) are
arranged as follows:

![EC of 011](http://danielprofili.github.io/resources/ec11.png)

Make sure that Maple's working directory is the same directory as the
`ECSymmetry.mpl` script. In a Maple worksheet, enter the following.

```
restart:
with(LinearAlgebra):
read("ECSymmetry.mpl");
F := Matrix([[1,-2],[-2,1]]);
G := Matrix([[1, p^2 + p, p + 1], [p^2 + p, -p, p-1], [-p + 1, p - 1, p]]);
Csym, Asym := ECSym(F,G);
```

![Maple output](http://danielprofili.github.io/resources/sym-output.png)

The outputs are a numeric matrix `Csym` and a column vector `Asym` of
polynomials in `p`, such that if `c` is any valid eigenvalue
configuration, then 

```
c = Csym . v(Asym),
```

where `v(Asym)` means to take the sign variation count of the
coefficients of each entry in `Asym`.

## Reference
1. Hong, Hoon, **Daniel Profili**, and J. Rafael Sendra. ["Conditions
   for eigenvalue configurations of two real symmetric matrices
   (symmetric function approach)."](https://arxiv.org/pdf/2401.00089)
   Journal of Symbolic Computation (pending) (2026).

# GIGcpp

This package provides a generator for generalized inverse Gaussian distributions that adaptively controls the rejection rate given a user-specified upper bound.



## Installation
```
devtools::install_github(repo = "Xiaozhu-Zhang1998/GIGcpp")
library(GIGcpp)
```

## Function
```
rgig(N, lambda, psi, chi, eps = 0.5, K = 0L)
```

`N`: number of variates.

`lambda`:	parameter lambda.

`psi`: parameter psi.

`chi`: parameter chi.

`eps`: desired rejection rate between 0 and 1; active only when `K` is 0.

`K`: desired number of cutoff points.

This function returns a vector of `N` random variates of the specified GIG distribution.


## Details

See [Zhang and Reiter (2022)](https://arxiv.org/abs/2211.13049).

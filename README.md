
# CompEngineR

Calculate CompEngine’s time-series features quickly in R using a C++
backend

## Installation

You can install the development version of `CompEngineR` from GitHub:

``` r
devtools::install_github("hendersontrent/CompEngineR")
```

## General purpose

[`tsfeatures`](https://github.com/robjhyndman/tsfeatures) is an
incredibly powerful and valued R package that provides access to 63
time-series features. Among these are a selection of features from the
self-organising time-series analysis database
[CompEngine](https://www.comp-engine.org/), which stem from the
extensive Matlab library of $>7700$ features
[`hctsa`](https://github.com/benfulcher/hctsa). The methods to produce
these features werte originally presented in [Fulcher, Little, and Jones
(2013)](https://royalsocietypublishing.org/doi/10.1098/rsif.2013.0048)
and [Fulcher and Jones
(2017)](https://www.sciencedirect.com/science/article/pii/S2405471217304386).
However, recent work by [Henderson and Fulcher
(2021)](https://ieeexplore.ieee.org/abstract/document/9679937) showed
that the computation time of `tsfeatures` was substantially slower than
other open-source time-series feature sets. Further investigation
revealed that the `"CompEngine"` features were the slowest to evaluate
by several orders of magnitude. As such `CompEngineR` aims to improve
this by coding these slower features in C++.

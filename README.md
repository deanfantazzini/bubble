## bubble

[![Travis build status](https://travis-ci.org/deanfantazzini/bubble.svg?branch=master)](https://travis-ci.org/deanfantazzini/bubble)

The `bubble` package and the `bitcoinFinance` package will be the companion materials for a textbook I am working on. This forthcoming book wants to
organize the materials that I developed over time for the Perm Summer/Winter School on Blockchain and CryptoEconomics, as well as for my course in Numerical
methods for finance using R.

The `bubble` package implements several models proposed to test for financial bubbles and explosive behavior. These tests span several fields of research, 
  from economics to econophysics, from statistics to computer science.
  
Note that the functions computing the DS LPPLS Confidence and Trust indicators and the critical values for the GSADF test require a lot of time, even with parallel computing. An optimized version would required coding these functions using C++: given that this goes beyond the scope
of my courses, this was not done, but it is left to the interested reader as a (hopefully) interesting exercise. 

You can install the package using the commands:
``` {.r}
if (!require("devtools")) install.packages("devtools")
devtools::install_github("deanfantazzini/bubble")
```

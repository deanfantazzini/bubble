## bubble

[![Travis build status](https://travis-ci.org/deanfantazzini/bubble.svg?branch=master)](https://travis-ci.org/deanfantazzini/bubble)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/deanfantazzini/bubble?branch=master&svg=true)](https://ci.appveyor.com/project/deanfantazzini/bubble)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

The `bubble` package and the `bitcoinFinance` package are the companion materials for my textbook titled [Quantitative Finance with R and Cryptocurrencies](https://www.amazon.com/dp/1090685319?ref_=pe_3052080_397514860). This book organizes the materials that I developed over time for the Perm Summer/Winter School on Blockchain and CryptoEconomics, as well as for my course in Numerical methods for finance using R.

The `bubble` package implements several models proposed to test for financial bubbles and explosive behaviour. These tests span several fields of research, 
  from economics to econophysics, from statistics to computer science.
  
Note that the functions computing the DS LPPLS Confidence and Trust indicators and the critical values for the GSADF test require a lot of time, even with parallel computing. An optimized version would require coding these functions using C++: given that this goes beyond the scope
of my courses, this was not done, but it is left to the interested reader as a (hopefully) interesting exercise. 

You can install the package using the commands:
``` {.r}
if (!require("devtools")) install.packages("devtools")
devtools::install_github("deanfantazzini/bubble")
```

---
title: NEWS/Changelog for package pvars
subtitle: pvars - VAR Modeling for Heterogeneous Panels
          | CRAN <https://cran.r-project.org/package=pvars> 
          | development repository <https://github.com/Lenni89/pvars>
---
***
# pvars 1.2.0

### Features: 
* `pid.iv` and `id.iv` now allow to impose all listed `S2` restrictions on SVD 
results automatically.

### Removed: 
* Unused argument `R0` with misleading documentation has been removed from 
`pid.iv` and `id.iv`.

### Changes or Improvements: 
* Argument `MG_IRF` in `irf.pvarx` and `sboot.pmb` is now set to the recommended 
default `FALSE`.

***
# pvars 1.1.0

### Documentation: 
* A vignette has been added to the package.

### Features: 
* `rboot.normality` now summarizes the function calls for bootstrapping the test 
distribution, calculating empirical test statistics, and retrieving p-values.

***
# pvars 1.0.0

Initial version for new submission to CRAN. 

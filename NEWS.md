# eventglm 1.2.0

- Remove tidyverse lifecycle badge

# eventglm 1.2.0

New features:
  - Time vector allowed in cumincglm. Use this to model multiple timepoints simultaneously, and allow for time varying covariate effects using tve() (for time varying effect) in the right side of the formula. 
  - Inverse probability of censoring weights are returned by the fitting functions in the list element called "ipcw.weights".
  - Small edits and improvements to the documentation.

# eventglm 1.1.1


Bug fixes:
  - Check and fix name clashes for reserved variables pseudo.vals, .Tci, and .Ci
  - Minor efficiency updates in corrected covariance estimation
  - Improve documentation 

# eventglm 1.1.0

New features and vignette: 
  - Methods for computing pseudo observations are now in modules. 
  - Users can define their own modules for computing pseudo observations
  - Option to use survival instead of cumulative incidence in standard models

# eventglm 1.0.2

* Update tests/jackknife-agree.R to do tolerance based comparison rather than floating point comparison

# eventglm 1.0.1

* fix typos and add reference to DESCRIPTION

# eventglm 1.0.0

* Initial release

# MEMWAS 0.1.1 ----------

* Initial CRAN submission (again).

## Methodology advances
### Multiple serial processes

* Each component may use AR(1), continuous-time exponential/OU, AR(p), ARMA(1,1), compound symmetry, Toeplitz, or unstructured covariance and its own one-column design.
* Execution of predictor-modulated latent processes of the form `diag(h_k) u_k`, simultaneous covariance assembly, identifiability checks, component-specific latent modes, conditional prediction, and simulation.
* Enabling tests covering simultaneous AR(1), CS, and UN components, multiple outcome families, covariance assembly, prediction, validation filtering, and confounding safeguards.


### Compatibility

* Create helper functions for extending family likelihood to non-Gaussian distributions.
* Update function `set_MEMWAS()` for extending family likelihood to non-Gaussian distributions.
* Update function `fit_MEMWAS()` for extending family likelihood to non-Gaussian distributions.
* The original Gaussian identity-link method remains the default for those models. Set `control = list(force_repaired = TRUE)` to use the joint latent engine for Gaussian models too. Non-Gaussian calls are routed through the corrected engine.


## Functions
* Add function `serial_component()` for fitting several independently parameterized serial processes in one mixed-effects likelihood.
* Update functino `set_MEMWAS()` with a `serial` model argument for fitting several independently parameterized serial processes in one mixed-effects likelihood.
* Update function `fit_MEMWAS()` with a `serial` model argument for fitting several independently parameterized serial processes in one mixed-effects likelihood.
* Preserved component definitions through `set_MEMWAS()`, grouped-fold `tune_MEMWAS()`, and the final refit. The legacy scalar `autocor` interface remains available.
* Add function `tune_MEMWAS()` and its print, summary, and prediction functions.


## Note
Typhoon Bavi (July 10–12) gave me the perfect uninterrupted time to get MEMWAS updated.


# MEMWAS 0.1.0 ----------

* Initial CRAN submission (failed).
* Document the package (revise for several times).


# MEMWAS 0.0.9 ----------

* Create helper functions for simulating panel data and lagging data.


# MEMWAS 0.0.8 ----------

* Add predicting function for the output of MEMWAS_fit class object using S3Method.
* Update function `fit_MEMWAS()` for making prediction and tuning.


# MEMWAS 0.0.7 ----------

* Add printing and summarizing functions for the output of MEMWAS_fit class
object using S3Method.
* Update function `fit_MEMWAS()` for printing and summarizing output.


# MEMWAS 0.0.6 ----------

* Add printing and summarizing functions for the output of MEMWAS class object
using S3Method.
* Update function `set_MEMWAS()` for printing and summarizing output.


# MEMWAS 0.0.5 ----------

* Create helper functions of assumption checks for function `set_MEMWAS()`.
* Update function `set_MEMWAS()` for checking assumptions.


# MEMWAS 0.0.4 ----------

* Create helper functions of approximation for function `set_MEMWAS` and `fit_MEMWAS()`.
* Update function `set_MEMWAS()` for separating initial approximation method from
the final approximation method.
* Update function `fit_MEMWAS()` for separating initial approximation method from
the final approximation method.


# MEMWAS 0.0.3 ----------

* Create common helper functions for both function `set_MEMWAS()` and `fit_MEMWAS()`.
* Update function `set_MEMWAS()` for separating helper functions.
* Update function `fit_MEMWAS()` for separating helper functions.


# MEMWAS 0.0.2 ----------

* Add prototype of function `fit_MEMWAS()`.


# MEMWAS 0.0.1 ----------

* Initial the MEMWAS package with prototype of function `set_MEMWAS()`.

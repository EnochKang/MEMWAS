# MEMWAS

 <table style="border: none; border-collapse: collapse; margin-left: 0; margin-right: auto;">
   <tr style="border: none;">
     <td style="border: none; vertical-align: middle;">
      <img src = "https://github.com/EnochKang/MEMWAS/blob/main/vignettes/MEMWAS_Logo.png?raw=true" align = "left" width = "120" style="margin-left: 20px;" />
     </td>
     <td style="border: none; vertical-align: middle; padding-left: 20px; text-align: left; line-height: 1.5;">
       <span style="color: #00008B; font-size: 1.2em; font-weight: bold;"> M</span>ixed<br>
       <span style="color: #00008B; font-size: 1.2em; font-weight: bold;">E</span>ffects<br>
       <span style="color: #00008B; font-size: 1.2em; font-weight: bold;">M</span>odels<br>
       <span style="color: #00008B; font-size: 1.2em; font-weight: bold;">W</span>ith<br>
       <span style="color: #00008B; font-size: 1.2em; font-weight: bold;">A</span>utocorrelation<br>
       <span style="color: #00008B; font-size: 1.2em; font-weight: bold;">S</span>tructures
     </td>
   </tr>
      <tr style="border: none;">
     <td style="border: none; vertical-align: middle; padding-left: 20px; text-align: center;">
       Version 0.1.1 <br>
       July 13, 2026
     </td>
     <td style="border: none; vertical-align: middle; padding-left: 20px; text-align: left; line-height: 1.5;">
     </td>
   </tr>
 </table>

&nbsp;

MEMWAS is a base-R package for fitting longitudinal mixed-effects models for Gaussian, 
Bernoulli/grouped-binomial, Poisson, negative-binomial, Gamma, and exponential outcomes. 
Models may combine fixed effects, Gaussian random effects, offsets, fixed-effect 
penalties, and multiple independent named serial processes in one likelihood. 
Each serial process has a one-column design and its own AR(1), exponential/Ornstein-Uhlenbeck, 
AR(p), ARMA(1,1), compound-symmetry, Toeplitz, or unstructured covariance, 
permitting residual and predictor-modulated autocorrelation structures to be 
estimated simultaneously. Grouped-binomial totals, family links, conditional predictions, 
and serial-component contributions are retained consistently through fitting and 
prediction. Besides, the MEMWAS package offers an integrated transparency advantage 
by embedding configurable assumption screening and structured diagnostic reporting 
within the model- analysis workflow. This design supports research integrity by 
making diagnostic choices, findings, limitations, and unavailable tests more visible 
and auditable.


## Environment

- R version R >= 4.6.0
- RTools version 4.5.6768
- C++ 11

## Installation


```r
# From a local source package archive
install.packages("MEMWAS_0.1.1.tar.gz", repos = NULL, type = "source")

# Or from GitHub
remotes::install_github("EnochKang/MEMWAS")
```

## Basic use

```r
library(MEMWAS)

settings <- set_MEMWAS(
  formula = y ~ x1 + x2,
  family = "gaussian",
  data = train_data,
  id = "subject_id",
  time = "visit",
  random = ~ 1,
  autocor = "AR(1)")

fit <- fit_MEMWAS(settings)
summary(fit)
predict(fit, newdata = test_data, type = "response")
```

## Expanded families and approximations

For non-Gaussian mixed models, MEMWAS uses a hierarchical strategy: variational inference initializes the model, and the selected final native-C++ marginal approximation is then optimized. The default final approximation is Laplace. Users can request adaptive Gauss-Hermite quadrature for small random-effect dimensions, saddlepoint approximation for tail-sensitive outcomes, skew-corrected Laplace for asymmetric conditional posteriors, final variational inference for fast screening, or legacy PQL.

```r
fit_nb <- fit_MEMWAS(
  emergency_visits ~ care_plan + baseline_risk,
  family = "negative_binomial",
  approximation = "laplace",
  init_approximation = "variational_inference",
  control = list(negative_binomial_theta = 1.8),
  data = clinic_panel,
  id = "patient_id",
  time = "month",
  random = ~ 1,
  autocor = "AR(1)"
)

fit_gamma <- fit_MEMWAS(
  care_cost ~ treatment_program + baseline_severity,
  family = "gamma",
  approximation = "laplace",
  control = list(gamma_shape = 2.5),
  data = cost_panel,
  id = "patient_id",
  time = "quarter",
  random = ~ 1,
  autocor = "TOEP"
)
```

Supported approximation values are `"variational_inference"`, `"laplace"`, `"saddlepoint"`, `"skew_corrected_laplace"`, `"adaptive_gaussian_quadrature"`, `"adaptive_gauss_hermite_quadrature"`, and legacy `"pql"`. Gaussian models use exact marginal Gaussian ML/REML regardless of the approximation argument.

Approximation quality can be inspected with `diagnose_approximation(fit)`. Sensitivity across final approximation methods can be checked with `compare_approximations(settings, approximations = c("laplace", "adaptive_gauss_hermite_quadrature", "saddlepoint"))`.


## Engine selection

The main functions accept `engine = "R"` or `engine = "cpp"`:

```r
setup <- set_MEMWAS(
  y ~ x1 + x2,
  data = dat,
  id = "id",
  time = "time",
  engine = "cpp")
fit <- fit_MEMWAS(setup)
```

The C++ engine uses `.Call` helpers in `src/memwas_native.cpp` for numerical kernels such as safe exponentials, diagonal matrices, linear solves, Cholesky factorizations, coordinate descent, spline bases, autocorrelation helper calculations, kernels, and metrics. Formula handling, data validation, and S3 object assembly remain in R.

See the vignettes for the model framework, methodological workflow, approximation details, and function arguments.

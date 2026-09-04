# Changelog

Development log for the `Stem` package: what changed, when, and why.

Commit messages in this repository are deliberately terse (`Update <date>`), so
that the file listing on GitHub stays readable. The reasoning behind each change
lives here instead. Entries are newest first.

For the user-facing summary of a release, see [`NEWS.md`](NEWS.md). For the
detailed account of the defects fixed in 2.0.0, see
[`dev/code-changes-report.html`](dev/code-changes-report.html). For the mapping
between the reference papers and the code, see
`inst/extdata/STEM_model_verification.html`.

---

## Unreleased

### 2026-09-04

**The conceptual map is now generated, in both formats.**
`inst/scripts/make-function-map.R` holds the map as tables of nodes and edges
and emits both `inst/extdata/STEM_function_map.svg` and the new
`inst/extdata/STEM_function_map.pdf`, so the two cannot drift apart. It ships
with the package, so a user can rerun it. The layout stays explicit rather than
computed by a graph-drawing algorithm, because the point of the map is the
reading imposed on the package - four bands, from the data up to the clustered
layer - which an automatic layout would not reproduce.

`STEM_function_map_original.pdf` has been removed. It was the map of version
1.0: dot-separated names (`Stem.Model`, `Stem.Estimation`) that no longer exist,
and no SC-STEM layer at all. It documented a package that is not this one.

The regeneration also caught a name the diagram had missed: it still said
`scstem_neighbours`, whereas the function was renamed `scstem_neighbors` in the
American-English pass.

**Documented the statistical features and the scope of the model.**
A new section of the README, and a matching block in `?SCSTEM_Estim`, state what
the family covers and what it does not: Gaussian response only; univariate
response, with the multivariate part being the latent state; exponential spatial
covariance, hence isotropic and stationary within a regime, the partition being
the only source of non-stationarity across the domain; point-referenced
locations, fixed over time; a discrete and regularly spaced time index; a known
loading matrix, common across regimes; `C0` fixed while `m0` is estimated.

The entry that matters most is that **missing values are not supported**, in the
response no more than in the covariates or the coordinates: `STEM_Model()`
rejects them and the Kalman recursion has no partial-observation branch. This is
a restriction of the implementation, not of the model, and the two cases are not
equally hard: for the response the standard state-space treatment applies
directly, whereas gaps in the covariates would require a stochastic E-step,
since the design matrix enters the closed-form M-step updates.

**Package logo.**
The hexagon now carries the three ideas the package is about: a relief
silhouette for space, a node-and-edge network split by color into two regimes,
and one time series per regime along the base. The faint gray edges are the ones
the partition cuts, which is what the spatial penalty pays for.

The wordmark was the hard part, and the difficulty was structural rather than
typographic. In a pointy-top hexagon the width collapses below the widest band -
at the baseline where the word had been sitting the shape is only about 24 units
across - so the word could not grow without being clipped. Moving it up into the
full-width band took it from 15 units to 26. Two further adjustments came out of
looking at the render: the series were made angular rather than wavy, because
smooth curves under a mountain range read as water, and the ridge was broadened,
because sharp peaks and sharp series were competing for the same reading.

`man/figures/logo.svg` is the scalable master and `man/figures/logo.png` the
rendered artifact, both produced from one set of coordinates. The README points
at the PNG: the SVG carries live text, so its wordmark depends on whichever
serif the viewer has, and GitHub sanitizes SVG attributes before serving them,
possibly including the `textLength` that guarantees the word fits.
`dev/make-logo.R` redraws the PNG at any size. It uses base graphics, which can
only clip to a rectangle, so every element of the composition is laid out to
fall inside the hexagon and no clipping is needed.

**American English throughout, and ASCII only.**
Every source, vignette and document converted to American spelling: `-ize`
rather than `-ise`, `neighbor` rather than `neighbour`, `modeling`, `center`,
`meter`, `analyze`. Non-ASCII characters removed as well - em and en dashes
folded to plain hyphens, curly quotes to straight ones, accented letters in the
references to their ASCII forms. Beyond consistency this keeps `R CMD check`
quiet about non-ASCII characters in the sources.

**Commit messages simplified.**
The narrative that used to live in commit subjects moved here. GitHub shows the
last commit message as the caption of every file it touched, and long subjects
made that listing unreadable.

**Fixed the README maths error.**
GitHub processes the markdown before handing the maths to MathJax, and its HTML
sanitizer read `<j}` in `\sum_{i<j}` as the start of a tag. What reached MathJax
was truncated, with a brace opened and never closed - exactly the reported
*"Extra open brace or missing close brace"*. This is why the error survived
three rewrites: every version carried `\sum_{i<j}`. It also explains why a
MathJax probe cleared the same expression: the probe fed the LaTeX straight to
MathJax and skipped the very stage that breaks it.

The penalty is now written as a sum over the edge set of the neighborhood graph,
which carries no character hostile to HTML, is the more natural way to state a
Potts term, and renders in every engine. The MathJax macro `\lt` would also have
worked on GitHub but is not standard LaTeX and would break the PDF manual, so
the edge-set notation is used in the README, in the SC-STEM vignette (GitHub
renders `.Rmd` files too) and in the `\deqn` of `?SCSTEM_Estim`.

**README: parameter roles, cluster-wise equations, more references.**
Equations returned to LaTeX, in fenced math blocks. A new table gives the role
and the reading of every parameter: `beta`, `y_t`, `K`, `G`, `Sigma_eta`,
`sigma2eps` as the nugget, `sigma2omega` as the partial sill, `theta` as the
range with `1/theta` its characteristic length, `Sigma_e`, and the initial
state. The SC-STEM part is now written cluster-wise - a measurement equation, a
state equation and a covariance per regime, then the full parameter set
`Psi_k` - making explicit that a regime differs from another not only in level
but in how it behaves in time and space. It is motivated by spatial
heterogeneity as the second law of geography, against Tobler's first law which
`Sigma_e` already encodes. Added Goodchild (2004), Zhu and Turner (2022), and
Maranzano, Mattera and Sugasawa (2026, arXiv:2608.13638).

**Authorship.**
Francesco Caccia restored among the package authors: he contributed to the
package, which is a separate matter from the paper. Michela Cameletti's ORCID
added, and the author list in the README turned into a table with ORCID and
GitHub handles.

**Documentation split in two.**
`inst/extdata/STEM_model_verification.html` now does one thing only: map the
equations of the three reference papers onto the functions that implement them.
It no longer mentions defects or changes, so it stays valid as documentation of
the package rather than of one release. The development record moved to
`dev/code-changes-report.html`, which `.Rbuildignore` excludes, so a
changelog-shaped document does not ship to CRAN.

**Verified the implementation against the three reference papers.**
Read Fasso, Cameletti and Nicolis (2007), the GRASPA technical report of Fasso
and Cameletti (2007) and Fasso and Cameletti (2010), and checked the estimation
code against them equation by equation. Verified as matching: the scaled
covariance and the reparameterization to `log(gamma)`; the likelihood by
prediction-error decomposition; `Qtilde` and `W`; the closed-form updates of
`beta`, `sigma2omega`, `G`, `Sigma_eta` and `m0`; the Newton-Raphson step,
whose gradient and Hessian were re-derived by hand, including the mixed term
that vanishes because the two covariance parameters enter additively; the
kriging predictor; the parametric bootstrap; and the two joint convergence
criteria.

One defect found: in `filterstep` the univariate branch computed the predictive
mean as `f - X beta` where Equation (9) gives `mu_t = X_t beta + K y_t`, so the
regression term must be added. The branch is reached only when a fit has a
single location - essentially never for a pooled model over a network, but
reachable cluster by cluster. It produced a silently wrong log-likelihood, hence
wrong information criteria and a wrong assignment, rather than an error.

The check also confirmed the SC-STEM assignment score: its per-location variance
`sigma2eps + sigma2omega` is exactly `Sigma_e[i,i]` under Equation (7) of the
technical report.

**Google Drive mirror.**
`dev/sync-gdrive.R` mirrors the git-tracked files onto the Google Drive working
folder and deletes what is no longer tracked, so the two copies cannot drift. It
refuses to run on a dirty tree, so the mirror always corresponds to a real
commit, and writes a `MIRROR.txt` recording which one.

**Restored the introductory vignette.**
It had been written as `vignettes/Stem.Rmd` and then deleted by an
`rm -f vignettes/STEM.Rmd` aimed at the older placeholder: the Windows
filesystem is case-insensitive. Rewritten as `vignettes/getting-started.Rmd`, a
name that cannot collide with the package name under any case folding.

---

## 2.0.0 - 2026-09-03

The full release notes are in [`NEWS.md`](NEWS.md). In summary:

* **New SC-STEM layer**: `SCSTEM_Estim()`, `SCSTEM_Infocrit()`,
  `SCSTEM_Select()`, `SCSTEM_Bootstrap()` and `SCSTEM_BootInference()`, with
  their classes and `print()` methods.
* **Algorithmic rewrite**: ICM label update, explicit penalized objective with
  cycle detection, minimum-size constraint, size-preserving swap pass,
  symmetrized neighborhood graph, data-driven penalty scale, final refit.
* **Four corrections to the cluster assignment score**, none of which raised an
  error and all of which changed the estimated partition.
* **Numerical robustness** of `kalman()` and `STEM_Estimation()`: pre-existing
  defects, harmless for a pooled fit and fatal for a clusterwise algorithm. On
  the `pm10` example they took the bootstrap from 4 usable draws out of 12 to
  12 out of 12.
* **Information criteria** on the exact log-likelihood of the final refit, with
  the correct parameter count and sample size.
* **CRAN compliance**: all five review comments addressed, plus the `pm10`
  dataset converted to `.rda`, the non-CRAN `SCDA` dependency removed, and
  `dplyr` and `sf` dropped as unused.
* **New material**: the `povalley` dataset, three vignettes, a `testthat` suite
  of 43 assertions, and a conceptual map of the package.

`R CMD check --as-cran`: `Status: OK`, no errors, warnings or notes.

---

## 1.0 - 2010

Original CRAN release by Michela Cameletti, subsequently archived.

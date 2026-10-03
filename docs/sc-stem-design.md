# Design notes

[Back to the README](../README.md) | [All pages](README.md)

The implementation choices that make SC-STEM behave sensibly in practice. The estimation algorithms, the stopping rule and the limits of the range are on the page [Computational aspects](computational-aspects.md).


Six points make the implementation behave sensibly in practice. The first five
concern SC-STEM and are documented in detail in `?SCSTEM_Estimation`; the sixth is
what makes a large network fittable at all.

**An unpenalized fit starts from the departures from the pooled model.** The
default `init_method = "departures"` fits the pooled model and clusters the
locations on how they depart from it: the mean of the residual, its slopes on
the covariates, and the autocorrelation and variance of what is left. Locations
of one regime share their departures; the covariate means, the earlier default
and still an option, say nothing about the regimes when the covariates are
exogenous to them.

**A penalized fit starts from the unpenalized one.** A fit with $\phi > 0$
starts from the partition of the fit with $\phi = 0$ at the same $K$, and the
automatic scale of the penalty is measured there, on regimes that are already
fitted. A penalty that is strong from the first sweep would freeze whatever
partition it is given, and a scale measured on the starting partition would
make the same $\phi$ mean different things for different initializations. On a
grid the rule costs nothing: `SCSTEM_Infocrit()` passes the partition of its
$\phi = 0$ fit on. Measured on fitted regimes the scale is of the order of the
whole gain of the right regime over the wrong ones, so the useful values of
$\phi$ are small: the default grid is
$\phi \in \{0, 0.025, 0.05, 0.1, 0.2, 0.5, 1\}$. `SCSTEM_Select()` chooses $K$
by the modal BIC over the moderate values $[0.025, 0.2]$, the pooled model
competing, and then $\phi$ by the smallest BIC at that $K$ over the whole grid:
since every penalized fit starts from the unpenalized one, a penalty is chosen
only when it leads to a partition with a higher likelihood, and the strong
values $0.5$ and $1$ are kept on the grid for data on which they do.

**Labels are updated sequentially (ICM).** Each location maximizes its own
penalized contribution given the current labels of all the others, so a sweep
cannot decrease the objective at fixed parameters and cannot cycle - unlike a
simultaneous update, which remains available as an option.

**Degeneracy is controlled.** Within-cluster homogeneity is exactly what the
assignment step pursues, so the regime-specific variances shrink and, left
unconstrained, the regime with the smallest residual variance attracts every
location. A minimum-size constraint prevents the collapse, and a size-preserving
swap pass restores mobility without breaking feasibility.

**Uncertainty includes the partition.** `SCSTEM_Bootstrap()` re-runs the whole
procedure, clustering included, on every draw, so the reported intervals are not
conditional on a partition that is itself estimated.

**The filter never inverts a $d \times d$ matrix.** The predictive covariance
$Q_t = A P_t A' + \Sigma_\varepsilon$ has $\Sigma_\varepsilon$ constant in $t$
and a rank-$p$ update on top of it, so the Woodbury identity and the matrix
determinant lemma reduce each step to $p \times p$ algebra against a Cholesky
factor of $\Sigma_\varepsilon$ that is computed once per pass. The cost of a
forward pass is $O(d^3 + T d^2)$ instead of $O(T d^3)$, and the log density is
evaluated in closed form rather than as the logarithm of a density, which on a
few hundred locations underflows to `-Inf`. On a network of 200 stations over
365 days one EM iteration is about ten times faster than the direct form, and a
400-station network, which the direct form cannot fit at all, takes about a
second per iteration. The two forms are algebraic identities of each other and
agree to about $10^{-13}$ in relative terms.

---

[Back to the README](../README.md) | [All pages](README.md)

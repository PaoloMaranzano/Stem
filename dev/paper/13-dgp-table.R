## ---------------------------------------------------------------------------
## The design in numbers: what every scenario actually sets, and how large the
## resulting cross is. Prints to the console and writes the LaTeX table used in
## the paper.
##
##   Rscript dev/paper/13-dgp-table.R
## ---------------------------------------------------------------------------

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  h <- if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  source(file.path(h, "00-setup.R"), chdir = TRUE)
}, envir = globalenv())
source(file.path(stem_paper_dir(), "design.R"), chdir = TRUE)
source(file.path(stem_paper_dir(), "06-dgp.R"), chdir = TRUE)
OUT <- stem_fig_dir()

b    <- dgp_base()
scen <- dgp_scenarios()
dims <- stem_design_dims()
tot  <- b$sigma2eps + b$sigma2omega

cat("BASELINE (regime 1 in every scenario), from the pooled fit on the network\n")
cat(sprintf("  beta            = (%.2f, %.2f)   intercept and one standardised covariate\n",
            b$beta[1], b$beta[2]))
cat(sprintf("  sigma2eps       = %.2f          nugget, share %.3f of the total\n",
            b$sigma2eps, b$sigma2eps / tot))
cat(sprintf("  sigma2omega     = %.2f          partial sill\n", b$sigma2omega))
cat(sprintf("  theta           = %.2e      range 1/theta = %.0f km\n",
            b$theta, 1 / b$theta / 1000))
cat(sprintf("  G               = %.2f          persistence of the latent process\n", b$G))
cat(sprintf("  var_y           = %.2f          stationary variance, so sigma2eta = %.2f\n",
            b$var_y, sigma_eta_of(b$G, b$var_y)))
cat(sprintf("  residual sd     = %.2f          sqrt(sigma2eps + sigma2omega)\n\n",
            sqrt(tot)))

cat("WHAT EACH SCENARIO SETS IN THE LAST REGIME (regime K; with K = 3 the\n")
cat("middle regime sits halfway, and with K = 1 every scenario is the baseline)\n\n")

rows <- do.call(rbind, lapply(seq_len(nrow(scen)), function(i) {
  s <- scen[i, ]
  p <- dgp_psi(s$scenario, s$level, K = 2L)[[2]]
  data.frame(id = s$id,
             beta1 = sprintf("%.2f", p$beta[2]),
             d_beta_sd = sprintf("%.2f", (p$beta[2] - b$beta[2]) / sqrt(tot)),
             G = sprintf("%.2f", p$G),
             s2eps = sprintf("%.1f", p$sigma2eps),
             nugget_share = sprintf("%.2f", p$sigma2eps / tot),
             range_km = sprintf("%.0f", 1 / p$theta / 1000),
             rho = sprintf("%.0f", s$rho),
             block = ifelse(s$force_block, "yes", "no"),
             label = s$label, stringsAsFactors = FALSE)
}))
print(rows, row.names = FALSE, right = FALSE)

cat("\nOVERLAP: standardised separation 2*omega/sqrt(nu_sp), the same for every K\n")
print(round(stats::setNames(2 * dims$omega / sqrt(NU_SP),
                            format(dims$omega)), 2))
cat(sprintf("nu_sp = %.1f, one abstract unit = %d km, so a cluster has sd %.0f km\n",
            NU_SP, UNIT_KM, sqrt(NU_SP) * UNIT_KM))
cat(sprintf("and at omega = 1 the centres are %d km apart\n", 2 * UNIT_KM))

cat("\nREGIME SIZES\n")
for (bal in dims$balance) {
  cat("  ", bal, "\n", sep = "")
  for (K in setdiff(dims$K, 1L)) {
    line <- vapply(dims$n, function(n) {
      s <- dgp_sizes(n, K, bal)
      ok <- dgp_feasible(n, K, bal)
      paste0(paste(s, collapse = "/"), if (ok) "" else " (X)")
    }, character(1))
    cat(sprintf("    K=%d  %s\n", K,
                paste(sprintf("n=%d: %s", dims$n, line), collapse = "   ")))
  }
}
cat("  (X) marks a cell that cannot carry regimes of at least ", N_MIN,
    " units, or whose imbalance does not survive that, and is dropped\n", sep = "")

## ---------------------------------------------------------------------------
## The size of the experiment
## ---------------------------------------------------------------------------
cells <- stem_design_cells()
keep  <- mapply(dgp_feasible, cells$n, cells$K, cells$balance)
cost  <- stem_design_cost(cells[keep, ], nrep = 100L)
cat(sprintf("\nDESIGN: %d cells, of which %d feasible; blocks %s\n",
            nrow(cells), sum(keep),
            paste(stem_design_blocks(), collapse = ", ")))
cat(sprintf("  %.0f core-hours at M = 100, %.1f h on %d cores\n",
            cost$core_hours, cost$wall_hours, cost$cores))
for (b in stem_design_blocks())
  cat(sprintf("    %-10s %4d cells\n", b, nrow(stem_design_cells(blocks = b))))

## ---------------------------------------------------------------------------
## LaTeX
## ---------------------------------------------------------------------------
tex <- c(
  "\\begin{tabular}{llrrrrrr}",
  "\\toprule",
  "& & \\multicolumn{6}{c}{parameters of the last regime} \\\\",
  "\\cmidrule(l){3-8}",
  "id & separation & $\\beta_1$ & $\\Delta\\beta_1/\\sigma$ & $G$ & nugget share & range (km) & $\\rho$ \\\\",
  "\\midrule",
  sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\",
          rows$id, rows$label, rows$beta1, rows$d_beta_sd, rows$G,
          rows$nugget_share, rows$range_km, rows$rho),
  "\\bottomrule",
  "\\end{tabular}")
writeLines(tex, file.path(OUT, "tab_dgp.tex"))
cat("wrote", file.path(OUT, "tab_dgp.tex"), "\n")

## ---------------------------------------------------------------------------
## The design in numbers: what every scenario actually sets, and how large the
## resulting cross is. Prints to the console and writes the LaTeX table used in
## the paper.
##
##   Rscript dev/paper/13-dgp-table.R
## ---------------------------------------------------------------------------

HERE <- "C:/Users/paulm/OneDrive/Documenti/GitHub/Stem"
OUT  <- "C:/Users/paulm/Dropbox/Applicazioni/Overleaf/SC-STEM package paper/Figures"
source(file.path(HERE, "dev", "paper", "06-dgp.R"))

b    <- dgp_base()
scen <- dgp_scenarios()
dims <- dgp_dims()
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
print(round(setNames(2 * dims$d / sqrt(NU_SP), c("0", "1/3", "2/3", "1")), 2))
cat(sprintf("nu_sp = %.1f, one abstract unit = %d km, so a cluster has sd %.0f km\n",
            NU_SP, UNIT_KM, sqrt(NU_SP) * UNIT_KM))
cat(sprintf("and at omega = 1 the centres are %d km apart\n", 2 * UNIT_KM))

cat("\nREGIME SIZES\n")
for (bal in dims$balance) {
  cat("  ", bal, "\n", sep = "")
  for (K in c(2L, 3L)) {
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
cells <- rbind(
  expand.grid(n = dims$n, TN = dims$TN, K = 1L, d = 2/3, id = "S0",
              balance = "balanced", stringsAsFactors = FALSE),
  expand.grid(n = dims$n, TN = dims$TN, K = c(2L, 3L), d = dims$d,
              id = scen$id, balance = dims$balance, stringsAsFactors = FALSE))
keep <- mapply(dgp_feasible, cells$n, cells$K, cells$balance)
cat(sprintf("\nFULL DESIGN: %d cells, of which %d feasible\n",
            nrow(cells), sum(keep)))

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

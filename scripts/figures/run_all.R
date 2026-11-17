# scripts/figures/run_all.R
#
# Regenerates every manuscript panel produced from the scRNA-seq data:
#   Figure 5c-g, Supplementary Figure 2a-h, Supplementary Figure 5b-e.
# Each figure script runs in a fresh R process so that one failure (e.g. no
# internet for Enrichr) does not stop the others.
#
# Run from the repository root:  Rscript scripts/figures/run_all.R

this.dir <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else file.path(getwd(), "scripts", "figures") })
rscript <- file.path(R.home("bin"), "Rscript")
for (s in c("Figure_5.R", "Supplementary_Figure_2.R", "Supplementary_Figure_5.R")) {
  cat("==== ", s, " ====\n")
  status <- system2(rscript, file.path(this.dir, s))
  if (status != 0) cat("!! ", s, " exited with status ", status, "\n")
}

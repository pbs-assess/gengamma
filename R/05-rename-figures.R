f <- c(
  "figure-1-probability-densities.pdf",
  "figure-2-cross-sim-qq.pdf",
  "figure-3-cross-sim-RE-AIC-CI-weight.pdf",
  "figure-4-50-fold.pdf",
  "figure-5-index-rqr.pdf"
  )
fclean <- c(
  "fig1.pdf",
  "fig2.pdf",
  "fig3.pdf",
  "fig4.pdf",
  "fig5.pdf"
)

if (Sys.info()[["user"]] == "jilliandunic") {
  for (i in seq_along(f)) {
    file.copy(paste0("figures/", f[i]), paste0("figures/submission/", fclean[i]), overwrite = TRUE)
  }
}
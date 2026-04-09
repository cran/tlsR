## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 6,
  fig.height = 5,
  eval = TRUE
)

## ----load-data----------------------------------------------------------------
library(tlsR)

data(toy_ldata)

# Structure of the built-in example dataset
str(toy_ldata)
table(toy_ldata[["ToySample"]]$phenotype)

## ----detect-tls---------------------------------------------------------------
# Ensure toy data has expected columns for the new validation
data(toy_ldata)
if (!"phenotype" %in% names(toy_ldata[["ToySample"]])) {
  toy_ldata[["ToySample"]]$phenotype <- toy_ldata[["ToySample"]]$coarse_phen_vec   # or whatever the correct mapping is
}
ldata <- detect_TLS(
  LSP                     = "ToySample",
  k                       = 30,     # neighbours for density estimation
  bcell_density_threshold = 15,     # min avg 1/k-distance (um)
  min_B_cells             = 50,     # min B cells per candidate TLS
  min_T_cells_nearby      = 30,     # min T cells within max_distance_T
  max_distance_T          = 50,     # search radius (um)
  ldata                   = toy_ldata
)

table(ldata[["ToySample"]]$tls_id_knn)

## ----base-plot, fig.alt="Scatter plot of ToySample cells coloured by TLS membership"----
df  <- ldata[["ToySample"]]
col <- ifelse(df$tls_id_knn == 0, "grey80",
              c("#0072B2", "#009E73", "#CC79A7")[df$tls_id_knn])
plot(df$x, df$y,
     col  = col, pch = 19, cex = 0.3,
     xlab = "x (um)", ylab = "y (um)",
     main = "Detected TLS — ToySample")
legend("topright",
       legend = c("Background", paste0("TLS ", sort(unique(df$tls_id_knn[df$tls_id_knn > 0])))),
       col    = c("grey80", "#0072B2", "#009E73", "#CC79A7"),
       pch    = 19, pt.cex = 1.2, bty = "n")

## ----scan, eval = FALSE-------------------------------------------------------
# # eval=FALSE because this step can take ~10–30 s on real data
# windows <- scan_clustering(
#   ws        = 500,          # window side (um)
#   sample    = "ToySample",
#   phenotype = "B cells",
#   nsim      = 39,           # Monte Carlo simulations (39 → p < 0.05)
#   plot      = FALSE,
#   ldata     = ldata
# )
# 
# cat("Significant windows:", length(windows), "\n")
# # Access the first window's centre and cell count:
# if (length(windows) > 0) {
#   cat("Centre:", windows[[1]]$window_center, "\n")
#   cat("Cells: ", windows[[1]]$n_cells, "\n")
# }

## ----icat---------------------------------------------------------------------
n_tls <- max(ldata[["ToySample"]]$tls_id_knn, na.rm = TRUE)

if (n_tls >= 1) {
  icat_scores <- vapply(
    seq_len(n_tls),
    function(id) calc_icat("ToySample", tlsID = id, ldata = ldata),
    numeric(1)
  )
  names(icat_scores) <- paste0("TLS", seq_len(n_tls))
  print(icat_scores)
}

## ----detect-tic---------------------------------------------------------------
ldata <- detect_tic(
  sample           = "ToySample",
  min_pts          = 10,   # HDBSCAN minPts
  min_cluster_size = 10,   # drop clusters smaller than this
  ldata            = ldata
)

table(ldata[["ToySample"]]$tcell_cluster_hdbscan, useNA = "ifany")

## ----summary------------------------------------------------------------------
sumtbl <- summarize_TLS(ldata, calc_icat_scores = FALSE)
print(sumtbl)

## ----plot-tls, fig.alt="ggplot2 spatial map of ToySample with TLS and TIC highlighted"----
p <- plot_TLS(
  sample     = "ToySample",
  ldata      = ldata,
  show_tic   = TRUE,
  point_size = 0.5,
  alpha      = 0.7
)

## ----plot-custom, fig.alt="Customised TLS plot with dark theme"---------------
library(ggplot2)
p + theme_dark() + labs(title = "ToySample — dark theme")

## ----multi-sample, eval = FALSE-----------------------------------------------
# samples <- names(ldata)
# 
# ldata <- Reduce(function(ld, s) detect_TLS(s, ldata = ld), samples, ldata)
# ldata <- Reduce(function(ld, s) detect_tic(s,  ldata = ld), samples, ldata)
# 
# summary_all <- summarize_TLS(ldata)
# print(summary_all)

## ----session------------------------------------------------------------------
sessionInfo()


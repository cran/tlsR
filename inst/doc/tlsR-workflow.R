## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 7,
  fig.height = 6,
  eval = TRUE
)

## ----load-data----------------------------------------------------------------
library(tlsR)

data(toy_ldata)

# Structure of the built-in example dataset
str(toy_ldata)
table(toy_ldata[["ToySample"]]$phenotype)

## ----detect-tls---------------------------------------------------------------
data(toy_ldata)

ldata <- detect_TLS(
  LSP                     = "ToySample",
  k                       = 10,     # neighbours for density estimation
  bcell_density_threshold = 17,     # min avg 1/k-distance (um)
  min_B_cells             = 100,    # min B cells per candidate TLS
  min_T_cells_nearby      = 5,      # min T cells within max_distance_T
  max_distance_T          = 50,     # search radius (um)
  expand_distance         = 100,    # expanding radius
  ldata                   = toy_ldata
)

table(ldata[["ToySample"]]$tls_id_knn)

## ----base-plot, fig.alt="Scatter plot of ToySample cells coloured by TLS membership"----
df <- ldata[["ToySample"]]

plot(df$x[df$tls_id_knn == 0],
     df$y[df$tls_id_knn == 0],
     col  = "grey80", pch = 19, cex = 0.3,
     xlab = "x (um)", ylab = "y (um)",
     main = "Detected TLS -- ToySample")

points(df$x[df$tls_id_knn > 0],
       df$y[df$tls_id_knn > 0],
       col = "#0072B2", pch = 19, cex = 0.4)

legend("bottomright",
       legend = c("Background", "TLS"),
       col    = c("grey80", "#0072B2"),
       pch    = 19, pt.cex = 1.2, bty = "n")

## ----scan-B, eval = FALSE-----------------------------------------------------
# # eval=FALSE because this can take ~10--30 s on real data
# L_B <- scan_clustering(
#   ws             = 1000,        # window side (um)
#   sample         = "ToySample",
#   phenotype      = "B cells",
#   plot           = TRUE,
#   creep          = 1L,
#   min_cells      = 10L,
#   min_phen_cells = 5L,
#   label_cex      = 1.1,        # increase if CI labels look small
#   ldata          = ldata
# )
# 
# cat("B-cell windows analysed:", length(L_B$B), "\n")

## ----scan-T, eval = FALSE-----------------------------------------------------
# L_T <- scan_clustering(
#   ws        = 500,
#   sample    = "ToySample",
#   phenotype = "T cells",
#   plot      = TRUE,
#   ldata     = ldata
# )
# 
# cat("T-cell windows analysed:", length(L_T$T), "\n")

## ----scan-both, eval = FALSE--------------------------------------------------
# L_both <- scan_clustering(
#   ws        = 3000,
#   sample    = "ToySample",
#   phenotype = "Both",
#   plot      = TRUE,
#   ldata     = ldata
# )
# 
# cat("B windows:", length(L_both$B), " | T windows:", length(L_both$T), "\n")

## ----icat---------------------------------------------------------------------
n_tls <- max(ldata[["ToySample"]]$tls_id_knn, na.rm = TRUE)

if (n_tls >= 1L) {
  icat_scores <- vapply(
    seq_len(n_tls),
    function(id) calc_icat("ToySample", tlsID = id, ldata = ldata),
    numeric(1L)
  )
  names(icat_scores) <- paste0("TLS", seq_len(n_tls))
  print(icat_scores)
}

## ----detect-tic---------------------------------------------------------------
ldata <- detect_tic(
  sample           = "ToySample",
  min_pts          = 20,    # HDBSCAN minPts
  min_cluster_size = 100,   # drop clusters smaller than this
  ldata            = ldata
)

table(
  ldata[["ToySample"]]$tcell_cluster_hdbscan[
    ldata[["ToySample"]]$tcell_cluster_hdbscan != 0
  ],
  useNA = "ifany"
)

## ----summary------------------------------------------------------------------
sumtbl <- summarize_TLS(ldata, calc_icat_scores = FALSE)
print(sumtbl)

## ----plot-tls, fig.alt="ggplot2 spatial map of ToySample with TLS and TIC highlighted"----
p <- plot_TLS(
  sample        = "ToySample",
  ldata         = ldata,
  show_tic      = TRUE,
  point_size    = 0.5,
  alpha         = 0.7,     # TLS / TIC cells
  bg_alpha      = 0.25,    # background cells (more transparent)
  tic_size_mult = 0.8      # TIC cells drawn 1.8x larger
)

## ----plot-custom, fig.alt="Customised TLS plot with additional title"---------
library(ggplot2)
p + labs(title = "ToySample -- Your custom title")

## ----multi-sample, eval = FALSE-----------------------------------------------
# samples <- names(ldata)
# 
# ldata <- Reduce(function(ld, s) detect_TLS(s, ldata = ld), samples, ldata)
# ldata <- Reduce(function(ld, s) detect_tic(s,  ldata = ld), samples, ldata)
# 
# summary_all <- summarize_TLS(ldata)
# print(summary_all)

## ----multi-scan, eval = FALSE-------------------------------------------------
# # Generate one spatial map per sample (side-by-side B and T panels)
# for (s in names(ldata)) {
#   scan_clustering(
#     ws        = 500,
#     sample    = s,
#     phenotype = "Both",    # two-panel plot: B cells | T cells
#     plot      = TRUE,
#     label_cex = 1.2,       # slightly larger CI labels for presentation
#     ldata     = ldata
#   )
# }

## ----session------------------------------------------------------------------
sessionInfo()


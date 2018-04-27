## ---- packages ----

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pacman::p_load_current_gh("sahirbhatnagar/sail")
pacman::p_load(ggplot2)
pacman::p_load(doMC)
registerDoMC(cores = 8)
pacman::p_load(latex2exp)
# pacman::p_load(multipanelfigure)

## ---- globals ----

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
trop <- RSkittleBrewer::RSkittleBrewer("trop")
gg_sy <- theme(legend.position = "bottom", axis.text = element_text(size = 20),
               axis.title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 20))

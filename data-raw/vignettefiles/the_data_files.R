## ---- include = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE,
  message = FALSE,
  warning = FALSE
)

options(knitr.kable.NA = "",
        knitr.table.format = "pandoc")

options("show.signif.stars"=FALSE,
        "stringsAsFactors"=FALSE,
        "max.print"=50000,
        "width"=240)

library(aMSE)
library(codeutils)
library(hplot)
library(knitr)
library(captioner)

tab_nums <- captioner(prefix = "Table")
fig_nums <- captioner(prefix = "Figure")

if (dir.exists("C:/Users/Malcolm/Dropbox")) {
  ddir <- "C:/Users/Malcolm/Dropbox/A_CodeUse/aMSEUse/documentation/"
} else {
  ddir <- "C:/Users/User/Dropbox/A_CodeUse/aMSEUse/documentation/"
}


## ----FIG1, echo=FALSE, out.width = '65%', fig.align='center'----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
filen <- paste0(ddir,"figures/AlternativeMSE.png")
suppressWarnings(knitr::include_graphics(filen))


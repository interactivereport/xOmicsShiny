
## pkgs from BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biomaRt")
BiocManager::install("pathview")
BiocManager::install("ComplexHeatmap")
BiocManager::install("Mfuzz")



filename <- "H:/Rcode/QuickomicsModule/expressionmodule.R"

library(NCmisc)
library(stringr)
library(dplyr)
setwd("H:/Rcode/QuickomicsModule")
NCmisc::list.functions.in.file(filename)


devtools::install_github("brshallo/funspotr")
library(funspotr)
library(dplyr)
filename <- "H:/Rcode/QuickomicsModule/expressionmodule.R"
funspotr::spot_pkgs(filename)

moduleFilelist <- c("qcplotmodule.R", "degmodule.R", "heatmapmodule.R", "expressionmodule.R",
                    "genesetmodule.R", "patternmodule.R", "networkmodule.R", "vennmodule.R", "pcsfmodule.R", "wgcnamodule.R", "drcmodule.R",
                    "dromicsmodule.R", "monotonicmodule.R", "mergedatamodule.R", "app.R","inputdata.R", "process_uploaded_files.R",
                     "outputmodule.R", "util-basicandfitfunc.R")

pkgslist <- c()
functionslist <- list()
for (filename in moduleFilelist){
   print(filename)
  pkgs <- funspotr::spot_pkgs(filename)
  pkgslist <- c(pkgslist, pkgs)
  print(pkgs)
  #functions <- funspotr::spot_funs(filename)
  #functions <- functions %>%
  #  as.data.frame() %>%
  #  dplyr::mutate(filename = filename)
  #print(functions)
  #functionslist[[filename]] <-   functions
}

pkgslist <- sort(unique(pkgslist))

func_df <- functionslist %>%
  map(as.data.frame) %>%
  bind_rows() %>%
  dplyr::pull(pkgs) %>%
  unique() %>% sort()


filename <- "qcplotmodule.R"

for (filename in moduleFilelist){
  print(filename)
  res <- NCmisc::list.functions.in.file(filename)
  print(names(res))
}




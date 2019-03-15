## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----metatable, echo = FALSE---------------------------------------------
metatable = read.csv("GSE86831_meta.csv")
knitr::kable(metatable)

## ----typical, eval=FALSE-------------------------------------------------
#  input_dir = "C:/path/to/IDAT_files"
#  sample_sheet = paste(input_dir, "meta.csv", sep = "/")
#  output_dir = "C:/path/for/output_files"
#  runPipeline(datadir = input_dir,
#              array = "450K",
#              metafile = sample_sheet,
#              outdir = output_dir)

## ----step01, eval=FALSE--------------------------------------------------
#  runPipeline(datadir = input_dir,
#              metafile = sample_sheet,
#              outdir = output_dir,
#              threads = 16,
#              array = "EPIC",
#              begin = 1, end = 1)

## ----step02, eval=FALSE--------------------------------------------------
#  runPipeline(datadir = input_dir,
#              metafile = sample_sheet,
#              outdir = output_dir,
#              array = "EPIC",
#              begin = 2, end = 2)

## ---- out.width="0.3\\linewidth", include=TRUE, fig.align="center", fig.cap=c("your caption"), echo=FALSE----
knitr::include_graphics("./QC02_plot_control_probes.pdf")

## ----step03, eval=FALSE--------------------------------------------------
#  runPipeline(datadir = input_dir,
#              metafile = sample_sheet,
#              outdir = output_dir,
#              begin = 3, end = 3)

## ----step04, eval=FALSE--------------------------------------------------
#  runPipeline(datadir = input_dir,
#              metafile = sample_sheet,
#              outdir = output_dir,
#              begin = 4, end = 4)

## ----freqtable, echo = FALSE---------------------------------------------
freqtable = read.csv("QC04_Detection_Failure_Table_By_Probe.csv")
knitr::kable(freqtable)

## ----step05, eval=FALSE--------------------------------------------------
#  runPipeline(datadir = input_dir,
#              metafile = sample_sheet,
#              outdir = output_dir,
#              begin = 5, end = 5)

## ---- out.width="0.3\\linewidth", include=TRUE, fig.align="center", fig.cap=c("your caption"), echo=FALSE----
knitr::include_graphics("./QC05_outliers.pdf")

## ---- out.width="0.3\\linewidth", include=TRUE, fig.align="center", fig.cap=c("your caption"), echo=FALSE----
knitr::include_graphics("./QC05_methy_colorBias.pdf")

## ---- out.width="0.3\\linewidth", include=TRUE, fig.align="center", fig.cap=c("your caption"), echo=FALSE----
knitr::include_graphics("./QC05_density.pdf")

## ----step06, eval=FALSE--------------------------------------------------
#  runPipeline(datadir = input_dir,
#              metafile = sample_sheet,
#              outdir = output_dir,
#              begin = 6, end = 6)

## ---- out.width="0.3\\linewidth", include=TRUE, fig.align="center", fig.cap=c("your caption"), echo=FALSE----
knitr::include_graphics("./QC06_afterColorAdj_sum_colorBias.pdf")

## ----step07, eval=FALSE--------------------------------------------------
#  runPipeline(datadir = input_dir,
#              metafile = sample_sheet,
#              outdir = output_dir,
#              begin = 7, end = 7)

## ---- out.width="0.3\\linewidth", include=TRUE, fig.align="center", fig.cap=c("your caption"), echo=FALSE----
knitr::include_graphics("./QC07_afterNorm_density.pdf")

## ----step08, eval=FALSE--------------------------------------------------
#  runPipeline(datadir = input_dir,
#              metafile = sample_sheet,
#              outdir = output_dir,
#              threads = 16,
#              begin = 8, end = 8)

## ---- out.width="0.3\\linewidth", include=TRUE, fig.align="center", fig.cap=c("your caption"), echo=FALSE----
knitr::include_graphics("./CheckBMIQ-GSM2309184_200134080019_R08C01.pdf")

## ----step09, eval=FALSE--------------------------------------------------
#  runPipeline(datadir = input_dir,
#              metafile = sample_sheet,
#              outdir = output_dir,
#              begin = 9, end = 9)

## ----step10, eval=FALSE--------------------------------------------------
#  runPipeline(datadir = input_dir,
#              metafile = sample_sheet,
#              outdir = output_dir,
#              begin = 10, end = 10)

## ----svatable, echo = FALSE----------------------------------------------
svatable = read.csv("QC04_Detection_Failure_Table_By_Probe.csv")
knitr::kable(svatable)

## ----cite----------------------------------------------------------------
MethyLiution::cite_me()


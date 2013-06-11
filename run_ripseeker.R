cat("Hello, World!\n")

args <- commandArgs(trailingOnly=T)

input <- args[1]
output.file.prefix <- args[2]
library(RIPSeeker)

outDir <- file.path(getwd(), "ripseeker")
seekOut  <- ripSeek(bamPath = input, strandType="+", outDir = outDir, multicore=TRUE)
export(do.call(c, unname(as.list(seekOut$RIPGRList))), paste(output.file.prefix, "_pos.bed", sep=""))

outDir <- file.path(getwd(), "ripseeker")
seekOut  <- ripSeek(bamPath = input, strandType="-", outDir = outDir, multicore=TRUE)
export(do.call(c, unname(as.list(seekOut$RIPGRList))), paste(output.file.prefix, "_neg.bed", sep=""))
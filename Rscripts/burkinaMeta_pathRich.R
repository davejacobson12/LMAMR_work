args <- commandArgs(TRUE)

inFile <- read.table(args[1], sep = "\t", header = T)
metadata <- read.table(args[1], sep = "\t", header = T)

suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(gridExtra))


cbind(

#!/usr/bin Rscript
## Made by Seongmun Jeong
## lovemun@kribb.re.kr

source("https://gist.githubusercontent.com/SiYangming/83a14a28994d54789b19bff923c53938/raw/33af46fa76217493626ea8fb3c4db7987c0f8f5d/pkgs_in.R")
pkgs_in(c("argparse", "data.table"))
source("https://gist.githubusercontent.com/SiYangming/7eb5c15a44e726839d383973559f6591/raw/81570aa4437af96fbe29c4dfa82953c71326eac0/thisPath.R")

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("data.table"))

parser <- ArgumentParser()
parser$add_argument("file", nargs=1, help="Input txt file")
parser$add_argument("-p", "--predefined", default = "NN")
parser$add_argument("-cv", "--coverage", default=99,
    help = "User defined coverage")
parser$add_argument("-d", "--delta", default=0.001,
    help = "user defined delta")
parser$add_argument("-o", "--output", default="Run",
    help = "output name")
parser$add_argument("-m", "--maf", default = 0, help = "Remove Minor allele frequency")

args <- parser$parse_args()
file <- args$file
cv <- args$coverage
delta <- args$delta
pfile <- args$predefined
output <- args$output
maf <- args$maf
source(file.path(thisPath, "genocore.R"))


if (grepl("csv", file)){
    tdata <- fread(file, header=TRUE, sep=",", check.names=FALSE, data.table=FALSE)
} else {
    tdata <- fread(file, header=TRUE, sep="\t", check.names=FALSE, data.table=FALSE)
}
data.set <- tdata[,-1]
rownames(data.set) <- tdata[,1]
rm(tdata)
data.set[data.set == -1] <- NA
if (pfile != "NN"){
    preset <- scan(pfile, what = "character")
} else {
    preset <- NULL
}
if (maf != 0){
    source(file.path(thisPath, "calc_maf.R"))
    cm <- calc_maf(data.set)
    rm.ix <- which(cm$maf < maf)
    write(rownames(data.set)[rm.ix], file = paste0(output, "Remove_markers.txt"))
    data.set <- data.set[-rm.ix,]
    cat("Remove ", length(rm.ix), " markers from dataset", "\n")
}

Temp_file = paste0(output, "_Temp.csv")
Cover_file = paste0(output, "_Coverage.csv")
Coreset_file = paste0(output, "_Coreset.csv")
core.set(data.set, preset = preset, coverage = cv, delta = delta, Temp_file = Temp_file, coverage_filename = Cover_file, Coreset = Coreset_file)


# Tool that generates chromosome diagrams of human genomic variants listed in a VCF file.
# https://lakshay-anand.github.io/chromoMap/docs.html

# check for option table package
if (!require("optparse"))
  install.packages("optparse")

library(optparse)

# !/usr/bin/env Rscript
option_list = list(
  make_option(
    c("-f", "--filter"),
    action = "store_true",
    default = NULL,
    help = "Show only variants with value “PASS” in FILTER field",
    metavar = "character"
  ),
  make_option(
    c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "Vcf file containing SNPs",
    metavar = "<file>.vcf"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "result.html",
    help = "Output html file visualising chromosome regions",
    metavar = "<file>.html"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# check if input file is included
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call. =
         FALSE)
}

# defining files from command
input_file = opt$input
output_file = opt$output
dir_name = normalizePath(dirname(output_file))

# check if the other dependencies are installed
if (!require("vcfR"))
  install.packages("vcfR")
if (!require("chromoMap"))
  install.packages("chromoMap")
if (!require("htmltools"))
  install.packages("htmltools")
if (!require("gtools"))
  install.packages("gtools")

library(vcfR)
library(chromoMap)
library(htmltools)
library(gtools)

# load of vcf file
vcf <-
  read.vcfR(input_file)

# extracting chromosome lengths by the contig tag in the vcf file
chrom_list <- queryMETA(vcf, element = "contig")
text <- c()
chrom_matrix <- matrix(, nrow = 1, ncol = 2)
for (chrom in chrom_list) {
  chrom_name <- substring(chrom[1], 11)
  chrom_end <- substring(chrom[2], 8)
  chrom_matrix <- rbind(chrom_matrix, c(chrom_name, chrom_end))
  
  text <-
    c(text, paste(c(chrom_name, "1", chrom_end), collapse = "\t"))
}

# exporting chromosome lengths to txt file
write(mixedsort(text),
      paste(dir_name, "/chromFile.txt", sep = ""))

# check for --filter option
ifelse(is.null(opt$filter),
       records <-
         getFIX(vcf),
       records <- getFIX(vcf)[getFIX(vcf)[, 7] == "PASS", ])

# extracting the annotation data containing id, chr, positions and adding link to existing reference SNPs
colnames(records) <- NULL
text <- c()
for (row_count in 1:nrow(records)) {
  elem_name <- records[row_count, 3]
  chrom_name <- records[row_count, 1]
  elem_start <- records[row_count, 2]
  ifelse(
    row_count < nrow(records) &&
      chrom_name == records[row_count + 1, 1],
    elem_end <-
      as.numeric(records[row_count + 1, 2]) - 1,
    elem_end <-
      as.numeric(chrom_matrix[which(chrom_matrix == chrom_name, arr.ind = TRUE), 2][1]) - as.numeric(elem_start)
  )
  ifelse(
    is.na(elem_name),
    data <-
      NA,
    data <-
      paste("https://www.ncbi.nlm.nih.gov/snp/", elem_name, sep = "")
  )
  
  text <-
    c(text, paste(
      c(elem_name, chrom_name, elem_start, elem_end, data),
      collapse = "\t"
    ))
}

# exporting annotation data to txt file
write(text,
      paste(dir_name, "/annoFile.txt", sep = ""))

# generating chromosome visual graph
chrom_map <- chromoMap(
  paste(dir_name, "/chromFile.txt", sep = ""),
  paste(dir_name, "/annoFile.txt", sep = ""),
  hlinks = T
)

# exporting the graph to a html file
save_html(
  chrom_map,
  output_file,
  background = "white",
  libdir = dir_name,
  lang = "en"
)
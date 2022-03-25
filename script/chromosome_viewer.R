
# Tool that generates chromosome diagrams of human genomic variants listed in a VCF file.
# https://lakshay-anand.github.io/chromoMap/docs.html

library(vcfR)
library(chromoMap)
library(htmltools)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
input_file = args[1]
output_file = args[2]
dir_name = normalizePath(dirname(output_file))

# check for argument input
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).n", call. = FALSE)
} else if (length(args) == 1) {
  # default output file
  args[2] = "result.html"
}

#load of vcf file
vcf <-
  read.vcfR(input_file)

#extracting chromosome lengths by the contig tag in the vcf file
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

#exporting chromosome lengths to txt file
write(text,
      paste(dir_name, "/chromFile.txt", sep = ""))

#extracting the annotation data containing id, chr, positions and adding link to existing reference SNPs
records <- (getFIX(vcf))
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

#exporting annotation data to txt file
write(text,
      paste(dir_name, "/annoFile.txt", sep = ""))

#generating chromosome visual graph
chrom_map <- chromoMap(
  paste(dir_name, "/chromFile.txt", sep = ""),
  paste(dir_name, "/annoFile.txt", sep = ""),
  hlinks = T
)

#exporting the graph to a html file
save_html(
  chrom_map,
  output_file,
  background = "white",
  libdir = dir_name,
  lang = "en"
)
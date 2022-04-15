
# Tool that generates chromosome diagrams of human genomic variants listed in a VCF file.
# https://lakshay-anand.github.io/chromoMap/docs.html

#automatic install of packages if they are not installed already
list.of.packages <- c(
  "vcfR",
  "chromoMap",
  "htmltools",
  "gtools",
  "readr",
  "optparse",
  "foreach",
  "doParallel",
  "filelock"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}

#loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}

#get number of cores for parallel package
n.cores <- parallel::detectCores() - 1

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

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
  ),
  make_option(
    c("-r", "--reference"),
    type = "character",
    default = "hg19",
    help = "If contig tags aren't provided you can choose chromosome lengths between hg19 and hg38 reference genomes",
    metavar = "<hg19/hg38>"
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

# load of vcf file
vcf <-
  read.vcfR(input_file)

# extracting chromosome lengths by the contig tag in the vcf file
# https://www.ncbi.nlm.nih.gov/grc/human/data
chrom_list <- queryMETA(vcf, element = "contig")

#check if contig is present
if (length(chrom_list) == 0) {
  ifelse(
    opt$reference == "hg19",
    file <-
      read_lines(file = "data/hg19"),
    file <-
      read_lines(file = "data/GRCh38")
  )
  vcf_chroms <- unique(getCHROM(vcf))
  
  text <- c()
  for (chrom in vcf_chroms) {
    for (line in file) {
      if (grepl(chrom, line, fixed = TRUE)) {
        chrom_name <- substring(line, first = 0, last = 4)
        chrom_end <- substring(line, 8)
        text <-
          c(text, paste(c(line), collapse = ""))
      }
    }
  }
} else {
  text <- c()
  for (chrom in chrom_list) {
    chrom_name <- substring(chrom[1], 11)
    chrom_end <- substring(chrom[2], 8)
    
    text <-
      c(text, paste(c(chrom_name, "1", chrom_end), collapse = "\t"))
  }
}

# exporting chromosome lengths to txt file
write(mixedsort(unique(text)),
      paste(dir_name, "/chromFile.txt", sep = ""))

#check if theres a need to change chromosome notation (notation must be chr*)
if (!grepl("chr", unique(getCHROM(vcf))[1], fixed = TRUE)) {
  changed_notation_file <-
    system(
      paste(
        'awk \'{gsub(/^chr/,""); print}\'',
        paste(dir_name, "/chromFile.txt", sep = ""),
        sep = " "
      ),
      intern = TRUE
    )

  write(changed_notation_file,  paste(dir_name, "/chromFile.txt", sep = ""))
}

# check for --filter option
ifelse(is.null(opt$filter),
       records <-
         getFIX(vcf),
       records <- getFIX(vcf)[getFIX(vcf)[, 7] == "PASS",])

colnames(records) <- NULL
anno_file <-paste(dir_name, "/annoFile.txt", sep = "")
#invisible(file.remove(anno_file))

chrom_matrix <- matrix(, nrow = 1, ncol = 2)
chrom_matrix <- read.table(paste(dir_name, "/chromFile.txt", sep = ""))
chrom_matrix <- unique(chrom_matrix)
chrom_matrix$V2 <- NULL

# extracting the annotation data containing id, chr, positions and adding link to existing reference SNPs
foreach (row_count=1:nrow(records), .packages='filelock') %dopar% {
  line <- c()

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
  
  line <- paste(
    c(elem_name, chrom_name, elem_start, elem_end, data),
    collapse = "\t"
  )
  
  lck <- lock("/tmp/anno_file.lock")
  write(line,
        anno_file, append=TRUE)
  unlock(lck)
}

#stop cluster
parallel::stopCluster(cl = my.cluster)

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
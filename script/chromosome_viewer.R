
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
  "filelock",
  "tools"
)

new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]

if (length(new.packages) > 0) {
  install.packages(new.packages, dep = TRUE)
}

#loading packages
for (package.i in list.of.packages) {
  suppressPackageStartupMessages(library(package.i,
                                         character.only = TRUE))
}

#get number of cores for parallel package
n.cores <- parallel::detectCores() - 1

#create the cluster
my.cluster <- parallel::makeCluster(n.cores,
                                    type = "PSOCK")

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
    default = NULL,
    help = "If contig tags aren't provided you can choose chromosome lengths between hg19 and hg38 reference genomes",
    metavar = "<hg19/hg38>"
  ),
  make_option(
    c("-c", "--chromosome"),
    type = "character",
    default = NULL,
    help = "Filter by chromosomes, ex. <chr#:from-to,chr#:from-to,...>",
    metavar = ""
  ),
  make_option(
    c("-p", "--pathogen"),
    type = "character",
    action = "store_true",
    default = FALSE,
    help = "Filter by pathogenicity",
    metavar = ""
  ),
  make_option(
    c("-g", "--clnsig"),
    type = "character",
    default = NULL,
    help = "Filter by clinical significance, USE NUMBERS ONLY, ex. Uncertain - 0, Not provided - 1, Benign - 2, Likely benign - 3, Likely pathogenic - 4, Pathogenic - 5,
Drug-response related - 6, Histocompatibility-related - 7, Other - 255",
    metavar = ""
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

# check for filtering options issues
if (!is.null(opt$clnsig) && isTRUE(opt$pathogen))
  stop("You can't filter by pathogenicity and by other clinical significance at once",
       call. =
         FALSE)

# defining files from command
input_file = opt$input
output_file = opt$output
dir_name = normalizePath(dirname(output_file))

#TODO filtered files should be deleted after completion
# filtering options
# filtering by chromosomes and lengths
if (!is.null(opt$chromosome)) {
  if (file_ext(input_file) != "gz") {
    system(paste(
      c("bcftools view ", input_file, " -Oz -o ", input_file, ".gz"),
      collapse = ""
    ))
    system(paste(c("bcftools index ", input_file, ".gz"), collapse = ""))
  }
  
  file_name <- (tools::file_path_sans_ext(input_file))
  
  system(paste(
    c(
      "bcftools filter ",
      input_file,
      ".gz -r ",
      opt$chromosome,
      " > ",
      file_name,
      "_filtered.vcf"
    ),
    collapse = ""
  ))
  input_file <- paste(c(file_name, "_filtered.vcf"), collapse = "")
}

#TODO what if the vcf is not annotated?
#filtering by pathogenicity
if (isTRUE(opt$pathogen)) {
  file_name <- (tools::file_path_sans_ext(input_file))
  system(paste(
    c(
      "bcftools view -i 'INFO/CLNSIG==\"Pathogenic\"' ",
      input_file,
      " > ",
      file_name,
      "_path.vcf"
    ),
    collapse = ""
  ))
  input_file <- paste(c(file_name, "_path.vcf"), collapse = "")
}

#filtering by clinical significance
if (!is.null(opt$clnsig)) {
  file_name <- (tools::file_path_sans_ext(input_file))
  
  clnsig_options <- c()
  for (option in as.list(unlist(strsplit(opt$clnsig, ",")))) {
    switch(
      option,
      "0" = {
        clnsig_options <- paste(c(
          clnsig_options,
          "INFO/CLNSIG~\"Uncertain_significance\""
        ),
        collapse = " | ")
      },
      "1" = {
        clnsig_options <-
          paste(c(clnsig_options, "INFO/CLNSIG~\"Uncertain\""),
                collapse = " | ")
      },
      "2" = {
        clnsig_options <-
          paste(c(clnsig_options, "INFO/CLNSIG~\"Benign\""), collapse = " | ")
      },
      "3" = {
        clnsig_options <-
          paste(c(clnsig_options, "INFO/CLNSIG~\"Likely_benign\""),
                collapse = " | ")
      },
      "4" = {
        clnsig_options <-
          paste(c(clnsig_options, "INFO/CLNSIG~\"Likely_pathogenic\""),
                collapse = " | ")
      },
      "5" = {
        clnsig_options <-
          paste(c(clnsig_options, "INFO/CLNSIG~\"Pathogenic\""),
                collapse = " | ")
      },
      "6" = {
        clnsig_options <-
          paste(c(clnsig_options, "INFO/CLNSIG~\"Drug_response\""),
                collapse = " | ")
      },
      "7" = {
        clnsig_options <-
          paste(c(clnsig_options, "INFO/CLNSIG~\"Histocompatibility\""),
                collapse = " | ")
      },
      "255" = {
        clnsig_options <-
          paste(c(clnsig_options, "INFO/CLNSIG~\"Other\""),
                collapse = " | ")
      },
      {
        print('default')
      }
    )
  }
  
  system(paste(
    c(
      "bcftools view -i '",
      clnsig_options,
      "' ",
      input_file,
      " > ",
      file_name,
      "_clnsig.vcf"
    ),
    collapse = ""
  ))
  
  input_file <- paste(c(file_name, "_clnsig.vcf"), collapse = "")
}

# load of vcf file
vcf <-
  read.vcfR(input_file)

# extracting chromosome lengths by the contig tag in the vcf file
# https://www.ncbi.nlm.nih.gov/grc/human/data
chrom_list <- queryMETA(vcf, element = "contig")
variant_chroms <- unique(getCHROM(vcf))

#TODO optimization
#check if contig is present
if (length(chrom_list) == 0 ||
    !grepl("length", chrom_list, fixed = TRUE)) {
  if (is.null(opt$reference))
    stop("Contig is missing in header. You must specify reference lengths with -r",
         call. =
           FALSE)
  
  ifelse(
    opt$reference == "hg19",
    file <-
      read_lines(file = "data/hg19"),
    file <-
      read_lines(file = "data/GRCh38")
  )
  vcf_chroms <- unique(getCHROM(vcf))
  
  if (!grepl("chr", unique(getCHROM(vcf))[1], fixed = TRUE)) {
    file <- gsub(".*chr", "", file)
  }
  
  text <- c()
  for (chrom in vcf_chroms) {
    for (line in file) {
      if (chrom == sub("\\\t1.*", "", line)) {
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
    
    if (chrom_name %in% variant_chroms)
      text <-
      c(text, paste(c(chrom_name, "1", chrom_end), collapse = "\t"))
  }
}

# exporting chromosome lengths to txt file
write(mixedsort(unique(text)),
      paste(dir_name, "/chromFile.txt", sep = ""))

# check for --filter option
ifelse(
  is.null(opt$filter),
  records <-
    getFIX(vcf, getINFO = TRUE),
  records <-
    getFIX(vcf, getINFO = TRUE)[getFIX(vcf)[, 7] == "PASS", ]
)

colnames(records) <- NULL

anno_file <- paste(dir_name, "/annoFile.txt", sep = "")
if (file.exists(anno_file)) {
  file.remove(anno_file)
}

chrom_matrix <- matrix(, nrow = 1, ncol = 2)
chrom_matrix <-
  read.table(paste(dir_name, "/chromFile.txt", sep = ""))
chrom_matrix <- unique(chrom_matrix)
chrom_matrix$V2 <- NULL

# extracting the annotation data containing id, chr, positions and adding link to existing reference SNPs
foreach (row_count = 1:nrow(records), .packages = 'filelock') %dopar% {
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
  
  ifelse(is.na(elem_name),
         data <-
           NA,
         data <- ifelse(
           !grepl("rs", elem_name, fixed = TRUE),
           paste(
             "https://www.ncbi.nlm.nih.gov/clinvar/variation/",
             elem_name,
             sep = ""
           ),
           paste("https://www.ncbi.nlm.nih.gov/snp/", elem_name, sep = "")
         ))
  
  #declaring pathogenicity for chromoMap to categorize
  ifelse(
    grepl("Pathogenic", records[row_count, 8], fixed = TRUE),
    pathogenic <- "pathogenic",
    pathogenic <- "nonpathogen"
  )
  
  line <-
    paste(c(elem_name, chrom_name, elem_start, elem_end, pathogenic, data),
          collapse = "\t")
  
  invisible(lck <- lock("/tmp/anno_file.lock"))
  write(line,
        anno_file, append = TRUE)
  invisible(unlock(lck))
}

#stop cluster
parallel::stopCluster(cl = my.cluster)

if (!all(sapply(
  c("nonpathogen", "pathogenic"),
  grepl,
  readChar(anno_file, file.info(anno_file)$size)
)))
  chrom_map <- chromoMap(
    title = basename(input_file),
    paste(dir_name, "/chromFile.txt", sep = ""),
    paste(dir_name, "/annoFile.txt", sep = ""),
    hlinks = T
  ) else
  chrom_map <- chromoMap(
    title = basename(input_file),
    paste(dir_name, "/chromFile.txt", sep = ""),
    paste(dir_name, "/annoFile.txt", sep = ""),
    hlinks = T,
    data_based_color_map = TRUE,
    data_type = c("categorical"),
    legend = TRUE
  )

# exporting the graph to a html file
save_html(
  chrom_map,
  output_file,
  background = "white",
  libdir = dir_name,
  lang = "en"
)
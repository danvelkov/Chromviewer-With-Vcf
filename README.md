# Chromviewer-With-Vcf
Tool that generates chromosome diagrams of human genomic variants listed in a VCF file.

Requires [R version 4.1.3 (2022-03-10)](https://www.r-project.org/) or newer version
and [bcftools-1.15.1](https://samtools.github.io/bcftools/) or newer version

## How to run
`Rscript ./script/chromosome_viewer.R -i <input.VCF> [-o <output.HTML>] [options]`

Options include

	-f, --filter
		Use only variants with value “PASS” in FILTER field
  
	-r <hg19|hg38>, --reference
		If contig tags with lenghts aren't provided you can choose chromosome lengths between hg19 and hg38 reference genomes

	-c , --chromosome
		Filter by chromosomes, ex. <chr#:from-to,chr#:from-to,...>, IMPORTNANT: notation should be same with vcf file

	-p, --pathogen
		Filter by pathogenicity

	-g , --clnsig
		Filter by clinical significance, USE NUMBERS ONLY, ex. Uncertain - 0, Not provided - 1, Benign - 2,
    Likely benign - 3, Likely pathogenic - 4, Pathogenic - 5, Drug-response related - 6,
    Histocompatibility-related - 7, Other - 255, ex. <-g 0,3,4>


There are sample vcf files in the data folder for test purposes

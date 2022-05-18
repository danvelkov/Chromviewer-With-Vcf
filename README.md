# Chromviewer-With-Vcf
Tool that generates chromosome diagrams of human genomic variants listed in a VCF file.

Requires [R version 4.1.3 (2022-03-10)](https://www.r-project.org/) or newer version
and [bcftools-1.15.1](https://samtools.github.io/bcftools/) or newer version

## How to run
**Important**: Firstly you have to run to install default askpass, incase your Linux distro doesn't provide the ssh default one

`sudo apt-get install ssh-askpass-gnome ssh-askpass`

To run script, just run the bash file in ChromViewer directory

`cd ChromViewer`

`./chromosome_viewer.sh -i <input.VCF> [-o <output.HTML>] [options]`

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

---

### Examples
- Visualise chromosome map with filter for chromosomes 1, 4 and X from the input file *test.vcf* in */data*

**Important**: when filtering by chromosmes ***-c*** flag should be called with the chromosome notation from the input file (ex. 1,2,3.. if #CHROM is 1,2,3 or chr1.. if #CHROM is chr1.. etc.)

`./chromosome_viewer.sh -i ./data/test.vcf -o ../output/result.html -c chr1,chr4,chrX`

Output:

[![non-anno-output.png](https://i.postimg.cc/4NvRrSnn/non-anno-output.png)](https://postimg.cc/3dwqG9XH)

- Visualise chromosome map from the input file *annotated_data.vcf* in */data*

*Note*: When contig tags are missing ***-r*** flag should be used to specify reference genome lengths

`./chromosome_viewer.sh -i ./data/annotated_data.vcf -o ../output/result.html -r hg38`

Output:

[![anno-output.png](https://i.postimg.cc/28k0Sc21/anno-output.png)](https://postimg.cc/8sYRyH8G)

- Visualise chromosome map from input file *annotated_data.vcf* in */data* only with variants which clinical significance (CLNSIG) is Bening and Likely Bening

*Note*: You can't use ***-g*** and ***-p*** flags simultaneously
 
`./chromosome_viewer.sh -i ./data/annotated_data.vcf -o ../output/result.html -r hg38 -g 2,3`

Output:

[![clnsig-output.png](https://i.postimg.cc/FF8X7CPs/clnsig-output.png)](https://postimg.cc/JtN2FxV9)

When you hover over a painted region it will show variant count and display their Id's
Depending on wheter it's annotated for clinical significance or not it will lead to NIH's ClinVar database or NIH's dbSNP database

In the output directory except exported diagram files there will be files with chromosome lengths and variant info required for <a href="https://lakshay-anand.github.io/chromoMap/docs.html" target="_blank">chromoMap</a>. They can be used for additional data presentation and analisys

---

If any issues or questions arise while using this script please feel free to mention them in Issues tab



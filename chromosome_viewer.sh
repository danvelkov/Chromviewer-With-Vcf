#!/bin/bash
#check if bcftools is installed
if ! command -v bcftools &> /dev/null
then
  echo "Installing bcftools"
  apt-get update;
  apt-get install gcc;
  apt-get install make;
  apt-get install libbz2-dev;
  apt-get install zlib1g-dev;
  apt-get install libncurses5-dev;
  apt-get install libncursesw5-dev;
  apt-get install liblzma-dev;
  cd ~
  wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
  tar -vxjf bcftools-1.9.tar.bz2
  cd bcftools-1.9
  make
  export PATH="$PATH:/usr/bin/bcftools-1.9"
  source ~/.profile
  
  apt install bcftools
fi

#check if r-base is installed
if [ $(dpkg-query -W -f='${Status}' r-base 2>/dev/null | grep -c "ok installed") -eq 0 ];
then
  echo "Installing R-base"
  apt install dirmngr gnupg apt-transport-https ca-certificates software-properties-common;
  apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9;
  add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/';
  apt install r-base;
  apt install build-essential;
fi

#when running Rscript for the first time packages will be installed
#TODO add variables from bash script to be added to Rscript call
Rscript ./script/chromosome_viewer.R




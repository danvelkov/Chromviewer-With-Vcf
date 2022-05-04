#!/bin/bash
#check if bcftools is installed
if ! command -v bcftools &> /dev/null
then
  echo "Installing bcftools"
  apt-get -y update;
  apt-get -y install gcc;
  apt-get -y install make;
  apt-get -y install libbz2-dev;
  apt-get -y install zlib1g-dev;
  apt-get -y install libncurses5-dev;
  apt-get -y install libncursesw5-dev;
  apt-get -y install liblzma-dev;
  cd ~
  wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
  tar -vxjf bcftools-1.9.tar.bz2
  cd bcftools-1.9
  make
  export PATH="$PATH:/usr/bin/bcftools-1.9"
  source ~/.profile
  
  apt -y install bcftools
fi

#check if r-base is installed
if [ $(dpkg-query -W -f='${Status}' r-base 2>/dev/null | grep -c "ok installed") -eq 0 ];
then
  echo "Installing R-base"
  apt -y install dirmngr gnupg apt-transport-https ca-certificates software-properties-common;
  apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9;
  add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/';
  apt -y install r-base;
  apt -y install build-essential;
fi

#call Rscript with input flags
Rscript ./script/chromosome_viewer.R "$@"

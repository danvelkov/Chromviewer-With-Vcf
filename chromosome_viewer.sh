#!/bin/bash
#check if bcftools is installed
if ! command -v bcftools &> /dev/null
then
  echo "Installing bcftools"
  sudo apt-get -y update;
  sudo apt-get -y install gcc;
  sudo apt-get -y install make;
  sudo apt-get -y install libbz2-dev;
  sudo apt-get -y install zlib1g-dev;
  sudo apt-get -y install libncurses5-dev;
  sudo apt-get -y install libncursesw5-dev;
  sudo apt-get -y install liblzma-dev;
  cd ~
  wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
  tar -vxjf bcftools-1.9.tar.bz2
  cd bcftools-1.9
  make
  export PATH="$PATH:/usr/bin/bcftools-1.9"
  source ~/.profile
  
  sudo apt -y install bcftools
fi

#check if r-base is installed
if [ $(dpkg-query -W -f='${Status}' r-base 2>/dev/null | grep -c "ok installed") -eq 0 ];
then
  echo "Installing R-base"
  sudo apt -y install dirmngr gnupg apt-transport-https ca-certificates software-properties-common;
  sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9;
  sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/';
  sudo apt -y install r-base;
  sudo apt -y install build-essential;
fi

#call to install Rscript
Rscript --version

#call Rscript with input flags
Rscript ./script/chromosome_viewer.R "$@"

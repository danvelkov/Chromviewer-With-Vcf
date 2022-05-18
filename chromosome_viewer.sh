#!/bin/bash
#askpass configuration
#if error occur, please read README.md
export SUDO_ASKPASS=/usr/bin/ssh-askpass

#start directory variable
STDIR=`pwd`

#check if bcftools is installed
if ! command -v bcftools &> /dev/null
then
  echo "Installing bcftools"
  sudo -A apt-get -y update;
  sudo -A apt-get -y install gcc;
  sudo -A apt-get -y install make;
  sudo -A apt-get -y install libbz2-dev;
  sudo -A apt-get -y install zlib1g-dev;
  sudo -A apt-get -y install libncurses5-dev;
  sudo -A apt-get -y install libncursesw5-dev;
  sudo -A apt-get -y install liblzma-dev;
  cd ~
  wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
  tar -vxjf bcftools-1.9.tar.bz2
  cd bcftools-1.9
  make
  export PATH="$PATH:/usr/bin/bcftools-1.9"
  source ~/.profile
  
  sudo -A apt -y install bcftools
fi

#check if r-base is installed
if [ $(dpkg-query -W -f='${Status}' r-base 2>/dev/null | grep -c "ok installed") -eq 0 ];
then
  echo "Installing R-base"
  sudo -A apt -y install dirmngr gnupg apt-transport-https ca-certificates software-properties-common;
  sudo -A apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9;
  sudo -A add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/';
  sudo -A apt -y install r-base;
  sudo -A apt -y install build-essential;
fi

#call to check Rscript
Rscript --version

#changing back to start directory
cd $STDIR

#call Rscript with input flags
Rscript ./script/chromosome_viewer.R "$@"

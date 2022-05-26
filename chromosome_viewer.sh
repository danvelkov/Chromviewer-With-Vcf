#!/bin/bash
if [[ "$OSTYPE" =~ ^linux ]]; then
  #askpass configuration
  #if error occur, please read README.md
  export SUDO_ASKPASS=/usr/bin/ssh-askpass

  #start directory variable
  STDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
  
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

  #changing back to start directory
  cd $STDIR
fi

#call to check Rscript
Rscript --version



sleep 5

#call Rscript with input flags
Rscript ./script/chromosome_viewer.R "$@"

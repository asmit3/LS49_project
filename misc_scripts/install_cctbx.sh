#!/bin/bash

# This file will install the necessary files for running psana-cctbx.xfel, build the required binaries, and set up the paths for their use.

#Fetch miniconda installer
function fetchConda(){
  if [ ! -e Miniconda2-latest-Linux-x86_64.sh ]; then
    echo "Conda installer not found. Acquiring."
    wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh --no-check-certificate
  fi
}

#install miniconda and necessary packages
function condaEnvSetup(){
  chmod +x ./Miniconda2-latest-Linux-x86_64.sh
  ./Miniconda2-latest-Linux-x86_64.sh -b -p ${PWD}/newconda 
  #conda update -y conda
  source ${PWD}/newconda/etc/profile.d/conda.sh
  conda create -y -n $1 #Pass argument for environment name
  conda activate $1 #Activate said environment
  conda install -y --channel lcls-rhel${2} psana-conda #Pass in RHEL version number: {5,6,7}
  conda install -y future h5py mpich2 wxpython pillow libtiff mock pytest jinja2 scikit-learn tabulate tqdm
  conda install -y --channel conda-forge orderedset
  python -m pip install procrunner msgpack blosc
}

#set the necessary data directory environemtn variables and paths
#function dataDirSetup(){
  #export PERM=$1; mkdir $PERM
  #mkdir -p $PERM/psdm/data/ExpNameDb;
  #rsync -t psexport.slac.stanford.edu:/reg/g/psdm/data/ExpNameDb/experiment-db.dat $PERM/psdm/data/ExpNameDb/
#}
function psanaEnvVar(){ #May not be necessary
  export SIT_DATA=$PERM/psdm/data
  export SIT_ROOT=$PERM/psdm/data
  export SIT_PSDM_DATA=$SIT_ROOT
}

#pass cci username, github username and number of compilation cores as args
function setupcctbx(){
  #if [ ! -e $PERM/cctbx.xfel ]; then
  #  echo "cctbx installer not found. Acquiring."
  #  mkdir $PERM/cctbx.xfel;
  #  cd $PERM/cctbx.xfel
  wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py --no-check-certificate --no-check-certificate
  #fi
  python bootstrap.py hot update --builder=dials
  python bootstrap.py build --builder=dials --with-python=`which python` --nproc=$1
}

##########################################
#                 main
##########################################

if [ $# -lt 3 ]
  then
    echo "Not enough arguments supplied"
    echo "Please specify the following arguments to correctly run the installation:
    {environmentName} {ELversion} {cciusername} {githubusername} {compilerCores}"
    echo "For an installation to ELversion=6, with  environmentName=myEnv, cciusername=me githubusername=me compilerCores=32 the command will be run as"
    echo "./install.sh myEnv 7 32"
    exit
fi
fetchConda
condaEnvSetup $1 $2 #for Centos 6 and a simple env name
#dataDirSetup $(pwd) #use current directory as data dir
#exit()
#psanaEnvVar
setupcctbx $3

Steps to install docker image for LS49 project

## On your Mac
Create a new folder on your mac - I call it dials and cd into it

1. First pull down miniconda & bootstrap.py
  a. curl -L -O https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py
  b. curl -L -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

2. Run hot, update for --builder=dials
  a. python bootstrap.py hot update --builder=dials

3. Clone exafel_project, LS49_project and LS49 in the modules directory
  a. cd modules
  b. git clone https://github.com/ExaFEL/exafel_project.git
  c. git clone https://github.com/asmit3/LS49_project.git
  d. git clone https://github.com/nksauter/LS49.git

4. Make sure you have checked out the right branches on dials, dxtbx and cctbx_project. These are specific for LS49 work
  a. cd dials; git checkout XTC_timestamp_identifier; cd ..
  b. cd dxtbx; git checkout XTC_timestamp_identifier; cd ..
  c. cd cctbx_project; git checkout LS49_exp; cd ..

5. Compile the docker image 
  a. docker build . [Note the .]

6. Upload the image to docker hub. My docker image is at asmit3/ls49:latest
  a. docker tag <imageID> asmit3/ls49:latest [imageID is the ID produced at the end of docker build] 
  b. docker login
  c. docker push asmit3/ls49:latest

## Now on Cori

1. Now you need to pull down the docker image from docker hub
  a. shifterimg pull asmit3/ls49:latest

2. Show that the pulled and converted docker image is installed to shifter
  a. shifterimg images | grep ls49

3. Acquire 30 minute KNL node and load the docker image into shifter
   Might have to use module=mpich-cle6
  a. salloc -N 1 -p debug --image=docker:asmit3/ls49:latest -t 00:30:00 -C knl -A m2859

4. All files from the docker build are placed at /   (e.g. /modules /build)
   Source the file paths to access the commands
  a. source /build/setpaths.sh


  

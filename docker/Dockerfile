# Designed for use with XTC data streams on NERSC's Cori
#
# Run bootstrap.py with hot and update first to get the modules directory.
# Next, add the Python2 Miniconda installer into the directory conttainer this file.
# Acquire the sources for MPICH3.2.1 and MPI4PY 3.0, untar them and place them here also.
# Run `docker build .` and wait for the image to build (may take some time)

FROM centos:centos7
MAINTAINER Asmit Bhowmick <abhowmick@lbl.gov>

RUN yum clean all &&  yum -y install bzip2.x86_64 libgomp.x86_64 telnet.x86_64 gcc-c++ gcc-gfortran && yum groupinstall "Development Tools" -y

# https://repo.continuum.io/miniconda/
ADD Miniconda3-latest-Linux-x86_64.sh miniconda.sh
ADD bootstrap.py bootstrap.py
RUN chmod +x miniconda.sh && echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && /bin/bash miniconda.sh -b -p /opt/conda && rm miniconda.sh
ENV PATH /opt/conda/bin:$PATH
ENV LD_LIBRARY_PATH /opt/conda/lib:$LD_LIBRARY_PATH

# psana-conda
ADD manifest_conda manifest_conda
RUN conda update -y conda && conda create -n psana_env --file manifest_conda && rm -rf /opt/conda/lib/python2.7/site-packages/numexpr-2.6.2-py2.7.egg-info

# cctbx
ADD modules modules
RUN python bootstrap.py build --builder=dials --with-python=/opt/conda/envs/psana_env/bin/python --nproc=4

RUN cd /build/ && source /build/setpaths.sh && libtbx.configure exafel_project LS49_project && make
#Build Mpi4py from scratch against mpich
ADD mpich-3.2.1 /mpich_src
ADD mpi4py-3.0.0 /mpi4py_src
#
RUN cd /mpich_src && ./configure && /usr/bin/make -j4 && /usr/bin/make install && /usr/bin/make clean && rm -rf /mpich_src && cd /mpi4py_src && /opt/conda/envs/psana_env/bin/python setup.py build && /opt/conda/envs/psana_env/bin/python setup.py install && rm -rf /mpi4py_src && /sbin/ldconfig

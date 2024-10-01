git clone https://github.com/epfl-ecps/channelflow.git
sudo apt-get update
sudo apt-get install libopenmpi-dev
sudo apt install libeigen3-dev
sudo apt-get install libfftw3-dev
sudo apt-get install libfftw3-mpi-dev
sudo apt install netcdf-bin libnetcdff-dev
# rm -r build 
mkdir build
cd build
cmake ../channelflow -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx -DCMAKE_BUILD_TYPE=release -DCMAKE_CXX_FLAGS_RELEASE:STRING=" -lfftw3 -lm " -DCMAKE_INSTALL_PREFIX=~/channelflow-2.0 -DWITH_HDF5CXX=On
# -DCMAKE_CXX_FLAGS_RELEASE:STRING=" -fPIC " (or -lm -lfftw3 if get errors regarding "undefined reference to "fftw_..."")
make -j4
sudo make install
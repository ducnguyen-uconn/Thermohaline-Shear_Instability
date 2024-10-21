CXX = /usr/bin/mpicxx
TYPE = release
PREFIX = /user/local/
CFLAGS = " -fPIC -lfftw3 -lm "

DDC = ON
NSOLVER = ON
HDF5 = ON

.PHONY: clean build
update:
	sudo apt-get update
	sudo apt-get install libopenmpi-dev
	sudo apt install libeigen3-dev
	sudo apt-get install libfftw3-dev
	sudo apt-get install libfftw3-mpi-dev
	sudo apt install netcdf-bin libnetcdff-dev

clone:
	rm -rf channelflow
	git clone https://github.com/epfl-ecps/channelflow.git

cloneilc:
	rm -rf channelflow
	git clone --single-branch --branch module/ilc https://github.com/epfl-ecps/channelflow.git

build:
	mkdir -p build
	cd build;\
	cmake ../channelflow -DWITH_DDC=$(DDC) -DWITH_NSOLVER=$(NSOLVER) -DCMAKE_CXX_COMPILER=$(CXX) -DCMAKE_BUILD_TYPE=$(TYPE) -DCMAKE_INSTALL_PREFIX=$(PREFIX) -DCMAKE_CXX_FLAGS_RELEASE:STRING=$(CFLAGS) -DWITH_HDF5CXX=$(HDF5);\
	make -j24

run:
	cp ./CMakeLists.txt ./channelflow/CMakeLists.txt
	mkdir -p ./channelflow/modules/
	rm -rf ./channelflow/modules/ddc
	cp -r ./ddc ./channelflow/modules/ddc
	make build
	rm -rf data
	mpiexec -n 16 ./build/modules/ddc/examples/2d_diffusive_convection


clean:
	rm -rf build
	rm -rf data

allclean:
	rm -rf channelflow
	rm -rf build
	rm -rf data

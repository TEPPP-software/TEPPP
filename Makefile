#################################################################
########################   MakeVars   ###########################
#################################################################

CXX=g++

MPICXX=mpic++

CXX_OPT= -g -I "./include" -Wall -Wextra -O2 -std=c++17

LD_LIB=

LD_OPT=

MKDIR=mkdir -p ./obj

MKDIR_BUILD=mkdir -p ./build

TARGET=convertor ./build/jones ./build/jones_scan ./build/lk ./build/lk_scan ./build/periodic_lk ./build/periodic_wr ./build/wr ./build/wr_scan ./build/jones_mpi ./build/lk_mpi ./build/lk_scan_mpi ./build/periodic_lk_mpi ./build/periodic_wr_mpi ./build/wr_mpi ./build/wr_scan_mpi

DCD_TARGET=convertor

TEMP_TARGET=convertor ./build/jones

SERIAL_TARGET=convertor ./build/jones ./build/jones_scan ./build/lk ./build/lk_scan ./build/periodic_lk ./build/periodic_wr ./build/wr ./build/wr_scan

MPI_TARGET=convertor ./build/jones_mpi ./build/lk_mpi ./build/lk_scan_mpi ./build/periodic_lk_mpi ./build/periodic_wr_mpi ./build/wr_mpi ./build/wr_scan_mpi


#################################################################
#################################################################

all:$(TARGET)
	@echo "=========================================================="
	@echo "Compilation Success"
	@echo "=========================================================="

serial:$(SERIAL_TARGET)
	@echo "=========================================================="
	@echo "Serial Compilation Success"
	@echo "=========================================================="

mpi:$(MPI_TARGET)
	@echo "=========================================================="
	@echo "Parallel Compilation Success"
	@echo "=========================================================="

./convertor:./main/convertor.cpp
	@echo "Convertor Compilation Success"
	$(CXX) $? -o $@

./build/jones:./main/jones.cpp
	@$(MKDIR_BUILD)
	$(CXX) $? -o $@
./build/jones_scan:./main/jones_scan.cpp
	@$(MKDIR_BUILD)
	$(CXX) $? -o $@
./build/lk:main/lk.cpp
	@$(MKDIR_BUILD)
	$(CXX) $? -o $@
./build/lk_scan:main/lk_scan.cpp
	@$(MKDIR_BUILD)
	$(CXX) $? -o $@
./build/periodic_lk:main/periodic_lk.cpp
	@$(MKDIR_BUILD)
	$(CXX) $? -o $@
./build/periodic_wr:main/periodic_wr.cpp
	@$(MKDIR_BUILD)
	$(CXX) $? -o $@
./build/wr:main/wr.cpp
	@$(MKDIR_BUILD)
	$(CXX) $? -o $@
./build/wr_scan:main/wr_scan.cpp
	@$(MKDIR_BUILD)
	$(CXX) $? -o $@
./build/jones_mpi:main/jones_mpi.cpp
	@$(MKDIR_BUILD)
	$(MPICXX) $? -o $@
./build/lk_mpi:main/lk_mpi.cpp
	@$(MKDIR_BUILD)
	$(MPICXX) $? -o $@
./build/lk_scan_mpi:main/lk_scan_mpi.cpp
	@$(MKDIR_BUILD)
	$(MPICXX) $? -o $@
./build/periodic_lk_mpi:main/periodic_lk_mpi.cpp
	@$(MKDIR_BUILD)
	$(MPICXX) $? -o $@
./build/periodic_wr_mpi:main/periodic_wr_mpi.cpp
	@$(MKDIR_BUILD)
	$(MPICXX) $? -o $@
./build/wr_mpi:main/wr_mpi.cpp
	@$(MKDIR_BUILD)
	$(MPICXX) $? -o $@
./build/wr_scan_mpi:main/wr_scan_mpi.cpp
	@$(MKDIR_BUILD)
	$(MPICXX) $? -o $@

# Remove all object files, executables, and the converted directory
# See the following stackoverflow posts for more information on the loops below:
# 1: https://stackoverflow.com/questions/1490949/how-to-write-loop-in-a-makefile
# 2: https://stackoverflow.com/questions/26564825/what-is-the-meaning-of-a-double-dollar-sign-in-bash-makefile
clean:
	rm -rfv ./build ./converted ./convertor ./output;

format:
	python ./scripts/clang-format-diff.py

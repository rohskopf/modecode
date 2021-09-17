.SUFFIXES: .h .cpp
#-----------------------------------------------
# Makefile for the program GEO.
# Please modify the variables properly.
# We recommend to use gcc >= 4.8 with openMPI >= 1.8
#-----------------------------------------------

# Enter the path to your LAMMPS directory
#LAMMPSDIR = /home/drew/lammps-5Jun19/src
LAMMPSDIR = /home/drew/lammps/src
#LAMMPSDIR = /home/rohskopf/lammps/src
#USE DEEPMD LAMMPS DIR IF YOU WANNA DO FC ON NNP. 
#LAMMPSDIR = /home/rohskopf/deepmd_mod/lammps/src

# Enter path to the shared library (.so) file
#LAMMPSLIB = /home/drew/lammps-5Jun19/src/liblammps_g++_openmpi.so
LAMMPSLIB = /home/drew/lammps/src/liblammps_g++_openmpi.so
#LAMMPSLIB = /home/rohskopf/lammps/src/liblammps_g++_openmpi.so
#LAMMPSLIB = /home/rohskopf/deepmd_mod/lammps/src/liblammps_g++_openmpi.so

#-----------------------------------------------

CXX = mpicxx
CXXFLAGS = -std=c++11 -Wno-write-strings

INCLUDE = -I${LAMMPSDIR}

CXXL =
LDFLAGS = 

LAPACK =
LIBS =

#-----------------------------------------------
# General rules 
#-----------------------------------------------

.cpp.o:
	${CXX} ${CXXFLAGS} ${INCLUDE} -c $<

PROG = modecode

CXXSRC= mc.cpp main.cpp mem.cpp in.cpp verify.cpp ifc2mcc.cpp asr.cpp compute.cpp tic.cpp qhgk.cpp visualize.cpp

OBJS= ${CXXSRC:.cpp=.o}

default: modecode

all: ${PROG}

modecode: ${OBJS}
	${CXX} ${CXXFLAGS} -o $@ ${OBJS} ${LIBS} ${LAMMPSLIB}

clean:
	rm -f ${OBJS}

.PHONY: clean

depend:
	gcc -MM *.cpp > .depend

-include .depend

################################################################################
# $Id: Makefile.mrblast.app 146491 2008-11-26 13:57:23Z camacho $
################################################################################

APP = mrblast
SRC = mrblast


### MPI ########################################################################
CC = mpicc
## Open MPI
#CXX = mpic++
## MVAPICH
CXX = mpicxx

MPI_COMPILE_FLAGS = $(shell mpic++ --showme:compile)
MPI_LINK_FLAGS = $(shell mpic++ --showme:link)


### dir setting ################################################################
## For me
MRMPI_HOME = /home/ssul/work/distros2/ncbi_cxx/ncbi_cxx--Jun_15_2010/src/app/mrblast/mrmpi
BOOST_HOME = /home/ssul/work/packages2

## For Ranger
#MRMPI_HOME = /work/01471/ssul/work/distros3/ncbi_cxx/ncbi_cxx--Jun_15_2010/src/app/mr-mpi-blast/mrmpi
#BOOST_HOME = /opt/apps/gcc4_4/boost/1.39.0
#BOOST_HOME = /work/01471/ssul/work/packages3


### MR-MPI Lib #################################################################
#MRMPI_USRLIB = -lmrmpi
#MRMPI_USRLIB = -L/home/ssul/work/distros2/ncbi_cxx/ncbi_cxx--Jun_15_2010/src/app/mrblast/mrmpi -lmrmpi

## For Ranger
#MRMPI_USRLIB = -L$(MRMPI_HOME) -lmrmpi_mpicxx

## For me
MRMPI_USRLIB = -L$(MRMPI_HOME) -lmrmpi_mpicc


### Boost Lib #################################################################
## For Ranger
#BOOST_INCLUDE = -I$(BOOST_HOME)/include/boost-1_39
#BOOST_USRLIB = -L$(BOOST_HOME)/lib -lboost_program_options-gcc44-mt

#BOOST_INCLUDE = -I$(BOOST_HOME)/include
#BOOST_USRLIB = -L$(BOOST_HOME)/lib -lboost_program_options -lboost_iostreams -lboost_filesystem

## For me: boost_1_45_0
BOOST_INCLUDE = -I$(BOOST_HOME)/include
BOOST_USRLIB = -L$(BOOST_HOME)/lib -lboost_program_options -lboost_iostreams -lboost_filesystem
### Boost.Log
#BOOST_USRLIB = -L$(BOOST_HOME)/lib -lboost_program_options -lboost_log -lboost_log_setup






CXXFLAGS = $(ORIG_CXXFLAGS) $(MPI_COMPILE_FLAGS)
ORIG_LIBS = $(MRMPI_USRLIB) $(BOOST_USRLIB)
LIBS = $(MPI_LINK_FLAGS) $(ORIG_LIBS)
CPPFLAGS = $(ORIG_CPPFLAGS) $(BOOST_INCLUDE)

 
# new_project.sh will copy everything in the following block to any
# Makefile.*_app generated from this sample project.  Do not change
# the lines reading "### BEGIN/END COPIED SETTINGS" in any way.

###############################################################################
### BEGIN COPIED SETTINGS
LIB_ = $(BLAST_INPUT_LIBS) ncbi_xloader_blastdb_rmt $(BLAST_LIBS) $(OBJMGR_LIBS)
LIB = $(LIB_:%=%$(STATIC)) 
LIBS = $(CMPRS_LIBS) $(NETWORK_LIBS) $(PCRE_LIBS) $(DL_LIBS) $(ORIG_LIBS) 

# These settings are necessary for optimized WorkShop builds, due to
# BLAST's own use of them.
CXXFLAGS = $(FAST_CXXFLAGS) $(ORIG_CXXFLAGS) 

LDFLAGS = $(FAST_LDFLAGS)

REQUIRES = objects -Cygwin

## These settings are necessary for optimized WorkShop builds, due to
## BLAST's own use of them.
#CXXFLAGS = $(FAST_CXXFLAGS) $(ORIG_CXXFLAGS) 
#LDFLAGS = $(FAST_LDFLAGS)
#REQUIRES = objects -Cygwin LIBXML LIBXSLT
###############################################################################

### END COPIED SETTINGS



################################################################################
# $Id: Makefile.mrblast.app 2012-01-18 ssul $
################################################################################

APP = mrblast
SRC = mrblast

## -------------------------------------------------------------------
## To build mr-mpi-blast, a MPI stack should be installed, for instance, Open
## MPI or MVAPICH. Also the MapReduce MPI library link path should be set.
##
## If you use Open MPI,
CC = mpicc
CXX = mpic++
MRMPI_USRLIB = -L../../../../src/app/mr-mpi-blast/mrmpi -lmrmpi_mpicc
##
## If you use MVAPICH,
#CC = mpicxx
#CXX = mpicxx
#MRMPI_USRLIB = -L../../../../src/app/mr-mpi-blast/mrmpi -lmrmpi_mpicxx

MPI_COMPILE_FLAGS = $(shell mpic++ --showme:compile)
MPI_LINK_FLAGS = $(shell mpic++ --showme:link)

## -------------------------------------------------------------------
## Also mr-mpi-blast implementation uses the Boost library. Check 'BOOST_HOME'
## on the system you are using and set it appropriately.
## The 'TACC_BOOST_DIR' is 'BOOST_HOME' in XSEDE Ranger cluster system.
##
#BOOST_HOME = $(TACC_BOOST_DIR)
BOOST_HOME = /home/ssul/work/packages2
BOOST_INCLUDE = -I$(BOOST_HOME)/include
BOOST_USRLIB = -L$(BOOST_HOME)/lib -lboost_program_options -lboost_iostreams -lboost_filesystem

CXXFLAGS = $(ORIG_CXXFLAGS) $(MPI_COMPILE_FLAGS)
ORIG_LIBS = $(MRMPI_USRLIB) $(BOOST_USRLIB)
LIBS = $(MPI_LINK_FLAGS) $(ORIG_LIBS)
CPPFLAGS = $(ORIG_CPPFLAGS) $(BOOST_INCLUDE)
 
dummy := $(shell date > last_build_timestamp.log) 

 
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



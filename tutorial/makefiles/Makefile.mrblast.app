################################################################################
# $Id: Makefile.mrblast.app 146491 2008-11-26 13:57:23Z camacho $
################################################################################

APP = mrblast
SRC = mrblast blast_app_util


### MPI ########################################################################
# For MVAPICH MPI
CC = mpicxx
CXX = mpicxx

MPI_COMPILE_FLAGS = $(shell mpic++ --showme:compile)
MPI_LINK_FLAGS = $(shell mpic++ --showme:link)


### env setting ################################################################
MRMPI_HOME=$(WORKING_DIR)/ncbi_cxx--7_0_0/src/app/mr-mpi-blast/mrmpi
BOOST_HOME=$(TACC_BOOST_DIR)
MRMPI_USRLIB = -L$(MRMPI_HOME) -lmrmpi_mpicxx


### Boost Lib #################################################################
BOOST_INCLUDE=-I$(BOOST_HOME)/include
BOOST_USRLIB=-L$(BOOST_HOME)/lib -lboost_program_options -lboost_iostreams -lboost_filesystem


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



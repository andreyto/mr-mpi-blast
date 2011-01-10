#################################
# $Id: Makefile.blast_demo.app 146491 2008-11-26 13:57:23Z camacho $
#################################

APP = mrblast
SRC = mrblast 

### MPI and MR-MPI #####################################################
CC = mpicc

### OPENMPI
#CXX = mpic++
### MVAPICH
CXX = mpicxx

MPI_COMPILE_FLAGS = $(shell mpic++ --showme:compile)
MPI_LINK_FLAGS = $(shell mpic++ --showme:link)
MRMPI_USRLIB = -L/home/ssul/work/distros2/ncbi_cxx/ncbi_cxx--Jun_15_2010/src/app/mrblast/mrmpi -lmrmpi
#MRMPI_USRLIB = -lmrmpi
BOOST_USRLIB = -lboost_program_options

### RANGER
#/opt/apps/gcc4_4/boost/1.39.0/include/boost-1_39
 
CXXFLAGS = $(ORIG_CXXFLAGS) $(MPI_COMPILE_FLAGS) $(BOOST_INCLUDE)
ORIG_LIBS = $(MRMPI_USRLIB) $(BOOST_USRLIB)
LIBS = $(MPI_LINK_FLAGS) $(ORIG_LIBS) 
 
### XML parsing ########################################################
#CPPFLAGS = $(LIBXML_INCLUDE) $(LIBXSLT_INCLUDE) $(ORIG_CPPFLAGS) 
########################################################################

### Boost include 
### if Boost is not installed /usr/inlcude, change this location!!!
#BOOST_INCLUDE = -I/usr/include/boost
#CPPFLAGS = $(BOOST_INCLUDE) $(ORIG_CPPFLAGS) 
########################################################################

# new_project.sh will copy everything in the following block to any
# Makefile.*_app generated from this sample project.  Do not change
# the lines reading "### BEGIN/END COPIED SETTINGS" in any way.

### BEGIN COPIED SETTINGS
LIB_ = $(BLAST_INPUT_LIBS) ncbi_xloader_blastdb_rmt $(BLAST_LIBS) $(OBJMGR_LIBS)
LIB = $(LIB_:%=%$(STATIC)) 
LIBS = $(CMPRS_LIBS) $(NETWORK_LIBS) $(PCRE_LIBS) $(DL_LIBS) $(ORIG_LIBS) 

# These settings are necessary for optimized WorkShop builds, due to
# BLAST's own use of them.
CXXFLAGS = $(FAST_CXXFLAGS) $(ORIG_CXXFLAGS) 

LDFLAGS = $(FAST_LDFLAGS)

REQUIRES = objects -Cygwin

### XML parsing ########################################################
#LIB_ = $(BLAST_INPUT_LIBS) ncbi_xloader_blastdb_rmt $(BLAST_LIBS) $(OBJMGR_LIBS)
#LIB = $(LIB_:%=%$(STATIC)) xmlwrapp xncbi
#LIBS = $(CMPRS_LIBS) $(NETWORK_LIBS) $(PCRE_LIBS) $(DL_LIBS) $(LIBXML_LIBS) $(LIBXSLT_LIBS) $(ORIG_LIBS)

## These settings are necessary for optimized WorkShop builds, due to
## BLAST's own use of them.
#CXXFLAGS = $(FAST_CXXFLAGS) $(ORIG_CXXFLAGS) 
#LDFLAGS = $(FAST_LDFLAGS)
#REQUIRES = objects -Cygwin LIBXML LIBXSLT
########################################################################

### END COPIED SETTINGS



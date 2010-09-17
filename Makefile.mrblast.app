#################################
# $Id: Makefile.blast_demo.app 146491 2008-11-26 13:57:23Z camacho $
#################################

APP = mrblast
SRC = mrblast

## Ssj
CC = mpicc
CXX = mpic++

MPI_COMPILE_FLAGS = $(shell mpic++ --showme:compile)
MPI_LINK_FLAGS = $(shell mpic++ --showme:link)
#MRMPI_USRLIB = -L./mrmpi -lmrmpi
MRMPI_USRLIB = -lmrmpi
 
CXXFLAGS = $(ORIG_CXXFLAGS) $(MPI_COMPILE_FLAGS)   
ORIG_LIBS = $(MRMPI_USRLIB)
LIBS = $(MPI_LINK_FLAGS) $(ORIG_LIBS) 
 
#LIB = $(LIB) -lmpi -lmrmpi
#ORIG_LIBS = $(ORIG_LIBS) $(MPI_LINK_FLAGS)
#ORIG_CXXFLAGS = $(MPI_COMPILE_FLAGS)
#ORIG_CXXFLAGS = $(ORIG_CXXFLAGS) -I/home/ssul/work/packages2/include

#CCC = mpic++
#F77 = mpif77
#FC = mpi
##

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
### END COPIED SETTINGS



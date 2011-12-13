BUILD - How to build mr-mpi-blast



Updated: 12/12/2011
By Seung-Jin Sul (ssul@jcvi.org)



* Requirements

    - The implementation of mr-mpi-blast uses the unmodified high-level NCBI 
      C++ Toolkit API library calls.
    - The mr-mpi-blast source tree should be inside of the NCBI C++ Toolkit 
      source tree for build.
    - The NCBI C++ Toolkit version 7.0.0 is tested for the current version of
      mr-mpi-blast. If you have problem with the toolkit version > 7.0.0, 
      please download the version 7.0.0.
    - The "MapReduce-MPI Library" from Sandia National Laboratories is used for
      MapReduce framework. The version 20Jun2011 is tested for the current 
      mr-mpi-blast and is included in "mrmpi" subdirectory of the source package.
    - To build mr-mpi-blast, the Boost library and a MPI stack should be 
      installed. The current version is tested with the MVAPICH MPI stack on 
      XSEDE TACC Ranger and on a local machine with Open MPI installed. The 
      MVAPICH is the current recommended MPI stack for large scale analysis on 
      Ranger.
    

1. Set installation directory

    The following assumes that you specify an installation target directory  
    using "MRMPIBLAST_PREFIX". 
      
    ex) $ mkdir $HOME/mr-mpi-blast
        $ export MRMPIBLAST_PREFIX=$HOME/mr-mpi-blast
           

2. Download the NCBI C++ Toolkit

    The following assumes that you are working in a directory, for 
    instance, "$HOME/mr-mpi-blast-src" for downloading and extracting the 
    C++ Toolkit source archive.
    
    ex) $ mkdir $HOME/mr-mpi-blast-src
        $ cd $HOME/mr-mpi-blast-src

    Download and extract the latest NCBI C++ Toolkit. The "<version>" represents 
    the latest version number.
    
    $ wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools++/ARCHIVE/<version>/
      ncbi_cxx--<version>.tar.gz 
    $ tar zxvf ncbi_cxx--<version>.tar.gz
       
   
3. Checkout the latest mr-mpi-blast. 

    * NOTE: The mr-mpi-blast should be inside of the NCBI C++ Toolkit source 
      tree for build.
  
    $ cd ./ncbi_cxx--<version>/src/app 
    $ git clone git@github.com:andreyto/mr-mpi-blast 
    $ cd ../../..
   

4. Edit meta-Make files  
   
    (1) Makefile.mrblast.app
    
        Edit "CC", "CXX", and "BOOST_HOME" in 
        ./ncbi_cxx--<version>/src/app/mr-mpi-blast/Makefile.mrblast.app
   
    (2) Makefile.in   
    
        Add "mr-mpi-blast" at end of "SUB_PROJ" line in 
        ./ncbi_cxx--<version>/src/app/Makefile.in

   
5. Configure and build the MapReduce-MPI library 

    (1) Makefile.mpicxx
    
        Edit "BOOST_HOME" in
        ./ncbi_cxx--<version>/src/app/mr-mpi-blast/mrmpi/MAKE/Makefile.mpicc (for Open MPI)
        or
        ./ncbi_cxx--<version>/src/app/mr-mpi-blast/mrmpi/MAKE/Makefile.mpicxx (for MVAPICH)
   
    (2) Build MR-MPI library
      
        $ cd ./ncbi_cxx--<version>/src/app/mr-mpi-blast/mrmpi 
        $ make mpicc (for Open MPI)
          or 
        $ make mpicxx (for MVAPICH)
        $ cd ../../../../..

    * NOTE: For details, please refer to the MapReduce MPI web pages 
      (http://www.sandia.gov/~sjplimp/mapreduce.html).
    
   
6. Configure, build, and install NCBI C++ toolkit with mr-mpi-blast
   
    $ cd ./ncbi_cxx--<version> 
    $ ./configure --prefix=$MRMPIBLAST_PREFIX 
    $ make && make install
        
    You can specify configuration flags as you want. For details, please refer 
    to the NCBI C++ Toolkit installation manual (http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/).
    
    ex) $ ./configure --prefix=$MRMPIBLAST_PREFIX --with-64 --without-dll 
          --without-gbench --without-gui --without-dbapi --without-connext 
          --without-xpm --without-sss --without-gnutls --without-openssl --with-mt 
 
    * NOTE: Once all the toolkit package is compiled, you can build only 
      mr-mpi-blast program in the build directory, if you want. 
    
    $ cd .../ncbi_cxx--<version>/GCC*/build/app/mr-mpi-blast/ && make


7. Remove unnecessary files

    Building and installing mr-mpi-blast with the NCBI C++ Toolkit will 
    install the toolkit related files and directories to 
    "$MRMPIBLAST_PREFIX". To remove those unnecessary files and directories 
    installed, please run the below.
   
    $ sh ./src/app/mr-mpi-blast/tools/remove_unnecessary.sh


9. Add "$MRMPIBLAST_PREFIX/bin" in PATH 

    For your convenience, add "$MRMPIBLAST_PREFIX/bin" to the executable search 
    path.
    
    $ export PATH=$MRMPIBLAST_PREFIX/bin:$PATH


    



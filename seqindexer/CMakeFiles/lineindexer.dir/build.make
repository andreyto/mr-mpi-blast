# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/ssul/work/packages/bin/cmake

# The command to remove a file.
RM = /home/ssul/work/packages/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /home/ssul/work/packages/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ssul/work/distros2/ncbi_cxx/ncbi_cxx--Jun_15_2010/GCC412-Debug64/build/app/mrblast/tools/lineindex

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ssul/work/distros2/ncbi_cxx/ncbi_cxx--Jun_15_2010/GCC412-Debug64/build/app/mrblast/tools/lineindex

# Include any dependencies generated for this target.
include CMakeFiles/lineindexer.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/lineindexer.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/lineindexer.dir/flags.make

CMakeFiles/lineindexer.dir/lineindexer.c.o: CMakeFiles/lineindexer.dir/flags.make
CMakeFiles/lineindexer.dir/lineindexer.c.o: lineindexer.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ssul/work/distros2/ncbi_cxx/ncbi_cxx--Jun_15_2010/GCC412-Debug64/build/app/mrblast/tools/lineindex/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/lineindexer.dir/lineindexer.c.o"
	/usr/local/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/lineindexer.dir/lineindexer.c.o   -c /home/ssul/work/distros2/ncbi_cxx/ncbi_cxx--Jun_15_2010/GCC412-Debug64/build/app/mrblast/tools/lineindex/lineindexer.c

CMakeFiles/lineindexer.dir/lineindexer.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/lineindexer.dir/lineindexer.c.i"
	/usr/local/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/ssul/work/distros2/ncbi_cxx/ncbi_cxx--Jun_15_2010/GCC412-Debug64/build/app/mrblast/tools/lineindex/lineindexer.c > CMakeFiles/lineindexer.dir/lineindexer.c.i

CMakeFiles/lineindexer.dir/lineindexer.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/lineindexer.dir/lineindexer.c.s"
	/usr/local/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/ssul/work/distros2/ncbi_cxx/ncbi_cxx--Jun_15_2010/GCC412-Debug64/build/app/mrblast/tools/lineindex/lineindexer.c -o CMakeFiles/lineindexer.dir/lineindexer.c.s

CMakeFiles/lineindexer.dir/lineindexer.c.o.requires:
.PHONY : CMakeFiles/lineindexer.dir/lineindexer.c.o.requires

CMakeFiles/lineindexer.dir/lineindexer.c.o.provides: CMakeFiles/lineindexer.dir/lineindexer.c.o.requires
	$(MAKE) -f CMakeFiles/lineindexer.dir/build.make CMakeFiles/lineindexer.dir/lineindexer.c.o.provides.build
.PHONY : CMakeFiles/lineindexer.dir/lineindexer.c.o.provides

CMakeFiles/lineindexer.dir/lineindexer.c.o.provides.build: CMakeFiles/lineindexer.dir/lineindexer.c.o
.PHONY : CMakeFiles/lineindexer.dir/lineindexer.c.o.provides.build

# Object files for target lineindexer
lineindexer_OBJECTS = \
"CMakeFiles/lineindexer.dir/lineindexer.c.o"

# External object files for target lineindexer
lineindexer_EXTERNAL_OBJECTS =

lineindexer: CMakeFiles/lineindexer.dir/lineindexer.c.o
lineindexer: CMakeFiles/lineindexer.dir/build.make
lineindexer: CMakeFiles/lineindexer.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable lineindexer"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lineindexer.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/lineindexer.dir/build: lineindexer
.PHONY : CMakeFiles/lineindexer.dir/build

CMakeFiles/lineindexer.dir/requires: CMakeFiles/lineindexer.dir/lineindexer.c.o.requires
.PHONY : CMakeFiles/lineindexer.dir/requires

CMakeFiles/lineindexer.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/lineindexer.dir/cmake_clean.cmake
.PHONY : CMakeFiles/lineindexer.dir/clean

CMakeFiles/lineindexer.dir/depend:
	cd /home/ssul/work/distros2/ncbi_cxx/ncbi_cxx--Jun_15_2010/GCC412-Debug64/build/app/mrblast/tools/lineindex && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ssul/work/distros2/ncbi_cxx/ncbi_cxx--Jun_15_2010/GCC412-Debug64/build/app/mrblast/tools/lineindex /home/ssul/work/distros2/ncbi_cxx/ncbi_cxx--Jun_15_2010/GCC412-Debug64/build/app/mrblast/tools/lineindex /home/ssul/work/distros2/ncbi_cxx/ncbi_cxx--Jun_15_2010/GCC412-Debug64/build/app/mrblast/tools/lineindex /home/ssul/work/distros2/ncbi_cxx/ncbi_cxx--Jun_15_2010/GCC412-Debug64/build/app/mrblast/tools/lineindex /home/ssul/work/distros2/ncbi_cxx/ncbi_cxx--Jun_15_2010/GCC412-Debug64/build/app/mrblast/tools/lineindex/CMakeFiles/lineindexer.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/lineindexer.dir/depend


# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/mdk/inti_navi/st_ekf_proj

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mdk/inti_navi/st_ekf_proj/build

# Include any dependencies generated for this target.
include CMakeFiles/estimator.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/estimator.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/estimator.dir/flags.make

CMakeFiles/estimator.dir/estimator.cpp.o: CMakeFiles/estimator.dir/flags.make
CMakeFiles/estimator.dir/estimator.cpp.o: ../estimator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mdk/inti_navi/st_ekf_proj/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/estimator.dir/estimator.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/estimator.dir/estimator.cpp.o -c /home/mdk/inti_navi/st_ekf_proj/estimator.cpp

CMakeFiles/estimator.dir/estimator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/estimator.dir/estimator.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mdk/inti_navi/st_ekf_proj/estimator.cpp > CMakeFiles/estimator.dir/estimator.cpp.i

CMakeFiles/estimator.dir/estimator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/estimator.dir/estimator.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mdk/inti_navi/st_ekf_proj/estimator.cpp -o CMakeFiles/estimator.dir/estimator.cpp.s

CMakeFiles/estimator.dir/rotation.cpp.o: CMakeFiles/estimator.dir/flags.make
CMakeFiles/estimator.dir/rotation.cpp.o: ../rotation.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mdk/inti_navi/st_ekf_proj/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/estimator.dir/rotation.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/estimator.dir/rotation.cpp.o -c /home/mdk/inti_navi/st_ekf_proj/rotation.cpp

CMakeFiles/estimator.dir/rotation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/estimator.dir/rotation.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mdk/inti_navi/st_ekf_proj/rotation.cpp > CMakeFiles/estimator.dir/rotation.cpp.i

CMakeFiles/estimator.dir/rotation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/estimator.dir/rotation.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mdk/inti_navi/st_ekf_proj/rotation.cpp -o CMakeFiles/estimator.dir/rotation.cpp.s

CMakeFiles/estimator.dir/earth_util.cpp.o: CMakeFiles/estimator.dir/flags.make
CMakeFiles/estimator.dir/earth_util.cpp.o: ../earth_util.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mdk/inti_navi/st_ekf_proj/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/estimator.dir/earth_util.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/estimator.dir/earth_util.cpp.o -c /home/mdk/inti_navi/st_ekf_proj/earth_util.cpp

CMakeFiles/estimator.dir/earth_util.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/estimator.dir/earth_util.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mdk/inti_navi/st_ekf_proj/earth_util.cpp > CMakeFiles/estimator.dir/earth_util.cpp.i

CMakeFiles/estimator.dir/earth_util.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/estimator.dir/earth_util.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mdk/inti_navi/st_ekf_proj/earth_util.cpp -o CMakeFiles/estimator.dir/earth_util.cpp.s

# Object files for target estimator
estimator_OBJECTS = \
"CMakeFiles/estimator.dir/estimator.cpp.o" \
"CMakeFiles/estimator.dir/rotation.cpp.o" \
"CMakeFiles/estimator.dir/earth_util.cpp.o"

# External object files for target estimator
estimator_EXTERNAL_OBJECTS =

estimator: CMakeFiles/estimator.dir/estimator.cpp.o
estimator: CMakeFiles/estimator.dir/rotation.cpp.o
estimator: CMakeFiles/estimator.dir/earth_util.cpp.o
estimator: CMakeFiles/estimator.dir/build.make
estimator: CMakeFiles/estimator.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/mdk/inti_navi/st_ekf_proj/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable estimator"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/estimator.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/estimator.dir/build: estimator

.PHONY : CMakeFiles/estimator.dir/build

CMakeFiles/estimator.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/estimator.dir/cmake_clean.cmake
.PHONY : CMakeFiles/estimator.dir/clean

CMakeFiles/estimator.dir/depend:
	cd /home/mdk/inti_navi/st_ekf_proj/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mdk/inti_navi/st_ekf_proj /home/mdk/inti_navi/st_ekf_proj /home/mdk/inti_navi/st_ekf_proj/build /home/mdk/inti_navi/st_ekf_proj/build /home/mdk/inti_navi/st_ekf_proj/build/CMakeFiles/estimator.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/estimator.dir/depend


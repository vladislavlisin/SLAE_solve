# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.21

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "D:\JetBrains\CLion 2021.3.3\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "D:\JetBrains\CLion 2021.3.3\bin\cmake\win\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = D:\CLionProjects\C++\SLAE_solve\Min_residual_method

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = D:\CLionProjects\C++\SLAE_solve\Min_residual_method\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/Min_residual_method.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/Min_residual_method.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/Min_residual_method.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Min_residual_method.dir/flags.make

CMakeFiles/Min_residual_method.dir/main.cpp.obj: CMakeFiles/Min_residual_method.dir/flags.make
CMakeFiles/Min_residual_method.dir/main.cpp.obj: ../main.cpp
CMakeFiles/Min_residual_method.dir/main.cpp.obj: CMakeFiles/Min_residual_method.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\CLionProjects\C++\SLAE_solve\Min_residual_method\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Min_residual_method.dir/main.cpp.obj"
	"D:\JetBrains\CLion 2021.3.3\bin\mingw\bin\g++.exe" $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/Min_residual_method.dir/main.cpp.obj -MF CMakeFiles\Min_residual_method.dir\main.cpp.obj.d -o CMakeFiles\Min_residual_method.dir\main.cpp.obj -c D:\CLionProjects\C++\SLAE_solve\Min_residual_method\main.cpp

CMakeFiles/Min_residual_method.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Min_residual_method.dir/main.cpp.i"
	"D:\JetBrains\CLion 2021.3.3\bin\mingw\bin\g++.exe" $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\CLionProjects\C++\SLAE_solve\Min_residual_method\main.cpp > CMakeFiles\Min_residual_method.dir\main.cpp.i

CMakeFiles/Min_residual_method.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Min_residual_method.dir/main.cpp.s"
	"D:\JetBrains\CLion 2021.3.3\bin\mingw\bin\g++.exe" $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:\CLionProjects\C++\SLAE_solve\Min_residual_method\main.cpp -o CMakeFiles\Min_residual_method.dir\main.cpp.s

# Object files for target Min_residual_method
Min_residual_method_OBJECTS = \
"CMakeFiles/Min_residual_method.dir/main.cpp.obj"

# External object files for target Min_residual_method
Min_residual_method_EXTERNAL_OBJECTS =

Min_residual_method.exe: CMakeFiles/Min_residual_method.dir/main.cpp.obj
Min_residual_method.exe: CMakeFiles/Min_residual_method.dir/build.make
Min_residual_method.exe: CMakeFiles/Min_residual_method.dir/linklibs.rsp
Min_residual_method.exe: CMakeFiles/Min_residual_method.dir/objects1.rsp
Min_residual_method.exe: CMakeFiles/Min_residual_method.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=D:\CLionProjects\C++\SLAE_solve\Min_residual_method\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Min_residual_method.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\Min_residual_method.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Min_residual_method.dir/build: Min_residual_method.exe
.PHONY : CMakeFiles/Min_residual_method.dir/build

CMakeFiles/Min_residual_method.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\Min_residual_method.dir\cmake_clean.cmake
.PHONY : CMakeFiles/Min_residual_method.dir/clean

CMakeFiles/Min_residual_method.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" D:\CLionProjects\C++\SLAE_solve\Min_residual_method D:\CLionProjects\C++\SLAE_solve\Min_residual_method D:\CLionProjects\C++\SLAE_solve\Min_residual_method\cmake-build-debug D:\CLionProjects\C++\SLAE_solve\Min_residual_method\cmake-build-debug D:\CLionProjects\C++\SLAE_solve\Min_residual_method\cmake-build-debug\CMakeFiles\Min_residual_method.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Min_residual_method.dir/depend

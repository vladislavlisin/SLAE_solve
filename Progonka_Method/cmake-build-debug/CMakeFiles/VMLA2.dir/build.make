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
CMAKE_SOURCE_DIR = D:\CLionProjects\C++\VMLA2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = D:\CLionProjects\C++\VMLA2\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/VMLA2.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/VMLA2.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/VMLA2.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/VMLA2.dir/flags.make

CMakeFiles/VMLA2.dir/main.cpp.obj: CMakeFiles/VMLA2.dir/flags.make
CMakeFiles/VMLA2.dir/main.cpp.obj: ../main.cpp
CMakeFiles/VMLA2.dir/main.cpp.obj: CMakeFiles/VMLA2.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\CLionProjects\C++\VMLA2\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/VMLA2.dir/main.cpp.obj"
	"D:\JetBrains\CLion 2021.3.3\bin\mingw\bin\g++.exe" $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/VMLA2.dir/main.cpp.obj -MF CMakeFiles\VMLA2.dir\main.cpp.obj.d -o CMakeFiles\VMLA2.dir\main.cpp.obj -c D:\CLionProjects\C++\VMLA2\main.cpp

CMakeFiles/VMLA2.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/VMLA2.dir/main.cpp.i"
	"D:\JetBrains\CLion 2021.3.3\bin\mingw\bin\g++.exe" $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\CLionProjects\C++\VMLA2\main.cpp > CMakeFiles\VMLA2.dir\main.cpp.i

CMakeFiles/VMLA2.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/VMLA2.dir/main.cpp.s"
	"D:\JetBrains\CLion 2021.3.3\bin\mingw\bin\g++.exe" $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:\CLionProjects\C++\VMLA2\main.cpp -o CMakeFiles\VMLA2.dir\main.cpp.s

# Object files for target VMLA2
VMLA2_OBJECTS = \
"CMakeFiles/VMLA2.dir/main.cpp.obj"

# External object files for target VMLA2
VMLA2_EXTERNAL_OBJECTS =

VMLA2.exe: CMakeFiles/VMLA2.dir/main.cpp.obj
VMLA2.exe: CMakeFiles/VMLA2.dir/build.make
VMLA2.exe: CMakeFiles/VMLA2.dir/linklibs.rsp
VMLA2.exe: CMakeFiles/VMLA2.dir/objects1.rsp
VMLA2.exe: CMakeFiles/VMLA2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=D:\CLionProjects\C++\VMLA2\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable VMLA2.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\VMLA2.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/VMLA2.dir/build: VMLA2.exe
.PHONY : CMakeFiles/VMLA2.dir/build

CMakeFiles/VMLA2.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\VMLA2.dir\cmake_clean.cmake
.PHONY : CMakeFiles/VMLA2.dir/clean

CMakeFiles/VMLA2.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" D:\CLionProjects\C++\VMLA2 D:\CLionProjects\C++\VMLA2 D:\CLionProjects\C++\VMLA2\cmake-build-debug D:\CLionProjects\C++\VMLA2\cmake-build-debug D:\CLionProjects\C++\VMLA2\cmake-build-debug\CMakeFiles\VMLA2.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/VMLA2.dir/depend

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
CMAKE_SOURCE_DIR = D:\CLionProjects\C++\find_inverse_matrix

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = D:\CLionProjects\C++\find_inverse_matrix\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/444_DFA.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/444_DFA.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/444_DFA.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/444_DFA.dir/flags.make

CMakeFiles/444_DFA.dir/main.cpp.obj: CMakeFiles/444_DFA.dir/flags.make
CMakeFiles/444_DFA.dir/main.cpp.obj: ../main.cpp
CMakeFiles/444_DFA.dir/main.cpp.obj: CMakeFiles/444_DFA.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\CLionProjects\C++\find_inverse_matrix\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/444_DFA.dir/main.cpp.obj"
	"D:\JetBrains\CLion 2021.3.3\bin\mingw\bin\g++.exe" $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/444_DFA.dir/main.cpp.obj -MF CMakeFiles\444_DFA.dir\main.cpp.obj.d -o CMakeFiles\444_DFA.dir\main.cpp.obj -c D:\CLionProjects\C++\find_inverse_matrix\main.cpp

CMakeFiles/444_DFA.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/444_DFA.dir/main.cpp.i"
	"D:\JetBrains\CLion 2021.3.3\bin\mingw\bin\g++.exe" $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\CLionProjects\C++\find_inverse_matrix\main.cpp > CMakeFiles\444_DFA.dir\main.cpp.i

CMakeFiles/444_DFA.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/444_DFA.dir/main.cpp.s"
	"D:\JetBrains\CLion 2021.3.3\bin\mingw\bin\g++.exe" $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:\CLionProjects\C++\find_inverse_matrix\main.cpp -o CMakeFiles\444_DFA.dir\main.cpp.s

# Object files for target 444_DFA
444_DFA_OBJECTS = \
"CMakeFiles/444_DFA.dir/main.cpp.obj"

# External object files for target 444_DFA
444_DFA_EXTERNAL_OBJECTS =

444_DFA.exe: CMakeFiles/444_DFA.dir/main.cpp.obj
444_DFA.exe: CMakeFiles/444_DFA.dir/build.make
444_DFA.exe: CMakeFiles/444_DFA.dir/linklibs.rsp
444_DFA.exe: CMakeFiles/444_DFA.dir/objects1.rsp
444_DFA.exe: CMakeFiles/444_DFA.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=D:\CLionProjects\C++\find_inverse_matrix\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable 444_DFA.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\444_DFA.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/444_DFA.dir/build: 444_DFA.exe
.PHONY : CMakeFiles/444_DFA.dir/build

CMakeFiles/444_DFA.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\444_DFA.dir\cmake_clean.cmake
.PHONY : CMakeFiles/444_DFA.dir/clean

CMakeFiles/444_DFA.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" D:\CLionProjects\C++\find_inverse_matrix D:\CLionProjects\C++\find_inverse_matrix D:\CLionProjects\C++\find_inverse_matrix\cmake-build-debug D:\CLionProjects\C++\find_inverse_matrix\cmake-build-debug D:\CLionProjects\C++\find_inverse_matrix\cmake-build-debug\CMakeFiles\444_DFA.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/444_DFA.dir/depend


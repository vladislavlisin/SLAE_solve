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
CMAKE_SOURCE_DIR = "D:\CLionProjects\C++\Shifted inverse iteration"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "D:\CLionProjects\C++\Shifted inverse iteration\cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/Shifted_inverse_iteration.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/Shifted_inverse_iteration.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/Shifted_inverse_iteration.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Shifted_inverse_iteration.dir/flags.make

CMakeFiles/Shifted_inverse_iteration.dir/main.cpp.obj: CMakeFiles/Shifted_inverse_iteration.dir/flags.make
CMakeFiles/Shifted_inverse_iteration.dir/main.cpp.obj: ../main.cpp
CMakeFiles/Shifted_inverse_iteration.dir/main.cpp.obj: CMakeFiles/Shifted_inverse_iteration.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="D:\CLionProjects\C++\Shifted inverse iteration\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Shifted_inverse_iteration.dir/main.cpp.obj"
	"D:\JetBrains\CLion 2021.3.3\bin\mingw\bin\g++.exe" $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/Shifted_inverse_iteration.dir/main.cpp.obj -MF CMakeFiles\Shifted_inverse_iteration.dir\main.cpp.obj.d -o CMakeFiles\Shifted_inverse_iteration.dir\main.cpp.obj -c "D:\CLionProjects\C++\Shifted inverse iteration\main.cpp"

CMakeFiles/Shifted_inverse_iteration.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Shifted_inverse_iteration.dir/main.cpp.i"
	"D:\JetBrains\CLion 2021.3.3\bin\mingw\bin\g++.exe" $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "D:\CLionProjects\C++\Shifted inverse iteration\main.cpp" > CMakeFiles\Shifted_inverse_iteration.dir\main.cpp.i

CMakeFiles/Shifted_inverse_iteration.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Shifted_inverse_iteration.dir/main.cpp.s"
	"D:\JetBrains\CLion 2021.3.3\bin\mingw\bin\g++.exe" $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "D:\CLionProjects\C++\Shifted inverse iteration\main.cpp" -o CMakeFiles\Shifted_inverse_iteration.dir\main.cpp.s

# Object files for target Shifted_inverse_iteration
Shifted_inverse_iteration_OBJECTS = \
"CMakeFiles/Shifted_inverse_iteration.dir/main.cpp.obj"

# External object files for target Shifted_inverse_iteration
Shifted_inverse_iteration_EXTERNAL_OBJECTS =

Shifted_inverse_iteration.exe: CMakeFiles/Shifted_inverse_iteration.dir/main.cpp.obj
Shifted_inverse_iteration.exe: CMakeFiles/Shifted_inverse_iteration.dir/build.make
Shifted_inverse_iteration.exe: CMakeFiles/Shifted_inverse_iteration.dir/linklibs.rsp
Shifted_inverse_iteration.exe: CMakeFiles/Shifted_inverse_iteration.dir/objects1.rsp
Shifted_inverse_iteration.exe: CMakeFiles/Shifted_inverse_iteration.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="D:\CLionProjects\C++\Shifted inverse iteration\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Shifted_inverse_iteration.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\Shifted_inverse_iteration.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Shifted_inverse_iteration.dir/build: Shifted_inverse_iteration.exe
.PHONY : CMakeFiles/Shifted_inverse_iteration.dir/build

CMakeFiles/Shifted_inverse_iteration.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\Shifted_inverse_iteration.dir\cmake_clean.cmake
.PHONY : CMakeFiles/Shifted_inverse_iteration.dir/clean

CMakeFiles/Shifted_inverse_iteration.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "D:\CLionProjects\C++\Shifted inverse iteration" "D:\CLionProjects\C++\Shifted inverse iteration" "D:\CLionProjects\C++\Shifted inverse iteration\cmake-build-debug" "D:\CLionProjects\C++\Shifted inverse iteration\cmake-build-debug" "D:\CLionProjects\C++\Shifted inverse iteration\cmake-build-debug\CMakeFiles\Shifted_inverse_iteration.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/Shifted_inverse_iteration.dir/depend


# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++17 -O3 -g -Wall -Wextra -pthread

# Source files
SOURCES = gSoD.cpp

# Header files
HEADERS = structures.hpp unique_ensurer.hpp

# Get names of executables from source files
EXECS = $(SOURCES:.cpp=)

.PHONY: all clean

all: $(EXECS)

# Rule for building executables, each depends on all headers
%: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(EXECS)


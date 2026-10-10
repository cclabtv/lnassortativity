# Define the C++ compiler to use
CXX = g++

#Define flags
OPT_FLAGS = -O2
DEBUG_FLAGS = -ggdb -O0

# Check if DEBUG is set, otherwise use optimization flags
ifeq ($(DEBUG), 1)
    CXXFLAGS = $(DEBUG_FLAGS)	
else
    CXXFLAGS = $(OPT_FLAGS)
endif

# Define the target files
TARGET = bin/simulate.out

# Define the C++ source files
SRC = src/simulate.cpp


# The default rule, compile our program
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)


# Define a clean rule
clean:
	$(RM) $(TARGET)

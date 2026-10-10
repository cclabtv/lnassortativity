# Define the C++ compiler to use
CXX = g++

#Define flags
OPT_FLAGS = -O2 -masm=intel
DEBUG_FLAGS = -ggdb -O0

# Check if DEBUG is set, otherwise use optimization flags
ifeq ($(DEBUG), 1)
    CXXFLAGS = $(DEBUG_FLAGS)	
else
    CXXFLAGS = $(OPT_FLAGS)
endif

# Define the target files
TARGET1 = bin/simulate.out
TARGET2 = bin/simulateBA.out
TARGET3 = bin/simulateWCap.out

# Define the C++ source files
SRC1 = src/simulate.cpp
SRC2 = src/simulateBA.cpp
SRC3 = src/simulateWCap.cpp

# The default rule, compile our program
all: $(TARGET1) $(TARGET2) $(TARGET3)

$(TARGET1): $(SRC1)
	$(CXX) $(CXXFLAGS) -o $(TARGET1) $(SRC1)

$(TARGET2): $(SRC2)
	$(CXX) $(CXXFLAGS) -o $(TARGET2) $(SRC2)

$(TARGET3): $(SRC3)
	$(CXX) $(CXXFLAGS) -o $(TARGET3) $(SRC3)


# Define a clean rule
clean:
	$(RM) $(TARGET3)

CXX = g++
CC = gcc

CXXFLAGS = -g -O2 -fsanitize=address -std=c++11 -pthread -march=native -Wall -I./utils 
CFLAGS   = -O2 -std=c99 -Wall -I./utils

LDLIBS = -lntl -lgmp -lflint -lm

BIN_DIR = tests/bin
SRC_DIR = tests

# All test sources except utils_test.cpp, being used to 
# declare useful functions for testing
TEST_SRCS := $(filter-out $(SRC_DIR)/utils_test.cpp,$(wildcard $(SRC_DIR)/*.cpp))
TEST_BINS := $(patsubst $(SRC_DIR)/%.cpp,$(BIN_DIR)/%,$(TEST_SRCS))

# Common sources for all tests
COMMON_SRCS := tests/utils_test.cpp \
               src/packing/cipher_packing.cpp \
               $(wildcard src/utils/*.cpp) \
               $(wildcard src/schemes/*.cpp) \
               src/utils/cyclotomic.c

COMMON_OBJS := $(COMMON_SRCS:.cpp=.o)
COMMON_OBJS := $(COMMON_OBJS:.c=.o)

all: $(BIN_DIR) $(TEST_BINS)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Rule to build each test binary
$(BIN_DIR)/%: $(SRC_DIR)/%.cpp $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(COMMON_OBJS) $(LDLIBS)

# Object file rules
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(BIN_DIR) $(COMMON_OBJS)

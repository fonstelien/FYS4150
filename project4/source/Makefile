CC := g++ -std=c++11
BIN := project4 test_project4
OBJ := utils.o modes.o
TEST_OBJ :=

ifdef optimize
OPTIMIZEFLAGS := $(optimize)
else
OPTIMIZEFLAGS := -O3
endif

ifdef debug
DEBUGFLAGS := -g
OPTIMIZEFLAGS := -O0
else
DEBUGFLAGS := -Wno-unused -Wno-unused-parameter -Wno-unused-result
endif

LIBS := -fopenmp -larmadillo

CCFLAGS := -Wall -Wextra $(LIBS) $(OPTIMIZEFLAGS) $(DEBUGFLAGS)

.PHONY: clean rm

all: project4

project4: project4.cpp $(OBJ)
	$(CC) $^ -o $@ $(CCFLAGS)

%.o: %.cpp utils.hpp
	$(CC) -c $< -o $@ $(CCFLAGS) 

clean:
	rm -f $(OBJ) $(TEST_OBJ)

rm:
	rm -f $(BIN) $(OBJ) $(TEST_OBJ)

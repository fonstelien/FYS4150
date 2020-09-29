CC := g++ -std=c++11
BIN := project2 test_project2
OBJ := utils.o jacobi.o polynomial_expansion.o
TEST_OBJ := $(OBJ) test_utils.o test_define.o

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

LIBS := -larmadillo

CCFLAGS := -Wall -Wextra $(LIBS) $(OPTIMIZEFLAGS) $(DEBUGFLAGS)

.PHONY: clean rm

all: project2 test_project2

project2: project2.cpp $(OBJ)
	$(CC) $^ -o $@ $(CCFLAGS)

test_project2: test_project2.cpp $(TEST_OBJ)
	$(CC) $^ -o $@ $(CCFLAGS)

%.o: %.cpp utils.hpp test_utils.hpp
	$(CC) -c $< -o $@ $(CCFLAGS) 

clean:
	rm -f $(OBJ) $(TEST_OBJ)

rm:
	rm -f $(BIN) $(OBJ) $(TEST_OBJ)

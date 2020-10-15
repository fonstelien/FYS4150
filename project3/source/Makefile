CC := g++ -std=c++11
BIN := project3 test_project3
OBJ := utils.o non_OO_earth_sun.o
TEST_OBJ := $(OBJ) catch.o

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

all: project3 test_project3

project3: project3.cpp $(OBJ)
	$(CC) $^ -o $@ $(CCFLAGS)

test_project3: test_project3.cpp $(TEST_OBJ)
	$(CC) $^ -o $@ $(CCFLAGS)

%.o: %.cpp utils.hpp
	$(CC) -c $< -o $@ $(CCFLAGS) 

clean:
	rm -f $(OBJ) $(TEST_OBJ)

rm:
	rm -f $(BIN) $(OBJ) $(TEST_OBJ)

CC := g++ -std=c++11
BIN := project1
OBJ := 

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


CCFLAGS := -Wall -Wextra -larmadillo $(OPTIMIZEFLAGS) $(DEBUGFLAGS)

.PHONY: clean rm

all: project1

project1: project1.cpp
	$(CC) $< -o $@ $(CCFLAGS)

%.o: %.cpp
	$(CC) -c $< -o $@ $(CCFLAGS) 

clean:
	rm -f $(OBJ)

rm:
	rm -f $(BIN) $(OBJ)

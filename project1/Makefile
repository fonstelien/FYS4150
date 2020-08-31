CC := g++ -std=c++11
BIN := project1
OBJ := 

ifdef optimize
OPTIMIZEFLAGS := $(optimize)
else
OPTIMIZEFLAGS := -O3
endif

CCFLAGS := -g -Wall -Wextra -Wno-unused -Wno-unused-parameter -Wno-unused-result -larmadillo $(OPTIMIZEFLAGS)

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

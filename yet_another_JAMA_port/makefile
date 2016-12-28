SRC=main.cpp SquareMatrix.cpp EigenvalueCalculator.cpp
CC=g++
CXXFLAGS=-std=c++11
APP=eigcalc

all: $(APP)

$(APP): $(SRC)
	$(CC) $(SRC) $(CXXFLAGS) -o $(APP)

clean:
	rm -f $(APP)

SRC=main.cpp SquareMatrix.cpp QRIterations.cpp
CC=g++
CXXFLAGS=-std=c++11 -fopenmp
APP=qriterations

all: $(APP)

$(APP): $(SRC)
	$(CC) $(SRC) $(CXXFLAGS) -o $(APP)

clean:
	rm $(APP)

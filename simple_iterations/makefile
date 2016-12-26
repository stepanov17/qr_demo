SRC=main.cpp SquareMatrix.cpp QRIterations.cpp
CC=g++
CXXFLAGS=-std=c++11 -fopenmp
STD=-std=c++11
APP=qriterations
#no-omp version for comparison
APP1=qriterations1

all: $(APP) $(APP1)

$(APP): $(SRC)
	$(CC) $(SRC) $(CXXFLAGS) -o $(APP)
$(APP1): $(SRC)
	$(CC) $(SRC) $(STD) -o $(APP1)


clean:
	rm -f $(APP)
	rm -f $(APP1)

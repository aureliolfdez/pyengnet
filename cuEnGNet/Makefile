program_NAME := cuEnGNet
program_CPP_SRCS := $(wildcard src/*.cpp)
CC ?=

all: $(program_NAME)

$(program_NAME): src/main.o src/preprocess.o src/correlations.o
	nvcc --cudart static --relocatable-device-code=false -gencode arch=compute_$(CC),code=compute_$(CC) -link -o "cuEnGNet" ./src/main.o ./src/preprocess.o ./src/correlations.o

src/main.o: src/main.cpp -lpthread 
	g++ -c src/main.cpp -o src/main.o 

src/preprocess.o: src/preprocess.cpp
	g++ -c src/preprocess.cpp -o src/preprocess.o 

src/correlations.o: src/correlations.cu
	nvcc -G -g -O0 -std=c++11 -gencode arch=compute_$(CC),code=sm_$(CC)  -odir "src" -M -o "src/correlations.d" "src/correlations.cu"
	nvcc -G -g -O0 -std=c++11 --compile --relocatable-device-code=false -gencode arch=compute_$(CC),code=compute_$(CC) -gencode arch=compute_$(CC),code=sm_$(CC)  -x cu -o  "src/correlations.o" "src/correlations.cu"

clean:
	rm -R src/*.o src/*.d results/* $(program_NAME) ./results
# dichiarazioni variabili
CL= mpicxx -Ofast -march=native -flto -fopenmp -pthread
CC= mpicxx -c -Ofast -march=native -flto -fopenmp
#pulizia
#.PHONY: *.o d3d

#definizioni dei target

all: d3d

d3d: d3d.o matrix.o utilities.o polynomials.o computationalElement.o physicalElement.o vectors.o cases.o inputOutput.o boundaryCondition.o physicalRelations.o tensor.o LES.o
	$(CL) d3d.o matrix.o utilities.o polynomials.o computationalElement.o physicalElement.o vectors.o cases.o inputOutput.o boundaryCondition.o physicalRelations.o tensor.o LES.o -o d3d
d3d.o: d3d.cpp
	$(CC) d3d.cpp	
matrix.o: matrix.cpp
	$(CC) matrix.cpp
utilities.o: utilities.cpp
	$(CC) utilities.cpp
polynomials.o: polynomials.cpp
	$(CC) polynomials.cpp
computationalElement.o: computationalElement.cpp
	$(CC) computationalElement.cpp
physicalElement.o: physicalElement.cpp
	$(CC) physicalElement.cpp
vectors.o: vectors.cpp
	$(CC) vectors.cpp
cases.o: cases.cpp
	$(CC) cases.cpp
inputOutput.o: inputOutput.cpp
	$(CC) inputOutput.cpp
boundaryCondition.o: boundaryCondition.cpp
	$(CC) boundaryCondition.cpp
physicalRelations.o: physicalRelations.cpp
	$(CC) physicalRelations.cpp
tensor.o: tensor.cpp
	$(CC) tensor.cpp
LES.o: LES.cpp
	$(CC) LES.cpp
	
clean:
	rm -rf *.o d3d
	
	

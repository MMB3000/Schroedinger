
all: Schroedinger


Schroedinger: MatAnyD.o GaussianElimination.o Eigen.o SchroedingerSolver.o Main.o
		g++ -o Schroedinger Main.o SchroedingerSolver.o Eigen.o GaussianElimination.o MatAnyD.o

MatAnyD.o: MatAnyD.cpp MatAnyD.h
		g++ -o MatAnyD.o -c MatAnyD.cpp

GaussianElimination.o: GaussianElimination.cpp MatAnyD.h
		g++ -o GaussianElimination.o -c GaussianElimination.cpp

Eigen.o: Eigen.cpp MatAnyD.h
		g++ -o Eigen.o -c Eigen.cpp

SchroedingerSolver.o: SchroedingerSolver.cpp SchroedingerSolver.h
		g++ -o SchroedingerSolver.o -c SchroedingerSolver.cpp

Main.o: Main.cpp MatAnyD.h
		g++ -o Main.o -c Main.cpp

clean:
		rm -f Schroedinger Main.o SchroedingerSolver.o Eigen.o GaussianElimination.o MatAnyD.o


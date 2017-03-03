default: heat
	mpirun -n 8 heat
heat: main.cpp heat.c heat.h
	mpic++ -std=c++11 main.cpp -o heat
test: test.cpp Matrix.h
	g++ -std=c++11 test.cpp -o test
.PHONY: clean
clean:
	-rm heat
	-rm test
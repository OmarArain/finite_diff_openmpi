default: heat
heat: main.cpp heat.c heat.h
	mpic++ -std=c++11 main.cpp -o heat
.PHONY: clean
clean:
	-rm heat
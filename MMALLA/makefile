INPUT = matriz.c
OUTPUT = matriz
OPTS = mpicc -w

all: 
	${OPTS} ${INPUT} -o ${OUTPUT} -O3 -lm
clear:
	rm ${OUTPUT}
debug:
	${OPTS} -g ${INPUT} -o ${OUTPUT} -lm
run:
	mpirun -np 4 ${OUTPUT}
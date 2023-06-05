INPUT = primos_test.c
OUTPUT = primos_test
OPTS = mpicc -w

all:
	${OPTS} ${INPUT} -o ${OUTPUT} -O3 -lm


clear: 
	rm ${OUTPUT}

redo: 
	rm ${OUTPUT} && ${OPTS} ${INPUT} -o ${OUTPUT}
	
run:
	mpirun -np 3 ${OUTPUT}
INPUT = anillo.c
OUTPUT = anillo
OPTS = mpicc -w

all:
	${OPTS} ${INPUT} -o ${OUTPUT} -O3


clear: 
	rm ${OUTPUT}

redo: 
	rm ${OUTPUT} && ${OPTS} ${INPUT} -o ${OUTPUT}
	
run:
	mpirun -np 5 ${OUTPUT}
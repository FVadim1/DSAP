# include <math.h> 
# include <stdio.h> 
# include <stdlib.h> 
# include <time.h> 
#include <mpi.h>

int numero_primos_par(int n, int rank, int size);
double get_max_time(double *times, int num_times);
int esPrimo(int p);
int numero_primos_sec( int n );

int main ( int argc, char *argv[] ) 
{ 
	//paralelo
	int mpi_primos, total_mpi_primos = 0;
	int rank, size;

	int n_min = 500; 
	int n_max = 5000000; 
	int n_factor = 10; 
	int n = n_min; 

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

	while ( n <= n_max ) 
	{ 
		printf("[P:%d, n:%d]Empiezo la iteración del bucle. \n", rank,n);

        mpi_primos = numero_primos_par(n, rank, size);
		printf("[P:%d, n:%d]---He obtenido mis números primos.  \n", rank,n);

		MPI_Reduce(&mpi_primos, &total_mpi_primos, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		printf("[P:%d, n:%d]---He sobrepasado el \033[1mreduce\033[0m. \n", rank,n);


		if(rank == 0){
			printf("▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓Primos menores que %10d: %10d. \n", n, total_mpi_primos);
        }
		
		printf("[P:%d, n:%d]He terminado la iteración del bucle. \n", rank,n);

		total_mpi_primos = 0;
		n = n * n_factor; //le añado un 0 a n:  500->5000->50000->...->5000000
	} 

	MPI_Finalize();
	printf(">>>>Llego al final, proceso: %d \n", rank);


	return 0; 
} 

int numero_primos_par(int n, int rank, int size)
{

    int inicio = (rank * n) / size + 2;
    int fin = ((rank + 1) * n) / size;
    int total = 0;  //el total de primos que ha encontrado este proceso
    for (int i = inicio; i <= fin; i++)
    {
        if (esPrimo(i) == 1)
            total++;
    }
    return total;
}

int esPrimo(int p)
{
	for (int i=2;i<=sqrt(p);i++)
	{
		if (p%i == 0)
			return 0;
	}
	return 1;
}



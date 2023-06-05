# include <math.h> 
# include <stdio.h> 
# include <stdlib.h> 
# include <time.h> 
#include <mpi.h>

int numero_primos_sec( int n );
int numero_primos_par(int n, int rank, int size);
double get_max_time(double *times, int num_times);
int esPrimo(int p);
int numero_primos_sec( int n );

int main ( int argc, char *argv[] ) 
{ 
	//secuencial
	int n, n_factor, n_min, n_max; 
	int primos, primos_parte;
	unsigned t0,t1; 

	//paralelo
	int mpi_primos, total_mpi_primos = 0;
	int rank, size;
	unsigned mpi_t0,mpi_t1; 
	double tiempo_paralelo;

	double *tiempos = NULL;

	n_min = 500; 
	n_max = 5000000; 
	n_factor = 10; 
	printf ( "\n" ); 
	printf ( " Programa para contar el número de primos menores que un valor.\n" ); 
	printf ( "\n" ); 
	n = n_min; 

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

	while ( n <= n_max ) 
	{ 
		//SECUENCIAL
		if(rank == 0){
			t0 = clock();
			primos = numero_primos_sec(n);
			t1 = clock();
			printf("Primos menores que %10d: %10d. Tiempo secuencial: %5.2f s.\n", n, primos, (double)(t1-t0)/CLOCKS_PER_SEC);
		}
		
		//PARALELO
		if (rank == 0) { tiempos = malloc(size * sizeof(double));}

		mpi_t0 = MPI_Wtime();
        mpi_primos = numero_primos_par(n, rank, size);
        mpi_t1 = MPI_Wtime();
		tiempo_paralelo = (double)(mpi_t1 - mpi_t0);
		
		MPI_Reduce(&mpi_primos, &total_mpi_primos, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    	MPI_Gather(&tiempo_paralelo, 1, MPI_DOUBLE, tiempos, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		//MPI_Barrier(MPI_COMM_WORLD);

		if(rank == 0){

			double max_time = get_max_time(tiempos, size);

			printf("Primos menores que %10d: %10d. Tiempo paralelo: %5.2f s.", n, total_mpi_primos, max_time);

			printf(" Tiempos parciales: ");

			for (int i = 0; i < size; i++) {
            	printf("%5.2f",tiempos[i]);
				printf(" ");
			}

			 printf("\n");

			 free(tiempos);
        }
		
		total_mpi_primos = 0;
		n = n * n_factor; //le añado un 0 a n:  500->5000->50000->...->5000000

	} 
	//printf("llego aqui, proceso: %d", rank);

	MPI_Finalize();

	return 0; 
} 

//En este método reparto los números que le tocarían a cada proceso 
// Devuelve el numero de primos menores que n , repartidos segun cual es el rank del proceso
int numero_primos_par(int n, int rank, int size)
{
	//SEPARO LOS NÚMEROS QUE LE TOCARÍAN A CADA PROCESO SEGÚN EL NÚMERO DE SU RANK
	/* Si tenemos que hay 4 procesadores (size = 4), rank de 0 a 3, teniendo n=500
	Proceso 1: del 1 al 125, 
	Proceso 2: del 126 a 250, 
	etc
	*/
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

double get_max_time(double *times, int num_times)
{
    double max_time = times[0];
    for (int i = 1; i < num_times; i++) {
        if (times[i] > max_time) {
            max_time = times[i];
        }
    }
    return max_time;
}

// Devuelve el numero de primos menores que n
int numero_primos_sec( int n ) 
{ 
	int total=0; 
	for (int i = 2; i <= n; i++) 
	{ 
		if (esPrimo(i)==1)
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


/*int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
MPI_REDUCE es una función de la biblioteca MPI que se utiliza para combinar los resultados de todas las tareas que participan en una comunicación. 
La función MPI_REDUCE toma los datos enviados por cada tarea y los combina en una sola tarea (por lo general, la tarea raíz) usando una operación e
specífica.

Después de realizar el MPI_Reduce, el proceso con rango 0 tendrá el resultado final. 

sendbuf: Puntero al búfer que contiene los datos que se enviarán.
recvbuf: Puntero al búfer que recibirá los datos.
count: Número de elementos en los búferes sendbuf y recvbuf.
datatype: Tipo de dato de los elementos en los búferes sendbuf y recvbuf.
op: Operación de reducción a realizar.
root: Rango de la tarea raíz.
comm: Comunicador utilizado para la operación.

MPI_REDUCE se utiliza típicamente para recoger información de todas las tareas y realizar algún tipo de operación en ella, 
como calcular el promedio, la suma total, el mínimo o el máximo. El resultado final se guarda en el búfer recvbuf de la tarea raíz.
MPI_REDUCE se puede usar para realizar operaciones de reducción, como la suma, el producto, el mínimo, el máximo, etc., en los datos 
enviados por las tareas. La operación que se utiliza se especifica mediante el parámetro "op". Los valores aceptados son MPI_SUM, MPI_PROD, 
MPI_MAX, MPI_MIN, MPI_LAND, MPI_BAND, MPI_LOR, MPI_BOR, MPI_LXOR, MPI_BXOR.


*/


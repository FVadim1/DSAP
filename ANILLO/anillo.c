#include <stdio.h>
#include <mpi.h>


int main(int argc, char* argv[]) {
  int my_rank, cantidadProcesadores, n;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); //devuelve el rango del proceso (su posición)
  MPI_Comm_size(MPI_COMM_WORLD, &cantidadProcesadores); //cantidad del procesos que haran el programa
  //printf("%d\n", size);

  if (my_rank == 0) {
    // El proceso root lee el número del input estándar
    printf("Introduce número: "); fflush( stdout );
    scanf("%d", &n);
    n++;
    // Envía el número al proceso 1
    MPI_Send(
      &n, //buf: Variable que contiene la información a comunicar.
      1,  //int count: Cantidad de elementos contenidos en buf
      MPI_INT, //datatype: Tipo de variable del buf
      1, //int dest: Al proceso al cual transfiero el buf
      0,  //int tag: Identifica el envío, generalmente es 0 y solo cámbia cuando se ha de comunicar más de un envio
      MPI_COMM_WORLD //comunicador
      );

    MPI_Recv(
      &n, //buf: variable que contiene la información a comunicar
      1, //int count: cantidad de elementos contendos en buf.
      MPI_INT, //datatype: tipo de variable buf.
      cantidadProcesadores-1, //int source: de quien espero el buf
      0, //int tag: identifica el envio, generalmente es 0
      MPI_COMM_WORLD, //comunicador
      MPI_STATUS_IGNORE //estructura interna que contiene los campos MPI_SOURCE, MPI_TAG y count
      );

      printf("Respuesta: %d\n", n);
      
  } else {
    // Recibe el número del proceso anterior
    MPI_Recv(&n,1, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
    n++;
    // Envía el número al siguiente proceso si es que existe
    if (my_rank < cantidadProcesadores - 1) {
      MPI_Send(&n, 1, MPI_INT, my_rank+1, 0, MPI_COMM_WORLD);
    } else {
      // Si es el último proceso, envía el número al proceso root
      MPI_Send(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
  }

  MPI_Finalize();
  return 0;
}
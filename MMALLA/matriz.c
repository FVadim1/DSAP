//echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>


void mult(double a[], double b[], double *c, int m) {
    int i, j, k;
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            for (k = 0;k < m;k++) {
                c[i*m+j] = c[i*m+j]+a[i*m+k]*b[k*m+j];
            }
        }
    }
}

void cyan() {printf("\033[1;36m");}
void green() {printf("\033[1;32m");}
void yellow() {printf("\033[1;33m");}
void red() {printf("\033[1;31m");}
void reset() {printf("\033[0m");}

int main(int argc, char **argv) {

    int myrank, numprocs, bloqtam, r,  maxbloqtam=100, minbloqtam = 1;
    int fila, columna, arriba, abajo,error;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    r = sqrt(numprocs);        

    // El proceso root recibe bloqtam
    /*
    █ █▄ █ █▀█ █ █ ▀█▀   █▄█   █▀▀ █▀█ █▀▄▀█ █▀█ █▀█ █▀█ █▄▄ ▄▀█ █▀▀ █ █▀█ █▄ █
    █ █ ▀█ █▀▀ █▄█  █     █    █▄▄ █▄█ █ ▀ █ █▀▀ █▀▄ █▄█ █▄█ █▀█ █▄▄ █ █▄█ █ ▀█*/
    if(myrank == 0) {

        //cyan();printf("████████████ nproc = %d █████████████\n", numprocs);reset(); //DEBUG
        //cyan();printf("████████████████ r = %d █████████████\n", r);reset(); //DEBUG

        cyan();printf("Introduce el tamaño del bloque(bloqtam): \n");fflush( stdout );reset();
        scanf("%d", &bloqtam);
        //bloqtam = 5;

        //cyan();printf("███████████bloqtam = %d █████████████\n", bloqtam);reset(); //DEBUG

        if(r*r != numprocs || bloqtam < minbloqtam || bloqtam > maxbloqtam ) {
            printf("ERROR: Número de procesos no es un cuadrado perfecto o tamaño de bloque demasiado grande.\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

	/*█▄▄ █▀█ █▀█ ▄▀█ █▀▄ █▀▀ ▄▀█ █▀ ▀█▀
	  █▄█ █▀▄ █▄█ █▀█ █▄▀ █▄▄ █▀█ ▄█  █ */
    // Envía bloqtam a todos los otros procesos
    MPI_Bcast(&bloqtam, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /*
    █▀▀ █ █   ▄▀█   █▀▀ █▀█ █   █ █ █▀▄▀█ █▄ █ ▄▀█   ▄▀█ █▀█ █▀█ █ █▄▄ ▄▀█   ▄▀█ █▄▄ ▄▀█   █ █▀█
    █▀  █ █▄▄ █▀█   █▄▄ █▄█ █▄▄ █▄█ █ ▀ █ █ ▀█ █▀█   █▀█ █▀▄ █▀▄ █ █▄█ █▀█   █▀█ █▄█ █▀█ █▄█ █▄█
    FILA Y COLUMNA se utilizan para saber la posición de cada proceso p0[0,0], p1[0,1], p2[1,0], p3[1,1] (p0 y p1 están en la misma fila, p2 y p3 están en la misma fila)*/
    fila = myrank / r;
    columna = myrank % r;

    //ARRIBA Y ABAJO se utilizan para la ROTACIÓN (conocer los procesos arriba y abajo)
    arriba = (fila == 0) ? (r - 1) * r + columna : (fila - 1) * r + columna; //primera fila
    abajo = (fila == r - 1) ? columna : (fila + 1) * r + columna; //ultima fila

    //yellow();printf("proceso= %d, [fila = %d , columna= %d ] # [arriba = %d, abajo = %d] \n",myrank,fila,columna,arriba,abajo);reset();

    /* Para r= 2 , 4 procesos
    arriba de fila 0:   2   3
    arriba de fila 1: | 0 | 1 |
                      |---|---|
    abajo fila fila 0 | 2 | 3 |
    abajo de fila 1:    0   1   */

    /*
    █▀▄▀█ █ █▀▀ █ █   ▄▀█   █▀   ▀█
    █ ▀ █ █ █▀  █ █▄▄ █▀█   █▄   ▄█
    SACO LOS PROCESOS QUE HAY EN MI FILA -> CADA PROCESO TENDRÁ SUS PROPIOS VALORES EN MIFILA[]*/
    int mifila[r-1]; //con 4 procesos, son 2 filas, por lo que solo 1 proceso (diferente a mi) puede estar en mi fila
    for (int i = 0, j = 0; i < numprocs; i++){
        if ((fila == i/r) && (columna != i%r)) {// Si [i] es la misma fila y diferente columna.
            mifila[j] = i;  // obtengo todos los procesos situados en la fila menos el proceso en el que estoy
            j++;
        }
    }

    /*
    █ █▄ █ █ █▀▀ █ ▄▀█ █   █ ▀█ ▄▀█ █▀█   ▄▀█ █▀█ █▀█ ▄▀█ █▄█ █▀
    █ █ ▀█ █ █▄▄ █ █▀█ █▄▄ █ █▄ █▀█ █▀▄   █▀█ █▀▄ █▀▄ █▀█  █  ▄█*/
    double a[maxbloqtam*maxbloqtam],
           b[maxbloqtam*maxbloqtam],
           c[maxbloqtam*maxbloqtam],
           atmp[maxbloqtam*maxbloqtam];

    for(int i = 0; i < bloqtam * bloqtam; ++i) {
        //Inicializar array 'a'
        a[i] = i * (double) (fila * columna + 1) / bloqtam;
        
        /*Inicializar matriz 'b'.Se inicializa como una matriz identidad solo en aquellos procesos 
        donde el número de fila y columna son los mismos. Para todos los demás procesos, 'b' se inicializa como una matriz de ceros. */
        int row = i / bloqtam;
        int col = i % bloqtam;
        b[i] = (row == col && columna == fila) ? 1.0 : 0.0; 

        //Inicializar array 'c'
        c[i] = 0.0;
    }

    /*█▀▀ ▄▀█ █▄░█ █▄░█ █▀█ █▄░█
      █▄▄ █▀█ █░▀█ █░▀█ █▄█ █░▀█*/
    /* ＩＮＩＣＩＡＲ ＢＵＦＦＥＲ ＰＡＲＡ ＵＴＩＬＩＺＡＲ ＢＳＥＮＤ
    MPI_Pack_size calcula el tamaño requerido para empaquetar una cierta cantidad de datos, en este caso, bloques de una matriz de tipo double. 
    Este tamaño se ajusta con MPI_BSEND_OVERHEAD (espacio adicional requerido por MPI_Bsend) y se multiplica por el número de procesos (numprocs). 
    Se crea un buffer con ese tamaño y se adjunta a MPI para su uso con MPI_Bsend. */
    int tamBsendBuffer = 0; double *buffer;
    MPI_Pack_size(bloqtam*bloqtam, MPI_DOUBLE, MPI_COMM_WORLD, &tamBsendBuffer); //Calcula el tamaño en bytes del espacio requerido para empaquetar una cantidad de datos dada.
    tamBsendBuffer = numprocs * (tamBsendBuffer + MPI_BSEND_OVERHEAD);//MPI_BSEND_OVERHEAD: Es la cantidad de espacio adicional requerido por el sistema al usar MPI_Bsend. Se suma al tamaño del buffer.
    buffer = malloc(tamBsendBuffer * sizeof(double));
    MPI_Buffer_attach(buffer, tamBsendBuffer); //Asigna un buffer de envío a MPI para su uso con MPI_Bsend. 

    /*ＡＬＧＯＲＩＴＭＯ ＤＥ ＣＡＮＮＯＮ: 
    En el bucle principal, cada procesador intercambia datos con otros procesadores y realiza la multiplicación de matrices. 
    Las matrices a y b son bloques de las matrices más grandes que se están multiplicando, y c es la matriz de resultados.*/
    for (int i = 0; i < r; i++) { 
        if (columna == (fila+i) % r) {//if para determinar si un proceso debe realizar el cálculo localmente o tiene que recibir de otros.
            // Mandar el array 'a' a todos los procesadores mi fila
            for (int j = 0; j < r-1; j++) {
                MPI_Send(&a, bloqtam*bloqtam, MPI_DOUBLE, mifila[j], i, MPI_COMM_WORLD);
            }
            mult(a, b, c, bloqtam);
        }
        else {
            //Recibo 'a' y lo guardo en 'atmp'
            MPI_Recv(&atmp, bloqtam*bloqtam, MPI_DOUBLE, fila*r + ((fila+i) % r), i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            mult(atmp, b, c, bloqtam);  
           
            /*  fila*r + ((fila+i) % r)
            La parte 'fila*r' de la expresión asegura que se está buscando procesos en la misma fila que el proceso actual. Esto se debe a que todos 
            los procesos en una fila dada tienen rangos que son múltiplos de 'r' (por ejemplo, en una grilla con r=2, los procesos en la fila 0 tendrán rangos 0, 2, 4,..., 
            y los procesos en la fila 1 tendrán rangos 1, 3, 5,...).

            La parte '((fila+i) % r)' de la expresión calcula el desplazamiento desde el proceso actual al proceso del que se va a recibir. En la primera iteración (i=0), cada 
            proceso necesita recibir de su propio proceso (desplazamiento 0). En la segunda iteración (i=1), cada proceso necesita recibir del proceso a su derecha (desplazamiento 1), 
            y así sucesivamente. El operador % r se usa para garantizar que el desplazamiento se mantenga dentro de los límites de la grilla (es decir, que no exceda r-1).
                 
            */     
        }

        /*ＲＯＴＡＣＩＯＮ
        Después de cada multiplicación, los bloques de la matriz b se rotan a lo largo de las columnas, utilizando la función MPI_Bsend 
        para enviar el bloque al proceso 'arriba' y MPI_Recv para recibir un nuevo bloque del proceso 'abajo'.  */
        //NOTA: Se utiliza 'Bsend' en vez de 'Send' porque  a diferencia de Send , Bsend no es bloqueante (se pueden dar problemas de sincronia y quedarse bloqueado.)
        MPI_Bsend(&b, bloqtam*bloqtam, MPI_DOUBLE, arriba, 2*i, MPI_COMM_WORLD); // envia b al proceso arriba
        MPI_Recv(&b, bloqtam*bloqtam, MPI_DOUBLE, abajo, 2*i, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // recibir b del proceso abajo 
    }

    /* ＬＩＢＥＲＡＲ ＢＵＦＦＥＲ */
     MPI_Buffer_detach(&buffer, &tamBsendBuffer); //  el buffer adjunto a MPI se desacopla (liberar buffer)
    /*██████████████████████████████████████████████████████████████████████████████████████████*/

    /*
    █▀▀ █▀█ █▀█ █▀█ █▀█ █▀▀ █▀
    ██▄ █▀▄ █▀▄ █▄█ █▀▄ ██▄ ▄█*/
    for (int i = 0; i < bloqtam*bloqtam; i++) { 
        if(abs(a[i] - c[i]) > 0.0000001) { error++;} //Si A != C hay error.
    }

    if (error != 0) { red();printf("Proceso [%d]. Errores = %d\n", myrank, error);reset();}   
    else { green();printf("Proceso [%d]. Errores = %d\n", myrank, error);reset();}

    //printArray(a, bloqtam*bloqtam,myrank,bloqtam, "Array a" );
    //printArray(b, bloqtam*bloqtam, myrank,bloqtam, "Array b");

    free(buffer);
    MPI_Finalize();
}


void printArray(double matriz[], int len, int rank,int bloqtam, char nombrematriz[])
{
    cyan();
    printf("[%6s] rank = %d \n",nombrematriz,rank);
    for (int i=0;i<len;i++){
        printf("%.1f ",matriz[i]);
        if ((i+1)%bloqtam  == 0) printf("\n");
    }
    printf("\n"); 
    reset();
}


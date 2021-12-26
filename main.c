#include <mpi.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char** argv) {
    int rc, myrank, nproc, namelen;
    char name[MPI_MAX_PROCESSOR_NAME];
    
    rc = MPI_Init(&argc, &argv);

    if (rc != MPI_SUCCESS) {
        printf("Error starting process, Code: %d\n", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
        return rc;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);  // Gets rank of process
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);   // Gets number of processors

    if (myrank == 0) {
        printf("This is the main process\n");
    }

    namelen = MPI_MAX_PROCESSOR_NAME;
    MPI_Get_processor_name(name, &namelen);
    printf("hello world %d from ’%s’\n", myrank, name);

    int x = 0;
    int y = 0;
    float z = 0.0f;
    while(x < 100000000){
      while (y < 100000000)
      {
        z += sqrtf((float)x + (float)y);
        y++;
      }
      x++;
    }
    printf("Exit with z = %f\n", z);

    MPI_Finalize();

    return 0;
}

#include "main.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

// Handles for custom datatypes
MPI_Datatype initial_data_handle;

void create_custom_data_types() {
    // init_data
    {
        int lengths[3] = {sizeof(int), sizeof(int), sizeof(double)};
        struct initial_data dummy_struct;
        MPI_Aint base_address;
        MPI_Aint displacements[3];
        MPI_Get_address(&dummy_struct, &base_address);
        MPI_Get_address(&dummy_struct.width, &displacements[0]);
        MPI_Get_address(&dummy_struct.height, &displacements[1]);
        MPI_Get_address(&dummy_struct.precision, &displacements[2]);
        displacements[0] = MPI_Aint_diff(displacements[0], base_address);
        displacements[1] = MPI_Aint_diff(displacements[1], base_address);
        displacements[2] = MPI_Aint_diff(displacements[2], base_address);

        MPI_Datatype types[3] = {MPI_INT, MPI_DOUBLE, MPI_CHAR};
        MPI_Type_create_struct(3, lengths, displacements, types,
                               &initial_data_handle);
        MPI_Type_commit(&initial_data_handle);
    }
}

// initialises grid
// init decided what pattern to initialise the grid with
void generateGrid(Grid* g, int init) {
    unsigned int width = g->width;
    unsigned int height = g->height;
    // Allocate memory for each grid cell
    g->val = (double*)malloc(width * height * sizeof(double));
    if (g->val == NULL) {
        printf("Cannot allocate memory for new grid");
        exit(-1);
        return;
    }
    unsigned int x, y;
    for (unsigned int i = 0; i < width; i++) {
        x = i % width;
        y = i / width;

        // Initialise each memory location with a value
        // If init == True assign left and top edges to 1 and the other
        // edges to 0. The rest are set to random values
        if (init == PATTERN_RANDOM) {
            g->val[i] = (double)(rand() / 10) / (double)(RAND_MAX / 10);
        } else if (init == PATTERN_GRADIENT) {
            if (x == 0 || (y == 0 && x != width - 1)) {
                g->val[i] = 1;
            } else if (x == width - 1 || y == height - 1) {
                g->val[i] = -1;
            } else {
                g->val[i] = 0;  //
                (double)(rand() / 10) / (double)(RAND_MAX / 10);
            }
        } else if (init == PATTERN_ZERO) {
            g->val[i] = 0;
        }
    }
}

void manager(const unsigned int width, const unsigned int height,
             const unsigned int n_process) {
    unsigned int n_workers = n_process - 1;
    float gap = (float)width / (n_workers);
    {
        MPI_Request init_data_requests[n_workers - 1];
        struct initial_data init_data;

        for (unsigned int i = 0; i < n_workers; i++) {
            init_data.height = height;
            init_data.width = (unsigned int)(gap * i);
            init_data.precision = 0.1;

            MPI_Isend(&init_data, 1, initial_data_handle, i + 1, TAG_INIT_GRID,
                      MPI_COMM_WORLD, &init_data_requests[i]);
        }
    }
}

void worker(const unsigned int my_rank) {
    struct initial_data init_data;
    MPI_Status stat;
    MPI_Recv(&init_data, 1, initial_data_handle, 0, TAG_INIT_GRID,
             MPI_COMM_WORLD, &stat);
}

int main(int argc, char** argv) {
    const unsigned int width = 16, height = 16;
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

    // Get name of hardware running this process
    namelen = MPI_MAX_PROCESSOR_NAME;
    MPI_Get_processor_name(name, &namelen);
    printf("%s connected.\n", name);

    // Initalise custom data types
    create_custom_data_types();

    if (myrank == 0) {
        manager(width, height, nproc);
    } else {
        worker(myrank);
    }

    MPI_Finalize();

    return 0;
}

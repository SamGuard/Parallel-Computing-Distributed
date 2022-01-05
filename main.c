#include "main.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

// Handles for custom datatypes
MPI_Datatype initial_data_handle;

void create_custom_data_types() {
    // init_data
    {
        int lengths[3] = {sizeof(int), sizeof(int), sizeof(double)};
        initial_data dummy_struct;
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

void print_grid(Grid* g) {
    int x, y;
    for (int i = 0; i < g->width * g->height; i++) {
        printf("%f", g->val[i]);

        if ((i + 1) % g->width != 0) {
            printf(" ");
        } else {
            printf("\n");
        }
    }
}

void print_error(int err_code) {
    char* stringBuff = (char*)malloc(2048 * sizeof(char));
    int resLength = 0;
    MPI_Error_string(err_code, stringBuff, &resLength);
    printf("Error in communication, error code: %s\n", stringBuff);
}

// initialises grid
// init decided what pattern to initialise the grid with
void generate_grid(Grid* g, int init) {
    unsigned int width = g->width;
    unsigned int height = g->height;
    // Allocate memory for each grid cell
    g->val = (double*)malloc(width * height * sizeof(double));
    if (g->val == NULL) {
        printf("Cannot allocate memory for new grid of size %d, %d", width,
               height);
        exit(-1);
        return;
    }
    unsigned int x, y;
    for (unsigned int i = 0; i < width * height; i++) {
        x = i % width;
        y = i / width;

        // Initialise each memory location with a value
        // If init == True assign left and top edges to 1 and the other
        // edges to 0. The rest are set to random values
        if (init == PATTERN_RANDOM) {
            g->val[i] = (double)(rand() / 10) / (double)(RAND_MAX / 10);
        } else if (init == PATTERN_GRADIENT) {
            if (x == 0 || (y == 0 && x != width - 1)) {
                g->val[i] = 1.0;
            } else if (x == width - 1 || y == height - 1) {
                g->val[i] = -1.0;
            } else {
                g->val[i] = 0.0;
            }
        } else if (init == PATTERN_ZERO) {
            g->val[i] = 0;
        }
    }
}

// Send data to workers of the size of grid they will work on
grid_metadata* send_init_data_to_workers(unsigned int width,
                                         unsigned int height, float gap,
                                         unsigned int n_workers) {
    MPI_Request init_data_requests[n_workers - 1];
    initial_data init_data;
    // Holds all information about the position, size and length of grids in
    // workers
    grid_metadata *g_data_array =
                      (grid_metadata*)malloc(sizeof(grid_metadata) * n_workers),
                  *g_data;

    unsigned int start_row, end_row;
    for (unsigned int i = 0; i < n_workers; i++) {
        g_data = &g_data_array[i];
        init_data.width = g_data->width = width;

        // Partition the grid into n_worker amount of groups of rows and then
        // add 1 or 2 rows for the perimeter around it
        start_row = (unsigned int)(gap * i) - (i == 0 ? 0 : 1);
        end_row = (unsigned int)(gap * (i + 1) + (i == n_workers - 1 ? 0 : 1));
        g_data->start_index = start_row * width;
        init_data.height = g_data->height = end_row - start_row;
        g_data->length = g_data->height * width;

        init_data.precision = 0.1;
        // Send asynchronously so the grid can be generated in the meantime
        MPI_Isend(&init_data, 1, initial_data_handle, i + 1, TAG_INIT_GRID,
                  MPI_COMM_WORLD, &init_data_requests[i]);
    }

    MPI_Status stat;
    for (int i = 0; i < n_workers; i++) {
        if (MPI_Wait(&init_data_requests[i], &stat) != MPI_SUCCESS) {
            print_error(stat.MPI_ERROR);
            return NULL;
        }
    }
    return g_data_array;
}

// Sends each worker their section of the grid to work on
// Requests are asynchronous
MPI_Request* send_entire_grid_to_workers(const float gap,
                                         const unsigned int n_workers, Grid* g,
                                         grid_metadata* g_data_array) {
    // Returned from function so it can wait until all data has been sent before
    // checking if they have finished
    MPI_Request* requests =
        (MPI_Request*)malloc(sizeof(MPI_Request) * n_workers);
    for (int i = 0; i < n_workers; i++) {
        /*
        printf("Start Index: %d, Length: %d\n", g_data_array[i].start_index,
               g_data_array[i].length);
        */
        MPI_Isend(g->val + g_data_array[i].start_index, g_data_array[i].length,
                  MPI_DOUBLE, i + 1, TAG_GRID_DATA, MPI_COMM_WORLD,
                  &requests[i]);
    }
    return requests;
}

// Gets all data from all workers and aggregates it into one grid
// g.val must be allocated before calling function
void retrieve_entire_grid_from_workers(const float gap,
                                       const unsigned int n_workers, Grid* g,
                                       grid_metadata* g_data_array) {
    if (g->val == NULL) {
        printf("grid not initialised in retrieve_entire_grid\n");
        return;
    }

    double* valBuff =
        (double*)malloc(((int)ceil(gap) + 2) * g->width * sizeof(double));

    grid_metadata* g_data;
    unsigned int index;
    for (int i = 0; i < n_workers; i++) {
        g_data = &g_data_array[i];

        /*
        printf("Start: %d, Length: %d, Height: %d\n", g_data->start_index,
               g_data->length, g_data->height);
        */
        MPI_Recv(valBuff, g_data->length, MPI_DOUBLE, i + 1, TAG_GRID_DATA,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Copies data from workers grid to main grid
        // Cannot be directly copied as the perimiter of each of the grids will
        // mostly be incorrect
        for (int y = 1; y < g_data->height - 1; y++) {
            for (int x = 1; x < g->width - 1; x++) {
                index = x + y * g->width;
                /*
                printf("x: %d, y: %d, index: %d, val_index: %d\n", x, y, index,
                       g_data->start_index + index);
                */
                g->val[g_data->start_index + index] = valBuff[index];
            }
        }
    }

    free(valBuff);
}

void send_outer_perimiter() {}

void recv_inner_perimiter(Grid* g, grid_metadata* g_data, unsigned int rank) {
    unsigned int length = g_data->width * 2;
    double* buff = (double*)malloc(length * sizeof(double));

    MPI_Recv(buff, length, MPI_DOUBLE, TAG_INNER_PERIM, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void manager(const unsigned int width, const unsigned int height,
             const unsigned int n_process) {
    const unsigned int n_workers = n_process - 1;  // Number of workers
    const float gap =
        (float)height / (float)n_workers;  // How many columns per worker

    grid_metadata* g_data_array;  // Stores meta data about the workers grids
    g_data_array = send_init_data_to_workers(width, height, gap, n_workers);

    Grid g;
    g.width = width;
    g.height = height;
    generate_grid(&g, PATTERN_GRADIENT);

    MPI_Request* send_requests =
        send_entire_grid_to_workers(gap, n_workers, &g, g_data_array);
    {
        MPI_Status stat;
        for (int i = 0; i < n_workers; i++) {
            if (MPI_Wait(&send_requests[i], &stat) != MPI_SUCCESS) {
                print_error(stat.MPI_ERROR);
                return;
            }
        }
    }

    for (int i = 0; i < 100; i++) {
    }
    retrieve_entire_grid_from_workers(gap, n_workers, &g, g_data_array);
    print_grid(&g);
    printf("------\n");
    free(g_data_array);
}

void compute_step(Grid* in, Grid* out) {
    double v0, v1, v2, v3;
    for (int y = 1; y < in->height - 1; y++) {
        for (int x = 1; x < in->width - 1; x++) {
            v0 = in->val[x + 1 + y * in->width];
            v1 = in->val[x - 1 + y * in->width];
            v2 = in->val[x + (y + 1) * in->width];
            v3 = in->val[x + (y - 1) * in->width];
            out->val[x + y * in->width] = (v0 + v1 + v2 + v3) / 4.0;
        }
    }
}

// Sends the top and bottom row of the grid as this is all that is needed in
// other workers
void send_inner_perimeter(Grid* g) {
    unsigned int width = g->width, height = g->height;
    unsigned int length = width * 2;
    double* buff = (double*)malloc(length * sizeof(double));

    int start_index;  // Index of row in grid
    int index;        // Index in buffer
    for (int i = 0; i < 2; i++) {
        start_index =
            i * (width * (height - 1));  // Either index 0 or the last row
        for (int j = 0; j < width; j++) {
            buff[index] = g->val[start_index + j];
            index++;
        }
    }

    MPI_Send(buff, length, MPI_DOUBLE, 0, TAG_INNER_PERIM, MPI_COMM_WORLD);
    free(buff);
}

// Recives the rows that make up the bounderies of this workers grid
// (the two columns that make up the other bounderies are constant and are sent
// at the start)
void recv_outer_perimeter(Grid* g) {}

void worker(const unsigned int my_rank) {
    initial_data init_data;
    MPI_Status stat;
    Grid g0, g1, temp;
    MPI_Recv(&init_data, 1, initial_data_handle, 0, TAG_INIT_GRID,
             MPI_COMM_WORLD, &stat);
    g0.width = g1.width = init_data.width;
    g0.height = g1.height = init_data.height;
    // printf("Width: %d, Height: %d\n", g0.width, g0.height);
    generate_grid(&g0, PATTERN_ZERO);
    generate_grid(&g1, PATTERN_ZERO);

    const double precision = init_data.precision;

    // Get starting values of the grid
    MPI_Recv(g0.val, g0.width * g0.height, MPI_DOUBLE, 0, TAG_GRID_DATA,
             MPI_COMM_WORLD, &stat);
    for (int i = 0; i < 1; i++) {
        recv_outer_perimeter(&g0);
        compute_step(&g0, &g1);
        /*
        if (my_rank == 1) {
            print_grid(&g0);
            printf("----\n");
            print_grid(&g1);
        }
        */
        send_inner_perimeter(&g1);
        temp = g1;
        g1 = g0;
        g0 = temp;
    }
    MPI_Send(g1.val, g1.width * g1.height, MPI_DOUBLE, 0, TAG_GRID_DATA,
             MPI_COMM_WORLD);
}

int main(int argc, char** argv) {
    const unsigned int width = 1024, height = 1024;
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

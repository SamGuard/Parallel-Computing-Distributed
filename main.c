#include "main.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

// Handles for custom datatypes
MPI_Datatype initial_data_handle;

// Initially I used datatypes but when running on the hpc there were errors i
// could not solve
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
    for (unsigned int i = 0; i < g->width * g->height; i++) {
        printf("%f", g->val[i]);

        if ((i + 1) % g->width != 0) {
            printf(" ");
        } else {
            printf("\n");
        }
    }
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
            } else if (x == width - 1 || (x != 0 && y == height - 1)) {
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
grid_metadata* send_init_data_to_workers(unsigned int width, float gap,
                                         unsigned int n_workers) {
    unsigned int init_data[2];
    // Holds all information about the position, size and length of grids in
    // workers
    grid_metadata *g_data_array =
                      (grid_metadata*)malloc(sizeof(grid_metadata) * n_workers),
                  *g_data;

    if (g_data_array == NULL) {
        printf("Failed to allocate memory in send_init_data_to_workers");
    }

    unsigned int start_row, end_row;
    for (unsigned int i = 0; i < n_workers; i++) {
        g_data = &g_data_array[i];
        init_data[0] = g_data->width = width;

        // Partition the grid into n_worker amount of groups of rows and then
        // add 1 or 2 rows for the perimeter around it
        start_row = (unsigned int)(gap * i) - (i == 0 ? 0 : 1);
        end_row = (unsigned int)(gap * (i + 1) + (i == n_workers - 1 ? 0 : 1));
        g_data->start_index = start_row * width;
        init_data[1] = g_data->height = end_row - start_row;
        g_data->length = g_data->height * width;

        MPI_Send(init_data, 2, MPI_UINT32_T, i + 1, TAG_INIT_GRID,
                 MPI_COMM_WORLD);
    }

    return g_data_array;
}

// Sends each worker their section of the grid to work on
void send_entire_grid_to_workers(const unsigned int n_workers,
                                 Grid* g, grid_metadata* g_data_array) {
    // Returned from function so it can wait until all data has been sent before
    // checking if they have finished
    MPI_Request* requests =
        (MPI_Request*)malloc(sizeof(MPI_Request) * n_workers);
    if (requests == NULL) {
        printf("Failed to allocate memory in send_entire_grid_to_workers");
    }
    for (unsigned int i = 0; i < n_workers; i++) {
        /*
        printf("Start Index: %d, Length: %d\n", g_data_array[i].start_index,
               g_data_array[i].length);
        */
        MPI_Send(g->val + g_data_array[i].start_index, g_data_array[i].length,
                 MPI_DOUBLE, i + 1, TAG_GRID_DATA, MPI_COMM_WORLD);
    }
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
    // Temporary place to store data from workers
    double* valBuff =
        (double*)malloc(((int)ceil(gap) * 2) * g->width * sizeof(double));
    if (valBuff == NULL) {
        printf(
            "Failed to allocate memory in retrieve_entire_grid_from_workers\n");
    }

    grid_metadata* g_data;
    unsigned int index;
    for (unsigned int i = 0; i < n_workers; i++) {
        g_data = &g_data_array[i];
        // Get grid from worker
        MPI_Recv(valBuff, g_data->length, MPI_DOUBLE, i + 1, TAG_GRID_DATA,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Copies data from workers grid to main grid
        // Cannot be directly copied as the perimiter of each of the grids will
        // mostly be incorrect
        for (unsigned int y = 1; y < g_data->height - 1; y++) {
            for (unsigned int x = 1; x < g->width - 1; x++) {
                index = x + y * g->width;
                g->val[g_data->start_index + index] = valBuff[index];
            }
        }
    }

    free(valBuff);
}

// Sends the top and bottom row of the grid as this is all that is needed in
// other workers
void send_perimeter(Grid* g, unsigned int row1, unsigned int row2,
                    unsigned int dest_rank) {
    unsigned int width = g->width;
    // Size of the two rows combined minus the 4 values for the perimeter of the
    // entire grid which are constant
    unsigned int length = width * 2 - 4;
    double* buff = (double*)malloc(length * sizeof(double));
    if (buff == NULL) {
        printf("Failed to allocate in send_inner_perimeter\n");
    }

    unsigned int start_index;  // Index of row in grid
    unsigned int index = 0;    // Index in buffer
    for (unsigned int i = 0; i < 2; i++) {
        // Gived index of each row
        start_index = width * (i == 0 ? row1 : row2);
        for (unsigned int j = 1; j < width - 1; j++) {
            buff[index] = g->val[start_index + j];
            index++;
        }
    }

    MPI_Send(buff, length, MPI_DOUBLE, dest_rank, TAG_INNER_PERIM,
             MPI_COMM_WORLD);
    free(buff);
}

// Recieves two rows which are placed in the grid
void recv_perimiter(Grid* g, unsigned int row1, unsigned int row2,
                    unsigned int dest_rank) {
    unsigned int width = g->width;
    // Size of the two rows combined minus the 4 values for the
    // perimeter of the entire grid which are constant
    unsigned int length = g->width * 2 - 4;
    double* buff = (double*)malloc(length * sizeof(double));
    if (buff == NULL) {
        printf("Failed to allocate memory in recv_inner_perimiter");
    }
    MPI_Recv(buff, length, MPI_DOUBLE, dest_rank, TAG_INNER_PERIM,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    unsigned int start_index;
    unsigned int index = 0;  // Index in buffer
    for (unsigned int i = 0; i < 2; i++) {
        // Gives index of where to put the data
        start_index = width * (i == 0 ? row1 : row2);
        for (unsigned int j = 1; j < width - 1; j++) {
            g->val[start_index + j] = buff[index];
            // printf("%f ", buff[index]);
            index++;
        }
        // printf("\n");
    }
    free(buff);
}

// Sends if the worker is done
int is_worker_done(int rank) {
    int res;
    MPI_Recv(&res, 1, MPI_INT, rank, TAG_IS_WORKER_DONE, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    return res;
}

void manager(const unsigned int width, const unsigned int height,
             const unsigned int n_process, const double precision,
             struct timeval startTime) {
    // printf("Started manager\n");
    const unsigned int n_workers = n_process - 1;  // Number of workers
    const float gap =
        (float)height / (float)n_workers;  // How many columns per worker

    grid_metadata* g_data_array;  // Stores data about the workers grids
    g_data_array = send_init_data_to_workers(width, gap, n_workers);

    Grid g;
    g.width = width;
    g.height = height;
    generate_grid(&g, PATTERN_GRADIENT);

    // Initalise each worker with grid
    send_entire_grid_to_workers(n_workers, &g, g_data_array);

    // Main loop
    unsigned int row1, row2;     // Index of rows to send
    int is_finished = FALSE;     // Should the program exit
    unsigned int iteration = 0;  // Counts the number of iterations
    while (is_finished == FALSE) {
        for (unsigned int i = 0; i < n_workers; i++) {
            row1 = g_data_array[i].start_index / width;
            row2 = row1 + g_data_array[i].height - 1;
            // Send if the worker should exit or not
            MPI_Send(&is_finished, 1, MPI_INT, i + 1, TAG_IS_WORKER_DONE,
                     MPI_COMM_WORLD);
            // Send only the information needed for the next iteration
            send_perimeter(&g, row1, row2, i + 1);
        }
        is_finished = TRUE;
        for (unsigned int i = 0; i < n_workers; i++) {
            row1 = g_data_array[i].start_index / width + 1;
            row2 = row1 + g_data_array[i].height - 3;
            // Take in only the data that is needed in other workers' grids
            recv_perimiter(&g, row1, row2, i + 1);
            if (is_worker_done(i + 1) == FALSE) {
                is_finished = FALSE;
            }
        }
        iteration++;
    }
    // Tell all workers to exit
    for (unsigned int i = 0; i < n_workers; i++) {
        MPI_Send(&is_finished, 1, MPI_INT, i + 1, TAG_IS_WORKER_DONE,
                 MPI_COMM_WORLD);
    }

    // Collect all data from workers
    retrieve_entire_grid_from_workers(gap, n_workers, &g, g_data_array);

    // Calculate time difference
    struct timeval endTime;
    gettimeofday(&endTime, NULL);
    // Adds milliseconds to seconds * 1000000 then finds the difference between
    // start and finishe time
    double diffTime = ((endTime.tv_sec * 1000000 + endTime.tv_usec) -
                       (startTime.tv_sec * 1000000 + startTime.tv_usec));
    // Scales time to seconds
    diffTime = diffTime / 1000000.0;
    printf("%d,%u,%u,%f,%d,%f,", n_process, width, height, precision, iteration,
           diffTime);
    print_grid(&g);

    free(g.val);
    free(g_data_array);
}

double compute_step(Grid* in, Grid* out) {
    // v0-3 are the values of the surrounding cells, max_diff is the maximum
    // difference between values, diff is used in the calculation of max_diff
    double v0, v1, v2, v3, max_diff = 0.0, diff;
    unsigned int index;
    for (unsigned int y = 1; y < in->height - 1; y++) {
        for (unsigned int x = 1; x < in->width - 1; x++) {
            index = x + y * in->width;
            v0 = in->val[x + 1 + y * in->width];
            v1 = in->val[x - 1 + y * in->width];
            v2 = in->val[x + (y + 1) * in->width];
            v3 = in->val[x + (y - 1) * in->width];
            // Calc difference
            out->val[index] = (v0 + v1 + v2 + v3) / 4.0;
            diff = fabs(out->val[index] - in->val[index]);
            // Store if bigger
            max_diff = max_diff < diff ? diff : max_diff;
        }
    }
    return max_diff;
}

void worker(const unsigned int my_rank, const double precision) {
    // Get width and height
    unsigned int init_data[2];
    if (init_data == NULL) {
        printf("Failed to allocate memory\n");
    }

    Grid g0, g1, temp;
    MPI_Recv(init_data, 2, MPI_UINT32_T, 0, TAG_INIT_GRID, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    // Init grids
    g0.width = g1.width = init_data[0];
    g0.height = g1.height = init_data[1];
    generate_grid(&g0, PATTERN_ZERO);
    generate_grid(&g1, PATTERN_ZERO);

    // Get starting values of the grid
    MPI_Recv(g0.val, g0.width * g0.height, MPI_DOUBLE, 0, TAG_GRID_DATA,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // copy perimeter from g0 to g1 as these are constant
    for (unsigned int x = 0; x < g0.width; x++) {
        g1.val[x] = g0.val[x];
        g1.val[x + g0.width * (g0.height - 1)] =
            g0.val[x + g0.width * (g0.height - 1)];
    }
    for (unsigned int y = 0; y < g0.height; y++) {
        g1.val[y * g0.width] = g0.val[y * g0.width];
        g1.val[(y + 1) * g0.width - 1] = g0.val[(y + 1) * g0.width - 1];
    }

    int is_finished = FALSE, is_less_than_prec;
    while (TRUE) {
        // See if program is finished
        MPI_Recv(&is_finished, 1, MPI_INT, 0, TAG_IS_WORKER_DONE,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (is_finished == TRUE) {
            break;
        }
        // Get data from main grid
        recv_perimiter(&g0, 0, g0.height - 1, 0);
        is_less_than_prec = compute_step(&g0, &g1) < precision;
        // Send data needed in other workers
        send_perimeter(&g1, 1, g0.height - 2, 0);

        // Swap pointers around
        temp = g1;
        g1 = g0;
        g0 = temp;
        // Send if this worker's grid is under the precision required
        MPI_Send(&is_less_than_prec, 1, MPI_INT, 0, TAG_IS_WORKER_DONE,
                 MPI_COMM_WORLD);
    }
    // Send all data back to main grid
    MPI_Send(g0.val, g0.width * g0.height, MPI_DOUBLE, 0, TAG_GRID_DATA,
             MPI_COMM_WORLD);
    free(g0.val);
    free(g1.val);
}

int main(int argc, char** argv) {
    struct timeval startTime;
    gettimeofday(&startTime, NULL);

    unsigned int width = 64, height = 64;
    double precision = 0.001;
    if (argc == 4) {
        width = atoi(argv[1]);
        height = atoi(argv[2]);
        precision = atof(argv[3]);
    }
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
    

    // Initalise custom data types
    // create_custom_data_types();

    // The program can only handle at most 1 worker per row so nproc must be less than the height
    // Threads with rank above or equal to the height are stopped
    if (myrank == 0) {
        manager(width, height, (unsigned int)nproc >= height ? height : (unsigned int) nproc, precision, startTime);
    } else if((unsigned int)myrank < height) {
        worker(myrank, precision);
    }
    MPI_Finalize();

    return 0;
}

// Tag numbers
#define TAG_INIT_GRID 100  // Meta data for the workers grid
#define TAG_GRID_DATA 101  // Communicating values in the grid

// Types of grid arrangement
#define PATTERN_ZERO 0
#define PATTERN_GRADIENT 1
#define PATTERN_RANDOM 2

// Stores all the values in the grid as well as meta data
typedef struct grid {
    double* val;                // Array for each value in cell
    unsigned int width, height;  // Size of the grid
} Grid;

// These are constants to send each of the workers
typedef struct initial_data {
    unsigned int width, height;
    double precision;
} initial_data_struct;
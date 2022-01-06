#define TRUE 1
#define FALSE 0

// Tag numbers
#define TAG_INIT_GRID 100  // Meta data for the workers grid
#define TAG_GRID_DATA 101  // Communicating values in the grid
#define TAG_INNER_PERIM 102 // Send inner perimeter of a grid
#define TAG_OUTTER_PERIM 103 // Send outter perimiter of grid
#define TAG_IS_WORKER_DONE 104 // Checking if worker is finished

// Types of grid arrangement
#define PATTERN_ZERO 0
#define PATTERN_GRADIENT 1
#define PATTERN_RANDOM 2

// Stores all the values in the grid as well as meta data
typedef struct grid {
    double* val;                 // Array for each value in cell
    unsigned int width, height;  // Size of the grid
} Grid;

// These are constants to send each of the workers
typedef struct initial_data {
    unsigned int width, height;
    double precision;
} initial_data;

// Data about the dimensions, size and the position of the workers grid
typedef struct grid_metadata {
    unsigned int width, height,
        start_index,  // Start of the workers area grid in the main grid
        length;       // How big this area is
} grid_metadata;

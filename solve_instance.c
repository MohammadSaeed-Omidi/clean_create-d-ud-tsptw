#include "solve_instance.h"

int run_concorde(const char *tsp_filename) {
    char command[1024];
    char directory[256] = {0};
    
    // Extract the directory part from tsp_filename
    char *last_slash = strrchr(tsp_filename, '/');
    if (last_slash) {
        // Copy the directory part including the slash
        strncpy(directory, tsp_filename, last_slash - tsp_filename + 1);
        directory[last_slash - tsp_filename + 1] = '\0';
    }
    
    // Change to the directory where the TSP file is located
    char current_dir[1024];
    if (getcwd(current_dir, sizeof(current_dir)) == NULL) {
        perror("Failed to get current directory");
        return -1;
    }
    
    if (directory[0] != '\0' && chdir(directory) != 0) {
        perror("Failed to change directory");
        return -1;
    }
    
    // Extract just the filename part
    const char *base_filename = last_slash ? last_slash + 1 : tsp_filename;
    
    // Run Concorde in the directory
    snprintf(command, sizeof(command), "%s -x %s", CONCORDE_PATH, base_filename);
    
    // First check if Concorde exists and is executable
    if (access(CONCORDE_PATH, X_OK) != 0) {
        fprintf(stderr, "Concorde executable not found or not executable at %s\n", CONCORDE_PATH);
        return -1;
    }
    
    printf("Running in directory %s: %s\n", directory[0] != '\0' ? directory : ".", command);
    int ret = system(command);
    
    // Change back to the original directory
    if (chdir(current_dir) != 0) {
        perror("Failed to change back to original directory");
    }
    
    // printf("\nret: %d\n", ret);
    ret = 0;
    return ret == 0 ? 0 : -1;
}
// ... existing code ...

// Add these new functions before solve_instance():

// Read edge weights from TSP file into a 2D array
int** read_edge_weights(const char* tsp_filename, int* dimension) {
    FILE* fp = fopen(tsp_filename, "r");
    if (!fp) {
        perror("Failed to open TSP file");
        return NULL;
    }

    char line[256];
    *dimension = 0;

    // Find DIMENSION in the header
    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "DIMENSION:", 10) == 0) {
            sscanf(line, "DIMENSION: %d", dimension);
            break;
        } else if (strncmp(line, "DIMENSION ", 10) == 0) {
            sscanf(line, "DIMENSION %d", dimension);
            break;
        }
    }

    if (*dimension == 0) {
        fclose(fp);
        return NULL;
    }

    // Allocate memory for edge weights
    int** weights = malloc(*dimension * sizeof(int*));
    for (int i = 0; i < *dimension; i++) {
        weights[i] = malloc(*dimension * sizeof(int));
    }

    // Find EDGE_WEIGHT_SECTION
    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "EDGE_WEIGHT_SECTION", 18) == 0) {
            break;
        }
    }

    // Read the weight matrix
    for (int i = 0; i < *dimension; i++) {
        for (int j = 0; j < *dimension; j++) {
            if (fscanf(fp, "%d", &weights[i][j]) != 1) {
                // Handle error
                for (int k = 0; k <= i; k++) {
                    free(weights[k]);
                }
                free(weights);
                fclose(fp);
                return NULL;
            }
        }
    }

    fclose(fp);
    return weights;
}

// Calculate tour length from .sol file using edge weights
double calculate_tour_length(const char* sol_filename, const char* tsp_filename) {
    int dimension;
    int** weights = read_edge_weights(tsp_filename, &dimension);
    if (!weights) return -1.0;

    FILE* fp = fopen(sol_filename, "r");
    if (!fp) {
        perror("Failed to open solution file");
        return -1.0;
    }

    // Read number of nodes
    int n;
    fscanf(fp, "%d", &n);
    if (n != dimension) {
        fclose(fp);
        return -1.0;
    }

    // Read the tour
    int* tour = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        if (fscanf(fp, "%d", &tour[i]) != 1) {
            free(tour);
            fclose(fp);
            return -1.0;
        }
    }

    // Calculate tour length
    double length = 0;
    for (int i = 0; i < n; i++) {
        int from = tour[i];
        int to = tour[(i + 1) % n];  // Wrap around to start for last node
        length += weights[from][to];
    }

    // Clean up
    free(tour);
    for (int i = 0; i < dimension; i++) {
        free(weights[i]);
    }
    free(weights);
    fclose(fp);

    return length;
}

double solve_tsp_instance(const char* tsp_filename, const char* sol_filename) {
    // Run Concorde - it will create the .sol file in the same directory as the .tsp
    if (run_concorde(tsp_filename) != 0) {
        fprintf(stderr, "Failed to solve instance %s\n", tsp_filename);
        return;
    }
    
    // Calculate the optimal tour length
    double optimal_value = calculate_tour_length(sol_filename, tsp_filename);
    
    return optimal_value;
}

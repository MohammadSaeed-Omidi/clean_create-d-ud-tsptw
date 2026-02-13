// solve_instance.h
#ifndef SOLVE_INSTANCE_H
#define SOLVE_INSTANCE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>

#define MAX_PATH_LENGTH 512
#define CONCORDE_PATH "/home/lonely/Downloads/concorde/concorde-bin"

// Function declarations
int run_concorde(const char *tsp_filename);
int** read_edge_weights(const char* tsp_filename, int* dimension);
double calculate_tour_length(const char* sol_filename, const char* tsp_filename);
double solve_tsp_instance(const char* tsp_filename, const char* sol_filename);

#endif // SOLVE_INSTANCE_H
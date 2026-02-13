#ifndef CREATE_TSPTW_INSTANCES_H
#define CREATE_TSPTW_INSTANCES_H

#include <igraph/igraph.h>
#include <string.h>          
#include <unistd.h>   
#include <fcntl.h>    
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "config.h"

#define MAX_LINE_LENGTH 1000000
#define MAX_FILENAME_LENGTH 4096

typedef struct {
    igraph_integer_t from;
    igraph_integer_t to;
} Edge;

typedef struct {
    Edge* edges;
    size_t num_edges;
    igraph_integer_t num_vertices;
} GraphEdges;

// Function declarations
GraphEdges create_and_check_sparse_graph_instance(int n, float target_fraction);
void create_tsp_instance(const GraphEdges* graph_edges, long int num_vertices, int instance_num, float target_fraction, igraph_matrix_t* distances);
igraph_t read_tsp_file(const char* filename, igraph_vector_t* weights, char* problem_name);

void create_transformed_tsp_instance(long int num_vertices, int instance_num, float target_fraction, igraph_matrix_t* distances);

void create_tsptw_instance_file(igraph_matrix_t* distances, long int num_vertices, int instance_num, 
                               double* a_i, double* b_i, double T,
                               int width,
                               const char* instance_type);
                               
void create_tsptw_dataset_from_tsp(long int num_vertices, int instance_num, int* width_list, int num_widths, igraph_matrix_t* vanila_distances, igraph_matrix_t* shortest_distances);
int* find_second_nearest_neighbor_tour(igraph_matrix_t* distances, int n);
#endif

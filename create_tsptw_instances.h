#ifndef CREATE_OP_INSTANSES
#define CREATE_OP_INSTANSES

#include <igraph/igraph.h>
#include <string.h>          // For string functions
#include <unistd.h>   // For access() and getcwd()
#include <fcntl.h>    // For F_OK and R_OK constants
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX_LINE_LENGTH 1000000
#define MAX_FILENAME_LENGTH 512

typedef struct {
    igraph_integer_t from;
    igraph_integer_t to;
} Edge;

typedef struct {
    Edge* edges;
    size_t num_edges;
    igraph_integer_t num_vertices;
} GraphEdges;

void create_and_check_sparse_graph_instance(int n, float target_fraction, igraph_matrix_t* distances);

void create_directed_and_undirected_tsp_instance(int num_vertices, int instance_num, float target_fraction, igraph_matrix_t* distances, igraph_matrix_t* directed_distances, igraph_matrix_t* undirected_distances);

void create_directed_tsp_instance(int num_vertices, int instance_num, float target_fraction, igraph_matrix_t* directed_distances);
void create_undirected_tsp_instance(int num_vertices, int instance_num, float target_fraction, igraph_matrix_t* undirected_distances);
// create graph functions
void create_directed_weighted_graph(igraph_t* directed_graph,
                                    igraph_vector_int_t* directed_edges, 
                                    igraph_vector_t* directed_weights,
                                    igraph_matrix_t* directed_distances,  // Changed from int** to igraph_matrix_t*
                                    long int n_vertices);

void create_undirected_weighted_graph(igraph_t* undirected_graph,
                                    igraph_vector_int_t* undirected_edges, 
                                    igraph_vector_t* undirected_weights,
                                    igraph_matrix_t* undirected_distances,  // Changed from int** to igraph_matrix_t*
                                    long int n_vertices);

void transform_tsp_file(const char* output_filename, const char* problem_name, 
                       igraph_t* graph, igraph_vector_t* weights, 
                       igraph_matrix_t* distances);

void create_directed_transformed_tsp_instance(igraph_t* directed_graph, igraph_vector_t* directed_weights, igraph_matrix_t* directed_shortest_path_distances, int instance_num, float target_fraction);

void create_undirected_transformed_tsp_instance(igraph_t* undirected_graph, igraph_vector_t* undirected_weights, igraph_matrix_t* undirected_shortest_path_distances, int instance_num, float target_fraction);   


void create_undirected_predecessors_matrix(igraph_t* undirected_graph, igraph_vector_t* undirected_weights, int instance_num);

void create_tsptw_dataset_from_tsp(long int num_vertices, int instance_num, int* width_list, int num_widths, igraph_matrix_t* distances, igraph_matrix_t* shortest_distances, char* type);
int* find_second_nearest_neighbor_tour(igraph_matrix_t* distances, long int num_vertices);
double* calculate_visit_times(int* tour, igraph_matrix_t* distances, int num_vertices);
void generate_time_windows(double* visit_times, igraph_matrix_t* distances, int n, 
                          int width, double** a_i, double** b_i, double* T);

#endif

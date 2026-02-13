#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "config.h"
#include "create_tsptw_instances.h"
#include "solve_instance.h"

int main() {
    printf("Creating TSPTW instances using second nearest neighbor method\n");
    
    // Set random seed once at program start
    srand(time(NULL));
    
    // Time window widths to test (as in the paper)
    int width_list[] = {20*MAX_WEIGHT/100, 30*MAX_WEIGHT/100, 40*MAX_WEIGHT/100};
    int num_widths = sizeof(width_list) / sizeof(width_list[0]);
    
    int number_of_node_list[] = NUMBER_OF_NODES_LIST;
    int len_number_of_node_list = sizeof(number_of_node_list) / sizeof(number_of_node_list[0]);
    float target_fraction = 1.0 / 3.0;
    
    // Create directories if they don't exist
    system("mkdir -p instances");

    for (int i = 0; i < len_number_of_node_list; i++) {
        long int num_vertices = number_of_node_list[i];
        // Create directories
        char mkdir_command[512];
        snprintf(mkdir_command, sizeof(mkdir_command), "mkdir -p instances/%ld_nodes", 
                num_vertices, mkdir_command);
        system(mkdir_command);
        snprintf(mkdir_command, sizeof(mkdir_command), "mkdir -p instances/%ld_nodes/tsp", 
                num_vertices, mkdir_command);
        system(mkdir_command);
        snprintf(mkdir_command, sizeof(mkdir_command), "mkdir -p instances/%ld_nodes/tsp_tw", 
                num_vertices, mkdir_command);
        system(mkdir_command);
        snprintf(mkdir_command, sizeof(mkdir_command), "mkdir -p instances/%ld_nodes/tsp_tw/un", 
                num_vertices, mkdir_command);
        system(mkdir_command);
        snprintf(mkdir_command, sizeof(mkdir_command), "mkdir -p instances/%ld_nodes/tsp_tw/d", 
                num_vertices, mkdir_command);
        system(mkdir_command);


        for (int instance = 0; instance < NUM_INSTANCES; instance++) {
            printf("\n=== Creating instance %d with %d vertices ===\n", instance, number_of_node_list[i]);
            
            long int num_vertices = number_of_node_list[i];
            float target_fraction = 1.0 / 3.0;
            
            igraph_matrix_t distances;
            // Initialize and calculate shortest path distances
            igraph_matrix_init(&distances, num_vertices, num_vertices);

            // Fill the entire matrix with MAX_WEIGHT before passing to the function
            for (igraph_integer_t i = 0; i < num_vertices; i++) {
                for (igraph_integer_t j = 0; j < num_vertices; j++) {
                    MATRIX(distances, i, j) = MAX_WEIGHT;
                }
            }

            // Create graph and get edges
            create_and_check_sparse_graph_instance(num_vertices, target_fraction, &distances);

            printf("check matrix.\n");
            for (int i = 0; i < num_vertices; i++) {
                for (int j = 0; j < num_vertices; j++) {
                    double dist = MATRIX(distances, i, j);
                    printf("%-4.0f", dist);
                }
                printf("\n");
            }

            igraph_matrix_t directed_distances;
            // Initialize and calculate shortest path distances
            igraph_matrix_init(&directed_distances, num_vertices, num_vertices);
            
            igraph_matrix_t undirected_distances;
            // Initialize and calculate shortest path distances
            igraph_matrix_init(&undirected_distances, num_vertices, num_vertices);
            
            create_directed_and_undirected_tsp_instance(num_vertices, instance, target_fraction, &distances, &directed_distances, &undirected_distances);

            igraph_t directed_graph;
            igraph_vector_int_t directed_edges;
            igraph_vector_t directed_weights;
    
            create_directed_weighted_graph(&directed_graph, &directed_edges, &directed_weights, &directed_distances, num_vertices);
            // print_weighted_graph_info(&directed_graph, &directed_weights, "DIRECTED WEIGHTED GRAPH", 1);

            igraph_t undirected_graph;
            igraph_vector_int_t undirected_edges;
            igraph_vector_t undirected_weights;
            create_directed_weighted_graph(&undirected_graph, &undirected_edges, &undirected_weights, &undirected_distances, num_vertices);
            // print_weighted_graph_info(&undirected_graph, &undirected_weights, "UNDIRECTED WEIGHTED GRAPH", 0);
            
            // Create directed_transformed TSP instance
            igraph_matrix_t directed_shortest_path_distances;
            // Initialize and calculate shortest path distances
            igraph_matrix_init(&directed_shortest_path_distances, num_vertices, num_vertices);
            create_directed_transformed_tsp_instance(&directed_graph, &directed_weights, &directed_shortest_path_distances, instance, target_fraction);

            igraph_matrix_t undirected_shortest_path_distances;
            // Initialize and calculate shortest path distances
            igraph_matrix_init(&undirected_shortest_path_distances, num_vertices, num_vertices);
            create_undirected_transformed_tsp_instance(&undirected_graph, &undirected_weights, &undirected_shortest_path_distances, instance, target_fraction);

            create_undirected_predecessors_matrix(&undirected_graph, &undirected_weights, instance);
            create_directed_predecessors_matrix(&directed_graph, &directed_weights, instance);

            create_tsptw_dataset_from_tsp(num_vertices, instance, width_list, num_widths, &undirected_distances, &undirected_shortest_path_distances, "un");
            create_tsptw_dataset_from_tsp(num_vertices, instance, width_list, num_widths, &directed_distances, &directed_shortest_path_distances, "d");
            

        }
    }
    return 0;
}
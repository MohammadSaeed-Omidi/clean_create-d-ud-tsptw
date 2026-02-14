#include "create_tsptw_instances.h"
#include "config.h"

void create_and_check_sparse_graph_instance(int n, float target_fraction, igraph_matrix_t* distances) {
    igraph_t graph;
    igraph_vector_int_t edges_to_add;
    igraph_integer_t n_vertices = n;
    igraph_integer_t max_possible_edges, target_edges, edges_needed;
    
    igraph_rng_seed(igraph_rng_default(), time(NULL));
    
    // Step 1: Create a simple chain (1-2-3-...-n)
    igraph_vector_int_t chain_edges;
    igraph_vector_int_init(&chain_edges, (n_vertices - 1) * 2);
    
    for (igraph_integer_t i = 0; i < n_vertices - 1; i++) {
        VECTOR(chain_edges)[2 * i] = i;
        VECTOR(chain_edges)[2 * i + 1] = i + 1;
    }
    
    // Create graph from chain edges
    igraph_create(&graph, &chain_edges, n_vertices, IGRAPH_UNDIRECTED);
    printf("Chain edges: %ld\n", igraph_ecount(&graph));
    
    // Step 2: Calculate how many additional edges we need
    max_possible_edges = n_vertices * (n_vertices - 1) / 2;
    target_edges = (igraph_integer_t)(max_possible_edges * target_fraction);
    edges_needed = target_edges - igraph_ecount(&graph);
    
    printf("Target edges: %ld\n", target_edges);
    printf("Additional edges needed: %ld\n", edges_needed);
    
    // Step 3: Create list of all possible edges that are NOT in the chain
    igraph_vector_int_t possible_edges;
    igraph_vector_int_init(&possible_edges, 0);
    
    for (igraph_integer_t i = 0; i < n_vertices; i++) {
        for (igraph_integer_t j = i + 1; j < n_vertices; j++) {
            // Skip chain edges (i connected to i+1)
            if (j == i + 1) {
                continue;
            }
            igraph_vector_int_push_back(&possible_edges, i);
            igraph_vector_int_push_back(&possible_edges, j);
        }
    }
    
    printf("Possible non-chain edges available: %ld\n", igraph_vector_int_size(&possible_edges) / 2);
    
    // Step 4: Shuffle and add required number of edges
    if (edges_needed > 0) {
        // Check if we have enough possible edges
        igraph_integer_t possible_edge_pairs = igraph_vector_int_size(&possible_edges) / 2;
        if (edges_needed > possible_edge_pairs) {
            printf("Warning: Only %ld possible edges available, but %ld needed. Using all available.\n", 
                   possible_edge_pairs, edges_needed);
            edges_needed = possible_edge_pairs;
        }
        
        // Create an array of edge indices and shuffle those
        igraph_integer_t num_edge_pairs = igraph_vector_int_size(&possible_edges) / 2;
        igraph_vector_int_t edge_indices;
        igraph_vector_int_init_range(&edge_indices, 0, num_edge_pairs);

        // Shuffle the edge indices
        igraph_vector_int_shuffle(&edge_indices);
        
        // Create edges_to_add using shuffled indices
        igraph_vector_int_init(&edges_to_add, edges_needed * 2);
        for (igraph_integer_t i = 0; i < edges_needed; i++) {
            igraph_integer_t edge_idx = VECTOR(edge_indices)[i];
            VECTOR(edges_to_add)[i * 2] = VECTOR(possible_edges)[edge_idx * 2];
            VECTOR(edges_to_add)[i * 2 + 1] = VECTOR(possible_edges)[edge_idx * 2 + 1];
        }
        
        // Add the edges
        igraph_add_edges(&graph, &edges_to_add, NULL);
        
        printf("Added %ld random edges\n", edges_needed);
        
        // Cleanup
        igraph_vector_int_destroy(&edge_indices);
        igraph_vector_int_destroy(&edges_to_add);
    }
    
    // Step 5: Verify results
    igraph_bool_t connected;
    igraph_is_connected(&graph, &connected, IGRAPH_STRONG);
    
    printf("Final graph - Connected: %s, Edges: %ld (%.2f%% of possible)\n", 
           connected ? "YES" : "NO", 
           igraph_ecount(&graph),
           (double)igraph_ecount(&graph) / max_possible_edges * 100);
    
    // Step 6: Extract all edges and return them 

    // Get the number of edges in the graph
    igraph_integer_t edge_count = igraph_ecount(&graph);

    // Use a single random weight for each edge (consistent in both directions)
    srand(time(NULL));

    for (igraph_integer_t k = 0; k < edge_count; k++) {
        igraph_integer_t from, to;
        
        // Get the edge endpoints directly from the graph
        igraph_edge(&graph, k, &from, &to);
        
        // For directed graph
        MATRIX(*distances, from, to)= (rand() % MAX_WEIGHT) + 1;
        MATRIX(*distances, to, from) = (rand() % MAX_WEIGHT) + 1;  // Same weight for undirected edge
    }

   
    // Cleanup
    igraph_vector_int_destroy(&chain_edges);
    igraph_vector_int_destroy(&possible_edges);
    igraph_destroy(&graph);
    
}

void create_directed_and_undirected_tsp_instance(int num_vertices, int instance_num, float target_fraction, igraph_matrix_t* distances, igraph_matrix_t* directed_distances, igraph_matrix_t* undirected_distances) {

    // Fill the Directed matrix
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < num_vertices; j++) {
            MATRIX(*directed_distances, i, j) = MATRIX(*distances, i, j)*10;  // Diagonal elements
        }
    }

    // Fill the Undirected matrix
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < num_vertices; j++) {
            MATRIX(*undirected_distances, i, j) = (MATRIX(*distances, i, j) + MATRIX(*distances, j, i))*5;
        }
    }


    // Create vanila directed TSP instance
    create_directed_tsp_instance(num_vertices, instance_num, target_fraction, directed_distances);

    // Create vanila directed TSP instance
    create_undirected_tsp_instance(num_vertices, instance_num, target_fraction, undirected_distances);

}

void create_directed_tsp_instance(int num_vertices, int instance_num, float target_fraction, igraph_matrix_t* directed_distances) {
    char filename[100];
    snprintf(filename, sizeof(filename), "instances/%d_nodes/tsp/vanila_directed_tsp_instance_n_%d_i_%d.tsp", 
             num_vertices, num_vertices, instance_num);

    
    FILE* file = fopen(filename, "w");
    if (!file) {
        printf("Error: Could not create file %s\n", filename);
        return;
    }
    
    // Write TSP header in Concorde format
    fprintf(file, "NAME: vanila_directed_tsp_instance_%d\n", instance_num);
    fprintf(file, "TYPE: TSP\n");
    // fprintf(file, "COMMENT: Generated from sparse graph with %d vertices and %ld edges\n", 
    //         num_vertices, (long)(num_vertices * (num_vertices - 1) / 2 * target_fraction));
    fprintf(file, "DIMENSION: %d\n", num_vertices);
    fprintf(file, "EDGE_WEIGHT_TYPE: EXPLICIT\n");
    fprintf(file, "EDGE_WEIGHT_FORMAT: FULL_MATRIX\n");
    fprintf(file, "EDGE_WEIGHT_SECTION\n");
    
    // Write FULL_MATRIX format (n x n matrix)
    for (igraph_integer_t i = 0; i < num_vertices; i++) {
        for (igraph_integer_t j = 0; j < num_vertices; j++) {
            double dist = MATRIX(*directed_distances, i, j);
            fprintf(file, "%.0f ", dist);
            
            // Add space between numbers, but not after the last number in row
            if (j < num_vertices - 1) {
                fprintf(file, " ");
            }
        }
        
        // Add newline after each row
        fprintf(file, "\n");
    }
    
    fprintf(file, "EOF\n");

    fclose(file);
    printf("TSP instance saved as: %s\n", filename);
}

void create_undirected_tsp_instance(int num_vertices, int instance_num, float target_fraction, igraph_matrix_t* undirected_distances) {
    char filename[100];
    snprintf(filename, sizeof(filename), "instances/%d_nodes/tsp/vanila_undirected_tsp_instance_n_%d_i_%d.tsp", 
             num_vertices, num_vertices, instance_num);
    
    FILE* file = fopen(filename, "w");
    if (!file) {
        printf("Error: Could not create file %s\n", filename);
        return;
    }
    
    // Write TSP header in Concorde format
    fprintf(file, "NAME: vanila_undirected_tsp_instance_%d\n", instance_num);
    fprintf(file, "TYPE: TSP\n");
    fprintf(file, "COMMENT: Generated from sparse graph with %d vertices and %ld edges\n", 
            num_vertices, (long)(num_vertices * (num_vertices - 1) / 2 * target_fraction));
    fprintf(file, "DIMENSION: %d\n", num_vertices);
    fprintf(file, "EDGE_WEIGHT_TYPE: EXPLICIT\n");
    fprintf(file, "EDGE_WEIGHT_FORMAT: FULL_MATRIX\n");
    fprintf(file, "EDGE_WEIGHT_SECTION\n");
    
    // Write FULL_MATRIX format (n x n matrix)
    for (igraph_integer_t i = 0; i < num_vertices; i++) {
        for (igraph_integer_t j = 0; j < num_vertices; j++) {
            double dist = MATRIX(*undirected_distances, i, j);
            fprintf(file, "%.0f ", dist);
            
            // Add space between numbers, but not after the last number in row
            if (j < num_vertices - 1) {
                fprintf(file, " ");
            }
        }
        
        // Add newline after each row
        fprintf(file, "\n");
    }
    
    fprintf(file, "EOF\n");

    fclose(file);
    printf("TSP instance saved as: %s\n", filename);
}

void create_directed_weighted_graph(igraph_t* directed_graph,
                                    igraph_vector_int_t* directed_edges, 
                                    igraph_vector_t* directed_weights,
                                    igraph_matrix_t* directed_distances,  // Changed from int** to igraph_matrix_t*
                                    long int n_vertices) {
    
    
    
    // Initialize edges vector
    igraph_vector_int_init(directed_edges, 0);
    
    // Initialize weights vector
    igraph_vector_init(directed_weights, 0);
    
    // Initialize empty directed graph
    igraph_empty(directed_graph, n_vertices, IGRAPH_DIRECTED);
    
    // Add arcs from matrix (all non-MAX_WEIGHT entries)
    for (int i = 0; i < n_vertices; i++) {
        for (int j = 0; j < n_vertices; j++) {
            double weight = MATRIX(*directed_distances, i, j);
            // Add arc if weight is not MAX_WEIGHT and not self-loop
            if (i != j && weight != MAX_WEIGHT*10) {
                igraph_vector_int_push_back(directed_edges, i);
                igraph_vector_int_push_back(directed_edges, j);
                igraph_vector_push_back(directed_weights, weight);
            }
        }
    }
    
    // Add edges to graph
    if (igraph_vector_int_size(directed_edges) > 0) {
        igraph_add_edges(directed_graph, directed_edges, NULL);
    }
    
}

void create_undirected_weighted_graph(igraph_t* undirected_graph,
                                    igraph_vector_int_t* undirected_edges, 
                                    igraph_vector_t* undirected_weights,
                                    igraph_matrix_t* undirected_distances,  // Changed from int** to igraph_matrix_t*
                                    long int n_vertices){
    
    // Initialize edges vector
    igraph_vector_int_init(undirected_edges, 0);
    
    // Initialize weights vector
    igraph_vector_init(undirected_weights, 0);
    
    // Initialize empty directed graph
    igraph_empty(undirected_graph, n_vertices, IGRAPH_UNDIRECTED);
    
    // Add arcs from matrix (all non-MAX_WEIGHT entries)
    for (int i = 0; i < n_vertices; i++) {
        for (int j = 0; j < n_vertices; j++) {
            double weight = MATRIX(*undirected_distances, i, j);
            // Add arc if weight is not MAX_WEIGHT and not self-loop
            if (i != j && weight != MAX_WEIGHT*10) {
                igraph_vector_int_push_back(undirected_edges, i);
                igraph_vector_int_push_back(undirected_edges, j);
                igraph_vector_push_back(undirected_weights, weight);
            }
        }
    }
    
    // Add edges to graph
    if (igraph_vector_int_size(undirected_edges) > 0) {
        igraph_add_edges(undirected_graph, undirected_edges, NULL);
    }
    
}


void transform_tsp_file(const char* output_filename, const char* problem_name, 
                       igraph_t* graph, igraph_vector_t* weights, 
                       igraph_matrix_t* distances) {

    igraph_integer_t n_vertices = igraph_vcount(graph);
    
    FILE* outfile = fopen(output_filename, "w");
    if (!outfile) {
        fprintf(stderr, "Error creating output file %s\n", output_filename);
        perror("fopen");
        return;
    }
    
    // Write header information
    fprintf(outfile, "NAME: walk_%s\n", problem_name);
    fprintf(outfile, "TYPE: TSP\n");
    fprintf(outfile, "COMMENT: Transformed TSP instance with %ld nodes (using shortest path distances)\n", n_vertices);
    fprintf(outfile, "DIMENSION: %ld\n", n_vertices);
    fprintf(outfile, "EDGE_WEIGHT_TYPE: EXPLICIT\n");
    fprintf(outfile, "EDGE_WEIGHT_FORMAT: FULL_MATRIX\n");
    fprintf(outfile, "EDGE_WEIGHT_SECTION\n");
    
    // Write the distance matrix
    for (int i = 0; i < n_vertices; i++) {
        for (int j = 0; j < n_vertices; j++) {
            double dist = MATRIX(*distances, i, j);
            if (dist > 1e9 || i == j) {
                fprintf(outfile, "%d ", MAX_WEIGHT*10);
            } else {
                fprintf(outfile, "%.0f ", dist);
            }
        }
        fprintf(outfile, "\n");
    }
    
    fprintf(outfile, "EOF\n");
    fclose(outfile);
    printf("Created walk file: %s\n", output_filename);
}

void create_directed_transformed_tsp_instance(igraph_t* directed_graph, igraph_vector_t* directed_weights, igraph_matrix_t* directed_shortest_path_distances, int instance_num, float target_fraction) {
    igraph_integer_t n_vertices = igraph_vcount(directed_graph);
    
    igraph_shortest_paths_dijkstra(directed_graph, directed_shortest_path_distances, igraph_vss_all(), 
                                  igraph_vss_all(), directed_weights, IGRAPH_OUT);  // IGRAPH_ALL undircted    IGRAPH_OUT   directed
    
    printf("transformed matrix: \n");
    for (int i = 0; i < n_vertices; i++) {
        for (int j = 0; j < n_vertices; j++) {
            double dist = MATRIX(*directed_shortest_path_distances, i, j);
            printf("%-5d", (int)dist);
        }
        printf("\n");
    }
    // Create output filename (replace "vanila" with "transformed")
    char output_filename[MAX_FILENAME_LENGTH];
    snprintf(output_filename, sizeof(output_filename), "instances/%ld_nodes/tsp/transformed_directed_tsp_instance_n_%ld_i_%d.tsp", n_vertices, n_vertices, instance_num);
    
    char problem_name[100];
    snprintf(problem_name, sizeof(problem_name), "transformed_directed_tsp_instance_n_%ld_i_%d.tsp", n_vertices, instance_num);

    // Transform and save the file
    transform_tsp_file(output_filename, problem_name, directed_graph, directed_weights, directed_shortest_path_distances);
    
    printf("Transformed TSP instance saved as: %s\n", output_filename);
}

void create_undirected_transformed_tsp_instance(igraph_t* undirected_graph, igraph_vector_t* undirected_weights, igraph_matrix_t* undirected_shortest_path_distances, int instance_num, float target_fraction) {
    igraph_integer_t n_vertices = igraph_vcount(undirected_graph);
    
    igraph_shortest_paths_dijkstra(undirected_graph, undirected_shortest_path_distances, igraph_vss_all(), 
                                  igraph_vss_all(), undirected_weights, IGRAPH_ALL);  // IGRAPH_ALL undircted    IGRAPH_OUT   directed
    
    printf("transformed matrix: \n");
    for (int i = 0; i < n_vertices; i++) {
        for (int j = 0; j < n_vertices; j++) {
            double dist = MATRIX(*undirected_shortest_path_distances, i, j);
            printf("%-5d", (int)dist);
        }
        printf("\n");
    }
    // Create output filename (replace "vanila" with "transformed")
    char output_filename[MAX_FILENAME_LENGTH];
    snprintf(output_filename, sizeof(output_filename), "instances/%ld_nodes/tsp/transformed_undirected_tsp_instance_n_%ld_i_%d.tsp", n_vertices, n_vertices, instance_num);
    
    char problem_name[100];
    snprintf(problem_name, sizeof(problem_name), "transformed_undirected_tsp_instance_n_%ld_i_%d.tsp", n_vertices, instance_num);

    // Transform and save the file
    transform_tsp_file(output_filename, problem_name, undirected_graph, undirected_weights, undirected_shortest_path_distances);
    
    printf("Transformed TSP instance saved as: %s\n", output_filename);
}

void create_undirected_predecessors_matrix(igraph_t* undirected_graph, igraph_vector_t* undirected_weights, int instance_num) {
    igraph_integer_t num_vertices = igraph_vcount(undirected_graph);

    igraph_matrix_t undirected_predecessors;
    // Initialize and calculate shortest path distances
    igraph_matrix_init(&undirected_predecessors, num_vertices, num_vertices);
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < num_vertices; j++) {
            MATRIX(undirected_predecessors, i, j) = -1;
        }
    }
    
    // Now calculate paths for each source to fill the predecessor matrix
    for (int source = 0; source < num_vertices; source++) {
        igraph_vector_int_list_t vertices;
        igraph_vector_int_list_t edges;
        igraph_vector_int_t parents;
        igraph_vector_int_t inbound_edges;
        
        // Initialize the vectors with proper capacity
        igraph_vector_int_list_init(&vertices, num_vertices);
        igraph_vector_int_list_init(&edges, num_vertices);
        igraph_vector_int_init(&parents, num_vertices);  // Initialize with capacity for n vertices
        igraph_vector_int_init(&inbound_edges, num_vertices);
        
        // Get shortest paths from source to all vertices
        igraph_error_t result = igraph_get_shortest_paths_dijkstra(
            undirected_graph,           // graph
            &vertices,       // vertices (paths)
            &edges,          // edges
            source,          // from
            igraph_vss_all(), // to all vertices
            undirected_weights,         // edge weights
            IGRAPH_ALL,      // directedness
            &parents,        // parents (predecessors)
            &inbound_edges   // inbound edges
        );
        
        if (result != IGRAPH_SUCCESS) {
            fprintf(stderr, "Error calculating shortest paths from vertex %d\n", source);
            // Cleanup and continue
            igraph_vector_int_list_destroy(&vertices);
            igraph_vector_int_list_destroy(&edges);
            igraph_vector_int_destroy(&parents);
            igraph_vector_int_destroy(&inbound_edges);
            continue;
        }
        
        // Store parents in matrix - make sure we have the right number of elements
        if (igraph_vector_int_size(&parents) == num_vertices) {
            for (int target = 0; target < num_vertices; target++) {
                MATRIX(undirected_predecessors, source, target) = VECTOR(parents)[target];
            }
        } else {
            fprintf(stderr, "Warning: parents vector size %ld doesn't match n=%ld for source %d\n", 
                   igraph_vector_int_size(&parents), num_vertices, source);
        }
        
        // Cleanup
        igraph_vector_int_list_destroy(&vertices);
        igraph_vector_int_list_destroy(&edges);
        igraph_vector_int_destroy(&parents);
        igraph_vector_int_destroy(&inbound_edges);
    }

    char dir_path[100];
    // Create the same directory structure as TSP files
    snprintf(dir_path, sizeof(dir_path), "instances/%ld_nodes/tsp_tw/%s/%d", num_vertices, "un", instance_num);

    // Create directories if they don't exist
    char mkdir_command[200];
    snprintf(mkdir_command, sizeof(mkdir_command), "mkdir -p %s", dir_path);
    system(mkdir_command);

    // Also save the paths to a separate file
    char filename[MAX_FILENAME_LENGTH];
    
    // Create .tsptw filename
    snprintf(filename, sizeof(filename), "instances/%ld_nodes/tsp_tw/%s/%d/predecessors_matrix_n_%ld_i_%d.paths", 
             num_vertices, "un", instance_num, num_vertices, instance_num);
    
    FILE* paths_file = fopen(filename, "w");
    fprintf(paths_file, "SHORTEST_PATHS_MATRIX for n_%ld_i_%d\n", num_vertices, instance_num);
    fprintf(paths_file, "DIMENSION: %ld\n", num_vertices);
    fprintf(paths_file, "PREDECESSOR_MATRIX_SECTION\n");
    
    // Write predecessor matrix
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < num_vertices; j++) {
            fprintf(paths_file, "%ld ", (long)MATRIX(undirected_predecessors, i, j));
        }
        fprintf(paths_file, "\n");
    }
    
    fprintf(paths_file, "PATH_EXAMPLES_SECTION\n");
    
    fprintf(paths_file, "EOF\n");
    fclose(paths_file);
    printf("Created paths file: %s\n", filename);
}

void create_directed_predecessors_matrix(igraph_t* directed_graph, igraph_vector_t* directed_weights, int instance_num) {
    igraph_integer_t num_vertices = igraph_vcount(directed_graph);

    igraph_matrix_t directed_predecessors;
    // Initialize and calculate shortest path distances
    igraph_matrix_init(&directed_predecessors, num_vertices, num_vertices);
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < num_vertices; j++) {
            MATRIX(directed_predecessors, i, j) = -1;
        }
    }
    
    // Now calculate paths for each source to fill the predecessor matrix
    for (int source = 0; source < num_vertices; source++) {
        igraph_vector_int_list_t vertices;
        igraph_vector_int_list_t edges;
        igraph_vector_int_t parents;
        igraph_vector_int_t inbound_edges;
        
        // Initialize the vectors with proper capacity
        igraph_vector_int_list_init(&vertices, num_vertices);
        igraph_vector_int_list_init(&edges, num_vertices);
        igraph_vector_int_init(&parents, num_vertices);  // Initialize with capacity for n vertices
        igraph_vector_int_init(&inbound_edges, num_vertices);
        
        // Get shortest paths from source to all vertices
        igraph_error_t result = igraph_get_shortest_paths_dijkstra(
            directed_graph,           // graph
            &vertices,       // vertices (paths)
            &edges,          // edges
            source,          // from
            igraph_vss_all(), // to all vertices
            directed_weights,         // edge weights
            IGRAPH_ALL,      // directedness
            &parents,        // parents (predecessors)
            &inbound_edges   // inbound edges
        );
        
        if (result != IGRAPH_SUCCESS) {
            fprintf(stderr, "Error calculating shortest paths from vertex %d\n", source);
            // Cleanup and continue
            igraph_vector_int_list_destroy(&vertices);
            igraph_vector_int_list_destroy(&edges);
            igraph_vector_int_destroy(&parents);
            igraph_vector_int_destroy(&inbound_edges);
            continue;
        }
        
        // Store parents in matrix - make sure we have the right number of elements
        if (igraph_vector_int_size(&parents) == num_vertices) {
            for (int target = 0; target < num_vertices; target++) {
                MATRIX(directed_predecessors, source, target) = VECTOR(parents)[target];
            }
        } else {
            fprintf(stderr, "Warning: parents vector size %ld doesn't match n=%ld for source %d\n", 
                   igraph_vector_int_size(&parents), num_vertices, source);
        }
        
        // Cleanup
        igraph_vector_int_list_destroy(&vertices);
        igraph_vector_int_list_destroy(&edges);
        igraph_vector_int_destroy(&parents);
        igraph_vector_int_destroy(&inbound_edges);
    }

    char dir_path[100];
    // Create the same directory structure as TSP files
    snprintf(dir_path, sizeof(dir_path), "instances/%ld_nodes/tsp_tw/%s/%d", num_vertices, "d", instance_num);

    // Create directories if they don't exist
    char mkdir_command[200];
    snprintf(mkdir_command, sizeof(mkdir_command), "mkdir -p %s", dir_path);
    system(mkdir_command);

    // Also save the paths to a separate file
    char filename[MAX_FILENAME_LENGTH];
    
    // Create .tsptw filename
    snprintf(filename, sizeof(filename), "instances/%ld_nodes/tsp_tw/%s/%d/predecessors_matrix_n_%ld_i_%d.paths", 
             num_vertices, "un", instance_num, num_vertices, instance_num);
    
    FILE* paths_file = fopen(filename, "w");
    fprintf(paths_file, "SHORTEST_PATHS_MATRIX for n_%ld_i_%d\n", num_vertices, instance_num);
    fprintf(paths_file, "DIMENSION: %ld\n", num_vertices);
    fprintf(paths_file, "PREDECESSOR_MATRIX_SECTION\n");
    
    // Write predecessor matrix
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < num_vertices; j++) {
            fprintf(paths_file, "%ld ", (long)MATRIX(directed_predecessors, i, j));
        }
        fprintf(paths_file, "\n");
    }
    
    fprintf(paths_file, "PATH_EXAMPLES_SECTION\n");
    
    fprintf(paths_file, "EOF\n");
    fclose(paths_file);
    printf("Created paths file: %s\n", filename);
}

// Fixed function to find second nearest neighbor tour
int* find_second_nearest_neighbor_tour(igraph_matrix_t* distances, long int num_vertices) {
    int* tour = (int*)malloc(num_vertices * sizeof(int));
    int* visited = (int*)calloc(num_vertices, sizeof(int));
    for (int k1 = 0; k1 < num_vertices; k1++) {
        visited[k1] = 0;
    }
    
    // Start at depot (node 0)
    tour[0] = 0;
    visited[0] = 1;
    int current_node = 0;
    
    // printf("Starting tour construction with %d nodes\n", n);
    
    for (int step = 1; step < num_vertices; step++) {
        // Find distances from current node to all unvisited nodes
        double min_dist = MAX_WEIGHT*100;
        double second_min_dist = MAX_WEIGHT*100;
        int min_node = -1;
        int second_min_node = -1;
        int unvisited_count = 0;
        
        // printf("Step %d: Current node = %d, looking for unvisited nodes...\n", step, current_node);
        
        for (int j = 0; j < num_vertices; j++) {
            if (!visited[j] && j != current_node) {
                unvisited_count++;
                int dist =  MATRIX(*distances, current_node, j);
                // printf("  Unvisited node %d, distance: %.2f\n", j, dist);
                
                if (dist < min_dist) {
                    second_min_dist = min_dist;
                    second_min_node = min_node;
                    min_dist = dist;
                    min_node = j;
                } else if (dist < second_min_dist) {
                    second_min_dist = dist;
                    second_min_node = j;
                }
            }
        }
        
        // printf("  Found %d unvisited nodes\n", unvisited_count);
        // printf("  Nearest: node %d (dist: %.2f)\n", min_node, min_dist);
        // printf("  Second nearest: node %d (dist: %.2f)\n", second_min_node, second_min_dist);
        
        // Choose the next node
        int next_node;
        if (unvisited_count > 1 && second_min_node != -1) {
            // If we have at least 2 unvisited nodes, choose second nearest
            next_node = second_min_node;
            // printf("  -> Choosing SECOND nearest: node %d\n", next_node);
        } else if (min_node != -1) {
            // If only one unvisited node or no second nearest, choose nearest
            next_node = min_node;
            // printf("  -> Choosing nearest: node %d\n", next_node);
        } else {
            // This should not happen in a connected graph
            printf("  ERROR: No valid next node found!\n");
            // Find ANY unvisited node as fallback
            for (int j = 0; j < num_vertices; j++) {
                if (!visited[j] && j != current_node) {
                    next_node = j;
                    printf("  -> Fallback to node %d\n", next_node);
                    break;
                }
            }
        }
        
        // Safety check
        if (next_node < 0 || next_node >= num_vertices || visited[next_node]) {
            printf("  ERROR: Invalid next node %d!\n", next_node);
            // Emergency fallback - find first unvisited node
            for (int j = 0; j < num_vertices; j++) {
                if (!visited[j] && j != current_node) {
                    next_node = j;
                    printf("  -> Emergency fallback to node %d\n", next_node);
                    break;
                }
            }
        }
        
        tour[step] = next_node;
        visited[next_node] = 1;
        current_node = next_node;
        
        // // Print partial tour
        // printf("  Partial tour: ");
        // for (int i = 0; i <= step; i++) {
        //     printf("%d ", tour[i]);
        // }
        // printf("\n");
    }
    
    free(visited);
    
    // // Print final tour
    // printf("Final tour: ");
    // for (int i = 0; i < n; i++) {
    //     printf("%d ", tour[i]);
    // }
    // printf("\n");
    
    return tour;
}

// Function to calculate visit times along a tour
double* calculate_visit_times(int* tour, igraph_matrix_t* distances, int num_vertices) {
    double* visit_times = (double*)malloc(num_vertices * sizeof(double));
    
    // Initialize all visit times to -1 (unvisited)
    for (int i = 0; i < num_vertices; i++) {
        visit_times[i] = -1.0;
    }
    
    // Start at depot (first node in tour) at time 0
    visit_times[tour[0]] = 0.0;
    double current_time = 0.0;
    
    // Calculate arrival times at each node in the tour order
    for (int i = 1; i < num_vertices; i++) {
        int prev_node = tour[i-1];
        int curr_node = tour[i];
        double travel_time = MATRIX(*distances, prev_node, curr_node);
        current_time += travel_time;
        visit_times[curr_node] = current_time;
    }
    
    return visit_times;
}

void generate_time_windows(double* visit_times, igraph_matrix_t* distances, int n, 
                          int width, double** a_i, double** b_i, double* T) {
    *a_i = (double*)malloc(n * sizeof(double));
    *b_i = (double*)malloc(n * sizeof(double));
    
    // Generate time windows for all nodes except depot
    for (int i = 0; i < n; i++) {
        if (i == 0) { // Depot
            (*a_i)[0] = 0.0;
        } else {
            // Random offset within [0, width]
            double offset_lower = (double)(rand() % (width + 1));
            double offset_upper = (double)(rand() % (width + 1));
            
            (*a_i)[i] = visit_times[i] - offset_lower;
            if ((*a_i)[i] < 0){
                (*a_i)[i] = 0;
            }
            (*b_i)[i] = visit_times[i] + offset_upper;
            
            // Ensure a_i <= b_i
            if ((*a_i)[i] > (*b_i)[i]) {
                double temp = (*a_i)[i];
                (*a_i)[i] = (*b_i)[i];
                (*b_i)[i] = temp;
            }
        }
    }
    
    // Calculate depot's latest time based on return from all nodes
    double max_return_time = 0.0;
    for (int i = 1; i < n; i++) {
        double return_time = (*b_i)[i] + MATRIX(*distances, i, 0);
        if (return_time > max_return_time) {
            max_return_time = return_time;
        }
    }
    
    (*b_i)[0] = max_return_time;
    
    // Calculate total resource T
    *T = 0.0;
    for (int i = 1; i < n; i++) {
        double candidate_T = (*b_i)[i] + MATRIX(*distances, i, 0);
        if (candidate_T > *T) {
            *T = candidate_T;
        }
    }
}

void create_transformed_tsptw_instance_file(igraph_matrix_t* distances, long int num_vertices, int instance_num, 
                               double* a_i, double* b_i, double T,
                               int width,
                               const char* instance_type,  char* type) {
    char filename[256];
    char dir_path[100];
    
    // Create the same directory structure as TSP files
    snprintf(dir_path, sizeof(dir_path), "instances/%ld_nodes/tsp_tw/%s/%d", num_vertices, type, instance_num);
    
    // Create .tsptw filename
    snprintf(filename, sizeof(filename), "%s/%s_%s_tsptw_instance_n_%ld_w%d_i_%d.tsptw", 
             dir_path, instance_type, type, num_vertices, width, instance_num);
    
    // // Create directories if they don't exist
    // char mkdir_command[200];
    // snprintf(mkdir_command, sizeof(mkdir_command), "mkdir -p %s", dir_path);
    // system(mkdir_command);
    
    FILE* file = fopen(filename, "w");
    if (!file) {
        printf("Error: Could not create file %s\n", filename);
        return;
    }
    
    int n = (int)num_vertices;
    
    // Write header
    fprintf(file, "NAME: %s_tsptw_instance_%d\n", instance_type, instance_num);
    fprintf(file, "TYPE: TSPTW\n");
    fprintf(file, "INSTANCE_TYPE: %s\n", instance_type);
    fprintf(file, "DIMENSION: %d\n", n);
    fprintf(file, "TOTAL_RESOURCE: %.2f\n", T);
    fprintf(file, "TIME_WINDOW_WIDTH: %d\n", width);
    fprintf(file, "EDGE_WEIGHT_TYPE: EXPLICIT\n");
    fprintf(file, "EDGE_WEIGHT_FORMAT: FULL_MATRIX\n");
    fprintf(file, "NODE_DATA_SECTION\n");
    
    // Write node data: node_id ready_time due_date service_time
    for (int i = 0; i < n; i++) {
        fprintf(file, "%d %.2f %.2f %d\n", i, a_i[i], b_i[i], SERVICE_TIME);
    }
    
    fprintf(file, "EDGE_WEIGHT_SECTION\n");
    
    // Write distance matrix with MAX_WEIGHT for diagonal and actual distances for others
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // Write MAX_WEIGHT for diagonal, actual distance for others
            if (i == j) {
                fprintf(file, "%d", MAX_WEIGHT*10);
            } else {
                double dist = MATRIX(*distances, i, j);
                // Ensure we write the distance as integer (like TSP format)
                fprintf(file, "%.0f", dist);
            }
            
            // Add space between numbers, but not after the last number in row
            if (j < n - 1) {
                fprintf(file, " ");
            }
        }
        
        // Add newline after each row
        fprintf(file, "\n");
    }
    
    fprintf(file, "EOF\n");
    fclose(file);
    printf("%s TSPTW instance saved as: %s\n", instance_type, filename);
}



// Main function to create TSPTW dataset from any TSP file
void create_tsptw_dataset_from_tsp(long int num_vertices, int instance_num, int* width_list, int num_widths, igraph_matrix_t* distances, igraph_matrix_t* shortest_distances, char* type) {
    // Create transformed TSP instance
    char dir_path[100];
    snprintf(dir_path, sizeof(dir_path), "instances/%ld_nodes/%d", 
        num_vertices, instance_num);

    printf("\n=== Creating TSPTW dataset from: %s ===\n", dir_path);
    
    // Find second nearest neighbor tour
    int* tour = find_second_nearest_neighbor_tour(distances, num_vertices);

    // Print the tour for verification
    printf("Second nearest neighbor tour: ");
    for (int i = 0; i < num_vertices; i++) {
        printf("%d ", tour[i]);
    }
    printf("\n");

    // Calculate visit times along the tour
    double* visit_times = calculate_visit_times(tour, shortest_distances, num_vertices);
    
    // Print visit times for verification
    printf("#######Visit times: \n");
    for (int i = 0; i < num_vertices; i++) {
        printf("%d: %.1f\n", i, visit_times[i]);
    }
    printf("\n");
    
    // Generate instances for different time window widths
    for (int w = 0; w < num_widths; w++) {
        int width = width_list[w];
        
        double* a_i, * b_i;
        double T;
        
        // Generate time windows
        generate_time_windows(visit_times, shortest_distances, num_vertices, width, &a_i, &b_i, &T);

        // Print times window for verification
        printf("Time Window: \n");
        for (int i = 0; i < num_vertices; i++) {
            printf("%d: %.1f, %.1f\n", i, a_i[i], b_i[i]);
        }
        printf("\n");
        
        create_transformed_tsptw_instance_file(distances, num_vertices, instance_num, 
                                 a_i, b_i, T, width, "vanilla", type);

        create_transformed_tsptw_instance_file(shortest_distances, num_vertices, instance_num, 
                                 a_i, b_i, T, width, "transformed", type);

        free(a_i);
        free(b_i);
    }
    
    // // Clean up
    // free(tour);
    // free(visit_times);
    
    
    // printf("TSPTW dataset creation completed for %s instance %d\n", instance_type, instance_num);
}
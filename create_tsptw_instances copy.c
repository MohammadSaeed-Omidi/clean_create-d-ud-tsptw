#include "create_tsptw_instances.h"
#include "config.h"

GraphEdges create_and_check_sparse_graph_instance(int n, float target_fraction) {
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
    GraphEdges result;
    result.num_vertices = n_vertices;
    result.num_edges = igraph_ecount(&graph);
    result.edges = (Edge*)malloc(result.num_edges * sizeof(Edge));
    
    for (igraph_integer_t i = 0; i < result.num_edges; i++) {
        igraph_edge(&graph, i, &result.edges[i].from, &result.edges[i].to);
    }
    
    // // Print final edge list
    // printf("Final edge list:\n");
    // for (size_t i = 0; i < result.num_edges; i++) {
    //     printf("  %ld -- %ld\n", (long)result.edges[i].from, (long)result.edges[i].to);
    // }
    
    // Cleanup
    igraph_vector_int_destroy(&chain_edges);
    igraph_vector_int_destroy(&possible_edges);
    igraph_destroy(&graph);
    
    return result;
}

void create_tsp_instance(const GraphEdges* graph_edges, long int num_vertices, int instance_num, float target_fraction, igraph_matrix_t* distances) {
    char filename[200];
    char dir_path[100];

    // Create node-specific directory path
    snprintf(dir_path, sizeof(dir_path), "instances/%ld_nodes/%d", num_vertices, instance_num);
    snprintf(filename, sizeof(filename), "%s/vanila_tsp_instance_n_%ld_i_%d.tsp", 
             dir_path, num_vertices, instance_num);

    // Create directories if they don't exist
    char mkdir_command[200];
    snprintf(mkdir_command, sizeof(mkdir_command), "mkdir -p %s", dir_path);
    system(mkdir_command);
    
    FILE* file = fopen(filename, "w");
    if (!file) {
        printf("Error: Could not create file %s\n", filename);
        return;
    }
    
    // Write TSP header in Concorde format
    fprintf(file, "NAME: vanila_tsp_instance_%d\n", instance_num);
    fprintf(file, "TYPE: TSP\n");
    fprintf(file, "COMMENT: Generated from sparse graph with %ld vertices and %ld edges\n", 
            (long)graph_edges->num_vertices, (long)graph_edges->num_edges);
    fprintf(file, "DIMENSION: %ld\n", (long)graph_edges->num_vertices);
    fprintf(file, "EDGE_WEIGHT_TYPE: EXPLICIT\n");
    fprintf(file, "EDGE_WEIGHT_FORMAT: FULL_MATRIX\n");
    fprintf(file, "EDGE_WEIGHT_SECTION\n");
    
    igraph_integer_t n = graph_edges->num_vertices;
    
    // Precompute weights for all edges
    int* edge_weights = (int*)malloc(graph_edges->num_edges * sizeof(int));
    for (size_t i = 0; i < graph_edges->num_edges; i++) {
        edge_weights[i] = 1 + rand() % (MAX_WEIGHT - 1);
    }

    // for (size_t i = 0; i < graph_edges->num_edges; i++) {
    //     printf("edgeweitght[%ld]: %d\n", i, edge_weights[i]);
    // }

    // Write FULL_MATRIX format (n x n matrix)
    for (igraph_integer_t i = 0; i < n; i++) {
        for (igraph_integer_t j = 0; j < n; j++) {
            int weight = MAX_WEIGHT; // Default: no edge (or 0 for same vertex)
            
            if (i != j) { // Skip diagonal (usually 0 in TSP)
                // Search for edge i-j in our edge list
                for (size_t k = 0; k < graph_edges->num_edges; k++) {
                    igraph_integer_t from = graph_edges->edges[k].from;
                    igraph_integer_t to = graph_edges->edges[k].to;
                    
                    if ((from == i && to == j) || (from == j && to == i)) {
                        weight = edge_weights[k];
                        break;
                    }
                }
            }
            
            fprintf(file, "%d", weight);
            MATRIX(*distances, i, j) = weight;
            
            // Add space between numbers, but not after the last number in row
            if (j < n - 1) {
                fprintf(file, " ");
            }
        }
        
        // Add newline after each row
        fprintf(file, "\n");
    }
    
    fprintf(file, "EOF\n");
    
    free(edge_weights);
    fclose(file);
    printf("TSP instance saved as: %s\n", filename);
}

// Function to read TSP file and create graph
igraph_t read_tsp_file(const char* filename, igraph_vector_t* weights, char* problem_name) {
    printf("Attempting to open file: %s\n", filename);  // Added newline
    
    // Check if file exists first
    if (access(filename, F_OK) == -1) {
        fprintf(stderr, "File does not exist: %s\n", filename);  // Added newline
        perror("access");
        exit(1);
    }
    
    // Check if file is readable
    if (access(filename, R_OK) == -1) {
        fprintf(stderr, "File exists but not readable: %s\n", filename);  // Added newline
        perror("access");
        exit(1);
    }
    
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file %s\n", filename);  // Added newline
        perror("fopen");
        exit(1);
    }
    
    // Print first few bytes to check content
    char buffer[100];
    size_t bytes_read = fread(buffer, 1, sizeof(buffer) - 1, file);
    buffer[bytes_read] = '\0';
    printf("File content preview: %s\n", buffer);  // Added newline
    
    // Reset file position to beginning
    rewind(file);
    
    char line[MAX_LINE_LENGTH];
    int dimension = 0;
    int matrix_started = 0;
    igraph_t graph;
    igraph_vector_int_t edges;  // Changed to igraph_vector_int_t
    
    // Initialize problem_name as empty
    if (problem_name) {
        problem_name[0] = '\0';
    }
    
    // Read header information
    while (fgets(line, MAX_LINE_LENGTH, file)) {
        if (strstr(line, "NAME:")) {
            if (problem_name) {
                sscanf(line, "NAME: %255[^\n\r]", problem_name);
            }
        } else if (strstr(line, "DIMENSION:")) {
            sscanf(line, "DIMENSION: %d", &dimension);
            printf("Number of nodes: %d\n", dimension);
        } else if (strstr(line, "EDGE_WEIGHT_SECTION")) {
            matrix_started = 1;
            break;
        }
    }
    
    if (dimension <= 0) {
        fprintf(stderr, "Invalid dimension in TSP file\n");
        fclose(file);
        exit(1);
    }
    
    // Initialize empty graph with n vertices
    igraph_empty(&graph, dimension, IGRAPH_UNDIRECTED);
    
    // Initialize vectors for edges and weights
    igraph_vector_int_init(&edges, 0);  // Changed to igraph_vector_int_init
    igraph_vector_init(weights, 0);
    
    // Read adjacency matrix
    int row = 0;
    while (matrix_started && row < dimension && fgets(line, MAX_LINE_LENGTH, file)) {
        char* token = strtok(line, " \t\n");
        int col = 0;
        
        while (token != NULL && col < dimension) {
            long weight = atol(token);
            
            // Add edge if i < j (to avoid duplicates in undirected graph)
            // and if weight is not the "infinity" value (100000000 in your case)
            if (row < col && weight < MAX_WEIGHT) {
                igraph_vector_int_push_back(&edges, row);  // Changed to igraph_vector_int_push_back
                igraph_vector_int_push_back(&edges, col);  // Changed to igraph_vector_int_push_back
                igraph_vector_push_back(weights, weight);
            }
            
            token = strtok(NULL, " \t\n");
            col++;
        }
        
        row++;
    }
    
    fclose(file);
    
    // Add edges to graph
    igraph_add_edges(&graph, &edges, NULL);
    
    // Ensure weights vector has correct size
    if (igraph_vector_size(weights) != igraph_ecount(&graph)) {
        fprintf(stderr, "Error: Weight vector size doesn't match edge count\n");
        exit(1);
    }
    
    igraph_vector_int_destroy(&edges);  // Changed to igraph_vector_int_destroy
    
    // // Print graph edges and their weights
    // printf("Graph edges and weights:\n");
    // printf("Total edges: %ld\n", igraph_ecount(&graph));
    // for (int i = 0; i < igraph_ecount(&graph); i++) {
    //     igraph_integer_t from, to;
    //     igraph_edge(&graph, i, &from, &to);
    //     double weight = VECTOR(*weights)[i];
    //     printf("  Edge %d: %ld -- %ld, Weight: %.0f\n", i, (long)from, (long)to, weight);
    // }
    // printf("\n");

    return graph;
}

// Function to reconstruct path from predecessor matrix
void reconstruct_path(int source, int target, igraph_matrix_int_t* predecessors, 
                     igraph_vector_int_t* path) {
    igraph_vector_int_clear(path);
    
    if (MATRIX(*predecessors, source, target) == -1) {
        // No path exists
        return;
    }
    
    // Reconstruct path backwards from target to source
    igraph_stack_int_t stack;
    igraph_stack_int_init(&stack, 0);
    
    int current = target;
    while (current != -1 && current != source) {
        igraph_stack_int_push(&stack, current);
        current = MATRIX(*predecessors, source, current);
    }
    
    if (current == source) {
        igraph_stack_int_push(&stack, source);
        
        // Reverse the stack to get path from source to target
        while (!igraph_stack_int_empty(&stack)) {
            igraph_vector_int_push_back(path, igraph_stack_int_pop(&stack));
        }
    }
    
    igraph_stack_int_destroy(&stack);
}

// Fixed function to get the actual shortest paths between all pairs
void get_shortest_paths_matrix_optimized(igraph_t* graph, igraph_vector_t* weights, 
                                        igraph_matrix_t* distances, igraph_matrix_int_t* predecessors) {
    int n = igraph_vcount(graph);
    
    // // Initialize distance matrix
    // igraph_matrix_init(distances, n, n);
    
    // Calculate all-pairs shortest paths distances in one call
    igraph_error_t dist_result = igraph_distances_dijkstra(graph, distances, 
                            igraph_vss_all(), // from all vertices
                            igraph_vss_all(), // to all vertices  
                            weights, 
                            IGRAPH_ALL);
    
    if (dist_result != IGRAPH_SUCCESS) {
        fprintf(stderr, "Error calculating all-pairs distances\n");
        return;
    }
    
    // Initialize predecessor matrix with -1 (no predecessor)
    igraph_matrix_int_init(predecessors, n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            MATRIX(*predecessors, i, j) = -1;
        }
    }
    
    // Now calculate paths for each source to fill the predecessor matrix
    for (int source = 0; source < n; source++) {
        igraph_vector_int_list_t vertices;
        igraph_vector_int_list_t edges;
        igraph_vector_int_t parents;
        igraph_vector_int_t inbound_edges;
        
        // Initialize the vectors with proper capacity
        igraph_vector_int_list_init(&vertices, n);
        igraph_vector_int_list_init(&edges, n);
        igraph_vector_int_init(&parents, n);  // Initialize with capacity for n vertices
        igraph_vector_int_init(&inbound_edges, n);
        
        // Get shortest paths from source to all vertices
        igraph_error_t result = igraph_get_shortest_paths_dijkstra(
            graph,           // graph
            &vertices,       // vertices (paths)
            &edges,          // edges
            source,          // from
            igraph_vss_all(), // to all vertices
            weights,         // edge weights
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
        if (igraph_vector_int_size(&parents) == n) {
            for (int target = 0; target < n; target++) {
                MATRIX(*predecessors, source, target) = VECTOR(parents)[target];
            }
        } else {
            fprintf(stderr, "Warning: parents vector size %ld doesn't match n=%d for source %d\n", 
                   igraph_vector_int_size(&parents), n, source);
        }
        
        // Cleanup
        igraph_vector_int_list_destroy(&vertices);
        igraph_vector_int_list_destroy(&edges);
        igraph_vector_int_destroy(&parents);
        igraph_vector_int_destroy(&inbound_edges);
    }
}

// Alternative simpler version that only calculates distances (if paths are not critical)
void get_shortest_paths_distances_only(igraph_t* graph, igraph_vector_t* weights, 
                                      igraph_matrix_t* distances) {
    int n = igraph_vcount(graph);
    
    // Initialize distance matrix
    igraph_matrix_init(distances, n, n);
    
    // Calculate all-pairs shortest paths distances in one call
    igraph_error_t result = igraph_distances_dijkstra(graph, distances, 
                            igraph_vss_all(), // from all vertices
                            igraph_vss_all(), // to all vertices  
                            weights, 
                            IGRAPH_ALL);
    
    if (result != IGRAPH_SUCCESS) {
        fprintf(stderr, "Error calculating all-pairs distances\n");
    }
}

// Updated transform_tsp_file to handle cases where predecessors might not be available
void transform_tsp_file(const char* input_filename, const char* problem_name, 
                       igraph_t* graph, igraph_vector_t* weights, 
                       igraph_matrix_t* distances, igraph_matrix_int_t* predecessors,
                       const char* output_filename) {
    int n_vertices = igraph_vcount(graph);
    
    FILE* outfile = fopen(output_filename, "w");
    if (!outfile) {
        fprintf(stderr, "Error creating output file %s\n", output_filename);
        perror("fopen");
        return;
    }
    
    // Write header information
    fprintf(outfile, "NAME: walk_%s\n", problem_name);
    fprintf(outfile, "TYPE: TSP\n");
    fprintf(outfile, "COMMENT: Transformed TSP instance with %d nodes (using shortest path distances)\n", n_vertices);
    fprintf(outfile, "DIMENSION: %d\n", n_vertices);
    fprintf(outfile, "EDGE_WEIGHT_TYPE: EXPLICIT\n");
    fprintf(outfile, "EDGE_WEIGHT_FORMAT: FULL_MATRIX\n");
    fprintf(outfile, "EDGE_WEIGHT_SECTION\n");
    
    // Write the distance matrix
    for (int i = 0; i < n_vertices; i++) {
        for (int j = 0; j < n_vertices; j++) {
            double dist = MATRIX(*distances, i, j);
            if (dist > 1e9 || i == j) {
                fprintf(outfile, "%d ", MAX_WEIGHT);
            } else {
                fprintf(outfile, "%.0f ", dist);
            }
        }
        fprintf(outfile, "\n");
    }
    
    fprintf(outfile, "EOF\n");
    fclose(outfile);
    printf("Created walk file: %s\n", output_filename);
    
    // Only save paths if we have valid predecessor data
    if (predecessors != NULL) {
        // Also save the paths to a separate file
        char paths_filename[MAX_FILENAME_LENGTH];
        snprintf(paths_filename, sizeof(paths_filename), "%s.paths", output_filename);
        
        FILE* paths_file = fopen(paths_filename, "w");
        if (paths_file) {
            fprintf(paths_file, "SHORTEST_PATHS_MATRIX for %s\n", problem_name);
            fprintf(paths_file, "DIMENSION: %d\n", n_vertices);
            fprintf(paths_file, "PREDECESSOR_MATRIX_SECTION\n");
            
            // Write predecessor matrix
            for (int i = 0; i < n_vertices; i++) {
                for (int j = 0; j < n_vertices; j++) {
                    fprintf(paths_file, "%ld ", (long)MATRIX(*predecessors, i, j));
                }
                fprintf(paths_file, "\n");
            }
            
            fprintf(paths_file, "PATH_EXAMPLES_SECTION\n");
            
            // Write some example paths
            int path_count = 0;
            for (int i = 0; i < n_vertices && path_count < 20; i += n_vertices/10 + 1) {
                for (int j = 0; j < n_vertices && path_count < 20; j += n_vertices/10 + 1) {
                    if (i != j && MATRIX(*distances, i, j) < MAX_WEIGHT) {
                        igraph_vector_int_t path;
                        igraph_vector_int_init(&path, 0);
                        
                        reconstruct_path(i, j, predecessors, &path);
                        
                        if (igraph_vector_int_size(&path) > 0) {
                            fprintf(paths_file, "Path from %d to %d (distance: %.0f): ", 
                                   i, j, MATRIX(*distances, i, j));
                            for (int k = 0; k < igraph_vector_int_size(&path); k++) {
                                fprintf(paths_file, "%ld", (long)VECTOR(path)[k]);
                                if (k < igraph_vector_int_size(&path) - 1) {
                                    fprintf(paths_file, " -> ");
                                }
                            }
                            fprintf(paths_file, "\n");
                            path_count++;
                        }
                        
                        igraph_vector_int_destroy(&path);
                    }
                }
            }
            
            fprintf(paths_file, "EOF\n");
            fclose(paths_file);
            printf("Created paths file: %s\n", paths_filename);
        }
    }
}

// Updated create_transformed_tsp_instance function with better error handling
void create_transformed_tsp_instance(long int num_vertices, int instance_num, float target_fraction, igraph_matrix_t* distances) {
    // Create transformed TSP instance
    char tsp_filename[200];
    snprintf(tsp_filename, sizeof(tsp_filename), "instances/%ld_nodes/%d/vanila_tsp_instance_n_%ld_i_%d.tsp", 
        num_vertices, instance_num, num_vertices, instance_num);

    printf("\n=== Creating transformed instance for: %s ===\n", tsp_filename);
    
    igraph_t graph;
    igraph_vector_t weights;
    
    igraph_matrix_int_t predecessors; // to store path information
    char problem_name[256] = {0};
    
    // Read the vanilla TSP file
    graph = read_tsp_file(tsp_filename, &weights, problem_name);
    
    int n_vertices = igraph_vcount(&graph);
    printf("Graph created with %d vertices and %ld edges\n", 
           n_vertices, (long)igraph_ecount(&graph));
    
    // Initialize matrices
    igraph_matrix_int_init(&predecessors, 0, 0);
    
    // Calculate shortest path distances AND paths
    get_shortest_paths_matrix_optimized(&graph, &weights, distances, &predecessors);
    
    // Create output filename
    char output_filename[MAX_FILENAME_LENGTH];
    const char *base_name = strrchr(tsp_filename, '/');
    base_name = base_name ? base_name + 1 : tsp_filename;
    
    // Get directory path
    char dir_path[MAX_FILENAME_LENGTH] = "";
    const char *last_slash = strrchr(tsp_filename, '/');
    if (last_slash) {
        size_t dir_len = last_slash - tsp_filename + 1;
        strncpy(dir_path, tsp_filename, dir_len);
        dir_path[dir_len] = '\0';
    }
    
    // Generate transformed filename with safer snprintf
    if (n_vertices > 10000) {
        // For very large n, use a safer filename
        snprintf(output_filename, sizeof(output_filename), "%stransformed_%d_%d.tsp", 
                 dir_path, n_vertices, instance_num);
    } else {
        snprintf(output_filename, sizeof(output_filename), "%stransformed_tsp_instance_n_%d_i_%d.tsp", 
                 dir_path, n_vertices, instance_num);
    }
    
    // Transform and save the file
    transform_tsp_file(tsp_filename, problem_name, &graph, &weights, 
                      distances, &predecessors, output_filename);
    
    // Clean up
    igraph_matrix_int_destroy(&predecessors);
    // igraph_matrix_destroy(&distances);
    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);
    
    printf("Transformed TSP instance saved as: %s\n", output_filename);
}

// Function to get distance matrix from graph edge weights
void get_distance_matrix_from_graph(igraph_t* graph, igraph_vector_t* weights, igraph_matrix_t* distances) {
    int n = igraph_vcount(graph);
    igraph_matrix_init(distances, n, n);
    
    // Initialize with MAX_WEIGHT (no connection)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            MATRIX(*distances, i, j) = MAX_WEIGHT;
        }
    }
    
    // Set diagonal to 0
    for (int i = 0; i < n; i++) {
        MATRIX(*distances, i, i) = 0;
    }
    
    // Fill with actual edge weights
    for (int i = 0; i < igraph_ecount(graph); i++) {
        igraph_integer_t from, to;
        igraph_edge(graph, i, &from, &to);
        double weight = VECTOR(*weights)[i];
        MATRIX(*distances, from, to) = weight;
        MATRIX(*distances, to, from) = weight; // For undirected graph
    }
}


// Fixed function to find second nearest neighbor tour
int* find_second_nearest_neighbor_tour(igraph_matrix_t* distances, int n) {
    int* tour = (int*)malloc(n * sizeof(int));
    int* visited = (int*)calloc(n, sizeof(int));
    
    // Start at depot (node 0)
    tour[0] = 0;
    visited[0] = 1;
    int current_node = 0;
    
    // printf("Starting tour construction with %d nodes\n", n);
    
    for (int step = 1; step < n; step++) {
        // Find distances from current node to all unvisited nodes
        double min_dist = MAX_WEIGHT;
        double second_min_dist = MAX_WEIGHT;
        int min_node = -1;
        int second_min_node = -1;
        int unvisited_count = 0;
        
        // printf("Step %d: Current node = %d, looking for unvisited nodes...\n", step, current_node);
        
        for (int j = 0; j < n; j++) {
            if (!visited[j] && j != current_node) {
                unvisited_count++;
                double dist = MATRIX(*distances, current_node, j);
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
            for (int j = 0; j < n; j++) {
                if (!visited[j] && j != current_node) {
                    next_node = j;
                    printf("  -> Fallback to node %d\n", next_node);
                    break;
                }
            }
        }
        
        // Safety check
        if (next_node < 0 || next_node >= n || visited[next_node]) {
            printf("  ERROR: Invalid next node %d!\n", next_node);
            // Emergency fallback - find first unvisited node
            for (int j = 0; j < n; j++) {
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

// // Function to calculate visit times along a tour
// double* calculate_visit_times(int* tour, igraph_matrix_t* distances, int n) {
//     // TODO: its wrong
//     double* visit_times = (double*)malloc(n * sizeof(double));
    
//     // Start at depot at time 0
//     visit_times[0] = 0.0;
    
//     // Calculate arrival times at each node (no waiting in base calculation)
//     for (int i = 1; i < n; i++) {
//         int prev_node = tour[i-1];
//         int curr_node = tour[i];
//         double travel_time = MATRIX(*distances, prev_node, curr_node);
//         visit_times[i] = visit_times[i-1] + travel_time;
//     }
    
//     return visit_times;
// }

// Function to calculate visit times along a tour
double* calculate_visit_times(int* tour, igraph_matrix_t* distances, int n) {
    double* visit_times = (double*)malloc(n * sizeof(double));
    
    // Initialize all visit times to -1 (unvisited)
    for (int i = 0; i < n; i++) {
        visit_times[i] = -1.0;
    }
    
    // Start at depot (first node in tour) at time 0
    visit_times[tour[0]] = 0.0;
    double current_time = 0.0;
    
    // Calculate arrival times at each node in the tour order
    for (int i = 1; i < n; i++) {
        int prev_node = tour[i-1];
        int curr_node = tour[i];
        double travel_time = MATRIX(*distances, prev_node, curr_node);
        current_time += travel_time;
        visit_times[curr_node] = current_time;
    }
    
    return visit_times;
}

// Function to generate time windows based on visit times
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

// Function to create TSPTW instance file
void create_tsptw_instance_file(igraph_matrix_t* distances, long int num_vertices, int instance_num, 
                               double* a_i, double* b_i, double T,
                               int width,
                               const char* instance_type) {
    char filename[256];
    char dir_path[100];
    
    // Create the same directory structure as TSP files
    snprintf(dir_path, sizeof(dir_path), "instances/%ld_nodes/%d", num_vertices, instance_num);
    
    // Create .tsptw filename
    snprintf(filename, sizeof(filename), "%s/%s_tsptw_instance_n_%ld_w%d_i_%d.tsptw", 
             dir_path, instance_type, num_vertices, width, instance_num);
    
    // Create directories if they don't exist
    char mkdir_command[200];
    snprintf(mkdir_command, sizeof(mkdir_command), "mkdir -p %s", dir_path);
    system(mkdir_command);
    
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
                fprintf(file, "%d", MAX_WEIGHT);
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
void create_tsptw_dataset_from_tsp(long int num_vertices, int instance_num, int* width_list, int num_widths, igraph_matrix_t* vanila_distances, igraph_matrix_t* shortest_distances) {
    // Create transformed TSP instance
    char dir_path[100];
    snprintf(dir_path, sizeof(dir_path), "instances/%ld_nodes/%d", 
        num_vertices, instance_num);

    printf("\n=== Creating TSPTW dataset from: %s ===\n", dir_path);
    
    int n = (int)num_vertices;
    // Find second nearest neighbor tour
    int* tour = find_second_nearest_neighbor_tour(shortest_distances, n);

    // Print the tour for verification
    printf("Second nearest neighbor tour: ");
    for (int i = 0; i < num_vertices; i++) {
        printf("%d ", tour[i]);
    }
    printf("\n");
    
    // Calculate visit times along the tour
    double* visit_times = calculate_visit_times(tour, shortest_distances, n);
    
    // Print visit times for verification
    printf("#######Visit times: \n");
    for (int i = 0; i < n; i++) {
        printf("%d: %.1f\n", i, visit_times[i]);
    }
    printf("\n");
    
    // Generate instances for different time window widths
    for (int w = 0; w < num_widths; w++) {
        int width = width_list[w];
        
        double* a_i, * b_i;
        double T;
        
        // Generate time windows
        generate_time_windows(visit_times, shortest_distances, n, width, &a_i, &b_i, &T);

        // Print times window for verification
        printf("Time Window: \n");
        for (int i = 0; i < n; i++) {
            printf("%d: %.1f, %.1f\n", i, a_i[i], b_i[i]);
        }
        printf("\n");
        
        create_tsptw_instance_file(vanila_distances, num_vertices, instance_num, 
                                 a_i, b_i, T, width, "vanilla");

        create_tsptw_instance_file(shortest_distances, num_vertices, instance_num, 
                                 a_i, b_i, T, width, "transformed");

        free(a_i);
        free(b_i);
    }
    
    // Clean up
    free(tour);
    free(visit_times);
    
    
    // printf("TSPTW dataset creation completed for %s instance %d\n", instance_type, instance_num);
}
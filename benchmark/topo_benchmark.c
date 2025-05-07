#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <time.h>
#include <getopt.h>

// Atom types commonly found in SMILES
typedef enum {
    ATOM_C,
    ATOM_N,
    ATOM_O,
    ATOM_S,
    ATOM_P,
    ATOM_F,
    ATOM_CL,
    ATOM_BR,
    ATOM_I,
    ATOM_OTHER,
    ATOM_TYPES_COUNT
} AtomType;

// Atom representation in a molecular graph
typedef struct {
    AtomType type;
    int index;
    int aromatic;
    int charge;
    int bonds_count;
    int* bonded_atoms;
    int* bond_types;
} Atom;

// Bond types
typedef enum {
    BOND_SINGLE = 1,
    BOND_DOUBLE = 2,
    BOND_TRIPLE = 3,
    BOND_AROMATIC = 4
} BondType;

// Molecular graph representation
typedef struct {
    int atom_count;
    Atom* atoms;
    int rings_count;
    int* ring_sizes;
    int max_path_length;
} MolecularGraph;

// Helper functions for SMILES parsing
AtomType determine_atom_type(const char* smiles, int pos) {
    if (smiles[pos] == 'C') {
        if (smiles[pos+1] == 'l') return ATOM_CL;
        return ATOM_C;
    } else if (smiles[pos] == 'N') {
        return ATOM_N;
    } else if (smiles[pos] == 'O') {
        return ATOM_O;
    } else if (smiles[pos] == 'S') {
        return ATOM_S;
    } else if (smiles[pos] == 'P') {
        return ATOM_P;
    } else if (smiles[pos] == 'F') {
        return ATOM_F;
    } else if (smiles[pos] == 'B' && smiles[pos+1] == 'r') {
        return ATOM_BR;
    } else if (smiles[pos] == 'I') {
        return ATOM_I;
    }
    return ATOM_OTHER;
}

bool is_aromatic(char c) {
    return (c == 'c' || c == 'n' || c == 'o' || c == 's' || c == 'p');
}

// Parse a SMILES string into a molecular graph
MolecularGraph* parse_smiles(const char* smiles) {
    int len = strlen(smiles);
    MolecularGraph* graph = (MolecularGraph*)malloc(sizeof(MolecularGraph));
    
    // First pass: count atoms
    int atom_count = 0;
    for (int i = 0; i < len; i++) {
        if (isalpha(smiles[i]) || smiles[i] == '[') {
            atom_count++;
            // Skip multi-character atoms
            if (smiles[i] == '[') {
                while (i < len && smiles[i] != ']') i++;
            } else if ((smiles[i] == 'C' && smiles[i+1] == 'l') || 
                      (smiles[i] == 'B' && smiles[i+1] == 'r')) {
                i++;
            }
        }
    }
    
    // Allocate atoms
    graph->atom_count = atom_count;
    graph->atoms = (Atom*)malloc(atom_count * sizeof(Atom));
    memset(graph->atoms, 0, atom_count * sizeof(Atom));
    
    // Initialize rings count
    graph->rings_count = 0;
    graph->ring_sizes = NULL;
    graph->max_path_length = 0;
    
    // Track brackets for ring closures
    int* ring_openings = (int*)malloc(10 * sizeof(int));
    memset(ring_openings, -1, 10 * sizeof(int));
    
    // Second pass: create atoms and basic connectivity
    int atom_idx = 0;
    int prev_atom = -1;
    int branch_points[100];
    int branch_count = 0;
    
    for (int i = 0; i < len; i++) {
        if (isalpha(smiles[i]) || smiles[i] == '[') {
            // Process atom
            graph->atoms[atom_idx].index = atom_idx;
            
            if (smiles[i] == '[') {
                // Complex atom notation - skip for simplicity
                graph->atoms[atom_idx].type = ATOM_OTHER;
                while (i < len && smiles[i] != ']') i++;
            } else if (islower(smiles[i])) {
                // Aromatic atom
                graph->atoms[atom_idx].aromatic = 1;
                if (smiles[i] == 'c') graph->atoms[atom_idx].type = ATOM_C;
                else if (smiles[i] == 'n') graph->atoms[atom_idx].type = ATOM_N;
                else if (smiles[i] == 'o') graph->atoms[atom_idx].type = ATOM_O;
                else if (smiles[i] == 's') graph->atoms[atom_idx].type = ATOM_S;
                else if (smiles[i] == 'p') graph->atoms[atom_idx].type = ATOM_P;
                else graph->atoms[atom_idx].type = ATOM_OTHER;
            } else {
                // Regular atom
                graph->atoms[atom_idx].type = determine_atom_type(smiles, i);
                if ((smiles[i] == 'C' && smiles[i+1] == 'l') || 
                    (smiles[i] == 'B' && smiles[i+1] == 'r')) {
                    i++;
                }
            }
            
            // Connect to previous atom if applicable
            if (prev_atom != -1) {
                // Allocate bond arrays if needed
                if (graph->atoms[prev_atom].bonds_count == 0) {
                    graph->atoms[prev_atom].bonded_atoms = (int*)malloc(4 * sizeof(int));
                    graph->atoms[prev_atom].bond_types = (int*)malloc(4 * sizeof(int));
                }
                if (graph->atoms[atom_idx].bonds_count == 0) {
                    graph->atoms[atom_idx].bonded_atoms = (int*)malloc(4 * sizeof(int));
                    graph->atoms[atom_idx].bond_types = (int*)malloc(4 * sizeof(int));
                }
                
                // Add bond
                int bond_type = BOND_SINGLE;
                if (graph->atoms[prev_atom].aromatic && graph->atoms[atom_idx].aromatic) {
                    bond_type = BOND_AROMATIC;
                }
                
                graph->atoms[prev_atom].bonded_atoms[graph->atoms[prev_atom].bonds_count] = atom_idx;
                graph->atoms[prev_atom].bond_types[graph->atoms[prev_atom].bonds_count] = bond_type;
                graph->atoms[prev_atom].bonds_count++;
                
                graph->atoms[atom_idx].bonded_atoms[graph->atoms[atom_idx].bonds_count] = prev_atom;
                graph->atoms[atom_idx].bond_types[graph->atoms[atom_idx].bonds_count] = bond_type;
                graph->atoms[atom_idx].bonds_count++;
            }
            
            prev_atom = atom_idx;
            atom_idx++;
            
            // Check for ring closures
            if (i+1 < len && isdigit(smiles[i+1])) {
                int ring_num = smiles[i+1] - '0';
                i++;
                
                if (ring_openings[ring_num] == -1) {
                    // Start of ring
                    ring_openings[ring_num] = prev_atom;
                } else {
                    // End of ring, connect atoms
                    int start_atom = ring_openings[ring_num];
                    int end_atom = prev_atom;
                    
                    // Allocate bond arrays if needed
                    if (graph->atoms[start_atom].bonds_count == 0) {
                        graph->atoms[start_atom].bonded_atoms = (int*)malloc(4 * sizeof(int));
                        graph->atoms[start_atom].bond_types = (int*)malloc(4 * sizeof(int));
                    }
                    if (graph->atoms[end_atom].bonds_count == 0) {
                        graph->atoms[end_atom].bonded_atoms = (int*)malloc(4 * sizeof(int));
                        graph->atoms[end_atom].bond_types = (int*)malloc(4 * sizeof(int));
                    }
                    
                    // Add bond
                    int bond_type = BOND_SINGLE;
                    if (graph->atoms[start_atom].aromatic && graph->atoms[end_atom].aromatic) {
                        bond_type = BOND_AROMATIC;
                    }
                    
                    graph->atoms[start_atom].bonded_atoms[graph->atoms[start_atom].bonds_count] = end_atom;
                    graph->atoms[start_atom].bond_types[graph->atoms[start_atom].bonds_count] = bond_type;
                    graph->atoms[start_atom].bonds_count++;
                    
                    graph->atoms[end_atom].bonded_atoms[graph->atoms[end_atom].bonds_count] = start_atom;
                    graph->atoms[end_atom].bond_types[graph->atoms[end_atom].bonds_count] = bond_type;
                    graph->atoms[end_atom].bonds_count++;
                    
                    // Increment rings count
                    graph->rings_count++;
                    
                    // Reset ring opening
                    ring_openings[ring_num] = -1;
                }
            }
        } else if (smiles[i] == '(') {
            // Branch start
            branch_points[branch_count++] = prev_atom;
        } else if (smiles[i] == ')') {
            // Branch end
            if (branch_count > 0) {
                prev_atom = branch_points[--branch_count];
            }
        } else if (smiles[i] == '=') {
            // Double bond - update previous bond
            if (prev_atom >= 0 && atom_idx > 0) {
                int prev_bond_idx = graph->atoms[prev_atom].bonds_count - 1;
                if (prev_bond_idx >= 0) {
                    graph->atoms[prev_atom].bond_types[prev_bond_idx] = BOND_DOUBLE;
                    
                    int other_atom = graph->atoms[prev_atom].bonded_atoms[prev_bond_idx];
                    for (int j = 0; j < graph->atoms[other_atom].bonds_count; j++) {
                        if (graph->atoms[other_atom].bonded_atoms[j] == prev_atom) {
                            graph->atoms[other_atom].bond_types[j] = BOND_DOUBLE;
                            break;
                        }
                    }
                }
            }
        } else if (smiles[i] == '#') {
            // Triple bond - update previous bond
            if (prev_atom >= 0 && atom_idx > 0) {
                int prev_bond_idx = graph->atoms[prev_atom].bonds_count - 1;
                if (prev_bond_idx >= 0) {
                    graph->atoms[prev_atom].bond_types[prev_bond_idx] = BOND_TRIPLE;
                    
                    int other_atom = graph->atoms[prev_atom].bonded_atoms[prev_bond_idx];
                    for (int j = 0; j < graph->atoms[other_atom].bonds_count; j++) {
                        if (graph->atoms[other_atom].bonded_atoms[j] == prev_atom) {
                            graph->atoms[other_atom].bond_types[j] = BOND_TRIPLE;
                            break;
                        }
                    }
                }
            }
        }
    }
    
    free(ring_openings);
    return graph;
}

// Free the molecular graph
void free_molecular_graph(MolecularGraph* graph) {
    if (!graph) return;
    
    for (int i = 0; i < graph->atom_count; i++) {
        free(graph->atoms[i].bonded_atoms);
        free(graph->atoms[i].bond_types);
    }
    
    free(graph->atoms);
    free(graph->ring_sizes);
    free(graph);
}

// Floyd-Warshall algorithm to compute all-pairs shortest paths
int** compute_distance_matrix(MolecularGraph* graph) {
    int n = graph->atom_count;
    int** dist = (int**)malloc(n * sizeof(int*));
    
    for (int i = 0; i < n; i++) {
        dist[i] = (int*)malloc(n * sizeof(int));
        for (int j = 0; j < n; j++) {
            dist[i][j] = (i == j) ? 0 : 9999; // "Infinity"
        }
    }
    
    // Initialize with direct connections
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < graph->atoms[i].bonds_count; j++) {
            int neighbor = graph->atoms[i].bonded_atoms[j];
            dist[i][neighbor] = 1;
        }
    }
    
    // Floyd-Warshall algorithm
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (dist[i][k] + dist[k][j] < dist[i][j]) {
                    dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
    }
    
    // Find maximum path length
    int max_path = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (dist[i][j] != 9999 && dist[i][j] > max_path) {
                max_path = dist[i][j];
            }
        }
    }
    graph->max_path_length = max_path;
    
    return dist;
}

// Free the distance matrix
void free_distance_matrix(int** dist_matrix, int atom_count) {
    for (int i = 0; i < atom_count; i++) {
        free(dist_matrix[i]);
    }
    free(dist_matrix);
}

// Molecular Mesh Descriptors
double calculate_molecular_mesh(MolecularGraph* graph) {
    if (graph->atom_count == 0) return 0.0;
    
    int** dist_matrix = compute_distance_matrix(graph);
    
    // Calculate mesh complexity based on connectivity patterns
    double mesh_value = 0.0;
    
    // Sum of all path lengths
    double total_paths = 0.0;
    int path_count = 0;
    
    for (int i = 0; i < graph->atom_count; i++) {
        for (int j = i+1; j < graph->atom_count; j++) {
            if (dist_matrix[i][j] < 9999) {
                total_paths += dist_matrix[i][j];
                path_count++;
            }
        }
    }
    
    // Average path length
    double avg_path = (path_count > 0) ? total_paths / path_count : 0;
    
    // Mesh complexity is a function of atom count, ring count, and path distribution
    mesh_value = (graph->atom_count * 0.5) + (graph->rings_count * 1.5) + (avg_path * 0.3);
    
    // Normalize by dividing by number of atoms
    mesh_value /= graph->atom_count;
    
    free_distance_matrix(dist_matrix, graph->atom_count);
    return mesh_value;
}

// Eccentric Distance Index
double calculate_eccentric_distance(MolecularGraph* graph) {
    if (graph->atom_count <= 1) return 0.0;
    
    int** dist_matrix = compute_distance_matrix(graph);
    
    // Calculate eccentricity for each atom (maximum distance to any other atom)
    int* eccentricity = (int*)malloc(graph->atom_count * sizeof(int));
    
    for (int i = 0; i < graph->atom_count; i++) {
        eccentricity[i] = 0;
        for (int j = 0; j < graph->atom_count; j++) {
            if (dist_matrix[i][j] < 9999 && dist_matrix[i][j] > eccentricity[i]) {
                eccentricity[i] = dist_matrix[i][j];
            }
        }
    }
    
    // Calculate eccentric distance sum
    double ecc_dist_sum = 0.0;
    
    for (int i = 0; i < graph->atom_count; i++) {
        ecc_dist_sum += eccentricity[i];
    }
    
    free(eccentricity);
    free_distance_matrix(dist_matrix, graph->atom_count);
    
    return ecc_dist_sum;
}

// Vertex-Edge Path Descriptors
double calculate_vertex_edge_path(MolecularGraph* graph) {
    if (graph->atom_count <= 1) return 0.0;
    
    int** dist_matrix = compute_distance_matrix(graph);
    
    // Calculate vertex-edge path index
    double vep_index = 0.0;
    int total_bonds = 0;
    
    // Count total bonds
    for (int i = 0; i < graph->atom_count; i++) {
        total_bonds += graph->atoms[i].bonds_count;
    }
    total_bonds /= 2;  // Each bond is counted twice
    
    // Calculate vertex contribution
    double vertex_contrib = graph->atom_count * 0.5;
    
    // Calculate edge contribution
    double edge_contrib = total_bonds * 0.7;
    
    // Calculate path contribution
    double path_contrib = 0.0;
    for (int i = 0; i < graph->atom_count; i++) {
        for (int j = i+1; j < graph->atom_count; j++) {
            if (dist_matrix[i][j] < 9999) {
                path_contrib += 1.0 / dist_matrix[i][j];
            }
        }
    }
    
    // Combine all contributions
    vep_index = vertex_contrib + edge_contrib + (path_contrib * 0.3);
    
    free_distance_matrix(dist_matrix, graph->atom_count);
    return vep_index;
}

// Mock implementation of Context and GetSmilesFunc for standalone use
typedef struct {
    const char* smiles;
} MockContext;

const char* getMockSmiles(const MockContext* ctx) {
    return ctx->smiles;
}

void print_usage() {
    printf("Usage: topo_benchmark [options]\n");
    printf("Options:\n");
    printf("  -i, --input <file>   Input file with SMILES strings (one per line)\n");
    printf("  -o, --output <file>  Output file for results\n");
    printf("  -s, --smiles <str>   Direct SMILES string input\n");
    printf("  -h, --help           Display this help message\n");
}

int main(int argc, char *argv[]) {
    char* input_file = NULL;
    char* output_file = NULL;
    char* smiles_str = NULL;
    
    // Parse command line arguments
    static struct option long_options[] = {
        {"input", required_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
        {"smiles", required_argument, 0, 's'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    
    int option_index = 0;
    int c;
    
    while ((c = getopt_long(argc, argv, "i:o:s:h", long_options, &option_index)) != -1) {
        switch (c) {
            case 'i':
                input_file = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 's':
                smiles_str = optarg;
                break;
            case 'h':
                print_usage();
                return 0;
            case '?':
                print_usage();
                return 1;
            default:
                abort();
        }
    }
    
    // Check if we have either input file or direct SMILES
    if (!input_file && !smiles_str) {
        fprintf(stderr, "Error: Either input file (-i) or SMILES string (-s) is required\n");
        print_usage();
        return 1;
    }
    
    // Open output file if specified, otherwise use stdout
    FILE* out_fp = stdout;
    if (output_file) {
        out_fp = fopen(output_file, "w");
        if (!out_fp) {
            fprintf(stderr, "Error: Could not open output file %s\n", output_file);
            return 1;
        }
    }
    
    // Write header
    fprintf(out_fp, "SMILES,MolecularMesh,EccentricDistance,VertexEdgePath,Time(ms)\n");
    
    // Process SMILES from direct input
    if (smiles_str) {
        clock_t start = clock();
        
        // Create mock context
        MockContext ctx = {smiles_str};
        
        // Parse SMILES
        MolecularGraph* graph = parse_smiles(smiles_str);
        
        // Calculate descriptors
        double mesh = calculate_molecular_mesh(graph);
        double ecc_dist = calculate_eccentric_distance(graph);
        double vep = calculate_vertex_edge_path(graph);
        
        // Calculate time
        clock_t end = clock();
        double time_ms = (double)(end - start) * 1000.0 / CLOCKS_PER_SEC;
        
        // Output results
        fprintf(out_fp, "%s,%.6f,%.6f,%.6f,%.2f\n", 
                smiles_str, mesh, ecc_dist, vep, time_ms);
        
        // Free resources
        free_molecular_graph(graph);
    }
    
    // Process SMILES from input file
    if (input_file) {
        FILE* in_fp = fopen(input_file, "r");
        if (!in_fp) {
            fprintf(stderr, "Error: Could not open input file %s\n", input_file);
            if (output_file) fclose(out_fp);
            return 1;
        }
        
        char line[1024];
        while (fgets(line, sizeof(line), in_fp)) {
            // Remove newline
            line[strcspn(line, "\n")] = 0;
            
            // Skip empty lines
            if (strlen(line) == 0) continue;
            
            clock_t start = clock();
            
            // Create mock context
            MockContext ctx = {line};
            
            // Parse SMILES
            MolecularGraph* graph = parse_smiles(line);
            
            // Calculate descriptors
            double mesh = calculate_molecular_mesh(graph);
            double ecc_dist = calculate_eccentric_distance(graph);
            double vep = calculate_vertex_edge_path(graph);
            
            // Calculate time
            clock_t end = clock();
            double time_ms = (double)(end - start) * 1000.0 / CLOCKS_PER_SEC;
            
            // Output results
            fprintf(out_fp, "%s,%.6f,%.6f,%.6f,%.2f\n", 
                    line, mesh, ecc_dist, vep, time_ms);
            
            // Free resources
            free_molecular_graph(graph);
        }
        
        fclose(in_fp);
    }
    
    // Close output file if opened
    if (output_file) fclose(out_fp);
    
    return 0;
}

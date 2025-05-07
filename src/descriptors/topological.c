#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <limits.h> // For INT_MAX
#include "../cregistry.h"

// Forward declaration for the Context interface to match C++ code

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
                if (graph->atoms[atom_idx].type == ATOM_CL || graph->atoms[atom_idx].type == ATOM_BR) {
                    i++; // Skip the second letter
                }
            }
            
            // Connect to previous atom if exists
            if (prev_atom >= 0) {
                // Allocate space for the bond
                if (graph->atoms[prev_atom].bonds_count == 0) {
                    graph->atoms[prev_atom].bonded_atoms = (int*)malloc(sizeof(int));
                    graph->atoms[prev_atom].bond_types = (int*)malloc(sizeof(int));
                } else {
                    graph->atoms[prev_atom].bonded_atoms = (int*)realloc(
                        graph->atoms[prev_atom].bonded_atoms, 
                        (graph->atoms[prev_atom].bonds_count + 1) * sizeof(int));
                    graph->atoms[prev_atom].bond_types = (int*)realloc(
                        graph->atoms[prev_atom].bond_types,
                        (graph->atoms[prev_atom].bonds_count + 1) * sizeof(int));
                }
                
                // Set the bond type based on context
                int bond_type = BOND_SINGLE;
                if (i > 0) {
                    if (smiles[i-1] == '=') bond_type = BOND_DOUBLE;
                    else if (smiles[i-1] == '#') bond_type = BOND_TRIPLE;
                    else if (graph->atoms[prev_atom].aromatic && graph->atoms[atom_idx].aromatic) 
                        bond_type = BOND_AROMATIC;
                }
                
                // Add the bond
                graph->atoms[prev_atom].bonded_atoms[graph->atoms[prev_atom].bonds_count] = atom_idx;
                graph->atoms[prev_atom].bond_types[graph->atoms[prev_atom].bonds_count] = bond_type;
                graph->atoms[prev_atom].bonds_count++;
                
                // Add the reverse bond
                if (graph->atoms[atom_idx].bonds_count == 0) {
                    graph->atoms[atom_idx].bonded_atoms = (int*)malloc(sizeof(int));
                    graph->atoms[atom_idx].bond_types = (int*)malloc(sizeof(int));
                } else {
                    graph->atoms[atom_idx].bonded_atoms = (int*)realloc(
                        graph->atoms[atom_idx].bonded_atoms, 
                        (graph->atoms[atom_idx].bonds_count + 1) * sizeof(int));
                    graph->atoms[atom_idx].bond_types = (int*)realloc(
                        graph->atoms[atom_idx].bond_types,
                        (graph->atoms[atom_idx].bonds_count + 1) * sizeof(int));
                }
                
                graph->atoms[atom_idx].bonded_atoms[graph->atoms[atom_idx].bonds_count] = prev_atom;
                graph->atoms[atom_idx].bond_types[graph->atoms[atom_idx].bonds_count] = bond_type;
                graph->atoms[atom_idx].bonds_count++;
            }
            
            prev_atom = atom_idx;
            atom_idx++;
        } else if (smiles[i] == '(') {
            // Start a branch - push the current atom onto the branch stack
            branch_points[branch_count++] = prev_atom;
        } else if (smiles[i] == ')') {
            // End a branch - pop the last branch point
            if (branch_count > 0) {
                prev_atom = branch_points[--branch_count];
            }
        } else if (isdigit(smiles[i])) {
            // Ring closure
            int ring_number = smiles[i] - '0';
            if (ring_openings[ring_number] == -1) {
                // Opening a ring
                ring_openings[ring_number] = prev_atom;
            } else {
                // Closing a ring - add bonds between the atoms
                int open_atom = ring_openings[ring_number];
                
                // Add bond from open to current
                if (graph->atoms[open_atom].bonds_count == 0) {
                    graph->atoms[open_atom].bonded_atoms = (int*)malloc(sizeof(int));
                    graph->atoms[open_atom].bond_types = (int*)malloc(sizeof(int));
                } else {
                    graph->atoms[open_atom].bonded_atoms = (int*)realloc(
                        graph->atoms[open_atom].bonded_atoms, 
                        (graph->atoms[open_atom].bonds_count + 1) * sizeof(int));
                    graph->atoms[open_atom].bond_types = (int*)realloc(
                        graph->atoms[open_atom].bond_types,
                        (graph->atoms[open_atom].bonds_count + 1) * sizeof(int));
                }
                
                // Determine bond type for ring closure
                int bond_type = BOND_SINGLE;
                if (i > 0 && smiles[i-1] == '=') bond_type = BOND_DOUBLE;
                else if (i > 0 && smiles[i-1] == '#') bond_type = BOND_TRIPLE;
                else if (graph->atoms[open_atom].aromatic && graph->atoms[prev_atom].aromatic) 
                    bond_type = BOND_AROMATIC;
                
                graph->atoms[open_atom].bonded_atoms[graph->atoms[open_atom].bonds_count] = prev_atom;
                graph->atoms[open_atom].bond_types[graph->atoms[open_atom].bonds_count] = bond_type;
                graph->atoms[open_atom].bonds_count++;
                
                // Add bond from current to open
                if (graph->atoms[prev_atom].bonds_count == 0) {
                    graph->atoms[prev_atom].bonded_atoms = (int*)malloc(sizeof(int));
                    graph->atoms[prev_atom].bond_types = (int*)malloc(sizeof(int));
                } else {
                    graph->atoms[prev_atom].bonded_atoms = (int*)realloc(
                        graph->atoms[prev_atom].bonded_atoms, 
                        (graph->atoms[prev_atom].bonds_count + 1) * sizeof(int));
                    graph->atoms[prev_atom].bond_types = (int*)realloc(
                        graph->atoms[prev_atom].bond_types,
                        (graph->atoms[prev_atom].bonds_count + 1) * sizeof(int));
                }
                
                graph->atoms[prev_atom].bonded_atoms[graph->atoms[prev_atom].bonds_count] = open_atom;
                graph->atoms[prev_atom].bond_types[graph->atoms[prev_atom].bonds_count] = bond_type;
                graph->atoms[prev_atom].bonds_count++;
                
                // Reset ring opening
                ring_openings[ring_number] = -1;
            }
        }
    }
    
    // Count and identify cycles using DFS
    graph->rings_count = 0;
    graph->ring_sizes = NULL;
    graph->max_path_length = 0;
    
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

// Helper for BFS to compute distances from a source_atom_idx
void bfs_from_source(MolecularGraph* graph, int source_atom_idx, int* dist_row) {
    int n = graph->atom_count;
    for (int i = 0; i < n; ++i) {
        dist_row[i] = INT_MAX; // Use INT_MAX for "infinity"
    }
    dist_row[source_atom_idx] = 0;

    // Simple array-based queue
    int* queue = (int*)malloc(n * sizeof(int));
    if (!queue) { // Allocation check
        perror("Failed to allocate queue for BFS");
        return; // Or handle error appropriately
    }
    int head = 0, tail = 0;
    queue[tail++] = source_atom_idx;

    while (head < tail) {
        int u = queue[head++];
        Atom* atom_u = &graph->atoms[u];
        for (int i = 0; i < atom_u->bonds_count; ++i) {
            int v = atom_u->bonded_atoms[i];
            if (dist_row[v] == INT_MAX) { // If not visited
                dist_row[v] = dist_row[u] + 1;
                if (tail < n) { // Basic bounds check for queue
                    queue[tail++] = v;
                } else {
                    // This case should ideally not happen if n is graph->atom_count
                    // and graph is connected. Handle error or resize if necessary.
                }
            }
        }
    }
    free(queue);
}

// Compute all-pairs shortest paths using N BFS runs
int** compute_distance_matrix(MolecularGraph* graph) { // Renamed BFS version to be the default
    int n = graph->atom_count;
    if (n == 0) return NULL;

    int** dist = (int**)malloc(n * sizeof(int*));
    if (!dist) {
        perror("Failed to allocate distance matrix rows");
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        dist[i] = (int*)malloc(n * sizeof(int));
        if (!dist[i]) {
            perror("Failed to allocate distance matrix columns");
            // Free previously allocated rows
            for (int k = 0; k < i; ++k) free(dist[k]);
            free(dist);
            return NULL;
        }
        bfs_from_source(graph, i, dist[i]);
    }

    // Update graph->max_path_length and handle "infinity" (INT_MAX)
    int max_path = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (dist[i][j] == INT_MAX) {
                dist[i][j] = -1; // Consistent "no path" representation
            } else if (dist[i][j] > max_path) {
                max_path = dist[i][j];
            }
        }
    }
    graph->max_path_length = max_path;
    
    return dist;
}

// Compute atom eccentricity (maximum distance from an atom to any other atom)
int* compute_eccentricity(int** dist_matrix, int atom_count) {
    int* eccentricity = (int*)malloc(atom_count * sizeof(int));
    
    for (int i = 0; i < atom_count; i++) {
        eccentricity[i] = 0;
        for (int j = 0; j < atom_count; j++) {
            if (dist_matrix[i][j] != 9999 && dist_matrix[i][j] > eccentricity[i]) {
                eccentricity[i] = dist_matrix[i][j];
            }
        }
    }
    
    return eccentricity;
}

// Free the distance matrix
void free_distance_matrix(int** dist_matrix, int atom_count) {
    for (int i = 0; i < atom_count; i++) {
        free(dist_matrix[i]);
    }
    free(dist_matrix);
}

// Extended Topological Atom Pairs (ETAP)
double calculate_etap(MolecularGraph* graph) {
    if (graph->atom_count < 2) return 0.0;
    
    int** dist_matrix = compute_distance_matrix(graph);
    
    // Calculate ETAP score using atom type and distance
    double etap_sum = 0.0;
    
    for (int i = 0; i < graph->atom_count; i++) {
        for (int j = i+1; j < graph->atom_count; j++) {
            if (dist_matrix[i][j] != 9999) {
                double atom_weight_i = 0.0, atom_weight_j = 0.0;
                
                // Weight by atom type
                switch (graph->atoms[i].type) {
                    case ATOM_C: atom_weight_i = 1.0; break;
                    case ATOM_N: atom_weight_i = 2.0; break;
                    case ATOM_O: atom_weight_i = 2.0; break;
                    case ATOM_S: atom_weight_i = 3.0; break;
                    case ATOM_P: atom_weight_i = 3.0; break;
                    case ATOM_F: 
                    case ATOM_CL: 
                    case ATOM_BR: 
                    case ATOM_I: atom_weight_i = 4.0; break;
                    default: atom_weight_i = 1.5; break;
                }
                
                switch (graph->atoms[j].type) {
                    case ATOM_C: atom_weight_j = 1.0; break;
                    case ATOM_N: atom_weight_j = 2.0; break;
                    case ATOM_O: atom_weight_j = 2.0; break;
                    case ATOM_S: atom_weight_j = 3.0; break;
                    case ATOM_P: atom_weight_j = 3.0; break;
                    case ATOM_F: 
                    case ATOM_CL: 
                    case ATOM_BR: 
                    case ATOM_I: atom_weight_j = 4.0; break;
                    default: atom_weight_j = 1.5; break;
                }
                
                // Additional weight for aromaticity
                if (graph->atoms[i].aromatic) atom_weight_i *= 1.5;
                if (graph->atoms[j].aromatic) atom_weight_j *= 1.5;
                
                // Distance factor: closer atoms contribute more
                double distance_factor = 1.0 / (double)dist_matrix[i][j];
                
                etap_sum += (atom_weight_i * atom_weight_j * distance_factor);
            }
        }
    }
    
    free_distance_matrix(dist_matrix, graph->atom_count);
    return etap_sum / graph->atom_count;
}

// Molecular Mesh Descriptors
double calculate_molecular_mesh(MolecularGraph* graph) {
    if (graph->atom_count < 3) return 0.0;
    
    int** dist_matrix = compute_distance_matrix(graph);
    
    // Calculate mesh density
    int edge_count = 0;
    for (int i = 0; i < graph->atom_count; i++) {
        edge_count += graph->atoms[i].bonds_count;
    }
    edge_count /= 2; // Each bond is counted twice
    
    // Mesh connectivity: ratio of actual edges to maximum possible edges
    double max_edges = (graph->atom_count * (graph->atom_count - 1)) / 2.0;
    double mesh_connectivity = edge_count / max_edges;
    
    // Mesh complexity: weighted by atom types and distances
    double mesh_complexity = 0.0;
    for (int i = 0; i < graph->atom_count; i++) {
        double atom_weight = 0.0;
        
        // Weight by atom type
        switch (graph->atoms[i].type) {
            case ATOM_C: atom_weight = 1.0; break;
            case ATOM_N: atom_weight = 2.0; break;
            case ATOM_O: atom_weight = 2.0; break;
            case ATOM_S: atom_weight = 3.0; break;
            case ATOM_P: atom_weight = 3.0; break;
            case ATOM_F: 
            case ATOM_CL: 
            case ATOM_BR: 
            case ATOM_I: atom_weight = 4.0; break;
            default: atom_weight = 1.5; break;
        }
        
        // Weight by connectivity degree (number of bonds)
        mesh_complexity += atom_weight * graph->atoms[i].bonds_count;
    }
    
    // Normalize by atom count
    mesh_complexity /= graph->atom_count;
    
    // Combine connectivity and complexity for final score
    double mesh_score = (mesh_connectivity * 0.4) + (mesh_complexity * 0.6);
    
    free_distance_matrix(dist_matrix, graph->atom_count);
    return mesh_score;
}

// Maximum Common Subgraph Size Index (approximation using atom type distribution)
double calculate_mcs_index(MolecularGraph* graph) {
    if (graph->atom_count == 0) return 0.0;
    
    // Count atom types
    int atom_type_counts[ATOM_TYPES_COUNT] = {0};
    for (int i = 0; i < graph->atom_count; i++) {
        atom_type_counts[graph->atoms[i].type]++;
    }
    
    // Calculate a simple diversity index
    double sum_squared = 0.0;
    for (int i = 0; i < ATOM_TYPES_COUNT; i++) {
        double fraction = (double)atom_type_counts[i] / graph->atom_count;
        sum_squared += fraction * fraction;
    }
    
    // Diversity index (1 = all different, 0 = all same)
    double diversity = 1.0 - sum_squared;
    
    // Calculate connectivity index
    double connectivity_sum = 0.0;
    for (int i = 0; i < graph->atom_count; i++) {
        connectivity_sum += graph->atoms[i].bonds_count;
    }
    double avg_connectivity = connectivity_sum / (2.0 * graph->atom_count); // Divide by 2 since each bond is counted twice
    
    // Normalize by theoretical maximum
    double max_connectivity = (graph->atom_count > 1) ? 3.0 : 0.0; // Maximum avg connectivity in organic compounds
    double connectivity_factor = avg_connectivity / max_connectivity;
    
    // Combine factors - lower diversity and higher connectivity suggest larger common subgraphs
    double mcs_index = (1.0 - diversity) * 0.5 + connectivity_factor * 0.5;
    
    return mcs_index;
}

// Vertex-Edge Path Descriptors
double calculate_vertex_edge_path(MolecularGraph* graph) {
    if (graph->atom_count < 2) return 0.0;
    
    int** dist_matrix = compute_distance_matrix(graph);
    
    // Calculate path descriptor based on path distribution
    int path_counts[10] = {0}; // Count paths of length 1-10
    int max_path_len = (graph->max_path_length > 10) ? 10 : graph->max_path_length;
    
    for (int i = 0; i < graph->atom_count; i++) {
        for (int j = i+1; j < graph->atom_count; j++) {
            if (dist_matrix[i][j] > 0 && dist_matrix[i][j] <= max_path_len) {
                path_counts[dist_matrix[i][j] - 1]++;
            }
        }
    }
    
    // Weight paths by length (longer paths are more significant)
    double weighted_sum = 0.0;
    int total_paths = 0;
    
    for (int i = 0; i < max_path_len; i++) {
        weighted_sum += path_counts[i] * (i + 1);
        total_paths += path_counts[i];
    }
    
    double path_descriptor = (total_paths > 0) ? weighted_sum / total_paths : 0.0;
    
    free_distance_matrix(dist_matrix, graph->atom_count);
    return path_descriptor;
}

// McFarland Complexity Index
double calculate_mcfarland_complexity(MolecularGraph* graph) {
    if (graph->atom_count == 0) return 0.0;
    
    // Count different atom environments
    int count_aromatic = 0;
    int count_heteroatoms = 0;
    int count_ring_atoms = 0;
    int count_branch_points = 0;
    
    for (int i = 0; i < graph->atom_count; i++) {
        // Count aromatic atoms
        if (graph->atoms[i].aromatic) {
            count_aromatic++;
        }
        
        // Count heteroatoms (non-carbon)
        if (graph->atoms[i].type != ATOM_C) {
            count_heteroatoms++;
        }
        
        // Count branch points (atoms with 3+ connections)
        if (graph->atoms[i].bonds_count >= 3) {
            count_branch_points++;
        }
    }
    
    // Use Floyd-Warshall results to identify ring atoms
    int** dist_matrix = compute_distance_matrix(graph);
    
    for (int i = 0; i < graph->atom_count; i++) {
        for (int j = 0; j < graph->atoms[i].bonds_count; j++) {
            int neighbor = graph->atoms[i].bonded_atoms[j];
            
            // If there's an alternative path between bonded atoms, they're in a ring
            if (dist_matrix[i][neighbor] < 9999 && dist_matrix[i][neighbor] > 1) {
                count_ring_atoms++;
                break;
            }
        }
    }
    
    // Calculate complexity index
    double complexity = 0.0;
    
    // Weight factors
    complexity += count_aromatic * 1.5;
    complexity += count_heteroatoms * 2.0;
    complexity += count_ring_atoms * 1.0;
    complexity += count_branch_points * 2.5;
    
    // Normalize by atom count
    double normalized_complexity = complexity / graph->atom_count;
    
    free_distance_matrix(dist_matrix, graph->atom_count);
    return normalized_complexity;
}

// Cycle Connectivity Index
double calculate_cycle_connectivity(MolecularGraph* graph) {
    if (graph->atom_count < 3) return 0.0; // Need at least 3 atoms for a cycle
    
    int** dist_matrix = compute_distance_matrix(graph);
    
    // Identify rings using path analysis
    int cycle_count = 0;
    int* cycle_atoms = (int*)calloc(graph->atom_count, sizeof(int));
    
    for (int i = 0; i < graph->atom_count; i++) {
        for (int j = 0; j < graph->atoms[i].bonds_count; j++) {
            int neighbor = graph->atoms[i].bonded_atoms[j];
            
            // If there's an alternative path between bonded atoms, a cycle exists
            if (i < neighbor && dist_matrix[i][neighbor] < 9999 && dist_matrix[i][neighbor] > 1) {
                cycle_count++;
                cycle_atoms[i] = 1;
                cycle_atoms[neighbor] = 1;
            }
        }
    }
    
    // Count atoms in rings
    int atoms_in_rings = 0;
    for (int i = 0; i < graph->atom_count; i++) {
        if (cycle_atoms[i]) atoms_in_rings++;
    }
    
    // If no cycles found, return 0
    if (cycle_count == 0) {
        free(cycle_atoms);
        free_distance_matrix(dist_matrix, graph->atom_count);
        return 0.0;
    }
    
    // Calculate connectivity between cycles
    int cycle_connectivity = 0;
    for (int i = 0; i < graph->atom_count; i++) {
        if (cycle_atoms[i]) {
            for (int j = 0; j < graph->atoms[i].bonds_count; j++) {
                int neighbor = graph->atoms[i].bonded_atoms[j];
                if (cycle_atoms[neighbor] && i < neighbor) {
                    cycle_connectivity++;
                }
            }
        }
    }
    
    // Normalize by potential connections
    double max_connectivity = atoms_in_rings * (atoms_in_rings - 1) / 2.0;
    double connectivity_index = (max_connectivity > 0) ? cycle_connectivity / max_connectivity : 0.0;
    
    free(cycle_atoms);
    free_distance_matrix(dist_matrix, graph->atom_count);
    return connectivity_index;
}

// Fragment Complexity Score
double calculate_fragment_complexity(MolecularGraph* graph) {
    if (graph->atom_count == 0) return 0.0;
    
    // Count different fragment types
    int count_aromatic_rings = 0;
    int count_aliphatic_rings = 0;
    int count_heteroaromatic = 0;
    int count_amide_bonds = 0;
    int count_ether_links = 0;
    int count_carbonyl_groups = 0;
    
    // Use distance matrix to identify rings
    int** dist_matrix = compute_distance_matrix(graph);
    
    // First identify all rings
    bool* in_ring = (bool*)calloc(graph->atom_count, sizeof(bool));
    bool* in_aromatic_ring = (bool*)calloc(graph->atom_count, sizeof(bool));
    
    for (int i = 0; i < graph->atom_count; i++) {
        for (int j = 0; j < graph->atoms[i].bonds_count; j++) {
            int neighbor = graph->atoms[i].bonded_atoms[j];
            
            // If there's an alternative path between bonded atoms, a cycle exists
            if (i < neighbor && dist_matrix[i][neighbor] < 9999 && dist_matrix[i][neighbor] > 1) {
                in_ring[i] = true;
                in_ring[neighbor] = true;
                
                // Check if aromatic
                if (graph->atoms[i].aromatic && graph->atoms[neighbor].aromatic) {
                    in_aromatic_ring[i] = true;
                    in_aromatic_ring[neighbor] = true;
                }
            }
        }
    }
    
    // Count ring types
    bool* counted_ring = (bool*)calloc(graph->atom_count, sizeof(bool));
    for (int i = 0; i < graph->atom_count; i++) {
        if (in_aromatic_ring[i] && !counted_ring[i]) {
            count_aromatic_rings++;
            
            // Mark all atoms in this aromatic ring as counted
            for (int j = 0; j < graph->atom_count; j++) {
                if (in_aromatic_ring[j] && dist_matrix[i][j] < 9999) {
                    counted_ring[j] = true;
                }
            }
            
            // Check if it's heteroaromatic (contains N, O, S)
            bool has_heteroatom = false;
            for (int j = 0; j < graph->atom_count; j++) {
                if (in_aromatic_ring[j] && dist_matrix[i][j] < 9999) {
                    if (graph->atoms[j].type == ATOM_N || 
                        graph->atoms[j].type == ATOM_O || 
                        graph->atoms[j].type == ATOM_S) {
                        has_heteroatom = true;
                        break;
                    }
                }
            }
            if (has_heteroatom) count_heteroaromatic++;
        }
        else if (in_ring[i] && !counted_ring[i] && !in_aromatic_ring[i]) {
            count_aliphatic_rings++;
            
            // Mark all atoms in this aliphatic ring as counted
            for (int j = 0; j < graph->atom_count; j++) {
                if (in_ring[j] && !in_aromatic_ring[j] && dist_matrix[i][j] < 9999) {
                    counted_ring[j] = true;
                }
            }
        }
    }
    
    // Count functional groups (simple approximation based on atom patterns)
    for (int i = 0; i < graph->atom_count; i++) {
        // Carbonyl group detection (C=O pattern)
        if (graph->atoms[i].type == ATOM_C) {
            for (int j = 0; j < graph->atoms[i].bonds_count; j++) {
                int neighbor = graph->atoms[i].bonded_atoms[j];
                if (graph->atoms[neighbor].type == ATOM_O && 
                    graph->atoms[i].bond_types[j] == BOND_DOUBLE) {
                    count_carbonyl_groups++;
                    
                    // Amide detection (look for N connected to C=O)
                    for (int k = 0; k < graph->atoms[i].bonds_count; k++) {
                        int n_neighbor = graph->atoms[i].bonded_atoms[k];
                        if (graph->atoms[n_neighbor].type == ATOM_N) {
                            count_amide_bonds++;
                            break;
                        }
                    }
                }
            }
        }
        
        // Ether link detection (C-O-C pattern)
        if (graph->atoms[i].type == ATOM_O) {
            int c_connections = 0;
            for (int j = 0; j < graph->atoms[i].bonds_count; j++) {
                int neighbor = graph->atoms[i].bonded_atoms[j];
                if (graph->atoms[neighbor].type == ATOM_C) c_connections++;
            }
            if (c_connections >= 2) count_ether_links++;
        }
    }
    
    // Calculate complexity score
    double complexity = 0.0;
    
    // Weight factors
    complexity += count_aromatic_rings * 2.0;
    complexity += count_aliphatic_rings * 1.5;
    complexity += count_heteroaromatic * 2.5;
    complexity += count_amide_bonds * 2.0;
    complexity += count_ether_links * 1.0;
    complexity += count_carbonyl_groups * 1.5;
    
    // Normalize by atom count
    double normalized_complexity = complexity / graph->atom_count;
    
    free(in_ring);
    free(in_aromatic_ring);
    free(counted_ring);
    free_distance_matrix(dist_matrix, graph->atom_count);
    return normalized_complexity;
}

// Topological Polar Buried Surface Area (approximation)
double calculate_topological_polar_bsa(MolecularGraph* graph) {
    if (graph->atom_count == 0) return 0.0;
    
    // Identify polar atoms
    bool* is_polar = (bool*)calloc(graph->atom_count, sizeof(bool));
    for (int i = 0; i < graph->atom_count; i++) {
        // N, O, S, P and halogens are considered polar
        if (graph->atoms[i].type == ATOM_N || 
            graph->atoms[i].type == ATOM_O || 
            graph->atoms[i].type == ATOM_S || 
            graph->atoms[i].type == ATOM_P ||
            graph->atoms[i].type == ATOM_F ||
            graph->atoms[i].type == ATOM_CL ||
            graph->atoms[i].type == ATOM_BR ||
            graph->atoms[i].type == ATOM_I) {
            is_polar[i] = true;
        }
    }
    
    // Count polar atoms
    int total_polar_atoms = 0;
    for (int i = 0; i < graph->atom_count; i++) {
        if (is_polar[i]) total_polar_atoms++;
    }
    
    // If no polar atoms, return 0
    if (total_polar_atoms == 0) {
        free(is_polar);
        return 0.0;
    }
    
    // Calculate distance matrix
    int** dist_matrix = compute_distance_matrix(graph);
    
    // Calculate buried surface area by estimating how many polar atoms are "buried"
    // (surrounded by non-polar atoms)
    int buried_polar_atoms = 0;
    for (int i = 0; i < graph->atom_count; i++) {
        if (is_polar[i]) {
            int non_polar_neighbors = 0;
            int total_neighbors = 0;
            
            // Count immediate neighbors
            for (int j = 0; j < graph->atoms[i].bonds_count; j++) {
                int neighbor = graph->atoms[i].bonded_atoms[j];
                total_neighbors++;
                if (!is_polar[neighbor]) non_polar_neighbors++;
            }
            
            // If more than 50% of neighbors are non-polar, consider it buried
            if (total_neighbors > 0 && (double)non_polar_neighbors / total_neighbors > 0.5) {
                buried_polar_atoms++;
            }
        }
    }
    
    // Calculate buried surface area ratio
    double buried_ratio = (double)buried_polar_atoms / total_polar_atoms;
    
    free(is_polar);
    free_distance_matrix(dist_matrix, graph->atom_count);
    return buried_ratio;
}

// Eccentric Distance Index
double calculate_eccentric_distance(MolecularGraph* graph) {
    if (graph->atom_count < 2) return 0.0;
    
    int** dist_matrix = compute_distance_matrix(graph);
    int* eccentricity = compute_eccentricity(dist_matrix, graph->atom_count);
    
    // Calculate eccentric distance index
    double ecc_dist_sum = 0.0;
    
    for (int i = 0; i < graph->atom_count; i++) {
        for (int j = 0; j < graph->atom_count; j++) {
            if (i != j && dist_matrix[i][j] != 9999) {
                ecc_dist_sum += eccentricity[i] * dist_matrix[i][j];
            }
        }
    }
    
    double normalized_ecc_dist = ecc_dist_sum / (graph->atom_count * graph->atom_count);
    
    free(eccentricity);
    free_distance_matrix(dist_matrix, graph->atom_count);
    return normalized_ecc_dist;
}

// Molecular Scaffolding Index
double calculate_molecular_scaffolding(MolecularGraph* graph) {
    if (graph->atom_count == 0) return 0.0;
    
    // Calculate distance matrix
    int** dist_matrix = compute_distance_matrix(graph);
    
    // Identify core scaffolding atoms (atoms with high centrality)
    double* centrality = (double*)calloc(graph->atom_count, sizeof(double));
    
    // Calculate closeness centrality for each atom
    for (int i = 0; i < graph->atom_count; i++) {
        double sum_distances = 0.0;
        int reachable_atoms = 0;
        
        for (int j = 0; j < graph->atom_count; j++) {
            if (i != j && dist_matrix[i][j] != 9999) {
                sum_distances += dist_matrix[i][j];
                reachable_atoms++;
            }
        }
        
        if (reachable_atoms > 0) {
            centrality[i] = (double)reachable_atoms / sum_distances;
        }
    }
    
    // Calculate average centrality
    double avg_centrality = 0.0;
    for (int i = 0; i < graph->atom_count; i++) {
        avg_centrality += centrality[i];
    }
    avg_centrality /= graph->atom_count;
    
    // Count atoms in the scaffold (atoms with above-average centrality)
    int scaffold_atoms = 0;
    for (int i = 0; i < graph->atom_count; i++) {
        if (centrality[i] > avg_centrality) scaffold_atoms++;
    }
    
    // Calculate scaffold ratio
    double scaffold_ratio = (double)scaffold_atoms / graph->atom_count;
    
    // Weight by connectivity of scaffold atoms
    double connectivity_weight = 0.0;
    for (int i = 0; i < graph->atom_count; i++) {
        if (centrality[i] > avg_centrality) {
            connectivity_weight += graph->atoms[i].bonds_count;
        }
    }
    
    double avg_scaffold_connectivity = (scaffold_atoms > 0) ? 
                                       connectivity_weight / scaffold_atoms : 0.0;
    
    // Combine scaffold ratio and connectivity
    double scaffolding_index = scaffold_ratio * 0.6 + 
                              (avg_scaffold_connectivity / 4.0) * 0.4; // Normalize by max connectivity
    
    free(centrality);
    free_distance_matrix(dist_matrix, graph->atom_count);
    return scaffolding_index;
}


double calculateMolecularMesh(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc((const Context*)context);
    if (!smiles) return 0.0;
    
    MolecularGraph* graph = parse_smiles(smiles);
    double result = calculate_molecular_mesh(graph);
    free_molecular_graph(graph);
    return result;
}

double calculateMCSIndex(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc((const Context*)context);
    if (!smiles) return 0.0;
    
    MolecularGraph* graph = parse_smiles(smiles);
    double result = calculate_mcs_index(graph);
    free_molecular_graph(graph);
    return result;
}

double calculateVertexEdgePath(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc((const Context*)context);
    if (!smiles) return 0.0;
    
    MolecularGraph* graph = parse_smiles(smiles);
    double result = calculate_vertex_edge_path(graph);
    free_molecular_graph(graph);
    return result;
}

double calculateMcFarlandComplexity(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc((const Context*)context);
    if (!smiles) return 0.0;
    
    MolecularGraph* graph = parse_smiles(smiles);
    double result = calculate_mcfarland_complexity(graph);
    free_molecular_graph(graph);
    return result;
}

double calculateCycleConnectivity(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc((const Context*)context);
    if (!smiles) return 0.0;
    
    MolecularGraph* graph = parse_smiles(smiles);
    double result = calculate_cycle_connectivity(graph);
    free_molecular_graph(graph);
    return result;
}

double calculateFragmentComplexity(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc((const Context*)context);
    if (!smiles) return 0.0;
    
    MolecularGraph* graph = parse_smiles(smiles);
    double result = calculate_fragment_complexity(graph);
    free_molecular_graph(graph);
    return result;
}

double calculateTopologicalPolarBSA(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc((const Context*)context);
    if (!smiles) return 0.0;
    
    MolecularGraph* graph = parse_smiles(smiles);
    double result = calculate_topological_polar_bsa(graph);
    free_molecular_graph(graph);
    return result;
}

double calculateEccentricDistance(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc((const Context*)context);
    if (!smiles) return 0.0;
    
    MolecularGraph* graph = parse_smiles(smiles);
    double result = calculate_eccentric_distance(graph);
    free_molecular_graph(graph);
    return result;
}

double calculateMolecularScaffolding(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc((const Context*)context);
    if (!smiles) return 0.0;
    
    MolecularGraph* graph = parse_smiles(smiles);
    double result = calculate_molecular_scaffolding(graph);
    free_molecular_graph(graph);
    return result;
}

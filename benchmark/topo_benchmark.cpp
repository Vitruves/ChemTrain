#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <chrono>
#include <getopt.h>
#include <ctime>
#include <cctype>
#include <iomanip>

// Atom types commonly found in SMILES
enum AtomType {
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
};

// Bond types
enum BondType {
    BOND_SINGLE = 1,
    BOND_DOUBLE = 2,
    BOND_TRIPLE = 3,
    BOND_AROMATIC = 4
};

// Atom representation in a molecular graph
struct Atom {
    AtomType type;
    int index;
    int aromatic;
    int charge;
    int bonds_count;
    std::vector<int> bonded_atoms;
    std::vector<int> bond_types;
    
    Atom() : type(ATOM_C), index(0), aromatic(0), charge(0), bonds_count(0) {}
};

// Molecular graph representation
class MolecularGraph {
public:
    int atom_count;
    std::vector<Atom> atoms;
    int rings_count;
    std::vector<int> ring_sizes;
    int max_path_length;
    
    MolecularGraph() : atom_count(0), rings_count(0), max_path_length(0) {}
    
    // Parse a SMILES string into a molecular graph
    void parse_smiles(const std::string& smiles);
    
    // Compute distance matrix
    std::vector<std::vector<int>> compute_distance_matrix();
    
    // Descriptor calculations
    double calculate_molecular_mesh();
    double calculate_eccentric_distance();
    double calculate_vertex_edge_path();
    
private:
    // Helper functions for SMILES parsing
    AtomType determine_atom_type(const std::string& smiles, int pos);
    bool is_aromatic(char c);
};

AtomType MolecularGraph::determine_atom_type(const std::string& smiles, int pos) {
    if (smiles[pos] == 'C') {
        if (pos+1 < smiles.length() && smiles[pos+1] == 'l') return ATOM_CL;
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
    } else if (smiles[pos] == 'B' && pos+1 < smiles.length() && smiles[pos+1] == 'r') {
        return ATOM_BR;
    } else if (smiles[pos] == 'I') {
        return ATOM_I;
    }
    return ATOM_OTHER;
}

bool MolecularGraph::is_aromatic(char c) {
    return (c == 'c' || c == 'n' || c == 'o' || c == 's' || c == 'p');
}

void MolecularGraph::parse_smiles(const std::string& smiles) {
    int len = smiles.length();
    
    // First pass: count atoms
    int atom_count = 0;
    for (int i = 0; i < len; i++) {
        if (std::isalpha(smiles[i]) || smiles[i] == '[') {
            atom_count++;
            // Skip multi-character atoms
            if (smiles[i] == '[') {
                while (i < len && smiles[i] != ']') i++;
            } else if ((smiles[i] == 'C' && i+1 < len && smiles[i+1] == 'l') || 
                      (smiles[i] == 'B' && i+1 < len && smiles[i+1] == 'r')) {
                i++;
            }
        }
    }
    
    // Allocate atoms
    this->atom_count = atom_count;
    this->atoms.resize(atom_count);
    
    // Initialize rings count
    this->rings_count = 0;
    this->max_path_length = 0;
    
    // Track brackets for ring closures
    std::vector<int> ring_openings(10, -1);
    
    // Second pass: create atoms and basic connectivity
    int atom_idx = 0;
    int prev_atom = -1;
    std::vector<int> branch_points;
    
    for (int i = 0; i < len; i++) {
        if (std::isalpha(smiles[i]) || smiles[i] == '[') {
            // Process atom
            this->atoms[atom_idx].index = atom_idx;
            
            if (smiles[i] == '[') {
                // Complex atom notation - skip for simplicity
                this->atoms[atom_idx].type = ATOM_OTHER;
                while (i < len && smiles[i] != ']') i++;
            } else if (std::islower(smiles[i])) {
                // Aromatic atom
                this->atoms[atom_idx].aromatic = 1;
                if (smiles[i] == 'c') this->atoms[atom_idx].type = ATOM_C;
                else if (smiles[i] == 'n') this->atoms[atom_idx].type = ATOM_N;
                else if (smiles[i] == 'o') this->atoms[atom_idx].type = ATOM_O;
                else if (smiles[i] == 's') this->atoms[atom_idx].type = ATOM_S;
                else if (smiles[i] == 'p') this->atoms[atom_idx].type = ATOM_P;
                else this->atoms[atom_idx].type = ATOM_OTHER;
            } else {
                // Regular atom
                this->atoms[atom_idx].type = determine_atom_type(smiles, i);
                if ((smiles[i] == 'C' && i+1 < len && smiles[i+1] == 'l') || 
                    (smiles[i] == 'B' && i+1 < len && smiles[i+1] == 'r')) {
                    i++;
                }
            }
            
            // Connect to previous atom if applicable
            if (prev_atom != -1) {
                // Add bond
                int bond_type = BOND_SINGLE;
                if (this->atoms[prev_atom].aromatic && this->atoms[atom_idx].aromatic) {
                    bond_type = BOND_AROMATIC;
                }
                
                this->atoms[prev_atom].bonded_atoms.push_back(atom_idx);
                this->atoms[prev_atom].bond_types.push_back(bond_type);
                this->atoms[prev_atom].bonds_count++;
                
                this->atoms[atom_idx].bonded_atoms.push_back(prev_atom);
                this->atoms[atom_idx].bond_types.push_back(bond_type);
                this->atoms[atom_idx].bonds_count++;
            }
            
            prev_atom = atom_idx;
            atom_idx++;
            
            // Check for ring closures
            if (i+1 < len && std::isdigit(smiles[i+1])) {
                int ring_num = smiles[i+1] - '0';
                i++;
                
                if (ring_openings[ring_num] == -1) {
                    // Start of ring
                    ring_openings[ring_num] = prev_atom;
                } else {
                    // End of ring, connect atoms
                    int start_atom = ring_openings[ring_num];
                    int end_atom = prev_atom;
                    
                    // Add bond
                    int bond_type = BOND_SINGLE;
                    if (this->atoms[start_atom].aromatic && this->atoms[end_atom].aromatic) {
                        bond_type = BOND_AROMATIC;
                    }
                    
                    this->atoms[start_atom].bonded_atoms.push_back(end_atom);
                    this->atoms[start_atom].bond_types.push_back(bond_type);
                    this->atoms[start_atom].bonds_count++;
                    
                    this->atoms[end_atom].bonded_atoms.push_back(start_atom);
                    this->atoms[end_atom].bond_types.push_back(bond_type);
                    this->atoms[end_atom].bonds_count++;
                    
                    // Increment rings count
                    this->rings_count++;
                    
                    // Reset ring opening
                    ring_openings[ring_num] = -1;
                }
            }
        } else if (smiles[i] == '(') {
            // Branch start
            branch_points.push_back(prev_atom);
        } else if (smiles[i] == ')') {
            // Branch end
            if (!branch_points.empty()) {
                prev_atom = branch_points.back();
                branch_points.pop_back();
            }
        } else if (smiles[i] == '=') {
            // Double bond - update previous bond
            if (prev_atom >= 0 && atom_idx > 0) {
                int prev_bond_idx = this->atoms[prev_atom].bonds_count - 1;
                if (prev_bond_idx >= 0) {
                    this->atoms[prev_atom].bond_types[prev_bond_idx] = BOND_DOUBLE;
                    
                    int other_atom = this->atoms[prev_atom].bonded_atoms[prev_bond_idx];
                    for (int j = 0; j < this->atoms[other_atom].bonds_count; j++) {
                        if (this->atoms[other_atom].bonded_atoms[j] == prev_atom) {
                            this->atoms[other_atom].bond_types[j] = BOND_DOUBLE;
                            break;
                        }
                    }
                }
            }
        } else if (smiles[i] == '#') {
            // Triple bond - update previous bond
            if (prev_atom >= 0 && atom_idx > 0) {
                int prev_bond_idx = this->atoms[prev_atom].bonds_count - 1;
                if (prev_bond_idx >= 0) {
                    this->atoms[prev_atom].bond_types[prev_bond_idx] = BOND_TRIPLE;
                    
                    int other_atom = this->atoms[prev_atom].bonded_atoms[prev_bond_idx];
                    for (int j = 0; j < this->atoms[other_atom].bonds_count; j++) {
                        if (this->atoms[other_atom].bonded_atoms[j] == prev_atom) {
                            this->atoms[other_atom].bond_types[j] = BOND_TRIPLE;
                            break;
                        }
                    }
                }
            }
        }
    }
}

std::vector<std::vector<int>> MolecularGraph::compute_distance_matrix() {
    int n = this->atom_count;
    std::vector<std::vector<int>> dist(n, std::vector<int>(n, 9999)); // "Infinity"
    
    // Initialize with direct connections
    for (int i = 0; i < n; i++) {
        dist[i][i] = 0; // Distance to self is 0
        for (int j = 0; j < this->atoms[i].bonds_count; j++) {
            int neighbor = this->atoms[i].bonded_atoms[j];
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
    this->max_path_length = max_path;
    
    return dist;
}

double MolecularGraph::calculate_molecular_mesh() {
    if (this->atom_count == 0) return 0.0;
    
    auto dist_matrix = compute_distance_matrix();
    
    // Calculate mesh complexity based on connectivity patterns
    double mesh_value = 0.0;
    
    // Sum of all path lengths
    double total_paths = 0.0;
    int path_count = 0;
    
    for (int i = 0; i < this->atom_count; i++) {
        for (int j = i+1; j < this->atom_count; j++) {
            if (dist_matrix[i][j] < 9999) {
                total_paths += dist_matrix[i][j];
                path_count++;
            }
        }
    }
    
    // Average path length
    double avg_path = (path_count > 0) ? total_paths / path_count : 0;
    
    // Mesh complexity is a function of atom count, ring count, and path distribution
    mesh_value = (this->atom_count * 0.5) + (this->rings_count * 1.5) + (avg_path * 0.3);
    
    // Normalize by dividing by number of atoms
    mesh_value /= this->atom_count;
    
    return mesh_value;
}

double MolecularGraph::calculate_eccentric_distance() {
    if (this->atom_count <= 1) return 0.0;
    
    auto dist_matrix = compute_distance_matrix();
    
    // Calculate eccentricity for each atom (maximum distance to any other atom)
    std::vector<int> eccentricity(this->atom_count, 0);
    
    for (int i = 0; i < this->atom_count; i++) {
        for (int j = 0; j < this->atom_count; j++) {
            if (dist_matrix[i][j] < 9999 && dist_matrix[i][j] > eccentricity[i]) {
                eccentricity[i] = dist_matrix[i][j];
            }
        }
    }
    
    // Calculate eccentric distance sum
    double ecc_dist_sum = 0.0;
    
    for (int i = 0; i < this->atom_count; i++) {
        ecc_dist_sum += eccentricity[i];
    }
    
    return ecc_dist_sum;
}

double MolecularGraph::calculate_vertex_edge_path() {
    if (this->atom_count <= 1) return 0.0;
    
    auto dist_matrix = compute_distance_matrix();
    
    // Calculate vertex-edge path index
    double vep_index = 0.0;
    int total_bonds = 0;
    
    // Count total bonds
    for (int i = 0; i < this->atom_count; i++) {
        total_bonds += this->atoms[i].bonds_count;
    }
    total_bonds /= 2;  // Each bond is counted twice
    
    // Calculate vertex contribution
    double vertex_contrib = this->atom_count * 0.5;
    
    // Calculate edge contribution
    double edge_contrib = total_bonds * 0.7;
    
    // Calculate path contribution
    double path_contrib = 0.0;
    for (int i = 0; i < this->atom_count; i++) {
        for (int j = i+1; j < this->atom_count; j++) {
            if (dist_matrix[i][j] < 9999) {
                path_contrib += 1.0 / dist_matrix[i][j];
            }
        }
    }
    
    // Combine all contributions
    vep_index = vertex_contrib + edge_contrib + (path_contrib * 0.3);
    
    return vep_index;
}

void print_usage() {
    std::cout << "Usage: topo_benchmark_cpp [options]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -i, --input <file>   Input file with SMILES strings (one per line)" << std::endl;
    std::cout << "  -o, --output <file>  Output file for results" << std::endl;
    std::cout << "  -s, --smiles <str>   Direct SMILES string input" << std::endl;
    std::cout << "  -h, --help           Display this help message" << std::endl;
}

int main(int argc, char *argv[]) {
    std::string input_file;
    std::string output_file;
    std::string smiles_str;
    
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
    if (input_file.empty() && smiles_str.empty()) {
        std::cerr << "Error: Either input file (-i) or SMILES string (-s) is required" << std::endl;
        print_usage();
        return 1;
    }
    
    // Open output file if specified, otherwise use stdout
    std::ostream* out_stream = &std::cout;
    std::ofstream out_file;
    
    if (!output_file.empty()) {
        out_file.open(output_file);
        if (!out_file.is_open()) {
            std::cerr << "Error: Could not open output file " << output_file << std::endl;
            return 1;
        }
        out_stream = &out_file;
    }
    
    // Write header
    *out_stream << "SMILES,MolecularMesh,EccentricDistance,VertexEdgePath,Time(ms)" << std::endl;
    
    // Process SMILES from direct input
    if (!smiles_str.empty()) {
        auto start = std::chrono::high_resolution_clock::now();
        
        // Create and parse molecular graph
        MolecularGraph graph;
        graph.parse_smiles(smiles_str);
        
        // Calculate descriptors
        double mesh = graph.calculate_molecular_mesh();
        double ecc_dist = graph.calculate_eccentric_distance();
        double vep = graph.calculate_vertex_edge_path();
        
        // Calculate time
        auto end = std::chrono::high_resolution_clock::now();
        double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        
        // Output results
        *out_stream << smiles_str << "," 
                   << std::fixed << std::setprecision(6) << mesh << "," 
                   << std::fixed << std::setprecision(6) << ecc_dist << "," 
                   << std::fixed << std::setprecision(6) << vep << "," 
                   << std::fixed << std::setprecision(2) << time_ms << std::endl;
    }
    
    // Process SMILES from input file
    if (!input_file.empty()) {
        std::ifstream in_file(input_file);
        if (!in_file.is_open()) {
            std::cerr << "Error: Could not open input file " << input_file << std::endl;
            return 1;
        }
        
        std::string line;
        while (std::getline(in_file, line)) {
            // Skip empty lines
            if (line.empty()) continue;
            
            auto start = std::chrono::high_resolution_clock::now();
            
            // Create and parse molecular graph
            MolecularGraph graph;
            graph.parse_smiles(line);
            
            // Calculate descriptors
            double mesh = graph.calculate_molecular_mesh();
            double ecc_dist = graph.calculate_eccentric_distance();
            double vep = graph.calculate_vertex_edge_path();
            
            // Calculate time
            auto end = std::chrono::high_resolution_clock::now();
            double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
            
            // Output results
            *out_stream << line << "," 
                       << std::fixed << std::setprecision(6) << mesh << "," 
                       << std::fixed << std::setprecision(6) << ecc_dist << "," 
                       << std::fixed << std::setprecision(6) << vep << "," 
                       << std::fixed << std::setprecision(2) << time_ms << std::endl;
        }
        
        in_file.close();
    }
    
    // Close output file if opened
    if (out_file.is_open()) {
        out_file.close();
    }
    
    return 0;
}

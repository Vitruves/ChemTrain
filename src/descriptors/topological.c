#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <limits.h> // For INT_MAX

#define INITIAL_ATOM_BUFFER_SIZE 32 
#define MAX_BONDS_PER_ATOM_ESTIMATE 8

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
    ATOM_B, // Boron
    ATOM_SI, // Silicon
    ATOM_SE, // Selenium
    ATOM_AS, // Arsenic
    ATOM_UNKNOWN, // For elements not explicitly handled or complex bracket atoms
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
    int atomic_number; // New field
    double electronegativity; // New field
    double atomic_weight; // New field
    int group_number; // New field
    int valence_electrons; // New field
    int hybridization; // 0:unknown, 1:sp, 2:sp2, 3:sp3 (Estimated)
    int formal_charge; // New field
    int explicit_hydrogens; // New field
    bool is_in_bracket; // New field
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
    int* ring_sizes; // ring_sizes should probably be an array of structs/arrays for multiple rings
    int max_path_length; // Calculated property, not directly from SMILES usually
    // Add fields to manage ring closures during parsing if not already present
    // For example, a temporary structure to hold ring number and atom index
} MolecularGraph;

// --- Constants for Descriptors ---
static const double AVG_BOND_LENGTH_A = 1.5; // Angstroms
static const double TETRAHEDRAL_BOND_ANGLE_DEG = 109.5;
static const double RING_STRAIN_INCREMENT_KCAL_MOL = 1.2;
static const double AVG_ATOMIC_RADIUS_ORGANIC_A = 3.8; // This seems high, usually ~0.7-1.5 A. Using as specified.
static const double AVG_RESONANCE_ENERGY_DIFF_KCAL_MOL = 30.0; // Placeholder
static const double GAS_CONSTANT_J_MOL_K = 8.314;
static const double LN10_TO_LOG2 = 2.303; // ln(10) - but should be ln(2) for bits if shannon
static const double AVOGADRO_NUMBER_MOL_INV = 6.022e23;
static const double ROTAMER_ENERGY_DIFF_KCAL_MOL = 0.8;
static const double RYDBERG_CONSTANT_EV = 13.6;
static const double AVG_RING_CLOSURE_DISTANCE_A = 2.8;
static const double PAULING_EN_CARBON = 2.55; // More common value than 2.20
static const double THERMOCHEMICAL_CALORIE_TO_JOULE = 4.184;
static const double BOLTZMANN_CONSTANT_J_K = 1.38e-23;
static const double LN2_INFO_CONTENT_BIT = 0.693147;
static const double H_BOND_ENERGY_KCAL_MOL_N = 9.0; // Placeholder
static const double CARBONYL_STABILIZATION_KCAL_MOL = 23.4; // Placeholder
static const double AVG_COVALENT_RADIUS_A = 0.77;
static const double AVG_CC_BOND_ENERGY_KCAL_MOL = 83.0;
static const double CONFORMATIONAL_BARRIER_KCAL_MOL = 4.0;
static const double VDW_RADIUS_CARBON_A = 1.7;
static const double PAULING_EN_FLUORINE = 3.98; // Common value
static const double PAULING_BOND_ORDER_CONVERSION = 1.412; // Placeholder
static const double BENZENE_RESONANCE_ENERGY_KCAL_MOL = 36.0;
static const double VOLUME_INCREMENT_BRANCH_A3 = 3.5; // Angstrom^3
static const double ATOMIC_MASS_UNIT_U = 16.0; // This is for Oxygen, generic placeholder

// Helper to get atomic number (simplified)
int get_atomic_number_simple(AtomType type, const char* atom_label_in_bracket) {
    if (atom_label_in_bracket && atom_label_in_bracket[0] != '\0') {
        // Basic support for elements in brackets
        if (strncmp(atom_label_in_bracket, "Cl", 2) == 0) return 17;
        if (strncmp(atom_label_in_bracket, "Br", 2) == 0) return 35;
        if (strncmp(atom_label_in_bracket, "Si", 2) == 0) return 14;
        if (strncmp(atom_label_in_bracket, "As", 2) == 0) return 33;
        if (strncmp(atom_label_in_bracket, "Se", 2) == 0) return 34;
        // Single letter elements
        char first_char = toupper(atom_label_in_bracket[0]);
        if (strlen(atom_label_in_bracket) == 1 || !islower(atom_label_in_bracket[1])) {
            switch (first_char) {
                case 'H': return 1;
                case 'B': return 5;
                case 'C': return 6;
                case 'N': return 7;
                case 'O': return 8;
                case 'F': return 9;
                case 'P': return 15;
                case 'S': return 16;
                case 'I': return 53;
                // Add more common elements as needed
                default: return 0; // Unknown
            }
        }
    }

    switch (type) {
        case ATOM_C: return 6;
        case ATOM_N: return 7;
        case ATOM_O: return 8;
        case ATOM_S: return 16;
        case ATOM_P: return 15;
        case ATOM_F: return 9;
        case ATOM_CL: return 17;
        case ATOM_BR: return 35;
        case ATOM_I: return 53;
        case ATOM_B: return 5;
        case ATOM_SI: return 14;
        case ATOM_SE: return 34;
        case ATOM_AS: return 33;
        default: return 0; // Unknown or ATOM_OTHER
    }
}

// Placeholder for Pauling electronegativity
double get_pauling_en(int atomic_number) {
    switch (atomic_number) {
        case 1: return 2.20; // H
        case 5: return 2.04; // B
        case 6: return 2.55; // C
        case 7: return 3.04; // N
        case 8: return 3.44; // O
        case 9: return 3.98; // F
        case 14: return 1.90; // Si
        case 15: return 2.19; // P
        case 16: return 2.58; // S
        case 17: return 3.16; // Cl
        case 33: return 2.18; // As
        case 34: return 2.55; // Se
        case 35: return 2.96; // Br
        case 53: return 2.66; // I
        default: return 0.0; // Unknown
    }
}

// Placeholder for atomic weight
double get_atomic_weight(int atomic_number) {
    switch (atomic_number) {
        case 1: return 1.008;   // H
        case 5: return 10.81;   // B
        case 6: return 12.011;  // C
        case 7: return 14.007;  // N
        case 8: return 15.999;  // O
        case 9: return 18.998;  // F
        case 14: return 28.085; // Si
        case 15: return 30.974; // P
        case 16: return 32.06;  // S
        case 17: return 35.45;  // Cl
        case 33: return 74.922; // As
        case 34: return 78.971; // Se
        case 35: return 79.904; // Br
        case 53: return 126.90; // I
        default: return 0.0;    // Unknown
    }
}

// Helper functions for SMILES parsing
AtomType determine_atom_type(const char* smiles, int pos, char* bracket_content, size_t bracket_content_size) {
    if (bracket_content) bracket_content[0] = '\0';

    if (smiles[pos] == '[') {
        if (bracket_content) {
            int k = 0;
            for (int j = pos + 1; smiles[j] != '\0' && smiles[j] != ']' && k < bracket_content_size - 1; ++j) {
                 // Extract only the element symbol part, stop at charge, isotope, H count etc.
                if (isalpha(smiles[j])) {
                    bracket_content[k++] = smiles[j];
                } else if (k > 0 && (smiles[j] == '+' || smiles[j] == '-' || isdigit(smiles[j]) || smiles[j] == 'H')) {
                    break; // Stop if we hit charge, isotope, or H count after an element symbol
                } else if (k == 0 && (smiles[j] == '+' || smiles[j] == '-')) {
                    // Allow leading charge for things like [Fe++]
                } else if (k==0 && isdigit(smiles[j])) {
                    // Isotope before element symbol
                }
                else {
                    // Other characters inside bracket - for now, we are interested in element symbol
                }
            }
            bracket_content[k] = '\0';
        }
        // More robust bracket parsing would be needed here for all features
        if (bracket_content && bracket_content[0] != '\0') {
            if (strcmp(bracket_content, "C") == 0) return ATOM_C;
            if (strcmp(bracket_content, "N") == 0) return ATOM_N;
            if (strcmp(bracket_content, "O") == 0) return ATOM_O;
            // ... add more from bracket_content
            if (strncmp(bracket_content, "Cl", 2) == 0) return ATOM_CL;
            if (strncmp(bracket_content, "Br", 2) == 0) return ATOM_BR;
            if (strncmp(bracket_content, "Si", 2) == 0) return ATOM_SI;
            // ...
        }
        return ATOM_UNKNOWN; // Or parse more specifically
    }
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
    } else if (smiles[pos] == 'B') {
        if (smiles[pos+1] == 'r') return ATOM_BR;
        return ATOM_B;
    } else if (smiles[pos] == 'I') {
        return ATOM_I;
    } else if (smiles[pos] == 's') { // aromatic sulfur or silicon?
        if (smiles[pos+1] == 'i' || smiles[pos+1] == 'I' ) return ATOM_SI; // crude check for [si]
        return ATOM_S; // assume aromatic sulfur
    } else if (smiles[pos] == 'A') { // As or Ar?
        if (smiles[pos+1] == 's') return ATOM_AS;
    } else if (smiles[pos] == 'S' && smiles[pos+1] == 'e') {
        return ATOM_SE;
    }
     // aromatic lowercase
    if (islower(smiles[pos])) {
        if (smiles[pos] == 'c') return ATOM_C;
        if (smiles[pos] == 'n') return ATOM_N;
        if (smiles[pos] == 'o') return ATOM_O;
        if (smiles[pos] == 'p') return ATOM_P;
        if (smiles[pos] == 's') return ATOM_S; // aromatic sulfur
    }
    return ATOM_OTHER;
}

// Parse a SMILES string into a molecular graph
MolecularGraph* parse_smiles(const char* smiles) {
    if (!smiles) return NULL;
    
    size_t len = strlen(smiles); 
    
    MolecularGraph* graph = (MolecularGraph*)calloc(1, sizeof(MolecularGraph));
    if (!graph) { 
        perror("Failed to allocate graph struct"); 
        return NULL; 
    }

    if (len == 0) {
        return graph; // atom_count is 0, atoms is NULL
    }
    
    // Dynamic allocation for atoms
    int allocated_atoms_count = 0;
    graph->atoms = NULL;


    // Ring closure tracking: (Using a simple array, might need improvement for complex cases)
    // Max 100 distinct ring closure points (0-99 for digits, plus more for %XX)
    // Each entry stores the atom_idx for the first time a ring number is seen.
    // -1 indicates not seen yet.
    int ring_closure_points[100]; 
    BondType ring_closure_bond_types[100];
    for(int k=0; k<100; ++k) {
        ring_closure_points[k] = -1;
        ring_closure_bond_types[k] = BOND_SINGLE; // Default bond type for rings
    }

    int atom_idx = 0;
    int prev_atom_chain = -1; // Index of the previously processed atom in the current chain segment
    
    // Stack for branches
    int branch_stack[100]; // Assuming max 100 nested branches
    int branch_stack_top = -1;

    BondType next_bond_type = BOND_SINGLE; // Default, can be changed by explicit bond chars

    char current_atom_label_buffer[16]; // Buffer for atom labels from brackets

    size_t main_loop_idx = 0;
    while(main_loop_idx < len) {
        char current_char = smiles[main_loop_idx];
        Atom* current_atom_ptr = NULL;
        int atom_chars_to_consume = 1; // How many chars from smiles string this atom takes

        if (isspace(current_char)) {
            main_loop_idx++;
            continue;
        }

        bool is_atom_char = (current_char == '[' || isalpha(current_char) || current_char == '*');

        if (is_atom_char) {
            if (atom_idx >= allocated_atoms_count) {
                int new_allocated_count = (allocated_atoms_count == 0) ? INITIAL_ATOM_BUFFER_SIZE : allocated_atoms_count * 2;
                Atom* temp_atoms = (Atom*)realloc(graph->atoms, new_allocated_count * sizeof(Atom));
                if (!temp_atoms) {
                    perror("Failed to realloc graph atoms");
                    if (graph->atoms) { // Free previously allocated atoms' bond arrays
                        for (int k_free = 0; k_free < atom_idx; ++k_free) {
                            if (graph->atoms[k_free].bonded_atoms) free(graph->atoms[k_free].bonded_atoms);
                            if (graph->atoms[k_free].bond_types) free(graph->atoms[k_free].bond_types);
                        }
                        free(graph->atoms);
                    }
                    free(graph);
                    return NULL;
                }
                graph->atoms = temp_atoms;
                memset(graph->atoms + allocated_atoms_count, 0, (new_allocated_count - allocated_atoms_count) * sizeof(Atom));
                allocated_atoms_count = new_allocated_count;
            }
            
            current_atom_ptr = &graph->atoms[atom_idx];
            current_atom_ptr->index = atom_idx;
            current_atom_ptr->bonds_count = 0; // Initialize bonds_count

            current_atom_ptr->bonded_atoms = (int*)calloc(MAX_BONDS_PER_ATOM_ESTIMATE, sizeof(int));
            current_atom_ptr->bond_types = (int*)calloc(MAX_BONDS_PER_ATOM_ESTIMATE, sizeof(int)); // <-- Fix: use int* not BondType*

            if (!current_atom_ptr->bonded_atoms || !current_atom_ptr->bond_types) {
                perror("Failed to allocate bond arrays for atom");
                 for(int k_free=0; k_free < atom_idx; ++k_free) { 
                    if(graph->atoms[k_free].bonded_atoms) free(graph->atoms[k_free].bonded_atoms);
                    if(graph->atoms[k_free].bond_types) free(graph->atoms[k_free].bond_types);
                 }
                if (current_atom_ptr->bonded_atoms) free(current_atom_ptr->bonded_atoms); // Free partially allocated for current
                // graph->atoms might point to valid memory from realloc, or new memory
                if (graph->atoms) free(graph->atoms);
                free(graph);
                return NULL;
            }
            
            memset(current_atom_label_buffer, 0, sizeof(current_atom_label_buffer));
            if (current_char == '[') {
                current_atom_ptr->is_in_bracket = true;
                size_t end_bracket_idx = main_loop_idx + 1;
                size_t label_idx = 0;
                bool first_char_of_symbol = true;
                // Basic extraction of element symbol from bracket:
                // e.g., [CH3], [N+], [O-], [SiH4], [Fe++]
                // This needs to be more robust to handle isotopes, explicit H, charge, class, etc.
                for (size_t k_bracket = main_loop_idx + 1; k_bracket < len && smiles[k_bracket] != ']'; ++k_bracket) {
                    if (isalpha(smiles[k_bracket])) {
                        if (label_idx < sizeof(current_atom_label_buffer) - 1) {
                            current_atom_label_buffer[label_idx++] = smiles[k_bracket];
                        }
                    } else if (label_idx > 0 && (strchr("+-0123456789H@", smiles[k_bracket]))) {
                        // Heuristic: stop collecting element symbol if we hit charge, H count, isotope etc.
                        // Or if it's part of stereochemistry info
                    } else if (label_idx == 0 && isdigit(smiles[k_bracket])) {
                        // Isotope number before element symbol
                    }
                     // For now, we just try to get the element symbol for determine_atom_type
                    end_bracket_idx = k_bracket +1;
                }
                 atom_chars_to_consume = (end_bracket_idx < len && smiles[end_bracket_idx] == ']') ? (end_bracket_idx - main_loop_idx + 1) : (len - main_loop_idx) ;


                // Full bracket parsing logic for charge, explicit_hydrogens, aromaticity (lowercase in bracket) should go here
                // For now, use determine_atom_type with the extracted label
                current_atom_ptr->type = determine_atom_type(smiles, main_loop_idx, current_atom_label_buffer, sizeof(current_atom_label_buffer));
                if (current_atom_label_buffer[0] != '\0' && islower(current_atom_label_buffer[0])) { // e.g. [se]
                    bool all_lower = true;
                    for(size_t l_idx =0; current_atom_label_buffer[l_idx] != '\0'; ++l_idx) {
                        if (!islower(current_atom_label_buffer[l_idx])) {
                            all_lower = false; break;
                        }
                    }
                    if(all_lower) current_atom_ptr->aromatic = 1;
                }


            } else { // Single char or two-char element
                current_atom_ptr->type = determine_atom_type(smiles, main_loop_idx, NULL, 0);
                current_atom_ptr->aromatic = islower(current_char);
                if (isupper(current_char) && (main_loop_idx + 1 < len) && islower(smiles[main_loop_idx + 1])) {
                    char two_char_symbol_check[3] = {current_char, smiles[main_loop_idx+1], '\0'};
                    if (strcmp(two_char_symbol_check, "Cl") == 0 || strcmp(two_char_symbol_check, "Br") == 0 ||
                        strcmp(two_char_symbol_check, "Si") == 0 || strcmp(two_char_symbol_check, "Se") == 0 ||
                        strcmp(two_char_symbol_check, "As") == 0) {
                        atom_chars_to_consume = 2;
                    }
                }
            }
            current_atom_ptr->atomic_number = get_atomic_number_simple(current_atom_ptr->type, current_atom_label_buffer[0] != '\0' ? current_atom_label_buffer : NULL);
            current_atom_ptr->electronegativity = get_pauling_en(current_atom_ptr->atomic_number);
            current_atom_ptr->atomic_weight = get_atomic_weight(current_atom_ptr->atomic_number);
            // TODO: Populate group_number, valence_electrons, hybridization, formal_charge, explicit_hydrogens more accurately


            // Add bond to previous atom in chain (if any)
            if (prev_atom_chain != -1) {
                Atom* prev_atom_ptr = &graph->atoms[prev_atom_chain];
                if (prev_atom_ptr->bonds_count < MAX_BONDS_PER_ATOM_ESTIMATE && current_atom_ptr->bonds_count < MAX_BONDS_PER_ATOM_ESTIMATE) {
                    prev_atom_ptr->bonded_atoms[prev_atom_ptr->bonds_count] = current_atom_ptr->index;
                    prev_atom_ptr->bond_types[prev_atom_ptr->bonds_count] = next_bond_type;
                    prev_atom_ptr->bonds_count++;

                    current_atom_ptr->bonded_atoms[current_atom_ptr->bonds_count] = prev_atom_ptr->index;
                    current_atom_ptr->bond_types[current_atom_ptr->bonds_count] = next_bond_type;
                    current_atom_ptr->bonds_count++;
                } else {
                    // fprintf(stderr, "Warning: Max bonds per atom estimate (%d) reached for atom %d or %d. Bond ignored.\n", MAX_BONDS_PER_ATOM_ESTIMATE, prev_atom_ptr->index, current_atom_ptr->index);
                }
            }
            next_bond_type = BOND_SINGLE; // Reset for next bond, unless specified otherwise
            prev_atom_chain = atom_idx;
            atom_idx++;
            main_loop_idx += atom_chars_to_consume;

        } else if (current_char == '(') {
            if (branch_stack_top < 99) {
                branch_stack[++branch_stack_top] = prev_atom_chain;
            } else { /*fprintf(stderr, "Warning: Branch stack overflow.\n");*/ }
            main_loop_idx++;
        } else if (current_char == ')') {
            if (branch_stack_top >= 0) {
                prev_atom_chain = branch_stack[branch_stack_top--];
            } else { /*fprintf(stderr, "Warning: Branch stack underflow (unmatched ')' ).\n");*/ }
            main_loop_idx++;
        } else if (isdigit(current_char)) {
            int ring_num = current_char - '0';
            if (main_loop_idx + 1 < len && smiles[main_loop_idx+1] == '%') { // check for % before two digits
                // This case is unusual, usually % is first e.g. %10
            }


            if (ring_closure_points[ring_num] == -1) { // First time seeing this ring number
                ring_closure_points[ring_num] = prev_atom_chain; // Atom that this ring connects back to
                ring_closure_bond_types[ring_num] = next_bond_type;
            } else { // Closing the ring
                int first_atom_idx_in_ring = ring_closure_points[ring_num];
                BondType ring_bond_type = ring_closure_bond_types[ring_num]; // Use bond type from first point
                                                                            // or use current next_bond_type if specified before digit.
                                                                            // SMILES standard: bond type before second closure digit applies.
                if (next_bond_type != BOND_SINGLE) ring_bond_type = next_bond_type;


                if (first_atom_idx_in_ring != -1 && prev_atom_chain != -1 && first_atom_idx_in_ring < atom_idx && prev_atom_chain < atom_idx) {
                    Atom* atom1 = &graph->atoms[first_atom_idx_in_ring];
                    Atom* atom2 = &graph->atoms[prev_atom_chain];

                    // Check if bond already exists (can happen with complex ring closures or errors)
                    bool bond_exists = false;
                    for(int k_bond=0; k_bond < atom1->bonds_count; ++k_bond) {
                        if (atom1->bonded_atoms[k_bond] == atom2->index) { bond_exists = true; break;}
                    }

                    if (!bond_exists && atom1->bonds_count < MAX_BONDS_PER_ATOM_ESTIMATE && atom2->bonds_count < MAX_BONDS_PER_ATOM_ESTIMATE) {
                        atom1->bonded_atoms[atom1->bonds_count] = atom2->index;
                        atom1->bond_types[atom1->bonds_count] = ring_bond_type;
                        atom1->bonds_count++;

                        atom2->bonded_atoms[atom2->bonds_count] = atom1->index;
                        atom2->bond_types[atom2->bonds_count] = ring_bond_type;
                        atom2->bonds_count++;
                        
                        // Simple ring counting (counts pairs of closure digits processed)
                        graph->rings_count++; 
                    } else {
                        // fprintf(stderr, "Warning: Max bonds estimate reached or bond exists, ring closure %d ignored between %d and %d.\n", ring_num, atom1->index, atom2->index);
                    }
                }
                ring_closure_points[ring_num] = -1; // Reset for potential reuse of this ring number (though not standard SMILES for same number simultaneously)
                ring_closure_bond_types[ring_num] = BOND_SINGLE;
            }
            next_bond_type = BOND_SINGLE; // Bond type used up by ring or was default
            main_loop_idx++;
        } else if (current_char == '%') { // Ring closures > 9, e.g., %10, %11
             if (main_loop_idx + 2 < len && isdigit(smiles[main_loop_idx+1]) && isdigit(smiles[main_loop_idx+2])) {
                int ring_num = (smiles[main_loop_idx+1] - '0') * 10 + (smiles[main_loop_idx+2] - '0');
                if (ring_num < 100) { // Our array limit
                     if (ring_closure_points[ring_num] == -1) {
                        ring_closure_points[ring_num] = prev_atom_chain;
                        ring_closure_bond_types[ring_num] = next_bond_type;
                    } else {
                        // ... (similar ring closing logic as for single digit) ...
                        int first_atom_idx_in_ring = ring_closure_points[ring_num];
                        BondType ring_bond_type = (next_bond_type != BOND_SINGLE) ? next_bond_type : ring_closure_bond_types[ring_num];

                        if (first_atom_idx_in_ring != -1 && prev_atom_chain != -1 && first_atom_idx_in_ring < atom_idx && prev_atom_chain < atom_idx) {
                            Atom* atom1 = &graph->atoms[first_atom_idx_in_ring];
                            Atom* atom2 = &graph->atoms[prev_atom_chain];
                            // ... (add bond logic, checking MAX_BONDS_PER_ATOM_ESTIMATE) ...
                             bool bond_exists = false;
                            for(int k_bond=0; k_bond < atom1->bonds_count; ++k_bond) {
                                if (atom1->bonded_atoms[k_bond] == atom2->index) { bond_exists = true; break;}
                            }
                            if (!bond_exists && atom1->bonds_count < MAX_BONDS_PER_ATOM_ESTIMATE && atom2->bonds_count < MAX_BONDS_PER_ATOM_ESTIMATE) {
                                atom1->bonded_atoms[atom1->bonds_count] = atom2->index;
                                atom1->bond_types[atom1->bonds_count] = ring_bond_type;
                                atom1->bonds_count++;

                                atom2->bonded_atoms[atom2->bonds_count] = atom1->index;
                                atom2->bond_types[atom2->bonds_count] = ring_bond_type;
                                atom2->bonds_count++;
                                graph->rings_count++;
                            }
                        }
                        ring_closure_points[ring_num] = -1;
                        ring_closure_bond_types[ring_num] = BOND_SINGLE;
                    }
                }
                next_bond_type = BOND_SINGLE;
                main_loop_idx += 3;
             } else { main_loop_idx++; /* Malformed %xx */ }
        } else { // Bond characters or other
            switch (current_char) {
                case '-': next_bond_type = BOND_SINGLE; main_loop_idx++; break;
                case '=': next_bond_type = BOND_DOUBLE; main_loop_idx++; break;
                case '#': next_bond_type = BOND_TRIPLE; main_loop_idx++; break;
                case '$': next_bond_type = BOND_AROMATIC; main_loop_idx++; break; // Quadruple or Aromatic? Often Aromatic.
                case ':': next_bond_type = BOND_AROMATIC; main_loop_idx++; break;
                // case '/': case '\\': // Stereochemistry bonds, treat as single for now
                //    next_bond_type = BOND_SINGLE; main_loop_idx++; break;
                case '.': // Disconnected structures. Reset prev_atom_chain.
                    prev_atom_chain = -1;
                    next_bond_type = BOND_SINGLE;
                    main_loop_idx++;
                    break;
                default:
                    // fprintf(stderr, "Warning: Unknown SMILES character '%c' at position %zu. Skipping.\n", current_char, main_loop_idx);
                    main_loop_idx++; // Skip unknown characters
                    break;
            }
        }
    } // End of while loop over SMILES string

    graph->atom_count = atom_idx;

    // Optional: Shrink graph->atoms to actual size if significantly overallocated
    if (graph->atom_count > 0 && graph->atom_count < allocated_atoms_count) {
        Atom* final_atoms = (Atom*)realloc(graph->atoms, graph->atom_count * sizeof(Atom));
        if (final_atoms) { // If realloc succeeds
            graph->atoms = final_atoms;
            // allocated_atoms_count = graph->atom_count; // Update this if needed elsewhere
        } // If realloc fails, keep the larger buffer.
    } else if (graph->atom_count == 0 && graph->atoms != NULL) { 
        // This means space was allocated but no atoms were parsed.
        free(graph->atoms);
        graph->atoms = NULL;
    }
    // If len == 0 initially, graph->atoms is already NULL and atom_count is 0.

    return graph;
}

void free_molecular_graph(MolecularGraph* graph) {
    if (!graph) return;
    if (graph->atoms) {
        for (int i = 0; i < graph->atom_count; ++i) {
            // Check if bonded_atoms and bond_types were actually allocated for this atom
            // The memset in realloc should make them NULL if atom_idx didn't reach them
            // Or if calloc inside loop failed for some reason (though it would likely return NULL earlier)
            if (graph->atoms[i].bonded_atoms) {
                free(graph->atoms[i].bonded_atoms);
                graph->atoms[i].bonded_atoms = NULL; 
            }
            if (graph->atoms[i].bond_types) {
                free(graph->atoms[i].bond_types);
                graph->atoms[i].bond_types = NULL;
            }
        }
        free(graph->atoms);
        graph->atoms = NULL; 
    }
    if (graph->ring_sizes) { // ring_sizes parsing is not fully implemented above
        free(graph->ring_sizes);
        graph->ring_sizes = NULL; 
    }
    free(graph);
}

// Helper for BFS to compute distances from a source_atom_idx
void bfs_from_source(MolecularGraph* graph, int source_atom_idx, int* dist_row) {
    if (!graph || !dist_row || source_atom_idx < 0 || source_atom_idx >= graph->atom_count) {
        if (graph && dist_row && graph->atom_count > 0) { // Check graph->atom_count before using
            for (int i = 0; i < graph->atom_count; ++i) dist_row[i] = -1; 
        } else if (dist_row) { // If graph is NULL or atom_count is 0, but dist_row might be for a fixed size
            // This case is tricky, assume dist_row cannot be safely initialized if graph is invalid
        }
        return;
    }

    int n = graph->atom_count;
    for (int i = 0; i < n; ++i) {
        dist_row[i] = INT_MAX; 
    }

    if (n == 0) return; 
    if (source_atom_idx >= n) { // Should be caught by initial check, but safeguard
        for (int i = 0; i < n; ++i) dist_row[i] = -1;
        return;
    }

    dist_row[source_atom_idx] = 0;

    if (n > 100000) { // Safety check: if n is abnormally large, might indicate prior corruption
        fprintf(stderr, "Warning: BFS attempted with very large atom count: %d. Skipping BFS.\n", n);
        for (int i = 0; i < n; ++i) dist_row[i] = -1; // Mark as error
        return;
    }

    int* queue = (int*)malloc(n * sizeof(int));
    if (!queue) { 
        perror("Failed to allocate queue for BFS");
        for (int i = 0; i < n; ++i) dist_row[i] = -1;
        return;
    }
    int head = 0, tail = 0;
    queue[tail++] = source_atom_idx;

    while (head < tail) {
        int u = queue[head++];
        if (u < 0 || u >= n) {
            continue; 
        }
        Atom* atom_u = &graph->atoms[u];
        if (!atom_u->bonded_atoms) continue; // Safeguard if bonded_atoms is NULL

        for (int i = 0; i < atom_u->bonds_count; ++i) {
            int v = atom_u->bonded_atoms[i];
            if (v < 0 || v >= n) {
                continue;
            }
            if (dist_row[v] == INT_MAX) { 
                dist_row[v] = dist_row[u] + 1;
                if (tail < n) { 
                    queue[tail++] = v;
                } else {
                    fprintf(stderr, "BFS queue overflow. Graph size: %d, tail: %d. SMILES might be malformed.\n", n, tail);
                    // Free queue and return, as graph state is suspect or too large
                    free(queue);
                    // Mark remaining unvisited as error
                    for(int k_err=0; k_err<n; ++k_err) if(dist_row[k_err] == INT_MAX) dist_row[k_err] = -1;
                    return;
                }
            }
        }
    }
    free(queue);
}

// Compute all-pairs shortest paths using N BFS runs
int** compute_distance_matrix(MolecularGraph* graph) {
    if (!graph) return NULL; 

    int n = graph->atom_count;
    if (n == 0) return NULL;
    if (n > 100000) { // Safety check similar to BFS
        fprintf(stderr, "Warning: compute_distance_matrix with very large atom count: %d. Skipping.\n", n);
        return NULL;
    }


    int** dist = (int**)malloc(n * sizeof(int*));
    if (!dist) {
        perror("Failed to allocate distance matrix rows");
        return NULL;
    }
    for (int i = 0; i < n; ++i) {
        dist[i] = NULL;
    }

    for (int i = 0; i < n; i++) {
        dist[i] = (int*)malloc(n * sizeof(int));
        if (!dist[i]) {
            perror("Failed to allocate distance matrix columns");
            for (int k = 0; k < i; ++k) free(dist[k]); // Only free what was allocated
            free(dist);
            return NULL;
        }
        bfs_from_source(graph, i, dist[i]);
    }

    int max_path = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (dist[i][j] == INT_MAX || dist[i][j] == -1) { 
                dist[i][j] = -1; 
            } else if (dist[i][j] > max_path) {
                max_path = dist[i][j];
            }
        }
    }
    if (graph) { // Ensure graph is not NULL before assigning
        graph->max_path_length = max_path;
    }
    
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


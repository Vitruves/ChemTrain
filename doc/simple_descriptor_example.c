#include "../cregistry.h"
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>

// Include the shared structures from topological.c
#include "topological.c"

// Forward declarations for utility functions that we'll use
extern MolecularGraph* parse_smiles(const char* smiles);
extern void free_molecular_graph(MolecularGraph* graph);

/**
 * Calculate the total number of bonds in the molecule
 */
double calculateBondCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles || *smiles == '\0') return 0.0;
    
    MolecularGraph* graph = parse_smiles(smiles);
    if (!graph) return 0.0;
    
    // Count total bonds (each bond connects two atoms)
    int total_bonds = 0;
    for (int i = 0; i < graph->atom_count; i++) {
        total_bonds += graph->atoms[i].bonds_count;
    }
    
    // Divide by 2 since each bond is counted twice (once for each atom)
    double result = total_bonds / 2.0;
    
    free_molecular_graph(graph);
    return result;
}

/**
 * Calculate the average number of bonds per atom
 */
double calculateAverageBondCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles || *smiles == '\0') return 0.0;
    
    MolecularGraph* graph = parse_smiles(smiles);
    if (!graph) return 0.0;
    
    if (graph->atom_count == 0) {
        free_molecular_graph(graph);
        return 0.0;
    }
    
    // Count total bonds
    int total_bonds = 0;
    for (int i = 0; i < graph->atom_count; i++) {
        total_bonds += graph->atoms[i].bonds_count;
    }
    
    // Calculate average (bonds/2 because each bond is counted twice)
    double result = (double)(total_bonds / 2.0) / graph->atom_count;
    
    free_molecular_graph(graph);
    return result;
} 
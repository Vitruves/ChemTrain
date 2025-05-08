# C Descriptors Implementation Guide

This document outlines how to create and implement new molecular descriptors in C for the ChemTrain project.

## Overview

The ChemTrain system supports descriptors written in both C++ and C. C descriptors offer performance benefits for computationally intensive operations. This guide focuses on implementing new descriptors in C that will integrate with the existing registry system.

## File Structure

1. Create a new C file in the `src/descriptors/` directory with a descriptive name (e.g., `my_descriptors.c`).

## Implementation Requirements

### Required Includes

Your C descriptor file must include the central registry header:

```c
#include "../cregistry.h"
```

This header provides the necessary function signatures and type definitions.

### Function Signature

Each descriptor must be implemented as a function with the following signature:

```c
double calculate<DescriptorName>(const void* context, GetSmilesFunc getSmilesFunc)
```

Where:
- `<DescriptorName>` is the name of your descriptor (e.g., `MyNewDescriptor`)
- `context` is a pointer to the molecule context
- `getSmilesFunc` is a function pointer for getting the SMILES string

### Function Implementation Structure

A typical implementation follows this pattern:

```c
double calculateMyNewDescriptor(const void* context, GetSmilesFunc getSmilesFunc) {
    // Get the SMILES string for the molecule
    const char* smiles = getSmilesFunc(context);
    if (!smiles || *smiles == '\0') return 0.0;
    
    // Parse the SMILES string into a molecular graph
    MolecularGraph* graph = parse_smiles(smiles);
    if (!graph) return 0.0;
    
    // Perform descriptor calculation
    double result = 0.0;
    
    // ... your calculation logic here ...
    
    // Free memory and return result
    free_molecular_graph(graph);
    return result;
}
```

## Common Utilities

Several utility functions are available in the existing C files that you can leverage:

- `parse_smiles(const char* smiles)`: Parses a SMILES string into a molecular graph
- `free_molecular_graph(MolecularGraph* graph)`: Frees memory allocated for a molecular graph
- `compute_distance_matrix(MolecularGraph* graph)`: Computes all-pairs shortest paths
- `free_distance_matrix(int** dist_matrix, int atom_count)`: Frees a distance matrix

## Molecular Graph Structure

The `MolecularGraph` structure provides access to:

```c
typedef struct {
    int atom_count;
    Atom* atoms;
    int rings_count;
    int* ring_sizes;
    int max_path_length;
} MolecularGraph;
```

Each `Atom` has properties like:

```c
typedef struct {
    AtomType type;
    int index;
    int aromatic;
    int charge;
    int bonds_count;
    int* bonded_atoms;
    int* bond_types;
    int atomic_number;
    double electronegativity;
    double atomic_weight;
    // ... other properties ...
} Atom;
```

## Integration Process

After creating your C descriptor file:

1. **No need to manually modify the registry** - the `update_registry.py` script will automatically:
   - Detect your new descriptor functions
   - Register them in the C registry header
   - Create C++ wrapper classes

2. **Run the registry update script**:
   ```bash
   python3 scripts/update_registry.py
   ```

3. **Rebuild the project**:
   ```bash
   make
   ```

## Examples

### Simple Atom Count Descriptor

```c
double calculateAtomCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles || *smiles == '\0') return 0.0;
    
    MolecularGraph* graph = parse_smiles(smiles);
    if (!graph) return 0.0;
    
    double result = graph->atom_count;
    
    free_molecular_graph(graph);
    return result;
}
```

### Descriptor Using Distance Matrix

```c
double calculateMaxDistance(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles || *smiles == '\0') return 0.0;
    
    MolecularGraph* graph = parse_smiles(smiles);
    if (!graph) return 0.0;
    
    int** dist_matrix = compute_distance_matrix(graph);
    if (!dist_matrix) {
        free_molecular_graph(graph);
        return 0.0;
    }
    
    int max_dist = 0;
    for (int i = 0; i < graph->atom_count; i++) {
        for (int j = 0; j < graph->atom_count; j++) {
            if (dist_matrix[i][j] > max_dist && dist_matrix[i][j] != 9999) {
                max_dist = dist_matrix[i][j];
            }
        }
    }
    
    free_distance_matrix(dist_matrix, graph->atom_count);
    free_molecular_graph(graph);
    return (double)max_dist;
}
```

## Debugging Tips

1. Always check the return value of memory allocation calls
2. Ensure proper memory management (free everything you allocate)
3. Handle edge cases: empty SMILES strings, single atoms, disconnected structures
4. Use the verbose logging in the main application to debug issues

## Performance Considerations

1. Minimize string operations on the SMILES string
2. Cache intermediate results when applicable
3. Avoid unnecessary memory allocations
4. Consider using bitwise operations for efficiency where appropriate
5. For computationally intensive algorithms, consider adding early termination conditions

## Common Pitfalls

1. Using `const const char*` instead of just `const char*` (this causes compiler warnings)
2. Forgetting to free allocated memory
3. Not handling NULL pointers or empty strings
4. Accessing array elements out of bounds

By following this guide, you can create efficient C descriptor implementations that integrate seamlessly with the ChemTrain system. 
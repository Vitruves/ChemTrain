// Template for C descriptor implementation

#include "../cregistry.h"
#include <string.h>

// Example descriptor implementation
double calculateExampleDescriptor(const Context* context, GetSmilesFunc getSmilesFunc) {
    // Get SMILES string using the provided function
    const char* smiles = getSmilesFunc(context);
    
    // Process SMILES and calculate descriptor value
    double result = 0.0;
    
    // Simple example: return length of SMILES string
    if (smiles) {
        result = strlen(smiles);
    }
    
    return result;
}

// Manually created C registry header
#ifndef CREGISTRY_H
#define CREGISTRY_H

#ifdef __cplusplus
extern "C" {
#endif

// Forward declaration for the Context interface
struct Context;
typedef struct Context Context;

// Function pointer type for getting SMILES string
typedef const char* (*GetSmilesFunc)(const Context*);

// Helper function that will be provided by the C++ implementation
const char* getSmilesCFunc(const Context* ctx);

// Extended Topological Atom Pairs
double calculateETAP(const Context* context, GetSmilesFunc getSmilesFunc);

// Molecular Mesh Descriptors
double calculateMolecularMesh(const void* context, GetSmilesFunc getSmilesFunc);

// Maximum Common Subgraph Size Index
double calculateMCSIndex(const void* context, GetSmilesFunc getSmilesFunc);

// Vertex-Edge Path Descriptors
double calculateVertexEdgePath(const void* context, GetSmilesFunc getSmilesFunc);

// McFarland Complexity Index
double calculateMcFarlandComplexity(const void* context, GetSmilesFunc getSmilesFunc);

// Cycle Connectivity Index
double calculateCycleConnectivity(const void* context, GetSmilesFunc getSmilesFunc);

// Fragment Complexity Score
double calculateFragmentComplexity(const void* context, GetSmilesFunc getSmilesFunc);

// Topological Polar Buried Surface Area
double calculateTopologicalPolarBSA(const void* context, GetSmilesFunc getSmilesFunc);

// Eccentric Distance Index
double calculateEccentricDistance(const void* context, GetSmilesFunc getSmilesFunc);

// Molecular Scaffolding Index
double calculateMolecularScaffolding(const void* context, GetSmilesFunc getSmilesFunc);

#ifdef __cplusplus
}
#endif

#endif // CREGISTRY_H
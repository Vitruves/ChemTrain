// src/cregistry.h
#ifndef CREGISTRY_H
#define CREGISTRY_H

#ifdef __cplusplus
extern "C" {
#endif

// Forward declaration for the Context interface to match C++ code
typedef struct Context Context;
typedef const char* (*GetSmilesFunc)(const Context*);

// Extended Topological Atom Pairs (ETAP)
double calculateETAP(const Context* context, GetSmilesFunc getSmilesFunc);

// Molecular Mesh Descriptors
double calculateMolecularMesh(const Context* context, GetSmilesFunc getSmilesFunc);

// Maximum Common Subgraph Size Index
double calculateMCSIndex(const Context* context, GetSmilesFunc getSmilesFunc);

// Vertex-Edge Path Descriptors 
double calculateVertexEdgePath(const Context* context, GetSmilesFunc getSmilesFunc);

// McFarland Complexity Index
double calculateMcFarlandComplexity(const Context* context, GetSmilesFunc getSmilesFunc);

// Cycle Connectivity Index
double calculateCycleConnectivity(const Context* context, GetSmilesFunc getSmilesFunc);

// Fragment Complexity Score
double calculateFragmentComplexity(const Context* context, GetSmilesFunc getSmilesFunc);

// Topological Polar Buried Surface Area
double calculateTopologicalPolarBSA(const Context* context, GetSmilesFunc getSmilesFunc);

// Eccentric Distance Index
double calculateEccentricDistance(const Context* context, GetSmilesFunc getSmilesFunc);

// Molecular Scaffolding Index
double calculateMolecularScaffolding(const Context* context, GetSmilesFunc getSmilesFunc);

#ifdef __cplusplus
}
#endif

#endif // CREGISTRY_H
// Manual implementation of C descriptor wrappers

#include "../common.hpp"
#include "../cregistry.h"
#include <string>

// Thread-local storage for SMILES string to maintain lifetime
static thread_local std::string tls_smiles_buffer;

// Implementation of the C function for getting SMILES
extern "C" {
    const char* getSmilesCFunc(const Context* ctx) {
        // Use C-style cast to bypass the C++/C type system barrier
        const desfact::Context* cpp_ctx = (const desfact::Context*)(ctx);
        tls_smiles_buffer = cpp_ctx->getSmiles();
        return tls_smiles_buffer.c_str();
    }
}

namespace desfact {

// --- MolecularMesh Wrapper ---
DECLARE_DESCRIPTOR(MolecularMesh, C_Descriptor, "Encoding molecular topology as a meshed network")
DESCRIPTOR_DEPENDENCIES(MolecularMesh) { return {}; }
DescriptorResult MolecularMeshDescriptor::calculate(Context& context) const {
    // Cast C++ Context to C Context structure
    const Context* c_ctx = (const Context*)(&context);
    return ::calculateMolecularMesh(c_ctx, getSmilesCFunc);
}

// --- MCSIndex Wrapper ---
DECLARE_DESCRIPTOR(MCSIndex, C_Descriptor, "Maximum Common Subgraph Size Index - Normalized by molecular size")
DESCRIPTOR_DEPENDENCIES(MCSIndex) { return {}; }
DescriptorResult MCSIndexDescriptor::calculate(Context& context) const {
    return ::calculateMCSIndex((const void*)&context, getSmilesCFunc);
}

// --- VertexEdgePath Wrapper ---
DECLARE_DESCRIPTOR(VertexEdgePath, C_Descriptor, "Unique paths between vertices passing specific edges")
DESCRIPTOR_DEPENDENCIES(VertexEdgePath) { return {}; }
DescriptorResult VertexEdgePathDescriptor::calculate(Context& context) const {
    return ::calculateVertexEdgePath((const void*)&context, getSmilesCFunc);
}

// --- McFarlandComplexity Wrapper ---
DECLARE_DESCRIPTOR(McFarlandComplexity, C_Descriptor, "Topological complexity based on fragment types")
DESCRIPTOR_DEPENDENCIES(McFarlandComplexity) { return {}; }
DescriptorResult McFarlandComplexityDescriptor::calculate(Context& context) const {
    return ::calculateMcFarlandComplexity((const void*)&context, getSmilesCFunc);
}

// --- CycleConnectivity Wrapper ---
DECLARE_DESCRIPTOR(CycleConnectivity, C_Descriptor, "Connectivity patterns of cyclic structures")
DESCRIPTOR_DEPENDENCIES(CycleConnectivity) { return {}; }
DescriptorResult CycleConnectivityDescriptor::calculate(Context& context) const {
    return ::calculateCycleConnectivity((const void*)&context, getSmilesCFunc);
}

// --- FragmentComplexity Wrapper ---
DECLARE_DESCRIPTOR(FragmentComplexity, C_Descriptor, "Complexity based on unique fragments")
DESCRIPTOR_DEPENDENCIES(FragmentComplexity) { return {}; }
DescriptorResult FragmentComplexityDescriptor::calculate(Context& context) const {
    return ::calculateFragmentComplexity((const void*)&context, getSmilesCFunc);
}

// --- TopologicalPolarBSA Wrapper ---
DECLARE_DESCRIPTOR(TopologicalPolarBSA, C_Descriptor, "Surface area inaccessible to solvent")
DESCRIPTOR_DEPENDENCIES(TopologicalPolarBSA) { return {}; }
DescriptorResult TopologicalPolarBSADescriptor::calculate(Context& context) const {
    return ::calculateTopologicalPolarBSA((const void*)&context, getSmilesCFunc);
}

// --- EccentricDistance Wrapper ---
DECLARE_DESCRIPTOR(EccentricDistance, C_Descriptor, "Sum of eccentricity-weighted distances")
DESCRIPTOR_DEPENDENCIES(EccentricDistance) { return {}; }
DescriptorResult EccentricDistanceDescriptor::calculate(Context& context) const {
    return ::calculateEccentricDistance((const void*)&context, getSmilesCFunc);
}

// --- MolecularScaffolding Wrapper ---
DECLARE_DESCRIPTOR(MolecularScaffolding, C_Descriptor, "Core structure representation")
DESCRIPTOR_DEPENDENCIES(MolecularScaffolding) { return {}; }
DescriptorResult MolecularScaffoldingDescriptor::calculate(Context& context) const {
    return ::calculateMolecularScaffolding((const void*)&context, getSmilesCFunc);
}

// Registration functions
void register_MolecularMeshDescriptor() {
    auto descriptor = std::make_shared<MolecularMeshDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_MCSIndexDescriptor() {
    auto descriptor = std::make_shared<MCSIndexDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_VertexEdgePathDescriptor() {
    auto descriptor = std::make_shared<VertexEdgePathDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_McFarlandComplexityDescriptor() {
    auto descriptor = std::make_shared<McFarlandComplexityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_CycleConnectivityDescriptor() {
    auto descriptor = std::make_shared<CycleConnectivityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_FragmentComplexityDescriptor() {
    auto descriptor = std::make_shared<FragmentComplexityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_TopologicalPolarBSADescriptor() {
    auto descriptor = std::make_shared<TopologicalPolarBSADescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_EccentricDistanceDescriptor() {
    auto descriptor = std::make_shared<EccentricDistanceDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_MolecularScaffoldingDescriptor() {
    auto descriptor = std::make_shared<MolecularScaffoldingDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

} // namespace desfact

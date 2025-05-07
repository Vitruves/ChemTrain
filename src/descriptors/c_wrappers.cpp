// Manual implementation of C descriptor wrappers

#include "../common.hpp"
#include "../cregistry.h"
#include <string>

// Thread-local storage for SMILES string to maintain lifetime
static thread_local std::string tls_smiles_buffer;

// Implementation of the C function for getting SMILES
extern "C" {
    const char* getSmilesCFunc(const void* ctx_void) {
        // Cast void* back to the C++ Context type
        const desfact::Context* cpp_ctx = static_cast<const desfact::Context*>(ctx_void);
        if (!cpp_ctx) {
             // Handle null context if necessary, perhaps return an empty string or null
             return ""; 
        }
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

// --- NEW C Descriptor Wrappers for Bit Operations ---

// Macro to define C wrappers easily
#define DECLARE_C_DESCRIPTOR_WRAPPER(ClassName, C_CalcFunc, Category, Description) \
    DECLARE_DESCRIPTOR(ClassName, Category, Description) \
    DESCRIPTOR_DEPENDENCIES(ClassName) { return {}; } \
    DescriptorResult ClassName##Descriptor::calculate(Context& context) const { \
        /* Cast to const void* when passing to C function */ \
        const void* c_ctx_void = (const void*)(&context); \
        return ::C_CalcFunc(c_ctx_void, getSmilesCFunc); \
    } \
    void register_##ClassName##Descriptor() { \
        auto descriptor = std::make_shared<ClassName##Descriptor>(); \
        auto& registry = DescriptorRegistry::getInstance(); \
        registry.registerDescriptor(descriptor); \
    }

// Declare wrappers using the macro
DECLARE_C_DESCRIPTOR_WRAPPER(StringLength, calculateStringLength, C_Descriptor, "Raw length of the SMILES string")
DECLARE_C_DESCRIPTOR_WRAPPER(BitCount, calculateBitCount, C_Descriptor, "Total count of binary 1 bits in the SMILES string")
DECLARE_C_DESCRIPTOR_WRAPPER(BitDensity, calculateBitDensity, C_Descriptor, "Density of 1 bits per byte")
DECLARE_C_DESCRIPTOR_WRAPPER(AsciiSum, calculateAsciiSum, C_Descriptor, "Sum of ASCII values of all characters")
DECLARE_C_DESCRIPTOR_WRAPPER(LongestBitRun, calculateLongestBitRun, C_Descriptor, "Length of longest consecutive bit run (0s or 1s)")
DECLARE_C_DESCRIPTOR_WRAPPER(ByteXorChecksum, calculateByteXorChecksum, C_Descriptor, "XOR of all bytes in the string")
DECLARE_C_DESCRIPTOR_WRAPPER(ByteEntropy, calculateByteEntropy, C_Descriptor, "Shannon entropy of byte distribution")
DECLARE_C_DESCRIPTOR_WRAPPER(FirstByteValue, calculateFirstByteValue, C_Descriptor, "Value of the first byte in the string")
DECLARE_C_DESCRIPTOR_WRAPPER(LastByteValue, calculateLastByteValue, C_Descriptor, "Value of the last byte in the string")
DECLARE_C_DESCRIPTOR_WRAPPER(OddEvenBitRatio, calculateOddEvenBitRatio, C_Descriptor, "Ratio of 1 bits at odd positions to even positions")
DECLARE_C_DESCRIPTOR_WRAPPER(CharPalindromeCheck, calculateCharPalindromeCheck, C_Descriptor, "Boolean indicating if string is a character palindrome")
DECLARE_C_DESCRIPTOR_WRAPPER(BitTransitionCount, calculateBitTransitionCount, C_Descriptor, "Count of bit transitions (0→1, 1→0)")
DECLARE_C_DESCRIPTOR_WRAPPER(LeftRightBitBalance, calculateLeftRightBitBalance, C_Descriptor, "Difference between bit counts in left and right halves")
DECLARE_C_DESCRIPTOR_WRAPPER(HammingWeightVariance, calculateHammingWeightVariance, C_Descriptor, "Variance of population counts across bytes")
DECLARE_C_DESCRIPTOR_WRAPPER(FastHash, calculateFastHash, C_Descriptor, "FNV-1a hash value of the string")
DECLARE_C_DESCRIPTOR_WRAPPER(Crc8Checksum, calculateCrc8Checksum, C_Descriptor, "CRC8 checksum of the string")
DECLARE_C_DESCRIPTOR_WRAPPER(MostFrequentByte, calculateMostFrequentByte, C_Descriptor, "Value of the most frequently occurring byte")
DECLARE_C_DESCRIPTOR_WRAPPER(LeastFrequentByte, calculateLeastFrequentByte, C_Descriptor, "Value of the least frequently occurring byte (that appears at least once)")
DECLARE_C_DESCRIPTOR_WRAPPER(UniqueByteCount, calculateUniqueByteCount, C_Descriptor, "Count of unique byte values in the string")
DECLARE_C_DESCRIPTOR_WRAPPER(EvenIndexedBitCount, calculateEvenIndexedBitCount, C_Descriptor, "Count of 1 bits at even-indexed positions")
DECLARE_C_DESCRIPTOR_WRAPPER(OddIndexedBitCount, calculateOddIndexedBitCount, C_Descriptor, "Count of 1 bits at odd-indexed positions")
DECLARE_C_DESCRIPTOR_WRAPPER(ByteWiseMinimum, calculateByteWiseMinimum, C_Descriptor, "Minimum byte value in the string")
DECLARE_C_DESCRIPTOR_WRAPPER(ByteWiseMaximum, calculateByteWiseMaximum, C_Descriptor, "Maximum byte value in the string")
DECLARE_C_DESCRIPTOR_WRAPPER(ByteWiseMedian, calculateByteWiseMedian, C_Descriptor, "Median byte value in the string")
DECLARE_C_DESCRIPTOR_WRAPPER(HighLowBitRatio, calculateHighLowBitRatio, C_Descriptor, "Ratio of high bits (values 1) to low bits (values 0)")
DECLARE_C_DESCRIPTOR_WRAPPER(CharSequenceComplexity, calculateCharSequenceComplexity, C_Descriptor, "Ratio of unique characters to total length")
DECLARE_C_DESCRIPTOR_WRAPPER(AvgBitPopulationCount, calculateAvgBitPopulationCount, C_Descriptor, "Average number of set bits per byte")
DECLARE_C_DESCRIPTOR_WRAPPER(AvgLeadingZeroCount, calculateAvgLeadingZeroCount, C_Descriptor, "Average number of leading zeros per byte")
DECLARE_C_DESCRIPTOR_WRAPPER(AvgTrailingZeroCount, calculateAvgTrailingZeroCount, C_Descriptor, "Average number of trailing zeros per byte")
DECLARE_C_DESCRIPTOR_WRAPPER(PopcountVariance, calculatePopcountVariance, C_Descriptor, "Variance of the bit count across bytes")
DECLARE_C_DESCRIPTOR_WRAPPER(ByteAndReduction, calculateByteAndReduction, C_Descriptor, "Result of ANDing all bytes together")
DECLARE_C_DESCRIPTOR_WRAPPER(ByteOrReduction, calculateByteOrReduction, C_Descriptor, "Result of ORing all bytes together")
DECLARE_C_DESCRIPTOR_WRAPPER(ByteXorReduction, calculateByteXorReduction, C_Descriptor, "Result of XORing all bytes together")
DECLARE_C_DESCRIPTOR_WRAPPER(NonZeroByteDensity, calculateNonZeroByteDensity, C_Descriptor, "Fraction of bytes that are non-zero")
DECLARE_C_DESCRIPTOR_WRAPPER(AsciiProductMod256, calculateAsciiProductMod256, C_Descriptor, "Product of all ASCII values, modulo 256")
DECLARE_C_DESCRIPTOR_WRAPPER(BitAlternationFreq, calculateBitAlternationFreq, C_Descriptor, "Frequency of bit flips (transitions) per bit")
DECLARE_C_DESCRIPTOR_WRAPPER(CharDistributionVariance, calculateCharDistributionVariance, C_Descriptor, "Variance of character frequencies")
DECLARE_C_DESCRIPTOR_WRAPPER(LowHighNibbleRatio, calculateLowHighNibbleRatio, C_Descriptor, "Ratio of frequency of low nibbles (0-7) to high nibbles (8-15)")
DECLARE_C_DESCRIPTOR_WRAPPER(CharPosWeightedSum, calculateCharPosWeightedSum, C_Descriptor, "Sum of byte values weighted by their position")
DECLARE_C_DESCRIPTOR_WRAPPER(NibbleEntropy, calculateNibbleEntropy, C_Descriptor, "Shannon entropy of the nibble (4-bit) distribution")
DECLARE_C_DESCRIPTOR_WRAPPER(CharFrequencyEntropy, calculateCharFrequencyEntropy, C_Descriptor, "Shannon entropy of character frequencies")
DECLARE_C_DESCRIPTOR_WRAPPER(AsciiRangeCoverage, calculateAsciiRangeCoverage, C_Descriptor, "Ratio of ASCII range covered by the string")
DECLARE_C_DESCRIPTOR_WRAPPER(NibbleHistogramBalance, calculateNibbleHistogramBalance, C_Descriptor, "Balance of nibble frequencies (higher = more uniform)")
DECLARE_C_DESCRIPTOR_WRAPPER(BinaryCompressionRatio, calculateBinaryCompressionRatio, C_Descriptor, "Estimate of compressibility based on bit patterns")
DECLARE_C_DESCRIPTOR_WRAPPER(BitReversalDistance, calculateBitReversalDistance, C_Descriptor, "Hamming distance between bit pattern and its reversal")
DECLARE_C_DESCRIPTOR_WRAPPER(ByteReversalDistance, calculateByteReversalDistance, C_Descriptor, "Hamming distance between byte sequence and its reversal")
DECLARE_C_DESCRIPTOR_WRAPPER(MurmurHash3, calculateMurmurHash3, C_Descriptor, "MurmurHash3 hash value of the string")
DECLARE_C_DESCRIPTOR_WRAPPER(CharRunLength, calculateCharRunLength, C_Descriptor, "Average length of runs of identical characters")
DECLARE_C_DESCRIPTOR_WRAPPER(ByteRunEntropy, calculateByteRunEntropy, C_Descriptor, "Entropy of run lengths of bytes")
DECLARE_C_DESCRIPTOR_WRAPPER(SingleBitPatternMatch, calculateSingleBitPatternMatch, C_Descriptor, "Frequency of single bit patterns ('1' followed by '0')")
DECLARE_C_DESCRIPTOR_WRAPPER(BytePairTransitions, calculateBytePairTransitions, C_Descriptor, "Shannon entropy of byte pair transitions")
DECLARE_C_DESCRIPTOR_WRAPPER(BitDifferenceSum, calculateBitDifferenceSum, C_Descriptor, "Sum of bit differences between consecutive bytes")
DECLARE_C_DESCRIPTOR_WRAPPER(FirstLastByteXor, calculateFirstLastByteXor, C_Descriptor, "XOR of first and last bytes")
DECLARE_C_DESCRIPTOR_WRAPPER(ByteMarkovProperty, calculateByteMarkovProperty, C_Descriptor, "Measure of how well bytes predict subsequent bytes")
DECLARE_C_DESCRIPTOR_WRAPPER(SwarPopcount, calculateSwarPopcount, C_Descriptor, "SIMD-within-a-register population count")
DECLARE_C_DESCRIPTOR_WRAPPER(FastBitHistogram, calculateFastBitHistogram, C_Descriptor, "Histogram uniformity of bit patterns")
DECLARE_C_DESCRIPTOR_WRAPPER(BitDispersionIndex, calculateBitDispersionIndex, C_Descriptor, "Measure of how evenly distributed the bits are")
DECLARE_C_DESCRIPTOR_WRAPPER(HalfByteRotationHash, calculateHalfByteRotationHash, C_Descriptor, "Hash value based on rotating half-bytes")
DECLARE_C_DESCRIPTOR_WRAPPER(AvalanchePatternSensitivity, calculateAvalanchePatternSensitivity, C_Descriptor, "Measure of bit pattern change propagation")
DECLARE_C_DESCRIPTOR_WRAPPER(SequentialXorDiff, calculateSequentialXorDiff, C_Descriptor, "Sum of sequential XOR differences between bytes")
DECLARE_C_DESCRIPTOR_WRAPPER(ByteDeltaEncodingSize, calculateByteDeltaEncodingSize, C_Descriptor, "Estimate of delta encoding efficiency")
DECLARE_C_DESCRIPTOR_WRAPPER(BitCorrelationCoef, calculateBitCorrelationCoef, C_Descriptor, "Correlation coefficient between even and odd bits")
DECLARE_C_DESCRIPTOR_WRAPPER(BitRunLengthSize, calculateBitRunLengthSize, C_Descriptor, "Estimate of bit-level run length encoding efficiency")
DECLARE_C_DESCRIPTOR_WRAPPER(ByteMirroringDistance, calculateByteMirroringDistance, C_Descriptor, "Hamming distance of bytes from their bit-reversed versions")
DECLARE_C_DESCRIPTOR_WRAPPER(RollingByteChecksum, calculateRollingByteChecksum, C_Descriptor, "Adler-32 type rolling checksum")
DECLARE_C_DESCRIPTOR_WRAPPER(ZeroBitRatio, calculateZeroBitRatio, C_Descriptor, "Ratio of 0 bits to total bits")
DECLARE_C_DESCRIPTOR_WRAPPER(BitAutocorrelation, calculateBitAutocorrelation, C_Descriptor, "Autocorrelation of bit sequence with lag 1")
DECLARE_C_DESCRIPTOR_WRAPPER(DiagonalBitAccumulation, calculateDiagonalBitAccumulation, C_Descriptor, "Sum of bits where index equals value")
DECLARE_C_DESCRIPTOR_WRAPPER(BitwiseMajorityFunction, calculateBitwiseMajorityFunction, C_Descriptor, "Majority vote from adjacent bits")
DECLARE_C_DESCRIPTOR_WRAPPER(LzcntPattern, calculateLzcntPattern, C_Descriptor, "Pattern of leading zero counts")
DECLARE_C_DESCRIPTOR_WRAPPER(TzcntPattern, calculateTzcntPattern, C_Descriptor, "Pattern of trailing zero counts")
DECLARE_C_DESCRIPTOR_WRAPPER(ByteGeometricMean, calculateByteGeometricMean, C_Descriptor, "Geometric mean of byte values")

// --- NEW Advanced Topological/Physicochemical Descriptor Wrappers ---
DECLARE_C_DESCRIPTOR_WRAPPER(BracketComplexity, calculateBracketComplexity, C_AdvancedTopo, "Bracket Complexity")
DECLARE_C_DESCRIPTOR_WRAPPER(HeteroatomDistribution, calculateHeteroatomDistribution, C_AdvancedTopo, "Heteroatom Distribution")
DECLARE_C_DESCRIPTOR_WRAPPER(CharacterTransitionEntropy, calculateCharacterTransitionEntropy, C_AdvancedTopo, "Character Transition Entropy")
DECLARE_C_DESCRIPTOR_WRAPPER(SyntacticNestingLevel, calculateSyntacticNestingLevel, C_AdvancedTopo, "Syntactic Nesting Level")
DECLARE_C_DESCRIPTOR_WRAPPER(HalogenIndex, calculateHalogenIndex, C_AdvancedTopo, "Halogen Index")
DECLARE_C_DESCRIPTOR_WRAPPER(AromaticCharacterRatio, calculateAromaticCharacterRatio, C_AdvancedTopo, "Aromatic Character Ratio")
DECLARE_C_DESCRIPTOR_WRAPPER(BranchPointDensity, calculateBranchPointDensity, C_AdvancedTopo, "Branch Point Density")
DECLARE_C_DESCRIPTOR_WRAPPER(ElementDiversityMeasure, calculateElementDiversityMeasure, C_AdvancedTopo, "Element Diversity Measure")

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

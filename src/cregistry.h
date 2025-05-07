// Manually created C registry header
#ifndef CREGISTRY_H
#define CREGISTRY_H

#ifdef __cplusplus
extern "C" {
#endif

// Remove the forward declaration for the C struct Context, as we'll use void*
// struct Context;
// typedef struct Context Context;

// Function pointer type for getting SMILES string using void* context
typedef const char* (*GetSmilesFunc)(const void*);

// Helper function that will be provided by the C++ implementation
// Note: The implementation in c_wrappers.cpp needs adjustment if it relies on Context*
const char* getSmilesCFunc(const void* ctx); // Adjusted to use void*

// Extended Topological Atom Pairs
double calculateETAP(const void* context, GetSmilesFunc getSmilesFunc);

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

// --- NEW C Descriptor Declarations for Bit Operations ---
double calculateStringLength(const void* context, GetSmilesFunc getSmilesFunc);
double calculateBitCount(const void* context, GetSmilesFunc getSmilesFunc);
double calculateBitDensity(const void* context, GetSmilesFunc getSmilesFunc);
double calculateAsciiSum(const void* context, GetSmilesFunc getSmilesFunc);
double calculateLongestBitRun(const void* context, GetSmilesFunc getSmilesFunc);
double calculateByteXorChecksum(const void* context, GetSmilesFunc getSmilesFunc);
double calculateByteEntropy(const void* context, GetSmilesFunc getSmilesFunc);
double calculateFirstByteValue(const void* context, GetSmilesFunc getSmilesFunc);
double calculateLastByteValue(const void* context, GetSmilesFunc getSmilesFunc);
double calculateOddEvenBitRatio(const void* context, GetSmilesFunc getSmilesFunc);
double calculateCharPalindromeCheck(const void* context, GetSmilesFunc getSmilesFunc);
double calculateBitTransitionCount(const void* context, GetSmilesFunc getSmilesFunc);
double calculateLeftRightBitBalance(const void* context, GetSmilesFunc getSmilesFunc);
double calculateHammingWeightVariance(const void* context, GetSmilesFunc getSmilesFunc);
double calculateFastHash(const void* context, GetSmilesFunc getSmilesFunc);
double calculateCrc8Checksum(const void* context, GetSmilesFunc getSmilesFunc);
double calculateMostFrequentByte(const void* context, GetSmilesFunc getSmilesFunc);
double calculateLeastFrequentByte(const void* context, GetSmilesFunc getSmilesFunc);
double calculateUniqueByteCount(const void* context, GetSmilesFunc getSmilesFunc);
double calculateEvenIndexedBitCount(const void* context, GetSmilesFunc getSmilesFunc);
double calculateOddIndexedBitCount(const void* context, GetSmilesFunc getSmilesFunc);
double calculateByteWiseMinimum(const void* context, GetSmilesFunc getSmilesFunc);
double calculateByteWiseMaximum(const void* context, GetSmilesFunc getSmilesFunc);
double calculateByteWiseMedian(const void* context, GetSmilesFunc getSmilesFunc);
double calculateHighLowBitRatio(const void* context, GetSmilesFunc getSmilesFunc);
double calculateCharSequenceComplexity(const void* context, GetSmilesFunc getSmilesFunc);
double calculateAvgBitPopulationCount(const void* context, GetSmilesFunc getSmilesFunc);
double calculateAvgLeadingZeroCount(const void* context, GetSmilesFunc getSmilesFunc);
double calculateAvgTrailingZeroCount(const void* context, GetSmilesFunc getSmilesFunc);
double calculatePopcountVariance(const void* context, GetSmilesFunc getSmilesFunc);
double calculateByteAndReduction(const void* context, GetSmilesFunc getSmilesFunc);
double calculateByteOrReduction(const void* context, GetSmilesFunc getSmilesFunc);
double calculateByteXorReduction(const void* context, GetSmilesFunc getSmilesFunc);
double calculateNonZeroByteDensity(const void* context, GetSmilesFunc getSmilesFunc);
double calculateAsciiProductMod256(const void* context, GetSmilesFunc getSmilesFunc);
double calculateBitAlternationFreq(const void* context, GetSmilesFunc getSmilesFunc);
double calculateCharDistributionVariance(const void* context, GetSmilesFunc getSmilesFunc);
double calculateLowHighNibbleRatio(const void* context, GetSmilesFunc getSmilesFunc);
double calculateCharPosWeightedSum(const void* context, GetSmilesFunc getSmilesFunc);
double calculateNibbleEntropy(const void* context, GetSmilesFunc getSmilesFunc);
double calculateCharFrequencyEntropy(const void* context, GetSmilesFunc getSmilesFunc);
double calculateAsciiRangeCoverage(const void* context, GetSmilesFunc getSmilesFunc);
double calculateNibbleHistogramBalance(const void* context, GetSmilesFunc getSmilesFunc);
double calculateBinaryCompressionRatio(const void* context, GetSmilesFunc getSmilesFunc);
double calculateBitReversalDistance(const void* context, GetSmilesFunc getSmilesFunc);
double calculateByteReversalDistance(const void* context, GetSmilesFunc getSmilesFunc);
double calculateMurmurHash3(const void* context, GetSmilesFunc getSmilesFunc);
double calculateCharRunLength(const void* context, GetSmilesFunc getSmilesFunc);
double calculateByteRunEntropy(const void* context, GetSmilesFunc getSmilesFunc);
double calculateSingleBitPatternMatch(const void* context, GetSmilesFunc getSmilesFunc);
double calculateBytePairTransitions(const void* context, GetSmilesFunc getSmilesFunc);
double calculateBitDifferenceSum(const void* context, GetSmilesFunc getSmilesFunc);
double calculateFirstLastByteXor(const void* context, GetSmilesFunc getSmilesFunc);
double calculateByteMarkovProperty(const void* context, GetSmilesFunc getSmilesFunc);
double calculateSwarPopcount(const void* context, GetSmilesFunc getSmilesFunc);
double calculateFastBitHistogram(const void* context, GetSmilesFunc getSmilesFunc);
double calculateBitDispersionIndex(const void* context, GetSmilesFunc getSmilesFunc);
double calculateHalfByteRotationHash(const void* context, GetSmilesFunc getSmilesFunc);
double calculateAvalanchePatternSensitivity(const void* context, GetSmilesFunc getSmilesFunc);
double calculateSequentialXorDiff(const void* context, GetSmilesFunc getSmilesFunc);
double calculateByteDeltaEncodingSize(const void* context, GetSmilesFunc getSmilesFunc);
double calculateBitCorrelationCoef(const void* context, GetSmilesFunc getSmilesFunc);
double calculateBitRunLengthSize(const void* context, GetSmilesFunc getSmilesFunc);
double calculateByteMirroringDistance(const void* context, GetSmilesFunc getSmilesFunc);
double calculateRollingByteChecksum(const void* context, GetSmilesFunc getSmilesFunc);
double calculateZeroBitRatio(const void* context, GetSmilesFunc getSmilesFunc);
double calculateBitAutocorrelation(const void* context, GetSmilesFunc getSmilesFunc);
double calculateDiagonalBitAccumulation(const void* context, GetSmilesFunc getSmilesFunc);
double calculateBitwiseMajorityFunction(const void* context, GetSmilesFunc getSmilesFunc);
double calculateLzcntPattern(const void* context, GetSmilesFunc getSmilesFunc);
double calculateTzcntPattern(const void* context, GetSmilesFunc getSmilesFunc);
double calculateByteGeometricMean(const void* context, GetSmilesFunc getSmilesFunc);

// --- NEW Advanced Topological/Physicochemical Descriptors ---
double calculateSMILESPathLengthIndex(const void* context, GetSmilesFunc getSmilesFunc);
double calculateChainBranchFactor(const void* context, GetSmilesFunc getSmilesFunc);
double calculateAtomicNumberSum(const void* context, GetSmilesFunc getSmilesFunc);
double calculateRingSizeEnergy(const void* context, GetSmilesFunc getSmilesFunc);
double calculateElectronegativityPattern(const void* context, GetSmilesFunc getSmilesFunc);
double calculateBondCharacterIndex(const void* context, GetSmilesFunc getSmilesFunc);
double calculateBracketComplexity(const void* context, GetSmilesFunc getSmilesFunc);
double calculateHeteroatomDistribution(const void* context, GetSmilesFunc getSmilesFunc);
double calculateSMILESCompressionRatio(const void* context, GetSmilesFunc getSmilesFunc);
double calculateCharacterTransitionEntropy(const void* context, GetSmilesFunc getSmilesFunc);
double calculateAtomicWeightDispersion(const void* context, GetSmilesFunc getSmilesFunc);
double calculateBranchDepthEnergy(const void* context, GetSmilesFunc getSmilesFunc);
double calculateValenceElectronDistribution(const void* context, GetSmilesFunc getSmilesFunc);
double calculateRingClosureDistance(const void* context, GetSmilesFunc getSmilesFunc);
double calculateElementPeriodicityIndex(const void* context, GetSmilesFunc getSmilesFunc);
double calculateBondOrderSum(const void* context, GetSmilesFunc getSmilesFunc);
double calculateCharacterClassEntropy(const void* context, GetSmilesFunc getSmilesFunc);
double calculateStructuralInformationContent(const void* context, GetSmilesFunc getSmilesFunc);
double calculateNitrogenPositionEffect(const void* context, GetSmilesFunc getSmilesFunc);
double calculateOxygenClusteringFactor(const void* context, GetSmilesFunc getSmilesFunc);
double calculateSequentialElementContrast(const void* context, GetSmilesFunc getSmilesFunc);
double calculateCarbonChainStability(const void* context, GetSmilesFunc getSmilesFunc);
double calculateSyntacticNestingLevel(const void* context, GetSmilesFunc getSmilesFunc);
double calculateFunctionalGroupPosition(const void* context, GetSmilesFunc getSmilesFunc);
double calculateHalogenIndex(const void* context, GetSmilesFunc getSmilesFunc);
double calculateHybridizationPattern(const void* context, GetSmilesFunc getSmilesFunc);
double calculateAromaticCharacterRatio(const void* context, GetSmilesFunc getSmilesFunc);
double calculateBranchPointDensity(const void* context, GetSmilesFunc getSmilesFunc);
double calculateElementDiversityMeasure(const void* context, GetSmilesFunc getSmilesFunc);
double calculateFunctionalGroupConnectionIndex(const void* context, GetSmilesFunc getSmilesFunc);

#ifdef __cplusplus
}
#endif

#endif // CREGISTRY_H
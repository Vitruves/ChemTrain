#include "../common.hpp"
#include <GraphMol/ROMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/MolOps.h>        // For topological indices, distance matrix etc.
#include <GraphMol/AtomIterators.h> // For iterating atoms
#include <GraphMol/BondIterators.h> // For iterating bonds
#include <GraphMol/RingInfo.h>      // For ring information
#include <vector>
#include <cmath>
#include <limits>
#include <set>
#include <map>

namespace desfact {

// --- Thread-local cache for graph-based statistics ---
struct GraphStatsCache {
    // --- Basic Counts ---
    unsigned int numAtoms = 0;
    unsigned int numBonds = 0;
    unsigned int numHeavyAtoms = 0;

    // --- Distance Matrix & Derived ---
    std::vector<std::vector<int>> distMatrix;
    double wienerIndex = 0.0;
    double meanWienerIndex = 0.0;
    double hararyIndex = 0.0; // Reciprocal Wiener
    int graphDiameter = 0;
    int graphRadius = std::numeric_limits<int>::max();
    std::vector<int> eccentricities;
    double minEccentricity = 0.0;
    double maxEccentricity = 0.0;
    double meanEccentricity = 0.0;
    double varianceEccentricity = 0.0;
    double petitjeanIndex = 0.0;

    // --- Degree Statistics ---
    std::vector<int> degrees;
    int minDegree = 0;
    int maxDegree = 0;
    double meanDegree = 0.0;
    double varianceDegree = 0.0;
    std::vector<int> degreeCounts; // Counts for degree 0, 1, 2, 3, 4, 5+
    std::vector<double> degreeFractions;
    double plattIndex = 0.0;        // Sum of edge degrees
    double meanBondDegreeSum = 0.0; // Avg sigma value for bonds
    double varianceBondDegreeSum = 0.0;
    double zagrebIndexM1 = 0.0;
    double zagrebIndexM2 = 0.0;
    double meanZagrebIndexM1 = 0.0;
    double meanZagrebIndexM2 = 0.0;
    double quadraticIndex = 0.0;
    double graphDensity = 0.0;
    double randicIndex = 0.0;

    // --- Cycle Information ---
    const RDKit::RingInfo* ringInfo = nullptr;
    unsigned int cyclomaticNumber = 0;
    unsigned int numRings = 0;
    unsigned int numAliphaticCarbocycles = 0;
    unsigned int numAliphaticHeterocycles = 0;
    unsigned int numAliphaticRings = 0;
    unsigned int numAromaticCarbocycles = 0;
    unsigned int numAromaticHeterocycles = 0;
    unsigned int numAromaticRings = 0;
    unsigned int numSaturatedCarbocycles = 0;
    unsigned int numSaturatedHeterocycles = 0;
    unsigned int numSaturatedRings = 0;
    unsigned int numBridgeheadAtoms = 0;
    unsigned int numSpiroAtoms = 0;
    unsigned int numFusedRings = 0; // Estimate
    unsigned int numAtomsInRings = 0;
    unsigned int numBondsInRings = 0;
    double fractionAtomsInRings = 0.0;
    double fractionBondsInRings = 0.0;
    int minRingSize = 0;
    int maxRingSize = 0;
    double meanRingSize = 0.0;
    std::vector<int> ringSizeCounts; // Counts for size 3, 4, 5, 6, 7, 8+
    std::vector<double> ringSizeFractions;
    std::vector<std::vector<int>> atomRings;

    // --- Topological & Shape Indices ---
    double balabanJIndex = 0.0;
    double kierKappa1 = 0.0;
    double kierKappa2 = 0.0;
    double kierKappa3 = 0.0;
    double kierAlpha = 0.0;
    double kierPhi = 0.0;

    // --- Branching ---
    unsigned int numTerminalNodes = 0;
    unsigned int numBranchingNodes = 0;
    double fractionTerminalNodes = 0.0;
    double fractionBranchingNodes = 0.0;

    // --- Complexity ---
    double bertzCT = 0.0;

    const RDKit::ROMol* populatedForMol = nullptr;

    void reset() {
        numAtoms = 0;
        numBonds = 0;
        numHeavyAtoms = 0;
        distMatrix.clear();
        wienerIndex = 0.0;
        meanWienerIndex = 0.0;
        hararyIndex = 0.0;
        graphDiameter = 0;
        graphRadius = std::numeric_limits<int>::max();
        eccentricities.clear();
        minEccentricity = 0.0;
        maxEccentricity = 0.0;
        meanEccentricity = 0.0;
        varianceEccentricity = 0.0;
        petitjeanIndex = 0.0;
        degrees.clear();
        minDegree = 0;
        maxDegree = 0;
        meanDegree = 0.0;
        varianceDegree = 0.0;
        degreeCounts.assign(6, 0);
        degreeFractions.assign(6, 0.0);
        plattIndex = 0.0;
        meanBondDegreeSum = 0.0;
        varianceBondDegreeSum = 0.0;
        zagrebIndexM1 = 0.0;
        zagrebIndexM2 = 0.0;
        meanZagrebIndexM1 = 0.0;
        meanZagrebIndexM2 = 0.0;
        quadraticIndex = 0.0;
        graphDensity = 0.0;
        randicIndex = 0.0;
        ringInfo = nullptr;
        cyclomaticNumber = 0;
        numRings = 0;
        numAliphaticCarbocycles = 0;
        numAliphaticHeterocycles = 0;
        numAliphaticRings = 0;
        numAromaticCarbocycles = 0;
        numAromaticHeterocycles = 0;
        numAromaticRings = 0;
        numSaturatedCarbocycles = 0;
        numSaturatedHeterocycles = 0;
        numSaturatedRings = 0;
        numBridgeheadAtoms = 0;
        numSpiroAtoms = 0;
        numFusedRings = 0;
        numAtomsInRings = 0;
        numBondsInRings = 0;
        fractionAtomsInRings = 0.0;
        fractionBondsInRings = 0.0;
        minRingSize = 0;
        maxRingSize = 0;
        meanRingSize = 0.0;
        ringSizeCounts.assign(6, 0);
        ringSizeFractions.assign(6, 0.0);
        atomRings.clear();
        balabanJIndex = 0.0;
        kierKappa1 = 0.0;
        kierKappa2 = 0.0;
        kierKappa3 = 0.0;
        kierAlpha = 0.0;
        kierPhi = 0.0;
        numTerminalNodes = 0;
        numBranchingNodes = 0;
        fractionTerminalNodes = 0.0;
        fractionBranchingNodes = 0.0;
        bertzCT = 0.0;
        populatedForMol = nullptr;
    }

    void computeFromMol(const RDKit::ROMol* mol) {
        if (!mol) return;
        reset();
        populatedForMol = mol;

        numAtoms = mol->getNumAtoms();
        numBonds = mol->getNumBonds();
        if (numAtoms == 0) return;

        numHeavyAtoms = mol->getNumHeavyAtoms();
        if (numAtoms > 1) {
            graphDensity = 2.0 * numBonds / (numAtoms * (numAtoms - 1.0));
        }

        // Precompute degrees
        degrees.resize(numAtoms);
        double sumDegree = 0.0;
        double sumDegreeSq = 0.0;
        degreeCounts.assign(6, 0); // Index 0=deg0, 1=deg1, ..., 5=deg5+
        minDegree = numAtoms > 0 ? mol->getAtomWithIdx(0)->getDegree() : 0;
        maxDegree = minDegree;
        for (unsigned int i = 0; i < numAtoms; ++i) {
            const RDKit::Atom* atom = mol->getAtomWithIdx(i);
            int deg = atom->getDegree();
            degrees[i] = deg;
            sumDegree += deg;
            sumDegreeSq += deg * deg;
            minDegree = std::min(minDegree, deg);
            maxDegree = std::max(maxDegree, deg);
            if (deg >= 5) degreeCounts[5]++;
            else degreeCounts[deg]++;

            if (atom->getAtomicNum() > 1) {
                if (deg == 1) numTerminalNodes++;
                if (deg > 2) numBranchingNodes++;
            }
        }
        quadraticIndex = sumDegreeSq; // Same as Zagreb M1
        zagrebIndexM1 = sumDegreeSq;
        if (numAtoms > 0) {
            meanDegree = sumDegree / numAtoms;
            varianceDegree = (sumDegreeSq / numAtoms) - (meanDegree * meanDegree);
            meanZagrebIndexM1 = zagrebIndexM1 / numAtoms;
             for(size_t i=0; i<degreeCounts.size(); ++i) {
                degreeFractions[i] = static_cast<double>(degreeCounts[i]) / numAtoms;
            }
        }
         if (numHeavyAtoms > 0) {
             fractionTerminalNodes = static_cast<double>(numTerminalNodes) / numHeavyAtoms;
             fractionBranchingNodes = static_cast<double>(numBranchingNodes) / numHeavyAtoms;
         }


        // Bond-based degree calculations (Zagreb M2, Randic, Platt, AvgSigma)
        plattIndex = 0.0;
        double sumBondDegreeSum = 0.0;
        double sumBondDegreeSumSq = 0.0;
        for (const auto& bond : mol->bonds()) {
            int deg1 = bond->getBeginAtom()->getDegree();
            int deg2 = bond->getEndAtom()->getDegree();
            zagrebIndexM2 += deg1 * deg2;
             if (deg1 > 0 && deg2 > 0) {
                randicIndex += 1.0 / std::sqrt(static_cast<double>(deg1 * deg2));
            }
            // Platt = sum(deg(edge)) = sum(deg(atom_i) + deg(atom_j) - 2) for edge (i,j)
            // Faster: sum(deg_i * deg_i) for atoms = ZagrebM1
            // Platt = sum over edges (deg_i + deg_j) - 2*numBonds
            // Sum over edges (deg_i + deg_j) = sum over nodes (deg_i * deg_i) = ZagrebM1
            // So, Platt = ZagrebM1 - 2*numBonds
            // Let's calculate bond degree sum (sigma = deg_i + deg_j) directly for AvgSigma
            double bondSigma = deg1 + deg2;
            sumBondDegreeSum += bondSigma;
            sumBondDegreeSumSq += bondSigma * bondSigma;
        }
        plattIndex = zagrebIndexM1; // Platt index often defined as sum(deg_i*deg_i) over atoms in cheminformatics
                                    // Alternative definition: sum(deg_i + deg_j - 2) over bonds. Using atom-based version.

        if (numBonds > 0) {
            meanBondDegreeSum = sumBondDegreeSum / numBonds;
            varianceBondDegreeSum = (sumBondDegreeSumSq / numBonds) - (meanBondDegreeSum * meanBondDegreeSum);
        }
         if (numAtoms > 0) {
            meanZagrebIndexM2 = zagrebIndexM2 / numAtoms;
        }

        // Distance Matrix
        distMatrix = mol::calculateDistanceMatrix(mol);

        // Compute Wiener index manually
        wienerIndex = 0.0;
        if (numAtoms > 1) {
            for (unsigned int i = 0; i < numAtoms; ++i) {
                for (unsigned int j = i + 1; j < numAtoms; ++j) {
                    wienerIndex += distMatrix[i][j];
                }
            }
            meanWienerIndex = wienerIndex / (static_cast<double>(numAtoms) * (numAtoms - 1.0) / 2.0);
            double sumReciprocalDist = 0.0;
            for (unsigned int i = 0; i < numAtoms; ++i) {
                for (unsigned int j = i + 1; j < numAtoms; ++j) {
                    if (distMatrix[i][j] > 0) {
                        sumReciprocalDist += 1.0 / static_cast<double>(distMatrix[i][j]);
                    }
                }
            }
            hararyIndex = sumReciprocalDist;
        }


        // Eccentricity, Radius, Diameter
        eccentricities.resize(numAtoms, 0);
        graphRadius = std::numeric_limits<int>::max();
        graphDiameter = 0;
        double sumEccentricity = 0.0;
        double sumEccentricitySq = 0.0;
        for (unsigned int i = 0; i < numAtoms; ++i) {
            int maxDist = 0;
            for (unsigned int j = 0; j < numAtoms; ++j) {
                if (distMatrix[i][j] > maxDist) {
                    maxDist = distMatrix[i][j];
                }
            }
            eccentricities[i] = maxDist;
            sumEccentricity += maxDist;
            sumEccentricitySq += maxDist * maxDist;
            if (maxDist > 0) { // Avoid radius 0 for single atom case affecting min
                 graphRadius = std::min(graphRadius, maxDist);
            }
            graphDiameter = std::max(graphDiameter, maxDist);
        }

        if (numAtoms == 1) graphRadius = 0; // Correct radius for single atom
        if (numAtoms > 0) {
            minEccentricity = static_cast<double>(graphRadius);
            maxEccentricity = static_cast<double>(graphDiameter);
            meanEccentricity = sumEccentricity / numAtoms;
            varianceEccentricity = (sumEccentricitySq / numAtoms) - (meanEccentricity * meanEccentricity);
             if (graphDiameter > 0 && graphRadius > 0) { // Petitjean index
                petitjeanIndex = static_cast<double>(graphDiameter - graphRadius) / graphRadius;
            }
        } else {
            graphRadius = 0; // Handle zero atoms
        }


        // Kier Shape Indices & Alpha
        kierKappa1 = 0.0;
        kierKappa2 = 0.0;
        kierKappa3 = 0.0;
        kierAlpha = 0.0;
        kierPhi = 0.0;


        // Balaban J Index
        balabanJIndex = 0.0;

        // Bertz Complexity Index
        bertzCT = 0.0;

        // Cycle Information
        ringInfo = mol->getRingInfo();
        if (!ringInfo->isInitialized()) {
             // RDKit computes SSSR lazily. Explicit call usually not needed.
             // RDKit::MolOps::findSSSR(*mol); // This function variant might not exist or be deprecated.
             // Accessing ringInfo methods should trigger computation if needed.
             mol->getRingInfo()->numRings(); // Trigger computation if not done yet
        }
        numRings = ringInfo->numRings();
        cyclomaticNumber = numBonds - numAtoms + 1;

        numAliphaticCarbocycles = 0;
        numAliphaticHeterocycles = 0;
        numAliphaticRings = 0;
        numAromaticCarbocycles = 0;
        numAromaticHeterocycles = 0;
        numAromaticRings = 0;
        numSaturatedCarbocycles = 0;
        numSaturatedHeterocycles = 0;
        numSaturatedRings = 0;

        atomRings = ringInfo->atomRings();
        std::set<int> atoms_in_rings_set;
        std::set<int> bonds_in_rings_set;
        std::map<int, int> ring_bond_counts; // bond_idx -> count across rings
        int totalRingSizeSum = 0;
        minRingSize = (numRings > 0) ? std::numeric_limits<int>::max() : 0;
        maxRingSize = 0;
        ringSizeCounts.assign(6, 0); // Size 3, 4, 5, 6, 7, 8+

        for (const auto& ring : atomRings) {
            int currentSize = ring.size();
             if (currentSize > 0) {
                 totalRingSizeSum += currentSize;
                 minRingSize = std::min(minRingSize, currentSize);
                 maxRingSize = std::max(maxRingSize, currentSize);
                 if(currentSize >= 8) ringSizeCounts[5]++;
                 else if (currentSize >= 3) ringSizeCounts[currentSize-3]++; // Index 0=size3, 1=size4...
             }
            for (int atomIdx : ring) {
                atoms_in_rings_set.insert(atomIdx);
            }
        }
        numAtomsInRings = atoms_in_rings_set.size();

        const auto& bondRings = ringInfo->bondRings();
         for (const auto& ring : bondRings) {
             for (int bondIdx : ring) {
                 bonds_in_rings_set.insert(bondIdx);
                 ring_bond_counts[bondIdx]++;
             }
         }
        numBondsInRings = bonds_in_rings_set.size();

        // Count fused rings: bonds shared by >1 ring
        numFusedRings = 0;
        for(const auto& pair : ring_bond_counts) {
            if (pair.second > 1) {
                numFusedRings += (pair.second - 1); // Each count > 1 implies fusion points
            }
        }
        // Note: This is a simple heuristic for fused ring count. More rigorous definitions exist.


        // Bridgehead and Spiro Atoms
        numBridgeheadAtoms = 0;
        numSpiroAtoms = 0;


        if (numAtoms > 0) {
            fractionAtomsInRings = static_cast<double>(numAtomsInRings) / numAtoms;
        }
        if (numBonds > 0) {
            fractionBondsInRings = static_cast<double>(numBondsInRings) / numBonds;
        }
        if (numRings > 0) {
            meanRingSize = static_cast<double>(totalRingSizeSum) / numRings;
            for(size_t i=0; i<ringSizeCounts.size(); ++i) {
                ringSizeFractions[i] = static_cast<double>(ringSizeCounts[i]) / numRings;
            }
        } else {
             minRingSize = 0;
        }
    }
};

static thread_local GraphStatsCache tlsGraphCache;

// Helper to ensure cache is populated
inline void ensureCachePopulated(Context& context) {
    const RDKit::ROMol* mol = context.getMolecule();
    if (tlsGraphCache.populatedForMol != mol || !mol) {
        // std::cout << "Populating graph cache for molecule: " << context.getSmiles() << std::endl; // Debugging
        tlsGraphCache.computeFromMol(mol);
    }
}

// --- Base Class ---
class GraphStatsDescriptor : public Descriptor {
public:
    using Descriptor::Descriptor;
    std::string getCategory() const override { return "GraphStats"; }
    std::vector<std::string> getDependencies() const override { return {}; }
};

// --- Descriptor Implementations ---

// Basic Counts & Ratios (4)
DECLARE_DESCRIPTOR(NumAtoms, GraphStatsDescriptor, "Total number of atoms (incl. H)")
DESCRIPTOR_DEPENDENCIES(NumAtoms) { return {}; }
DescriptorResult NumAtomsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return static_cast<int>(tlsGraphCache.numAtoms); }

DECLARE_DESCRIPTOR(NumBonds, GraphStatsDescriptor, "Total number of bonds")
DESCRIPTOR_DEPENDENCIES(NumBonds) { return {}; }
DescriptorResult NumBondsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return static_cast<int>(tlsGraphCache.numBonds); }

DECLARE_DESCRIPTOR(NumHeavyAtomsGraph, GraphStatsDescriptor, "Number of heavy atoms")
DESCRIPTOR_DEPENDENCIES(NumHeavyAtomsGraph) { return {}; }
DescriptorResult NumHeavyAtomsGraphDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return static_cast<int>(tlsGraphCache.numHeavyAtoms); }

DECLARE_DESCRIPTOR(EdgeNodeRatio, GraphStatsDescriptor, "Ratio B/N")
DESCRIPTOR_DEPENDENCIES(EdgeNodeRatio) { return {}; }
DescriptorResult EdgeNodeRatioDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.numAtoms > 0 ? tlsGraphCache.numBonds / static_cast<double>(tlsGraphCache.numAtoms) : 0.0; }

// Connectivity & Distance-Based (11)
DECLARE_DESCRIPTOR(WienerIndex, GraphStatsDescriptor, "Wiener index (sum of topological distances)")
DESCRIPTOR_DEPENDENCIES(WienerIndex) { return {}; }
DescriptorResult WienerIndexDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.wienerIndex; }

DECLARE_DESCRIPTOR(MeanWienerIndex, GraphStatsDescriptor, "Mean Wiener index (Wiener / (N*(N-1)/2))")
DESCRIPTOR_DEPENDENCIES(MeanWienerIndex) { return {}; }
DescriptorResult MeanWienerIndexDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.meanWienerIndex; }

DECLARE_DESCRIPTOR(HararyIndex, GraphStatsDescriptor, "Harary index (sum of reciprocal topological distances)")
DESCRIPTOR_DEPENDENCIES(HararyIndex) { return {}; }
DescriptorResult HararyIndexDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.hararyIndex; }

DECLARE_DESCRIPTOR(BalabanJIndex, GraphStatsDescriptor, "Balaban J index")
DESCRIPTOR_DEPENDENCIES(BalabanJIndex) { return {}; }
DescriptorResult BalabanJIndexDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.balabanJIndex; }

DECLARE_DESCRIPTOR(GraphDiameter, GraphStatsDescriptor, "Graph diameter (longest shortest path)")
DESCRIPTOR_DEPENDENCIES(GraphDiameter) { return {}; }
DescriptorResult GraphDiameterDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.graphDiameter; }

DECLARE_DESCRIPTOR(GraphRadius, GraphStatsDescriptor, "Graph radius (minimum atomic eccentricity)")
DESCRIPTOR_DEPENDENCIES(GraphRadius) { return {}; }
DescriptorResult GraphRadiusDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.graphRadius; }

DECLARE_DESCRIPTOR(MinAtomEccentricity, GraphStatsDescriptor, "Minimum atomic eccentricity")
DESCRIPTOR_DEPENDENCIES(MinAtomEccentricity) { return {}; }
DescriptorResult MinAtomEccentricityDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.minEccentricity; }

DECLARE_DESCRIPTOR(MaxAtomEccentricity, GraphStatsDescriptor, "Maximum atomic eccentricity")
DESCRIPTOR_DEPENDENCIES(MaxAtomEccentricity) { return {}; }
DescriptorResult MaxAtomEccentricityDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.maxEccentricity; }

DECLARE_DESCRIPTOR(MeanAtomEccentricity, GraphStatsDescriptor, "Mean atomic eccentricity")
DESCRIPTOR_DEPENDENCIES(MeanAtomEccentricity) { return {}; }
DescriptorResult MeanAtomEccentricityDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.meanEccentricity; }

DECLARE_DESCRIPTOR(VarianceAtomEccentricity, GraphStatsDescriptor, "Variance of atomic eccentricities")
DESCRIPTOR_DEPENDENCIES(VarianceAtomEccentricity) { return {}; }
DescriptorResult VarianceAtomEccentricityDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.varianceEccentricity; }

DECLARE_DESCRIPTOR(PetitjeanIndex, GraphStatsDescriptor, "Petitjean graph index ((Diameter - Radius) / Radius)")
DESCRIPTOR_DEPENDENCIES(PetitjeanIndex) { return {}; }
DescriptorResult PetitjeanIndexDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.petitjeanIndex; }


// Degree-Based (18)
DECLARE_DESCRIPTOR(MinAtomDegree, GraphStatsDescriptor, "Minimum atomic degree")
DESCRIPTOR_DEPENDENCIES(MinAtomDegree) { return {}; }
DescriptorResult MinAtomDegreeDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.minDegree; }

DECLARE_DESCRIPTOR(MaxAtomDegree, GraphStatsDescriptor, "Maximum atomic degree")
DESCRIPTOR_DEPENDENCIES(MaxAtomDegree) { return {}; }
DescriptorResult MaxAtomDegreeDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.maxDegree; }

DECLARE_DESCRIPTOR(MeanAtomDegree, GraphStatsDescriptor, "Mean atomic degree")
DESCRIPTOR_DEPENDENCIES(MeanAtomDegree) { return {}; }
DescriptorResult MeanAtomDegreeDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.meanDegree; }

DECLARE_DESCRIPTOR(VarianceAtomDegree, GraphStatsDescriptor, "Variance of atomic degrees")
DESCRIPTOR_DEPENDENCIES(VarianceAtomDegree) { return {}; }
DescriptorResult VarianceAtomDegreeDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.varianceDegree; }

DECLARE_DESCRIPTOR(NumDegreeZeroAtoms, GraphStatsDescriptor, "Number of atoms with degree 0")
DESCRIPTOR_DEPENDENCIES(NumDegreeZeroAtoms) { return {}; }
DescriptorResult NumDegreeZeroAtomsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.degreeCounts[0]; }

DECLARE_DESCRIPTOR(NumDegreeOneAtoms, GraphStatsDescriptor, "Number of atoms with degree 1")
DESCRIPTOR_DEPENDENCIES(NumDegreeOneAtoms) { return {}; }
DescriptorResult NumDegreeOneAtomsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.degreeCounts[1]; }

DECLARE_DESCRIPTOR(NumDegreeTwoAtoms, GraphStatsDescriptor, "Number of atoms with degree 2")
DESCRIPTOR_DEPENDENCIES(NumDegreeTwoAtoms) { return {}; }
DescriptorResult NumDegreeTwoAtomsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.degreeCounts[2]; }

DECLARE_DESCRIPTOR(NumDegreeThreeAtoms, GraphStatsDescriptor, "Number of atoms with degree 3")
DESCRIPTOR_DEPENDENCIES(NumDegreeThreeAtoms) { return {}; }
DescriptorResult NumDegreeThreeAtomsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.degreeCounts[3]; }

DECLARE_DESCRIPTOR(NumDegreeFourAtoms, GraphStatsDescriptor, "Number of atoms with degree 4")
DESCRIPTOR_DEPENDENCIES(NumDegreeFourAtoms) { return {}; }
DescriptorResult NumDegreeFourAtomsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.degreeCounts[4]; }

DECLARE_DESCRIPTOR(NumDegreeFivePlusAtoms, GraphStatsDescriptor, "Number of atoms with degree 5 or more")
DESCRIPTOR_DEPENDENCIES(NumDegreeFivePlusAtoms) { return {}; }
DescriptorResult NumDegreeFivePlusAtomsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.degreeCounts[5]; }

DECLARE_DESCRIPTOR(FractionDegreeZeroAtoms, GraphStatsDescriptor, "Fraction of atoms with degree 0")
DESCRIPTOR_DEPENDENCIES(FractionDegreeZeroAtoms) { return {}; }
DescriptorResult FractionDegreeZeroAtomsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.degreeFractions[0]; }

DECLARE_DESCRIPTOR(FractionDegreeOneAtoms, GraphStatsDescriptor, "Fraction of atoms with degree 1")
DESCRIPTOR_DEPENDENCIES(FractionDegreeOneAtoms) { return {}; }
DescriptorResult FractionDegreeOneAtomsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.degreeFractions[1]; }

DECLARE_DESCRIPTOR(FractionDegreeTwoAtoms, GraphStatsDescriptor, "Fraction of atoms with degree 2")
DESCRIPTOR_DEPENDENCIES(FractionDegreeTwoAtoms) { return {}; }
DescriptorResult FractionDegreeTwoAtomsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.degreeFractions[2]; }

DECLARE_DESCRIPTOR(FractionDegreeThreeAtoms, GraphStatsDescriptor, "Fraction of atoms with degree 3")
DESCRIPTOR_DEPENDENCIES(FractionDegreeThreeAtoms) { return {}; }
DescriptorResult FractionDegreeThreeAtomsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.degreeFractions[3]; }

DECLARE_DESCRIPTOR(FractionDegreeFourAtoms, GraphStatsDescriptor, "Fraction of atoms with degree 4")
DESCRIPTOR_DEPENDENCIES(FractionDegreeFourAtoms) { return {}; }
DescriptorResult FractionDegreeFourAtomsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.degreeFractions[4]; }

DECLARE_DESCRIPTOR(FractionDegreeFivePlusAtoms, GraphStatsDescriptor, "Fraction of atoms with degree 5 or more")
DESCRIPTOR_DEPENDENCIES(FractionDegreeFivePlusAtoms) { return {}; }
DescriptorResult FractionDegreeFivePlusAtomsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.degreeFractions[5]; }

DECLARE_DESCRIPTOR(PlattIndex, GraphStatsDescriptor, "Platt index (sum of degree*degree over atoms)")
DESCRIPTOR_DEPENDENCIES(PlattIndex) { return {}; }
DescriptorResult PlattIndexDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.plattIndex; }

DECLARE_DESCRIPTOR(MeanBondDegreeSum, GraphStatsDescriptor, "Mean sum of degrees for atoms forming each bond")
DESCRIPTOR_DEPENDENCIES(MeanBondDegreeSum) { return {}; }
DescriptorResult MeanBondDegreeSumDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.meanBondDegreeSum; }

DECLARE_DESCRIPTOR(VarianceBondDegreeSum, GraphStatsDescriptor, "Variance of sum of degrees for atoms forming each bond")
DESCRIPTOR_DEPENDENCIES(VarianceBondDegreeSum) { return {}; }
DescriptorResult VarianceBondDegreeSumDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.varianceBondDegreeSum; }

DECLARE_DESCRIPTOR(ZagrebIndexM1, GraphStatsDescriptor, "First Zagreb index (sum of deg*deg)")
DESCRIPTOR_DEPENDENCIES(ZagrebIndexM1) { return {}; }
DescriptorResult ZagrebIndexM1Descriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.zagrebIndexM1; }

DECLARE_DESCRIPTOR(ZagrebIndexM2, GraphStatsDescriptor, "Second Zagreb index (sum of deg*deg over bonds)")
DESCRIPTOR_DEPENDENCIES(ZagrebIndexM2) { return {}; }
DescriptorResult ZagrebIndexM2Descriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.zagrebIndexM2; }

DECLARE_DESCRIPTOR(MeanZagrebIndexM1, GraphStatsDescriptor, "Mean First Zagreb index (M1/N)")
DESCRIPTOR_DEPENDENCIES(MeanZagrebIndexM1) { return {}; }
DescriptorResult MeanZagrebIndexM1Descriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.meanZagrebIndexM1; }

DECLARE_DESCRIPTOR(MeanZagrebIndexM2, GraphStatsDescriptor, "Mean Second Zagreb index (M2/N)")
DESCRIPTOR_DEPENDENCIES(MeanZagrebIndexM2) { return {}; }
DescriptorResult MeanZagrebIndexM2Descriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.meanZagrebIndexM2; }

DECLARE_DESCRIPTOR(QuadraticIndex, GraphStatsDescriptor, "Quadratic index (sum of deg*deg, same as M1)")
DESCRIPTOR_DEPENDENCIES(QuadraticIndex) { return {}; }
DescriptorResult QuadraticIndexDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.quadraticIndex; }

DECLARE_DESCRIPTOR(GraphDensity, GraphStatsDescriptor, "Graph density (2*B / (N*(N-1)))")
DESCRIPTOR_DEPENDENCIES(GraphDensity) { return {}; }
DescriptorResult GraphDensityDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.graphDensity; }

DECLARE_DESCRIPTOR(RandicIndex, GraphStatsDescriptor, "Randic connectivity index")
DESCRIPTOR_DEPENDENCIES(RandicIndex) { return {}; }
DescriptorResult RandicIndexDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.randicIndex; }


// Topological & Shape (5)
DECLARE_DESCRIPTOR(KierKappa1, GraphStatsDescriptor, "Kier kappa 1 shape index")
DESCRIPTOR_DEPENDENCIES(KierKappa1) { return {}; }
DescriptorResult KierKappa1Descriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.kierKappa1; }

DECLARE_DESCRIPTOR(KierKappa2, GraphStatsDescriptor, "Kier kappa 2 shape index")
DESCRIPTOR_DEPENDENCIES(KierKappa2) { return {}; }
DescriptorResult KierKappa2Descriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.kierKappa2; }

DECLARE_DESCRIPTOR(KierKappa3, GraphStatsDescriptor, "Kier kappa 3 shape index")
DESCRIPTOR_DEPENDENCIES(KierKappa3) { return {}; }
DescriptorResult KierKappa3Descriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.kierKappa3; }

DECLARE_DESCRIPTOR(KierAlpha, GraphStatsDescriptor, "Hall-Kier alpha value")
DESCRIPTOR_DEPENDENCIES(KierAlpha) { return {}; }
DescriptorResult KierAlphaDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.kierAlpha; }

DECLARE_DESCRIPTOR(KierPhi, GraphStatsDescriptor, "Kier flexibility index (K1*K3 / K2^2)")
DESCRIPTOR_DEPENDENCIES(KierPhi) { return {}; }
DescriptorResult KierPhiDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.kierPhi; }

// Cycle-Based (23)
DECLARE_DESCRIPTOR(CyclomaticNumber, GraphStatsDescriptor, "Cyclomatic number (B - N + C)")
DESCRIPTOR_DEPENDENCIES(CyclomaticNumber) { return {}; }
DescriptorResult CyclomaticNumberDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return static_cast<int>(tlsGraphCache.cyclomaticNumber); }

DECLARE_DESCRIPTOR(NumRings, GraphStatsDescriptor, "Number of SSSR rings")
DESCRIPTOR_DEPENDENCIES(NumRings) { return {}; }
DescriptorResult NumRingsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return static_cast<int>(tlsGraphCache.numRings); }

DECLARE_DESCRIPTOR(NumAliphaticCarbocycles, GraphStatsDescriptor, "Number of aliphatic carbocycles")
DESCRIPTOR_DEPENDENCIES(NumAliphaticCarbocycles) { return {}; }
DescriptorResult NumAliphaticCarbocyclesDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return static_cast<int>(tlsGraphCache.numAliphaticCarbocycles); }

DECLARE_DESCRIPTOR(NumAliphaticHeterocycles, GraphStatsDescriptor, "Number of aliphatic heterocycles")
DESCRIPTOR_DEPENDENCIES(NumAliphaticHeterocycles) { return {}; }
DescriptorResult NumAliphaticHeterocyclesDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return static_cast<int>(tlsGraphCache.numAliphaticHeterocycles); }

DECLARE_DESCRIPTOR(NumAliphaticRings, GraphStatsDescriptor, "Total number of aliphatic rings")
DESCRIPTOR_DEPENDENCIES(NumAliphaticRings) { return {}; }
DescriptorResult NumAliphaticRingsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return static_cast<int>(tlsGraphCache.numAliphaticRings); }

DECLARE_DESCRIPTOR(NumAromaticCarbocycles, GraphStatsDescriptor, "Number of aromatic carbocycles")
DESCRIPTOR_DEPENDENCIES(NumAromaticCarbocycles) { return {}; }
DescriptorResult NumAromaticCarbocyclesDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return static_cast<int>(tlsGraphCache.numAromaticCarbocycles); }

DECLARE_DESCRIPTOR(NumAromaticHeterocycles, GraphStatsDescriptor, "Number of aromatic heterocycles")
DESCRIPTOR_DEPENDENCIES(NumAromaticHeterocycles) { return {}; }
DescriptorResult NumAromaticHeterocyclesDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return static_cast<int>(tlsGraphCache.numAromaticHeterocycles); }

DECLARE_DESCRIPTOR(NumAromaticRings, GraphStatsDescriptor, "Total number of aromatic rings")
DESCRIPTOR_DEPENDENCIES(NumAromaticRings) { return {}; }
DescriptorResult NumAromaticRingsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return static_cast<int>(tlsGraphCache.numAromaticRings); }

DECLARE_DESCRIPTOR(NumSaturatedCarbocycles, GraphStatsDescriptor, "Number of saturated carbocycles")
DESCRIPTOR_DEPENDENCIES(NumSaturatedCarbocycles) { return {}; }
DescriptorResult NumSaturatedCarbocyclesDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return static_cast<int>(tlsGraphCache.numSaturatedCarbocycles); }

DECLARE_DESCRIPTOR(NumSaturatedHeterocycles, GraphStatsDescriptor, "Number of saturated heterocycles")
DESCRIPTOR_DEPENDENCIES(NumSaturatedHeterocycles) { return {}; }
DescriptorResult NumSaturatedHeterocyclesDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return static_cast<int>(tlsGraphCache.numSaturatedHeterocycles); }

DECLARE_DESCRIPTOR(NumSaturatedRings, GraphStatsDescriptor, "Total number of saturated rings")
DESCRIPTOR_DEPENDENCIES(NumSaturatedRings) { return {}; }
DescriptorResult NumSaturatedRingsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return static_cast<int>(tlsGraphCache.numSaturatedRings); }

DECLARE_DESCRIPTOR(NumBridgeheadAtoms, GraphStatsDescriptor, "Number of bridgehead atoms")
DESCRIPTOR_DEPENDENCIES(NumBridgeheadAtoms) { return {}; }
DescriptorResult NumBridgeheadAtomsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return static_cast<int>(tlsGraphCache.numBridgeheadAtoms); }

DECLARE_DESCRIPTOR(NumSpiroAtoms, GraphStatsDescriptor, "Number of spiro atoms")
DESCRIPTOR_DEPENDENCIES(NumSpiroAtoms) { return {}; }
DescriptorResult NumSpiroAtomsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return static_cast<int>(tlsGraphCache.numSpiroAtoms); }

DECLARE_DESCRIPTOR(NumFusedRingsEstimate, GraphStatsDescriptor, "Estimated number of fused rings (based on shared bonds)")
DESCRIPTOR_DEPENDENCIES(NumFusedRingsEstimate) { return {}; }
DescriptorResult NumFusedRingsEstimateDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return static_cast<int>(tlsGraphCache.numFusedRings); }


DECLARE_DESCRIPTOR(FractionAtomsInRings, GraphStatsDescriptor, "Fraction of atoms in rings")
DESCRIPTOR_DEPENDENCIES(FractionAtomsInRings) { return {}; }
DescriptorResult FractionAtomsInRingsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.fractionAtomsInRings; }

DECLARE_DESCRIPTOR(FractionBondsInRings, GraphStatsDescriptor, "Fraction of bonds in rings")
DESCRIPTOR_DEPENDENCIES(FractionBondsInRings) { return {}; }
DescriptorResult FractionBondsInRingsDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.fractionBondsInRings; }

DECLARE_DESCRIPTOR(MinRingSize, GraphStatsDescriptor, "Minimum ring size in SSSR")
DESCRIPTOR_DEPENDENCIES(MinRingSize) { return {}; }
DescriptorResult MinRingSizeDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.minRingSize; }

DECLARE_DESCRIPTOR(MaxRingSize, GraphStatsDescriptor, "Maximum ring size in SSSR")
DESCRIPTOR_DEPENDENCIES(MaxRingSize) { return {}; }
DescriptorResult MaxRingSizeDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.maxRingSize; }

DECLARE_DESCRIPTOR(MeanRingSize, GraphStatsDescriptor, "Mean ring size in SSSR")
DESCRIPTOR_DEPENDENCIES(MeanRingSize) { return {}; }
DescriptorResult MeanRingSizeDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.meanRingSize; }

DECLARE_DESCRIPTOR(NumRingSize3, GraphStatsDescriptor, "Number of 3-membered rings")
DESCRIPTOR_DEPENDENCIES(NumRingSize3) { return {}; }
DescriptorResult NumRingSize3Descriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.ringSizeCounts[0]; }

DECLARE_DESCRIPTOR(NumRingSize4, GraphStatsDescriptor, "Number of 4-membered rings")
DESCRIPTOR_DEPENDENCIES(NumRingSize4) { return {}; }
DescriptorResult NumRingSize4Descriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.ringSizeCounts[1]; }

DECLARE_DESCRIPTOR(NumRingSize5, GraphStatsDescriptor, "Number of 5-membered rings")
DESCRIPTOR_DEPENDENCIES(NumRingSize5) { return {}; }
DescriptorResult NumRingSize5Descriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.ringSizeCounts[2]; }

DECLARE_DESCRIPTOR(NumRingSize6, GraphStatsDescriptor, "Number of 6-membered rings")
DESCRIPTOR_DEPENDENCIES(NumRingSize6) { return {}; }
DescriptorResult NumRingSize6Descriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.ringSizeCounts[3]; }

DECLARE_DESCRIPTOR(NumRingSize7, GraphStatsDescriptor, "Number of 7-membered rings")
DESCRIPTOR_DEPENDENCIES(NumRingSize7) { return {}; }
DescriptorResult NumRingSize7Descriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.ringSizeCounts[4]; }

DECLARE_DESCRIPTOR(NumRingSize8Plus, GraphStatsDescriptor, "Number of rings with size 8 or more")
DESCRIPTOR_DEPENDENCIES(NumRingSize8Plus) { return {}; }
DescriptorResult NumRingSize8PlusDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.ringSizeCounts[5]; }

DECLARE_DESCRIPTOR(FractionRingSize3, GraphStatsDescriptor, "Fraction of rings that are 3-membered")
DESCRIPTOR_DEPENDENCIES(FractionRingSize3) { return {}; }
DescriptorResult FractionRingSize3Descriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.ringSizeFractions[0]; }

DECLARE_DESCRIPTOR(FractionRingSize4, GraphStatsDescriptor, "Fraction of rings that are 4-membered")
DESCRIPTOR_DEPENDENCIES(FractionRingSize4) { return {}; }
DescriptorResult FractionRingSize4Descriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.ringSizeFractions[1]; }

DECLARE_DESCRIPTOR(FractionRingSize5, GraphStatsDescriptor, "Fraction of rings that are 5-membered")
DESCRIPTOR_DEPENDENCIES(FractionRingSize5) { return {}; }
DescriptorResult FractionRingSize5Descriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.ringSizeFractions[2]; }

DECLARE_DESCRIPTOR(FractionRingSize6, GraphStatsDescriptor, "Fraction of rings that are 6-membered")
DESCRIPTOR_DEPENDENCIES(FractionRingSize6) { return {}; }
DescriptorResult FractionRingSize6Descriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.ringSizeFractions[3]; }

DECLARE_DESCRIPTOR(FractionRingSize7, GraphStatsDescriptor, "Fraction of rings that are 7-membered")
DESCRIPTOR_DEPENDENCIES(FractionRingSize7) { return {}; }
DescriptorResult FractionRingSize7Descriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.ringSizeFractions[4]; }

DECLARE_DESCRIPTOR(FractionRingSize8Plus, GraphStatsDescriptor, "Fraction of rings with size 8 or more")
DESCRIPTOR_DEPENDENCIES(FractionRingSize8Plus) { return {}; }
DescriptorResult FractionRingSize8PlusDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.ringSizeFractions[5]; }


// Branching (4)
DECLARE_DESCRIPTOR(NumTerminalNodes, GraphStatsDescriptor, "Number of terminal heavy atoms (degree 1)")
DESCRIPTOR_DEPENDENCIES(NumTerminalNodes) { return {}; }
DescriptorResult NumTerminalNodesDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return static_cast<int>(tlsGraphCache.numTerminalNodes); }

DECLARE_DESCRIPTOR(FractionTerminalNodes, GraphStatsDescriptor, "Fraction of heavy atoms that are terminal")
DESCRIPTOR_DEPENDENCIES(FractionTerminalNodes) { return {}; }
DescriptorResult FractionTerminalNodesDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.fractionTerminalNodes; }

DECLARE_DESCRIPTOR(NumBranchingNodes, GraphStatsDescriptor, "Number of branching heavy atoms (degree > 2)")
DESCRIPTOR_DEPENDENCIES(NumBranchingNodes) { return {}; }
DescriptorResult NumBranchingNodesDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return static_cast<int>(tlsGraphCache.numBranchingNodes); }

DECLARE_DESCRIPTOR(FractionBranchingNodes, GraphStatsDescriptor, "Fraction of heavy atoms that are branching")
DESCRIPTOR_DEPENDENCIES(FractionBranchingNodes) { return {}; }
DescriptorResult FractionBranchingNodesDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.fractionBranchingNodes; }

// Complexity (1)
DECLARE_DESCRIPTOR(BertzCT, GraphStatsDescriptor, "Bertz Complexity Index")
DESCRIPTOR_DEPENDENCIES(BertzCT) { return {}; }
DescriptorResult BertzCTDescriptor::calculate(Context& context) const { ensureCachePopulated(context); return tlsGraphCache.bertzCT; }


// Total: 4 + 11 + 18 + 5 + 23 + 4 + 1 = 66 descriptors

// --- Registration Functions ---
// Need to run ./scripts/update_registry.py after adding DECLARE_DESCRIPTOR macros


void register_NumAtomsDescriptor() {
    auto descriptor = std::make_shared<NumAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumBondsDescriptor() {
    auto descriptor = std::make_shared<NumBondsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumHeavyAtomsGraphDescriptor() {
    auto descriptor = std::make_shared<NumHeavyAtomsGraphDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_EdgeNodeRatioDescriptor() {
    auto descriptor = std::make_shared<EdgeNodeRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_WienerIndexDescriptor() {
    auto descriptor = std::make_shared<WienerIndexDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MeanWienerIndexDescriptor() {
    auto descriptor = std::make_shared<MeanWienerIndexDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_HararyIndexDescriptor() {
    auto descriptor = std::make_shared<HararyIndexDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_BalabanJIndexDescriptor() {
    auto descriptor = std::make_shared<BalabanJIndexDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_GraphDiameterDescriptor() {
    auto descriptor = std::make_shared<GraphDiameterDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_GraphRadiusDescriptor() {
    auto descriptor = std::make_shared<GraphRadiusDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MinAtomEccentricityDescriptor() {
    auto descriptor = std::make_shared<MinAtomEccentricityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MaxAtomEccentricityDescriptor() {
    auto descriptor = std::make_shared<MaxAtomEccentricityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MeanAtomEccentricityDescriptor() {
    auto descriptor = std::make_shared<MeanAtomEccentricityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_VarianceAtomEccentricityDescriptor() {
    auto descriptor = std::make_shared<VarianceAtomEccentricityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_PetitjeanIndexDescriptor() {
    auto descriptor = std::make_shared<PetitjeanIndexDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MinAtomDegreeDescriptor() {
    auto descriptor = std::make_shared<MinAtomDegreeDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MaxAtomDegreeDescriptor() {
    auto descriptor = std::make_shared<MaxAtomDegreeDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MeanAtomDegreeDescriptor() {
    auto descriptor = std::make_shared<MeanAtomDegreeDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_VarianceAtomDegreeDescriptor() {
    auto descriptor = std::make_shared<VarianceAtomDegreeDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumDegreeZeroAtomsDescriptor() {
    auto descriptor = std::make_shared<NumDegreeZeroAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumDegreeOneAtomsDescriptor() {
    auto descriptor = std::make_shared<NumDegreeOneAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumDegreeTwoAtomsDescriptor() {
    auto descriptor = std::make_shared<NumDegreeTwoAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumDegreeThreeAtomsDescriptor() {
    auto descriptor = std::make_shared<NumDegreeThreeAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumDegreeFourAtomsDescriptor() {
    auto descriptor = std::make_shared<NumDegreeFourAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumDegreeFivePlusAtomsDescriptor() {
    auto descriptor = std::make_shared<NumDegreeFivePlusAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionDegreeZeroAtomsDescriptor() {
    auto descriptor = std::make_shared<FractionDegreeZeroAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionDegreeOneAtomsDescriptor() {
    auto descriptor = std::make_shared<FractionDegreeOneAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionDegreeTwoAtomsDescriptor() {
    auto descriptor = std::make_shared<FractionDegreeTwoAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionDegreeThreeAtomsDescriptor() {
    auto descriptor = std::make_shared<FractionDegreeThreeAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionDegreeFourAtomsDescriptor() {
    auto descriptor = std::make_shared<FractionDegreeFourAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionDegreeFivePlusAtomsDescriptor() {
    auto descriptor = std::make_shared<FractionDegreeFivePlusAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_PlattIndexDescriptor() {
    auto descriptor = std::make_shared<PlattIndexDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MeanBondDegreeSumDescriptor() {
    auto descriptor = std::make_shared<MeanBondDegreeSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_VarianceBondDegreeSumDescriptor() {
    auto descriptor = std::make_shared<VarianceBondDegreeSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_ZagrebIndexM1Descriptor() {
    auto descriptor = std::make_shared<ZagrebIndexM1Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_ZagrebIndexM2Descriptor() {
    auto descriptor = std::make_shared<ZagrebIndexM2Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MeanZagrebIndexM1Descriptor() {
    auto descriptor = std::make_shared<MeanZagrebIndexM1Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MeanZagrebIndexM2Descriptor() {
    auto descriptor = std::make_shared<MeanZagrebIndexM2Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_QuadraticIndexDescriptor() {
    auto descriptor = std::make_shared<QuadraticIndexDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_GraphDensityDescriptor() {
    auto descriptor = std::make_shared<GraphDensityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_RandicIndexDescriptor() {
    auto descriptor = std::make_shared<RandicIndexDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_KierKappa1Descriptor() {
    auto descriptor = std::make_shared<KierKappa1Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_KierKappa2Descriptor() {
    auto descriptor = std::make_shared<KierKappa2Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_KierKappa3Descriptor() {
    auto descriptor = std::make_shared<KierKappa3Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_KierAlphaDescriptor() {
    auto descriptor = std::make_shared<KierAlphaDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_KierPhiDescriptor() {
    auto descriptor = std::make_shared<KierPhiDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_CyclomaticNumberDescriptor() {
    auto descriptor = std::make_shared<CyclomaticNumberDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumRingsDescriptor() {
    auto descriptor = std::make_shared<NumRingsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumAliphaticCarbocyclesDescriptor() {
    auto descriptor = std::make_shared<NumAliphaticCarbocyclesDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumAliphaticHeterocyclesDescriptor() {
    auto descriptor = std::make_shared<NumAliphaticHeterocyclesDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumAliphaticRingsDescriptor() {
    auto descriptor = std::make_shared<NumAliphaticRingsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumAromaticCarbocyclesDescriptor() {
    auto descriptor = std::make_shared<NumAromaticCarbocyclesDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumAromaticHeterocyclesDescriptor() {
    auto descriptor = std::make_shared<NumAromaticHeterocyclesDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumAromaticRingsDescriptor() {
    auto descriptor = std::make_shared<NumAromaticRingsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumSaturatedCarbocyclesDescriptor() {
    auto descriptor = std::make_shared<NumSaturatedCarbocyclesDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumSaturatedHeterocyclesDescriptor() {
    auto descriptor = std::make_shared<NumSaturatedHeterocyclesDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumSaturatedRingsDescriptor() {
    auto descriptor = std::make_shared<NumSaturatedRingsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumBridgeheadAtomsDescriptor() {
    auto descriptor = std::make_shared<NumBridgeheadAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumSpiroAtomsDescriptor() {
    auto descriptor = std::make_shared<NumSpiroAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumFusedRingsEstimateDescriptor() {
    auto descriptor = std::make_shared<NumFusedRingsEstimateDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionAtomsInRingsDescriptor() {
    auto descriptor = std::make_shared<FractionAtomsInRingsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionBondsInRingsDescriptor() {
    auto descriptor = std::make_shared<FractionBondsInRingsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MinRingSizeDescriptor() {
    auto descriptor = std::make_shared<MinRingSizeDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MaxRingSizeDescriptor() {
    auto descriptor = std::make_shared<MaxRingSizeDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MeanRingSizeDescriptor() {
    auto descriptor = std::make_shared<MeanRingSizeDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumRingSize3Descriptor() {
    auto descriptor = std::make_shared<NumRingSize3Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumRingSize4Descriptor() {
    auto descriptor = std::make_shared<NumRingSize4Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumRingSize5Descriptor() {
    auto descriptor = std::make_shared<NumRingSize5Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumRingSize6Descriptor() {
    auto descriptor = std::make_shared<NumRingSize6Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumRingSize7Descriptor() {
    auto descriptor = std::make_shared<NumRingSize7Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumRingSize8PlusDescriptor() {
    auto descriptor = std::make_shared<NumRingSize8PlusDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionRingSize3Descriptor() {
    auto descriptor = std::make_shared<FractionRingSize3Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionRingSize4Descriptor() {
    auto descriptor = std::make_shared<FractionRingSize4Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionRingSize5Descriptor() {
    auto descriptor = std::make_shared<FractionRingSize5Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionRingSize6Descriptor() {
    auto descriptor = std::make_shared<FractionRingSize6Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionRingSize7Descriptor() {
    auto descriptor = std::make_shared<FractionRingSize7Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionRingSize8PlusDescriptor() {
    auto descriptor = std::make_shared<FractionRingSize8PlusDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumTerminalNodesDescriptor() {
    auto descriptor = std::make_shared<NumTerminalNodesDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionTerminalNodesDescriptor() {
    auto descriptor = std::make_shared<FractionTerminalNodesDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NumBranchingNodesDescriptor() {
    auto descriptor = std::make_shared<NumBranchingNodesDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionBranchingNodesDescriptor() {
    auto descriptor = std::make_shared<FractionBranchingNodesDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_BertzCTDescriptor() {
    auto descriptor = std::make_shared<BertzCTDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

} // namespace desfact

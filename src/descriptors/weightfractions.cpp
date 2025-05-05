#include "../common.hpp"
#include <GraphMol/ROMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/Descriptors/Lipinski.h> // For HBA/HBD counts
#include <GraphMol/RingInfo.h> // Needed for atom->isInRing replacement
#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <array>
#include <functional>

namespace desfact {

// --- Base Classes for Fraction Descriptors ---

class WeightFractionDescriptor : public Descriptor {
public:
    std::string getCategory() const override { return "WeightFraction"; }
    std::vector<std::string> getDependencies() const override { return {}; }
};

class AtomCountFractionDescriptor : public Descriptor {
public:
    std::string getCategory() const override { return "AtomCountFraction"; }
    std::vector<std::string> getDependencies() const override { return {}; }
};

// --- Fast Weight Calculation ---

namespace {
// Single pass calculation with minimal branching and memory access

// Maximum atomic number we'll handle
constexpr int MAX_Z = 118;

// Pack boolean flags into a single byte per element to improve cache efficiency
struct ElementFlags {
    uint8_t flags = 0;
    enum Flag : uint8_t {
        IS_HALOGEN = 1,
        IS_METAL = 2,
        IS_TRANSITION_METAL = 4,
        IS_ALKALI = 8,
        IS_ALKALINE = 16,
        IS_NOBLE = 32,
        IS_MAIN_GROUP = 64,
        IS_HETEROATOM = 128
    };
};

// Lookup table structure with better memory layout
struct ElementData {
    float weight;          // Using float saves memory and is sufficient for this purpose
    ElementFlags flags;    // Bit-packed flags
    uint8_t enScaled;      // Electronegativity * 100 (scaled for compact storage)
};

// Initialize element data once at program startup - thread safe per C++11
const std::array<ElementData, MAX_Z+1>& getElementData() {
    static const std::array<ElementData, MAX_Z+1> data = []() {
        std::array<ElementData, MAX_Z+1> result = {};
        
        // Get data from RDKit once
        const RDKit::PeriodicTable* pt = RDKit::PeriodicTable::getTable();
        
        for (int z = 1; z <= MAX_Z; ++z) {
            ElementData& elem = result[z];
            
            // Set atomic weight
            elem.weight = static_cast<float>(pt->getAtomicWeight(z));
            
            // Set electronegativity (scaled by 100 to fit in a byte)
            float en = static_cast<float>(element::getElectronegativity(z));
            elem.enScaled = static_cast<uint8_t>(std::min(255.0f, en * 100.0f));
            
            // Set element type flags
            ElementFlags::Flag flags = static_cast<ElementFlags::Flag>(0);
            
            if (z == 9 || z == 17 || z == 35 || z == 53 || z == 85 || z == 117)
                flags = static_cast<ElementFlags::Flag>(flags | ElementFlags::IS_HALOGEN);
                
            if ((z >= 3 && z <= 12) || (z == 11 || z == 19 || z == 37 || z == 55 || z == 87) ||
                (z == 4 || z == 12 || z == 20 || z == 38 || z == 56 || z == 88))
                flags = static_cast<ElementFlags::Flag>(flags | ElementFlags::IS_METAL);
                
            if ((z >= 21 && z <= 30) || (z >= 39 && z <= 48) || (z >= 72 && z <= 80) || (z >= 104 && z <= 112))
                flags = static_cast<ElementFlags::Flag>(flags | ElementFlags::IS_TRANSITION_METAL);
                
            if (z == 3 || z == 11 || z == 19 || z == 37 || z == 55 || z == 87)
                flags = static_cast<ElementFlags::Flag>(flags | ElementFlags::IS_ALKALI);
                
            if (z == 4 || z == 12 || z == 20 || z == 38 || z == 56 || z == 88)
                flags = static_cast<ElementFlags::Flag>(flags | ElementFlags::IS_ALKALINE);
                
            if (z == 2 || z == 10 || z == 18 || z == 36 || z == 54 || z == 86 || z == 118)
                flags = static_cast<ElementFlags::Flag>(flags | ElementFlags::IS_NOBLE);
                
            if ((z >= 1 && z <= 2) || (z >= 3 && z <= 10) || (z >= 11 && z <= 18))
                flags = static_cast<ElementFlags::Flag>(flags | ElementFlags::IS_MAIN_GROUP);
                
            if (z != 6 && z != 1)
                flags = static_cast<ElementFlags::Flag>(flags | ElementFlags::IS_HETEROATOM);
                
            elem.flags.flags = flags;
        }
        
        return result;
    }();
    
    return data;
}

// Packed weight accumulators
struct __attribute__((aligned(64))) WeightSums {
    float total = 0.0f;                // Force cache line alignment
    float weights[34] = {0.0f};        // All weights in a single array for better cache behavior
    
    // Indices for the weights array
    enum Index {
        HEAVY = 0,
        CARBON, HYDROGEN, NITROGEN, OXYGEN, SULFUR, PHOSPHORUS, BORON, SILICON,
        HALOGEN, FLUORINE, CHLORINE, BROMINE, IODINE, METAL,
        AROMATIC, ALIPHATIC, RING, NON_RING, HETEROATOM, MAIN_GROUP, TRANSITION_METAL,
        SP3, SP2, SP, POS, NEG, POLAR, NON_POLAR, ALKALI, ALKALINE, NOBLE, BASIC_N, ACIDIC_O
    };
};

// Thread-local weight data for each thread
thread_local WeightSums tls_weights;

// Fast computation of all weights in a single pass
inline void computeAllWeights(const RDKit::ROMol* mol, WeightSums& sums) {
    if (!mol || mol->getNumAtoms() == 0) return;
    
    // Reset all weights
    sums = WeightSums();
    
    // Pre-fetch element data
    const auto& elemData = getElementData();
    
    // Ensure ring info is available
    mol::ensureRingInfo(mol);
    const auto* ringInfo = mol->getRingInfo();
    
    // Used to avoid branching in tight loops
    alignas(64) float atomWeights[34];  // Cache line aligned
    
    // Process all atoms in a single pass
    for (const auto atom : mol->atoms()) {
        int z = atom->getAtomicNum();
        if (z <= 0 || z > MAX_Z) continue;
        
        // Get cached element data
        const auto& elem = elemData[z];
        float weight = elem.weight;
        
        // Clear temporary array (faster than branching)
        for (int i = 0; i < 34; i++) atomWeights[i] = 0.0f;
        
        // Update basic counters
        sums.total += weight;
        if (z != 1) atomWeights[WeightSums::HEAVY] = weight;
        
        // Handle elements
        switch (z) {
            case 1: atomWeights[WeightSums::HYDROGEN] = weight; break;
            case 6: atomWeights[WeightSums::CARBON] = weight; break;
            case 7: atomWeights[WeightSums::NITROGEN] = weight; break;
            case 8: atomWeights[WeightSums::OXYGEN] = weight; break;
            case 15: atomWeights[WeightSums::PHOSPHORUS] = weight; break;
            case 16: atomWeights[WeightSums::SULFUR] = weight; break;
            case 5: atomWeights[WeightSums::BORON] = weight; break;
            case 14: atomWeights[WeightSums::SILICON] = weight; break;
            case 9: atomWeights[WeightSums::FLUORINE] = weight; atomWeights[WeightSums::HALOGEN] = weight; break;
            case 17: atomWeights[WeightSums::CHLORINE] = weight; atomWeights[WeightSums::HALOGEN] = weight; break;
            case 35: atomWeights[WeightSums::BROMINE] = weight; atomWeights[WeightSums::HALOGEN] = weight; break;
            case 53: atomWeights[WeightSums::IODINE] = weight; atomWeights[WeightSums::HALOGEN] = weight; break;
        }
        
        // Process element type flags (branchless using bit ops)
        uint8_t flags = elem.flags.flags;
        if (flags & ElementFlags::IS_HALOGEN) atomWeights[WeightSums::HALOGEN] = weight;
        if (flags & ElementFlags::IS_METAL) atomWeights[WeightSums::METAL] = weight;
        if (flags & ElementFlags::IS_TRANSITION_METAL) atomWeights[WeightSums::TRANSITION_METAL] = weight;
        if (flags & ElementFlags::IS_ALKALI) atomWeights[WeightSums::ALKALI] = weight;
        if (flags & ElementFlags::IS_ALKALINE) atomWeights[WeightSums::ALKALINE] = weight;
        if (flags & ElementFlags::IS_NOBLE) atomWeights[WeightSums::NOBLE] = weight;
        if (flags & ElementFlags::IS_MAIN_GROUP) atomWeights[WeightSums::MAIN_GROUP] = weight;
        if (flags & ElementFlags::IS_HETEROATOM) atomWeights[WeightSums::HETEROATOM] = weight;
        
        // Process atom properties
        if (atom->getIsAromatic()) {
            atomWeights[WeightSums::AROMATIC] = weight;
        } else {
            atomWeights[WeightSums::ALIPHATIC] = weight;
        }
        
        if (ringInfo->numAtomRings(atom->getIdx()) > 0) {
            atomWeights[WeightSums::RING] = weight;
        } else {
            atomWeights[WeightSums::NON_RING] = weight;
        }
        
        auto hyb = atom->getHybridization();
        if (hyb == RDKit::Atom::SP3) atomWeights[WeightSums::SP3] = weight;
        else if (hyb == RDKit::Atom::SP2) atomWeights[WeightSums::SP2] = weight;
        else if (hyb == RDKit::Atom::SP) atomWeights[WeightSums::SP] = weight;
        
        int charge = atom->getFormalCharge();
        if (charge > 0) atomWeights[WeightSums::POS] = weight;
        else if (charge < 0) atomWeights[WeightSums::NEG] = weight;
        
        if (elem.enScaled > 250) { // EN > 2.5
            atomWeights[WeightSums::POLAR] = weight;
        } else {
            atomWeights[WeightSums::NON_POLAR] = weight;
        }
        
        // These are more expensive checks, only do for relevant atoms
        if (z == 7 && mol::isBasicNitrogen(atom)) atomWeights[WeightSums::BASIC_N] = weight;
        if (z == 8 && mol::isAcidicOxygen(atom)) atomWeights[WeightSums::ACIDIC_O] = weight;
        
        // Vectorized accumulation of weights (compiler can use SIMD)
        for (int i = 0; i < 34; i++) {
            sums.weights[i] += atomWeights[i];
        }
    }
}

// Simple SIMD-friendly fraction calculation
inline void calculateFractions(const WeightSums& sums, float* fractions) {
    // Compute denominator safely
    float denom = (sums.total > 1e-6f) ? sums.total : 1.0f;
    float invDenom = 1.0f / denom;
    
    // Vectorizable loop
    for (int i = 0; i < 34; i++) {
        fractions[i] = sums.weights[i] * invDenom;
    }
}

// Thread-local array for fractions
thread_local float tls_fractions[34];

} // anonymous namespace

// Combined result accumulator for all fractions in a single pass
struct __attribute__((aligned(64))) MoleculeStats {
    // Weight data
    float totalWeight = 0.0f;
    float weights[34] = {0.0f};
    
    // Atom count data
    int totalAtoms = 0;
    int counts[34] = {0};
    
    // Indices (shared between weight and count arrays)
    enum Index {
        HEAVY = 0,
        CARBON, HYDROGEN, NITROGEN, OXYGEN, SULFUR, PHOSPHORUS, BORON, SILICON,
        HALOGEN, FLUORINE, CHLORINE, BROMINE, IODINE, METAL,
        AROMATIC, ALIPHATIC, RING, NON_RING, HETEROATOM, MAIN_GROUP, TRANSITION_METAL,
        SP3, SP2, SP, POS, NEG, POLAR, NON_POLAR, ALKALI, ALKALINE, NOBLE, BASIC_N, ACIDIC_O
    };
};

// Thread-local molecular stats for each thread
thread_local MoleculeStats tls_stats;
thread_local float tls_weightfractions[34];
thread_local float tls_atomfractions[34];

// Single-pass computation of both weights and counts
inline void computeAllStats(const RDKit::ROMol* mol, MoleculeStats& stats) {
    if (!mol || mol->getNumAtoms() == 0) return;
    
    // Reset all counters
    stats = MoleculeStats();
    
    // Pre-fetch element data
    const auto& elemData = getElementData();
    
    // Ensure ring info is available
    mol::ensureRingInfo(mol);
    const auto* ringInfo = mol->getRingInfo();
    
    // Temporary arrays to avoid branching in inner loop
    alignas(64) float atomWeights[34] = {0.0f};
    alignas(64) int atomCounts[34] = {0};
    
    // Process all atoms in a single pass
    for (const auto atom : mol->atoms()) {
        int z = atom->getAtomicNum();
        if (z <= 0 || z > MAX_Z) continue;
        
        // Get cached element data
        const auto& elem = elemData[z];
        float weight = elem.weight;
        
        // Clear temporary arrays (faster than conditional logic)
        for (int i = 0; i < 34; i++) {
            atomWeights[i] = 0.0f;
            atomCounts[i] = 0;
        }
        
        // Update basic counters
        stats.totalWeight += weight;
        stats.totalAtoms++;
        
        if (z != 1) {
            atomWeights[MoleculeStats::HEAVY] = weight;
            atomCounts[MoleculeStats::HEAVY] = 1;
        }
        
        // Handle elements
        switch (z) {
            case 1: 
                atomWeights[MoleculeStats::HYDROGEN] = weight; 
                atomCounts[MoleculeStats::HYDROGEN] = 1;
                break;
            case 6: 
                atomWeights[MoleculeStats::CARBON] = weight; 
                atomCounts[MoleculeStats::CARBON] = 1;
                break;
            case 7: 
                atomWeights[MoleculeStats::NITROGEN] = weight; 
                atomCounts[MoleculeStats::NITROGEN] = 1;
                break;
            case 8: 
                atomWeights[MoleculeStats::OXYGEN] = weight; 
                atomCounts[MoleculeStats::OXYGEN] = 1;
                break;
            case 15: 
                atomWeights[MoleculeStats::PHOSPHORUS] = weight; 
                atomCounts[MoleculeStats::PHOSPHORUS] = 1;
                break;
            case 16: 
                atomWeights[MoleculeStats::SULFUR] = weight; 
                atomCounts[MoleculeStats::SULFUR] = 1;
                break;
            case 5: 
                atomWeights[MoleculeStats::BORON] = weight; 
                atomCounts[MoleculeStats::BORON] = 1;
                break;
            case 14: 
                atomWeights[MoleculeStats::SILICON] = weight; 
                atomCounts[MoleculeStats::SILICON] = 1;
                break;
            case 9: 
                atomWeights[MoleculeStats::FLUORINE] = weight; 
                atomCounts[MoleculeStats::FLUORINE] = 1;
                atomWeights[MoleculeStats::HALOGEN] = weight; 
                atomCounts[MoleculeStats::HALOGEN] = 1;
                break;
            case 17: 
                atomWeights[MoleculeStats::CHLORINE] = weight; 
                atomCounts[MoleculeStats::CHLORINE] = 1;
                atomWeights[MoleculeStats::HALOGEN] = weight; 
                atomCounts[MoleculeStats::HALOGEN] = 1;
                break;
            case 35: 
                atomWeights[MoleculeStats::BROMINE] = weight; 
                atomCounts[MoleculeStats::BROMINE] = 1;
                atomWeights[MoleculeStats::HALOGEN] = weight; 
                atomCounts[MoleculeStats::HALOGEN] = 1;
                break;
            case 53: 
                atomWeights[MoleculeStats::IODINE] = weight; 
                atomCounts[MoleculeStats::IODINE] = 1;
                atomWeights[MoleculeStats::HALOGEN] = weight; 
                atomCounts[MoleculeStats::HALOGEN] = 1;
                break;
        }
        
        // Process element type flags (branchless using bit ops)
        uint8_t flags = elem.flags.flags;
        if (flags & ElementFlags::IS_HALOGEN) {
            atomWeights[MoleculeStats::HALOGEN] = weight;
            atomCounts[MoleculeStats::HALOGEN] = 1;
        }
        if (flags & ElementFlags::IS_METAL) {
            atomWeights[MoleculeStats::METAL] = weight;
            atomCounts[MoleculeStats::METAL] = 1;
        }
        if (flags & ElementFlags::IS_TRANSITION_METAL) {
            atomWeights[MoleculeStats::TRANSITION_METAL] = weight;
            atomCounts[MoleculeStats::TRANSITION_METAL] = 1;
        }
        if (flags & ElementFlags::IS_ALKALI) {
            atomWeights[MoleculeStats::ALKALI] = weight;
            atomCounts[MoleculeStats::ALKALI] = 1;
        }
        if (flags & ElementFlags::IS_ALKALINE) {
            atomWeights[MoleculeStats::ALKALINE] = weight;
            atomCounts[MoleculeStats::ALKALINE] = 1;
        }
        if (flags & ElementFlags::IS_NOBLE) {
            atomWeights[MoleculeStats::NOBLE] = weight;
            atomCounts[MoleculeStats::NOBLE] = 1;
        }
        if (flags & ElementFlags::IS_MAIN_GROUP) {
            atomWeights[MoleculeStats::MAIN_GROUP] = weight;
            atomCounts[MoleculeStats::MAIN_GROUP] = 1;
        }
        if (flags & ElementFlags::IS_HETEROATOM) {
            atomWeights[MoleculeStats::HETEROATOM] = weight;
            atomCounts[MoleculeStats::HETEROATOM] = 1;
        }
        
        // Process atom properties
        if (atom->getIsAromatic()) {
            atomWeights[MoleculeStats::AROMATIC] = weight;
            atomCounts[MoleculeStats::AROMATIC] = 1;
        } else {
            atomWeights[MoleculeStats::ALIPHATIC] = weight;
            atomCounts[MoleculeStats::ALIPHATIC] = 1;
        }
        
        if (ringInfo->numAtomRings(atom->getIdx()) > 0) {
            atomWeights[MoleculeStats::RING] = weight;
            atomCounts[MoleculeStats::RING] = 1;
        } else {
            atomWeights[MoleculeStats::NON_RING] = weight;
            atomCounts[MoleculeStats::NON_RING] = 1;
        }
        
        auto hyb = atom->getHybridization();
        if (hyb == RDKit::Atom::SP3) {
            atomWeights[MoleculeStats::SP3] = weight;
            atomCounts[MoleculeStats::SP3] = 1;
        } else if (hyb == RDKit::Atom::SP2) {
            atomWeights[MoleculeStats::SP2] = weight;
            atomCounts[MoleculeStats::SP2] = 1;
        } else if (hyb == RDKit::Atom::SP) {
            atomWeights[MoleculeStats::SP] = weight;
            atomCounts[MoleculeStats::SP] = 1;
        }
        
        int charge = atom->getFormalCharge();
        if (charge > 0) {
            atomWeights[MoleculeStats::POS] = weight;
            atomCounts[MoleculeStats::POS] = 1;
        } else if (charge < 0) {
            atomWeights[MoleculeStats::NEG] = weight;
            atomCounts[MoleculeStats::NEG] = 1;
        }
        
        if (elem.enScaled > 250) { // EN > 2.5
            atomWeights[MoleculeStats::POLAR] = weight;
            atomCounts[MoleculeStats::POLAR] = 1;
        } else {
            atomWeights[MoleculeStats::NON_POLAR] = weight;
            atomCounts[MoleculeStats::NON_POLAR] = 1;
        }
        
        // These are more expensive checks, only do for relevant atoms
        if (z == 7 && mol::isBasicNitrogen(atom)) {
            atomWeights[MoleculeStats::BASIC_N] = weight;
            atomCounts[MoleculeStats::BASIC_N] = 1;
        }
        if (z == 8 && mol::isAcidicOxygen(atom)) {
            atomWeights[MoleculeStats::ACIDIC_O] = weight;
            atomCounts[MoleculeStats::ACIDIC_O] = 1;
        }
        
        // Vectorized accumulation (compiler can use SIMD)
        for (int i = 0; i < 34; i++) {
            stats.weights[i] += atomWeights[i];
            stats.counts[i] += atomCounts[i];
        }
    }
}

// Calculate both weight and atom fractions in one pass
inline void calculateAllFractions(const MoleculeStats& stats, float* weightFractions, float* atomFractions) {
    float weightDenom = (stats.totalWeight > 1e-6f) ? stats.totalWeight : 1.0f;
    float weightInvDenom = 1.0f / weightDenom;
    
    float atomDenom = (stats.totalAtoms > 0) ? static_cast<float>(stats.totalAtoms) : 1.0f;
    float atomInvDenom = 1.0f / atomDenom;
    
    // SIMD-friendly loops
    for (int i = 0; i < 34; i++) {
        weightFractions[i] = stats.weights[i] * weightInvDenom;
        atomFractions[i] = stats.counts[i] * atomInvDenom;
    }
}

// Base class for all fraction descriptors - handles unified calculation
class FractionDescriptorBase {
public:
    static void ensureCalculated(const RDKit::ROMol* mol) {
        if (!mol) return;
        
        // Get thread-local storage
        MoleculeStats& stats = tls_stats;
        float* wfractions = tls_weightfractions;
        float* afractions = tls_atomfractions;
        
        // Calculate everything in one pass
        computeAllStats(mol, stats);
        calculateAllFractions(stats, wfractions, afractions);
    }
};

// Weight fraction descriptor implementation
class WeightFractionDescriptorFast : public WeightFractionDescriptor, private FractionDescriptorBase {
private:
    int fractionIndex;
    std::string name;
    std::string description;

public:
    WeightFractionDescriptorFast(int idx, const std::string& n, const std::string& desc)
        : fractionIndex(idx), name(n + "Descriptor"), description(desc) {}

    std::string getName() const override { return name; }
    std::string getDescription() const override { return description; }

    DescriptorResult calculate(Context& context) const override {
        const RDKit::ROMol* mol = context.getMolecule();
        if (!mol) return 0.0;
        
        // Calculate all stats and fractions at once (shared across descriptors)
        ensureCalculated(mol);
        
        // Return specific fraction
        return tls_weightfractions[fractionIndex];
    }
};

// Atom count fraction descriptor implementation
class AtomCountFractionDescriptorFast : public AtomCountFractionDescriptor, private FractionDescriptorBase {
private:
    int fractionIndex;
    std::string name;
    std::string description;
public:
    AtomCountFractionDescriptorFast(int idx, const std::string& n, const std::string& desc)
        : fractionIndex(idx), name(n + "Descriptor"), description(desc) {}
    std::string getName() const override { return name; }
    std::string getDescription() const override { return description; }
    DescriptorResult calculate(Context& context) const override {
        const RDKit::ROMol* mol = context.getMolecule();
        if (!mol) return 0.0;
        
        // Calculate all stats and fractions at once (shared across descriptors)
        ensureCalculated(mol);
        
        // Return specific fraction
        return tls_atomfractions[fractionIndex];
    }
};

// Macros for registering optimized descriptors
#define REGISTER_FAST_WEIGHT_FRACTION(NAME, INDEX, DESC) \
void register_##NAME##Descriptor() { \
    auto descriptor = std::make_shared<WeightFractionDescriptorFast>( \
        MoleculeStats::INDEX, #NAME, DESC); \
    auto& registry = DescriptorRegistry::getInstance(); \
    registry.registerDescriptor(std::static_pointer_cast<Descriptor>(descriptor)); \
}

#define REGISTER_FAST_ATOMCOUNT_FRACTION(NAME, INDEX, DESC) \
void register_##NAME##Descriptor() { \
    auto descriptor = std::make_shared<AtomCountFractionDescriptorFast>( \
        MoleculeStats::INDEX, #NAME, DESC); \
    auto& registry = DescriptorRegistry::getInstance(); \
    registry.registerDescriptor(std::static_pointer_cast<Descriptor>(descriptor)); \
}

// Register all descriptors with optimized implementation
REGISTER_FAST_WEIGHT_FRACTION(HeavyAtomWeightFraction, HEAVY, "Fraction of molecular weight from heavy atoms")
REGISTER_FAST_WEIGHT_FRACTION(CarbonWeightFraction, CARBON, "Fraction of molecular weight from carbon atoms")
REGISTER_FAST_WEIGHT_FRACTION(HydrogenWeightFraction, HYDROGEN, "Fraction of molecular weight from hydrogen atoms")
REGISTER_FAST_WEIGHT_FRACTION(NitrogenWeightFraction, NITROGEN, "Fraction of molecular weight from nitrogen atoms")
REGISTER_FAST_WEIGHT_FRACTION(OxygenWeightFraction, OXYGEN, "Fraction of molecular weight from oxygen atoms")
REGISTER_FAST_WEIGHT_FRACTION(SulfurWeightFraction, SULFUR, "Fraction of molecular weight from sulfur atoms")
REGISTER_FAST_WEIGHT_FRACTION(PhosphorusWeightFraction, PHOSPHORUS, "Fraction of molecular weight from phosphorus atoms")
REGISTER_FAST_WEIGHT_FRACTION(BoronWeightFraction, BORON, "Fraction of molecular weight from boron atoms")
REGISTER_FAST_WEIGHT_FRACTION(SiliconWeightFraction, SILICON, "Fraction of molecular weight from silicon atoms")
REGISTER_FAST_WEIGHT_FRACTION(HalogenWeightFraction, HALOGEN, "Fraction of molecular weight from halogen atoms")
REGISTER_FAST_WEIGHT_FRACTION(FluorineWeightFraction, FLUORINE, "Fraction of molecular weight from fluorine atoms")
REGISTER_FAST_WEIGHT_FRACTION(ChlorineWeightFraction, CHLORINE, "Fraction of molecular weight from chlorine atoms")
REGISTER_FAST_WEIGHT_FRACTION(BromineWeightFraction, BROMINE, "Fraction of molecular weight from bromine atoms")
REGISTER_FAST_WEIGHT_FRACTION(IodineWeightFraction, IODINE, "Fraction of molecular weight from iodine atoms")
REGISTER_FAST_WEIGHT_FRACTION(MetalWeightFraction, METAL, "Fraction of molecular weight from metal atoms")
REGISTER_FAST_WEIGHT_FRACTION(AromaticWeightFraction, AROMATIC, "Fraction of molecular weight from aromatic atoms")
REGISTER_FAST_WEIGHT_FRACTION(AliphaticWeightFraction, ALIPHATIC, "Fraction of molecular weight from aliphatic atoms")
REGISTER_FAST_WEIGHT_FRACTION(RingWeightFraction, RING, "Fraction of molecular weight from atoms in rings")
REGISTER_FAST_WEIGHT_FRACTION(NonRingWeightFraction, NON_RING, "Fraction of molecular weight from atoms not in rings")
REGISTER_FAST_WEIGHT_FRACTION(HeteroatomWeightFraction, HETEROATOM, "Fraction of molecular weight from heteroatoms (non C/H)")
REGISTER_FAST_WEIGHT_FRACTION(MainGroupWeightFraction, MAIN_GROUP, "Fraction of molecular weight from main group elements")
REGISTER_FAST_WEIGHT_FRACTION(TransitionMetalWeightFraction, TRANSITION_METAL, "Fraction of molecular weight from transition metals")
REGISTER_FAST_WEIGHT_FRACTION(Sp3WeightFraction, SP3, "Fraction of molecular weight from sp3 hybridized atoms")
REGISTER_FAST_WEIGHT_FRACTION(Sp2WeightFraction, SP2, "Fraction of molecular weight from sp2 hybridized atoms")
REGISTER_FAST_WEIGHT_FRACTION(SpWeightFraction, SP, "Fraction of molecular weight from sp hybridized atoms")
REGISTER_FAST_WEIGHT_FRACTION(PositiveChargeWeightFraction, POS, "Fraction of molecular weight from positively charged atoms")
REGISTER_FAST_WEIGHT_FRACTION(NegativeChargeWeightFraction, NEG, "Fraction of molecular weight from negatively charged atoms")
REGISTER_FAST_WEIGHT_FRACTION(PolarWeightFraction, POLAR, "Fraction of molecular weight from polar atoms (EN > 2.5)")
REGISTER_FAST_WEIGHT_FRACTION(NonPolarWeightFraction, NON_POLAR, "Fraction of molecular weight from non-polar atoms (EN <= 2.5)")
REGISTER_FAST_WEIGHT_FRACTION(AlkaliMetalWeightFraction, ALKALI, "Fraction of molecular weight from alkali metal atoms")
REGISTER_FAST_WEIGHT_FRACTION(AlkalineEarthWeightFraction, ALKALINE, "Fraction of molecular weight from alkaline earth atoms")
REGISTER_FAST_WEIGHT_FRACTION(NobleGasWeightFraction, NOBLE, "Fraction of molecular weight from noble gas atoms")
REGISTER_FAST_WEIGHT_FRACTION(BasicNitrogenWeightFraction, BASIC_N, "Fraction of molecular weight from basic nitrogen atoms")
REGISTER_FAST_WEIGHT_FRACTION(AcidicOxygenWeightFraction, ACIDIC_O, "Fraction of molecular weight from acidic oxygen atoms")

// ----- Begin Custom Advanced Fraction Descriptors -----

// Custom descriptor class for weight fractions using our unified one-pass computation:
class CustomWeightFractionDescriptor : public WeightFractionDescriptor, private FractionDescriptorBase {
public:
    CustomWeightFractionDescriptor(const std::string& n, const std::string& desc, std::function<double(const float*)> fn)
      : name(n), description(desc), computeFn(fn) {}
    std::string getName() const override { return name; }
    std::string getDescription() const override { return description; }
    DescriptorResult calculate(Context& context) const override {
         const RDKit::ROMol* mol = context.getMolecule();
         if (!mol) return 0.0;
         ensureCalculated(mol);
         return computeFn(tls_weightfractions);
    }
private:
    std::string name;
    std::string description;
    std::function<double(const float*)> computeFn;
};

// Custom descriptor class for atom count fractions:
class CustomAtomCountFractionDescriptor : public AtomCountFractionDescriptor, private FractionDescriptorBase {
public:
    CustomAtomCountFractionDescriptor(const std::string& n, const std::string& desc, std::function<double(const float*)> fn)
      : name(n), description(desc), computeFn(fn) {}
    std::string getName() const override { return name; }
    std::string getDescription() const override { return description; }
    DescriptorResult calculate(Context& context) const override {
         const RDKit::ROMol* mol = context.getMolecule();
         if (!mol) return 0.0;
         ensureCalculated(mol);
         return computeFn(tls_atomfractions);
    }
private:
    std::string name;
    std::string description;
    std::function<double(const float*)> computeFn;
};

// Helper lambda for computing standard deviation (population) for a subset of fractions:
inline double stddev(const float* arr, const std::vector<int>& indices) {
    double sum = 0.0, count = 0.0;
    for (int idx : indices) {
       sum += arr[idx];
       count += 1.0;
    }
    if (count == 0) return 0.0;
    double mean = sum / count;
    double sq_sum = 0.0;
    for (int idx : indices) {
       double diff = arr[idx] - mean;
       sq_sum += diff * diff;
    }
    return std::sqrt(sq_sum / count);
}

// Register new descriptors (using a single registration call per descriptor):
void registerCustomWeightAndAtomFractionDescriptors() {
    auto& registry = DescriptorRegistry::getInstance(); // Get registry instance once

    // (1) WeightFractionDiversity: stddev over {CARBON, HYDROGEN, NITROGEN, OXYGEN, SULFUR, PHOSPHORUS}
    {
        auto wfDiversityFn = [](const float* wf) -> double {
            std::vector<int> sel = {1, 2, 3, 4, 5, 6};
            return stddev(wf, sel);
        };
        // Explicitly create std::string objects
        std::string name = "WeightFractionDiversity";
        std::string desc_text = "Standard deviation of weight fractions for C, H, N, O, S, P";
        auto desc = std::make_shared<CustomWeightFractionDescriptor>(name, desc_text, wfDiversityFn);
        registry.registerDescriptor(std::static_pointer_cast<Descriptor>(desc)); // Use static_pointer_cast
    }
    // (2) HeavyToLightWeightRatio: wf[HEAVY] / (1 - wf[HEAVY] + 1e-6)
    {
        auto heavyToLightFn = [](const float* wf) -> double {
            return wf[0] / (1.0 - wf[0] + 1e-6);
        };
        // Explicitly create std::string objects
        std::string name = "HeavyToLightWeightRatio";
        std::string desc_text = "Ratio of heavy atom weight fraction to remainder";
        auto desc = std::make_shared<CustomWeightFractionDescriptor>(name, desc_text, heavyToLightFn);
        registry.registerDescriptor(std::static_pointer_cast<Descriptor>(desc));
    }
    // (3) AromaticToAliphaticWeightRatio: wf[AROMATIC] / (wf[ALIPHATIC] + 1e-6) [indices 15/16]
    {
        auto aromaticToAliphaticFn = [](const float* wf) -> double {
            return wf[15] / (wf[16] + 1e-6);
        };
        // Explicitly create std::string objects
        std::string name = "AromaticToAliphaticWeightRatio";
        std::string desc_text = "Ratio of aromatic to aliphatic weight fractions";
        auto desc = std::make_shared<CustomWeightFractionDescriptor>(name, desc_text, aromaticToAliphaticFn);
        registry.registerDescriptor(std::static_pointer_cast<Descriptor>(desc));
    }
    // (4) MetalToNonMetalWeightRatio: wf[METAL] / (1 - wf[METAL] + 1e-6) [index 14]
    {
        auto metalToNonMetalFn = [](const float* wf) -> double {
            return wf[14] / ((1.0 - wf[14]) + 1e-6);
        };
        // Explicitly create std::string objects
        std::string name = "MetalToNonMetalWeightRatio";
        std::string desc_text = "Ratio of metal weight fraction to non-metal";
        auto desc = std::make_shared<CustomWeightFractionDescriptor>(name, desc_text, metalToNonMetalFn);
        registry.registerDescriptor(std::static_pointer_cast<Descriptor>(desc));
    }
    // (5) TransitionToMainGroupWeightRatio - REMOVED (Zero Variance)
    // (6) PolarToNonPolarWeightRatio: wf[POLAR] / (wf[NON_POLAR] + 1e-6) [27/28]
    {
        auto polarToNonPolarFn = [](const float* wf) -> double {
            return wf[27] / (wf[28] + 1e-6);
        };
        std::string name = "PolarToNonPolarWeightRatio";
        std::string desc_text = "Ratio of polar to non-polar weight fractions";
        auto desc = std::make_shared<CustomWeightFractionDescriptor>(name, desc_text, polarToNonPolarFn);
        registry.registerDescriptor(std::static_pointer_cast<Descriptor>(desc));
    }
    // (7) SpHybridizationWeightRatio: wf[SP] / (wf[SP2] + wf[SP3] + 1e-6) [24,23,22]
    {
        auto spHybridFn = [](const float* wf) -> double {
            return wf[24] / ((wf[23] + wf[22]) + 1e-6);
        };
        std::string name = "SpHybridizationWeightRatio";
        std::string desc_text = "Ratio of sp to (sp2+sp3) weight fractions";
        auto desc = std::make_shared<CustomWeightFractionDescriptor>(name, desc_text, spHybridFn);
        registry.registerDescriptor(std::static_pointer_cast<Descriptor>(desc));
    }
    // (8) ChargeBalanceWeight - Removed
    // (9) NobleToAlkaliWeightRatio - REMOVED (Zero Variance)
    // (10) BasicToAcidicWeightRatio: wf[BASIC_N] / (wf[ACIDIC_O] + 1e-6) [32/33]
    {
        auto basicToAcidicFn = [](const float* wf) -> double {
            return wf[32] / (wf[33] + 1e-6);
        };
        // Explicitly create std::string objects
        std::string name = "BasicToAcidicWeightRatio";
        std::string desc_text = "Ratio of basic nitrogen to acidic oxygen weight fractions";
        auto desc = std::make_shared<CustomWeightFractionDescriptor>(name, desc_text, basicToAcidicFn);
        registry.registerDescriptor(std::static_pointer_cast<Descriptor>(desc));
    }
    // (11) HalogenWeightDiversity: stddev over {FLUORINE, CHLORINE, BROMINE, IODINE} [indices 10,11,12,13]
    {
        auto halogenDiversityFn = [](const float* wf) -> double {
            std::vector<int> sel = {10,11,12,13};
            return stddev(wf, sel);
        };
        // Explicitly create std::string objects
        std::string name = "HalogenWeightDiversity";
        std::string desc_text = "Stddev of individual halogen weight fractions";
        auto desc = std::make_shared<CustomWeightFractionDescriptor>(name, desc_text, halogenDiversityFn);
        registry.registerDescriptor(std::static_pointer_cast<Descriptor>(desc));
    }
    // (12) NonMetalToMetalWeightRatio: (1 - wf[METAL]) / (wf[METAL] + 1e-6) [index 14]
    {
        auto nonMetalToMetalFn = [](const float* wf) -> double {
            return (1.0 - wf[14]) / (wf[14] + 1e-6);
        };
        // Explicitly create std::string objects
        std::string name = "NonMetalToMetalWeightRatio";
        std::string desc_text = "Ratio of non-metal to metal weight fraction";
        auto desc = std::make_shared<CustomWeightFractionDescriptor>(name, desc_text, nonMetalToMetalFn);
        registry.registerDescriptor(std::static_pointer_cast<Descriptor>(desc));
    }
    // (13) SpHybridizationDominance: (wf[SP3] + wf[SP2]) / (wf[SP] + 1e-6) [indices 22,23,24]
    {
        auto spHybridDominanceFn = [](const float* wf) -> double {
            return (wf[22] + wf[23]) / (wf[24] + 1e-6);
        };
        std::string name = "SpHybridizationDominance";
        std::string desc_text = "Ratio of (sp3+sp2) to sp weight fractions";
        auto desc = std::make_shared<CustomWeightFractionDescriptor>(name, desc_text, spHybridDominanceFn);
        registry.registerDescriptor(std::static_pointer_cast<Descriptor>(desc));
    }
    // (14) RingAtomWeightRatio: wf[RING] / (wf[NON_RING] + 1e-6) [indices 17,18]
    {
        auto ringWeightRatioFn = [](const float* wf) -> double {
            return wf[17] / (wf[18] + 1e-6);
        };
        std::string name = "RingAtomWeightRatio";
        std::string desc_text = "Ratio of ring to non-ring weight fractions";
        auto desc = std::make_shared<CustomWeightFractionDescriptor>(name, desc_text, ringWeightRatioFn);
        registry.registerDescriptor(std::static_pointer_cast<Descriptor>(desc));
    }
    // (15) HeteroToCarbonWeightRatio: wf[HETEROATOM] / (wf[CARBON] + 1e-6) [indices 19,1]
    {
        auto heteroToCarbonFn = [](const float* wf) -> double {
            return wf[19] / (wf[1] + 1e-6);
        };
        // Explicitly create std::string objects
        std::string name = "HeteroToCarbonWeightRatio";
        std::string desc_text = "Ratio of heteroatom to carbon weight fractions";
        auto desc = std::make_shared<CustomWeightFractionDescriptor>(name, desc_text, heteroToCarbonFn);
        registry.registerDescriptor(std::static_pointer_cast<Descriptor>(desc));
    }
    // (16) MainGroupDominanceCount: af[MAIN_GROUP] / (af[TRANSITION_METAL] + 1e-6) [indices 20,21]
    {
        auto mgDominanceFn = [](const float* af) -> double {
            return af[20] / (af[21] + 1e-6);
        };
        // Explicitly create std::string objects
        std::string name = "MainGroupDominanceCount";
        std::string desc_text = "Ratio of main group to transition metal atom count fractions";
        auto desc = std::make_shared<CustomAtomCountFractionDescriptor>(name, desc_text, mgDominanceFn);
        registry.registerDescriptor(std::static_pointer_cast<Descriptor>(desc));
    }
    // (17) PolarAtomCountRatio: af[POLAR] / (af[NON_POLAR] + 1e-6) [indices 27,28]
    {
        auto polarAtomCountFn = [](const float* af) -> double {
            return af[27] / (af[28] + 1e-6);
        };
        // Explicitly create std::string objects
        std::string name = "PolarAtomCountRatio";
        std::string desc_text = "Ratio of polar to non-polar atom count fractions";
        auto desc = std::make_shared<CustomAtomCountFractionDescriptor>(name, desc_text, polarAtomCountFn);
        registry.registerDescriptor(std::static_pointer_cast<Descriptor>(desc));
    }
    // (18) ChargedAtomImbalance: af[POS] - af[NEG] [indices 25,26]
    {
        auto chargedAtomImbalanceFn = [](const float* af) -> double {
            return af[25] - af[26];
        };
        // Explicitly create std::string objects
        std::string name = "ChargedAtomImbalance";
        std::string desc_text = "Difference between positive and negative atom count fractions";
        auto desc = std::make_shared<CustomAtomCountFractionDescriptor>(name, desc_text, chargedAtomImbalanceFn);
        registry.registerDescriptor(std::static_pointer_cast<Descriptor>(desc));
    }
    // (19) SpHybridizationCountBalance: af[SP] / (af[SP2] + af[SP3] + 1e-6) [indices 24,23,22]
    {
        auto spHybridCountFn = [](const float* af) -> double {
            return af[24] / ((af[23] + af[22]) + 1e-6);
        };
        // Explicitly create std::string objects
        std::string name = "SpHybridizationCountBalance";
        std::string desc_text = "Ratio of sp to (sp2+sp3) atom count fractions";
        auto desc = std::make_shared<CustomAtomCountFractionDescriptor>(name, desc_text, spHybridCountFn);
        registry.registerDescriptor(std::static_pointer_cast<Descriptor>(desc));
    }
    // (20) ChargeImbalanceWeightRatio (weight-based): abs(wf[POS]-wf[NEG])/(wf[POS]+wf[NEG]+1e-6)
    {
        auto chargeImbalanceWeightFn = [](const float* wf) -> double {
            double pos = wf[25], neg = wf[26];
            return std::abs(pos - neg) / (pos + neg + 1e-6);
        };
        // Explicitly create std::string objects
        std::string name = "ChargeImbalanceWeightRatio";
        std::string desc_text = "Relative imbalance between positive and negative weight fractions";
        auto desc = std::make_shared<CustomWeightFractionDescriptor>(name, desc_text, chargeImbalanceWeightFn);
        registry.registerDescriptor(std::static_pointer_cast<Descriptor>(desc));
    }
}

// Call this function from a single translation unit (e.g., in a registration .cpp or at the end of this file)
struct CustomWeightFractionDescriptorAutoRegister {
    CustomWeightFractionDescriptorAutoRegister() { registerCustomWeightAndAtomFractionDescriptors(); }
};
// static CustomWeightFractionDescriptorAutoRegister _customWeightFractionDescriptorAutoRegister; // Moved registration call

// ----- End Custom Advanced Fraction Descriptors -----

} // namespace desfact

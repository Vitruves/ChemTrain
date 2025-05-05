#include "../common.hpp"
#include <GraphMol/ROMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/AtomIterators.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <map>
#include <limits>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>

namespace desfact {

// --- Thread-local cache for radius-based statistics ---
struct RadiusStatsCache {
    // --- Precomputed Atom Data ---
    std::vector<double> atomicRadii;
    std::vector<double> covalentRadii;
    std::vector<double> vanDerWaalsRadii;
    std::vector<double> bondingRadii;
    std::vector<int> atomicNumbers;
    std::vector<bool> isOrganic;
    std::vector<bool> isHetero;
    std::vector<bool> isAromatic;
    std::vector<bool> inRing;
    
    // --- Basic Stats ---
    uint32_t atomCount = 0;
    uint32_t heavyAtomCount = 0;
    uint32_t organicAtomCount = 0;
    uint32_t aromaticAtomCount = 0;
    uint32_t heteroAtomCount = 0;
    uint32_t ringAtomCount = 0;
    
    // --- Sums ---
    double sumAtomicRadii = 0.0;
    double sumCovalentRadii = 0.0;
    double sumVdwRadii = 0.0;
    double sumBondingRadii = 0.0;
    
    double sumHeavyAtomicRadii = 0.0;
    double sumHeavyCovalentRadii = 0.0;
    double sumHeavyVdwRadii = 0.0;
    double sumHeavyBondingRadii = 0.0;
    
    double sumOrganicAtomicRadii = 0.0;
    double sumOrganicCovalentRadii = 0.0;
    double sumOrganicVdwRadii = 0.0;
    double sumOrganicBondingRadii = 0.0;
    
    double sumHeteroAtomicRadii = 0.0;
    double sumHeteroCovalentRadii = 0.0;
    double sumHeteroVdwRadii = 0.0;
    double sumHeteroBondingRadii = 0.0;
    
    double sumAromaticAtomicRadii = 0.0;
    double sumAromaticCovalentRadii = 0.0;
    double sumAromaticVdwRadii = 0.0;
    double sumAromaticBondingRadii = 0.0;
    
    double sumRingAtomicRadii = 0.0;
    double sumRingCovalentRadii = 0.0;
    double sumRingVdwRadii = 0.0;
    double sumRingBondingRadii = 0.0;
    
    // --- Min/Max/Mean ---
    double minAtomicRadius = std::numeric_limits<double>::max();
    double maxAtomicRadius = 0.0;
    double meanAtomicRadius = 0.0;
    
    double minCovalentRadius = std::numeric_limits<double>::max();
    double maxCovalentRadius = 0.0;
    double meanCovalentRadius = 0.0;
    
    double minVdwRadius = std::numeric_limits<double>::max();
    double maxVdwRadius = 0.0;
    double meanVdwRadius = 0.0;
    
    double minBondingRadius = std::numeric_limits<double>::max();
    double maxBondingRadius = 0.0;
    double meanBondingRadius = 0.0;
    
    // --- Variance/Stddev ---
    double stddevAtomicRadius = 0.0;
    double stddevCovalentRadius = 0.0;
    double stddevVdwRadius = 0.0;
    double stddevBondingRadius = 0.0;
    
    // --- Ratios ---
    double vdwToCovalentRatio = 0.0;
    double atomicToCovalentRatio = 0.0;
    double bondingToAtomicRatio = 0.0;
    double bondingToVdwRatio = 0.0;
    
    // --- Detailed Distributions ---
    std::map<int, int> atomicRadiusBins; // Rounded to nearest 0.1Å
    std::map<int, int> covalentRadiusBins;
    std::map<int, int> vdwRadiusBins;
    
    // --- Most Common Values ---
    double mostCommonAtomicRadius = 0.0;
    double mostCommonCovalentRadius = 0.0;
    double mostCommonVdwRadius = 0.0;
    
    const RDKit::ROMol* populatedForMol = nullptr;
    
    void reset() {
        // Clear vectors
        atomicRadii.clear();
        covalentRadii.clear();
        vanDerWaalsRadii.clear();
        bondingRadii.clear();
        atomicNumbers.clear();
        isOrganic.clear();
        isHetero.clear();
        isAromatic.clear();
        inRing.clear();
        
        // Reset counts
        atomCount = 0;
        heavyAtomCount = 0;
        organicAtomCount = 0;
        aromaticAtomCount = 0;
        heteroAtomCount = 0;
        ringAtomCount = 0;
        
        // Reset sums
        sumAtomicRadii = 0.0;
        sumCovalentRadii = 0.0;
        sumVdwRadii = 0.0;
        sumBondingRadii = 0.0;
        
        sumHeavyAtomicRadii = 0.0;
        sumHeavyCovalentRadii = 0.0;
        sumHeavyVdwRadii = 0.0;
        sumHeavyBondingRadii = 0.0;
        
        sumOrganicAtomicRadii = 0.0;
        sumOrganicCovalentRadii = 0.0;
        sumOrganicVdwRadii = 0.0;
        sumOrganicBondingRadii = 0.0;
        
        sumHeteroAtomicRadii = 0.0;
        sumHeteroCovalentRadii = 0.0;
        sumHeteroVdwRadii = 0.0;
        sumHeteroBondingRadii = 0.0;
        
        sumAromaticAtomicRadii = 0.0;
        sumAromaticCovalentRadii = 0.0;
        sumAromaticVdwRadii = 0.0;
        sumAromaticBondingRadii = 0.0;
        
        sumRingAtomicRadii = 0.0;
        sumRingCovalentRadii = 0.0;
        sumRingVdwRadii = 0.0;
        sumRingBondingRadii = 0.0;
        
        // Reset min/max/mean
        minAtomicRadius = std::numeric_limits<double>::max();
        maxAtomicRadius = 0.0;
        meanAtomicRadius = 0.0;
        
        minCovalentRadius = std::numeric_limits<double>::max();
        maxCovalentRadius = 0.0;
        meanCovalentRadius = 0.0;
        
        minVdwRadius = std::numeric_limits<double>::max();
        maxVdwRadius = 0.0;
        meanVdwRadius = 0.0;
        
        minBondingRadius = std::numeric_limits<double>::max();
        maxBondingRadius = 0.0;
        meanBondingRadius = 0.0;
        
        // Reset variance/stddev
        stddevAtomicRadius = 0.0;
        stddevCovalentRadius = 0.0;
        stddevVdwRadius = 0.0;
        stddevBondingRadius = 0.0;
        
        // Reset ratios
        vdwToCovalentRatio = 0.0;
        atomicToCovalentRatio = 0.0;
        bondingToAtomicRatio = 0.0;
        bondingToVdwRatio = 0.0;
        
        // Reset distribution bins
        atomicRadiusBins.clear();
        covalentRadiusBins.clear();
        vdwRadiusBins.clear();
        
        // Reset most common
        mostCommonAtomicRadius = 0.0;
        mostCommonCovalentRadius = 0.0;
        mostCommonVdwRadius = 0.0;
        
        populatedForMol = nullptr;
    }
    
    void computeFromMol(const RDKit::ROMol* mol) {
        if (!mol) return;
        reset();
        populatedForMol = mol;
        
        atomCount = mol->getNumAtoms();
        if (atomCount == 0) return;
        
        // Pre-allocate vectors
        atomicRadii.reserve(atomCount);
        covalentRadii.reserve(atomCount);
        vanDerWaalsRadii.reserve(atomCount);
        bondingRadii.reserve(atomCount);
        atomicNumbers.reserve(atomCount);
        isOrganic.reserve(atomCount);
        isHetero.reserve(atomCount);
        isAromatic.reserve(atomCount);
        inRing.reserve(atomCount);
        
        // Process all atoms in parallel
        tbb::concurrent_vector<std::tuple<double, double, double, double, int, bool, bool, bool, bool>> atomData(atomCount);
        
        tbb::parallel_for(tbb::blocked_range<size_t>(0, atomCount),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    const RDKit::Atom* atom = mol->getAtomWithIdx(i);
                    int atomicNum = atom->getAtomicNum();
                    
                    // Get all radii in one pass - more efficient
                    double atomicRad = element::getAtomicMass(atomicNum) / 10.0; // Approximate atomic radius from mass
                    double covalentRad = element::getCovalentRadius(atomicNum);
                    double vdwRad = element::getCovalentRadius(atomicNum) * 1.4; // Approximate vdW radius
                    double bondingRad = covalentRad * 0.85; // Simplified bonding radius model
                    
                    bool organic = (atomicNum == 6 || atomicNum == 7 || atomicNum == 8 || 
                                    atomicNum == 9 || atomicNum == 15 || atomicNum == 16 ||
                                    atomicNum == 17 || atomicNum == 35 || atomicNum == 53);
                    bool hetero = (atomicNum != 6 && atomicNum != 1); // Heteroatoms are atoms that are not C or H
                    bool aromatic = atom->getIsAromatic();
                    bool ring = mol->getRingInfo()->numAtomRings(i) > 0;
                    
                    atomData[i] = std::make_tuple(atomicRad, covalentRad, vdwRad, bondingRad, 
                                                 atomicNum, organic, hetero, aromatic, ring);
                }
            });
        
        // Process collected data sequentially
        for (const auto& data : atomData) {
            double atomicRad, covalentRad, vdwRad, bondingRad;
            int atomicNum;
            bool organic, hetero, aromatic, ring;
            
            std::tie(atomicRad, covalentRad, vdwRad, bondingRad, atomicNum, 
                     organic, hetero, aromatic, ring) = data;
            
            // Store values
            atomicRadii.push_back(atomicRad);
            covalentRadii.push_back(covalentRad);
            vanDerWaalsRadii.push_back(vdwRad);
            bondingRadii.push_back(bondingRad);
            atomicNumbers.push_back(atomicNum);
            isOrganic.push_back(organic);
            isHetero.push_back(hetero);
            isAromatic.push_back(aromatic);
            inRing.push_back(ring);
            
            // Update sums
            sumAtomicRadii += atomicRad;
            sumCovalentRadii += covalentRad;
            sumVdwRadii += vdwRad;
            sumBondingRadii += bondingRad;
            
            // Update min/max
            minAtomicRadius = std::min(minAtomicRadius, atomicRad);
            maxAtomicRadius = std::max(maxAtomicRadius, atomicRad);
            
            minCovalentRadius = std::min(minCovalentRadius, covalentRad);
            maxCovalentRadius = std::max(maxCovalentRadius, covalentRad);
            
            minVdwRadius = std::min(minVdwRadius, vdwRad);
            maxVdwRadius = std::max(maxVdwRadius, vdwRad);
            
            minBondingRadius = std::min(minBondingRadius, bondingRad);
            maxBondingRadius = std::max(maxBondingRadius, bondingRad);
            
            // Update distribution bins (rounded to nearest 0.1Å)
            atomicRadiusBins[static_cast<int>(std::round(atomicRad * 10.0))]++;
            covalentRadiusBins[static_cast<int>(std::round(covalentRad * 10.0))]++;
            vdwRadiusBins[static_cast<int>(std::round(vdwRad * 10.0))]++;
            
            // Update counters and specialized sums
            if (atomicNum > 1) {
                heavyAtomCount++;
                sumHeavyAtomicRadii += atomicRad;
                sumHeavyCovalentRadii += covalentRad;
                sumHeavyVdwRadii += vdwRad;
                sumHeavyBondingRadii += bondingRad;
            }
            
            if (organic) {
                organicAtomCount++;
                sumOrganicAtomicRadii += atomicRad;
                sumOrganicCovalentRadii += covalentRad;
                sumOrganicVdwRadii += vdwRad;
                sumOrganicBondingRadii += bondingRad;
            }
            
            if (hetero) {
                heteroAtomCount++;
                sumHeteroAtomicRadii += atomicRad;
                sumHeteroCovalentRadii += covalentRad;
                sumHeteroVdwRadii += vdwRad;
                sumHeteroBondingRadii += bondingRad;
            }
            
            if (aromatic) {
                aromaticAtomCount++;
                sumAromaticAtomicRadii += atomicRad;
                sumAromaticCovalentRadii += covalentRad;
                sumAromaticVdwRadii += vdwRad;
                sumAromaticBondingRadii += bondingRad;
            }
            
            if (ring) {
                ringAtomCount++;
                sumRingAtomicRadii += atomicRad;
                sumRingCovalentRadii += covalentRad;
                sumRingVdwRadii += vdwRad;
                sumRingBondingRadii += bondingRad;
            }
        }
        
        // Calculate means
        meanAtomicRadius = atomCount > 0 ? sumAtomicRadii / atomCount : 0.0;
        meanCovalentRadius = atomCount > 0 ? sumCovalentRadii / atomCount : 0.0;
        meanVdwRadius = atomCount > 0 ? sumVdwRadii / atomCount : 0.0;
        meanBondingRadius = atomCount > 0 ? sumBondingRadii / atomCount : 0.0;
        
        // Calculate standard deviations
        double varAtomicRadius = 0.0;
        double varCovalentRadius = 0.0;
        double varVdwRadius = 0.0;
        double varBondingRadius = 0.0;
        
        for (size_t i = 0; i < atomCount; ++i) {
            varAtomicRadius += std::pow(atomicRadii[i] - meanAtomicRadius, 2);
            varCovalentRadius += std::pow(covalentRadii[i] - meanCovalentRadius, 2);
            varVdwRadius += std::pow(vanDerWaalsRadii[i] - meanVdwRadius, 2);
            varBondingRadius += std::pow(bondingRadii[i] - meanBondingRadius, 2);
        }
        
        stddevAtomicRadius = atomCount > 0 ? std::sqrt(varAtomicRadius / atomCount) : 0.0;
        stddevCovalentRadius = atomCount > 0 ? std::sqrt(varCovalentRadius / atomCount) : 0.0;
        stddevVdwRadius = atomCount > 0 ? std::sqrt(varVdwRadius / atomCount) : 0.0;
        stddevBondingRadius = atomCount > 0 ? std::sqrt(varBondingRadius / atomCount) : 0.0;
        
        // Calculate global ratios
        vdwToCovalentRatio = meanCovalentRadius > 0 ? meanVdwRadius / meanCovalentRadius : 0.0;
        atomicToCovalentRatio = meanCovalentRadius > 0 ? meanAtomicRadius / meanCovalentRadius : 0.0;
        bondingToAtomicRatio = meanAtomicRadius > 0 ? meanBondingRadius / meanAtomicRadius : 0.0;
        bondingToVdwRatio = meanVdwRadius > 0 ? meanBondingRadius / meanVdwRadius : 0.0;
        
        // Find most common radius values
        int maxAtomicCount = 0;
        int maxCovalentCount = 0;
        int maxVdwCount = 0;
        
        for (const auto& pair : atomicRadiusBins) {
            if (pair.second > maxAtomicCount) {
                maxAtomicCount = pair.second;
                mostCommonAtomicRadius = pair.first / 10.0;
            }
        }
        
        for (const auto& pair : covalentRadiusBins) {
            if (pair.second > maxCovalentCount) {
                maxCovalentCount = pair.second;
                mostCommonCovalentRadius = pair.first / 10.0;
            }
        }
        
        for (const auto& pair : vdwRadiusBins) {
            if (pair.second > maxVdwCount) {
                maxVdwCount = pair.second;
                mostCommonVdwRadius = pair.first / 10.0;
            }
        }
    }
};

static thread_local RadiusStatsCache tlsRadiusCache;

// Helper to ensure cache is populated
inline void ensureCachePopulated(Context& context) {
    const RDKit::ROMol* mol = context.getMolecule();
    if (tlsRadiusCache.populatedForMol != mol || !mol) {
        tlsRadiusCache.computeFromMol(mol);
    }
}

// --- Base Class ---
class RadiusStatsDescriptor : public Descriptor {
public:
    using Descriptor::Descriptor;
    std::string getCategory() const override { return "RadiusStats"; }
};

// --- Descriptor Implementations ---

// 1. Sum of atomic radii
// 2. Sum of covalent radii
// 3. Sum of van der Waals radii
// 4. Sum of bonding radii
// 5. Mean atomic radius
// 6. Mean covalent radius
// 7. Mean van der Waals radius
// 8. Mean bonding radius
// 9. Maximum atomic radius
// 10. Maximum covalent radius
// 11. Maximum van der Waals radius
// 12. Minimum atomic radius
// 13. Minimum covalent radius
// 14. Minimum van der Waals radius
// 15. Standard deviation of atomic radii
// 16. Standard deviation of covalent radii
// 17. Standard deviation of van der Waals radii
// 18. Ratio of vdW to covalent radius
// 19. Ratio of atomic to covalent radius
// 20. Sum of heavy atom covalent radii
// 21. Sum of heavy atom vdW radii
// 22. Mean heavy atom covalent radius
// 23. Mean heavy atom vdW radius

// 24. Sum of heteroatom covalent radii
DECLARE_DESCRIPTOR(SumHeteroCovalentRadii, RadiusStatsDescriptor, "Sum of heteroatom covalent radii")
DESCRIPTOR_DEPENDENCIES(SumHeteroCovalentRadii) { return {}; }
DescriptorResult SumHeteroCovalentRadiiDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsRadiusCache.sumHeteroCovalentRadii;
}

// 25. Sum of heteroatom vdW radii
// 26. Mean heteroatom covalent radius
// 27. Mean heteroatom vdW radius
// 28. Sum of aromatic atom covalent radii
// 29. Sum of aromatic atom vdW radii
// 30. Mean aromatic atom covalent radius
// 31. Mean aromatic atom vdW radius
// 32. Sum of ring atom covalent radii
// 33. Sum of ring atom vdW radii
// 34. Mean ring atom covalent radius
// 35. Mean ring atom vdW radius

// 36. Fraction of vdW radius sum from heteroatoms
DECLARE_DESCRIPTOR(FractionHeteroVdwRadii, RadiusStatsDescriptor, "Fraction of vdW radius sum from heteroatoms")
DESCRIPTOR_DEPENDENCIES(FractionHeteroVdwRadii) { return {}; }
DescriptorResult FractionHeteroVdwRadiiDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsRadiusCache.sumVdwRadii > 0 ? tlsRadiusCache.sumHeteroVdwRadii / tlsRadiusCache.sumVdwRadii : 0.0;
}

// 37. Fraction of covalent radius sum from heteroatoms
DECLARE_DESCRIPTOR(FractionHeteroCovalentRadii, RadiusStatsDescriptor, "Fraction of covalent radius sum from heteroatoms")
DESCRIPTOR_DEPENDENCIES(FractionHeteroCovalentRadii) { return {}; }
DescriptorResult FractionHeteroCovalentRadiiDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsRadiusCache.sumCovalentRadii > 0 ? tlsRadiusCache.sumHeteroCovalentRadii / tlsRadiusCache.sumCovalentRadii : 0.0;
}

// 38. Fraction of vdW radius sum from aromatic atoms
DECLARE_DESCRIPTOR(FractionAromaticVdwRadii, RadiusStatsDescriptor, "Fraction of vdW radius sum from aromatic atoms")
DESCRIPTOR_DEPENDENCIES(FractionAromaticVdwRadii) { return {}; }
DescriptorResult FractionAromaticVdwRadiiDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsRadiusCache.sumVdwRadii > 0 ? tlsRadiusCache.sumAromaticVdwRadii / tlsRadiusCache.sumVdwRadii : 0.0;
}

// 39. Fraction of covalent radius sum from aromatic atoms
DECLARE_DESCRIPTOR(FractionAromaticCovalentRadii, RadiusStatsDescriptor, "Fraction of covalent radius sum from aromatic atoms")
DESCRIPTOR_DEPENDENCIES(FractionAromaticCovalentRadii) { return {}; }
DescriptorResult FractionAromaticCovalentRadiiDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsRadiusCache.sumCovalentRadii > 0 ? tlsRadiusCache.sumAromaticCovalentRadii / tlsRadiusCache.sumCovalentRadii : 0.0;
}

// 40. Fraction of vdW radius sum from ring atoms
DECLARE_DESCRIPTOR(FractionRingVdwRadii, RadiusStatsDescriptor, "Fraction of vdW radius sum from ring atoms")
DESCRIPTOR_DEPENDENCIES(FractionRingVdwRadii) { return {}; }
DescriptorResult FractionRingVdwRadiiDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsRadiusCache.sumVdwRadii > 0 ? tlsRadiusCache.sumRingVdwRadii / tlsRadiusCache.sumVdwRadii : 0.0;
}

// 41. Fraction of covalent radius sum from ring atoms
DECLARE_DESCRIPTOR(FractionRingCovalentRadii, RadiusStatsDescriptor, "Fraction of covalent radius sum from ring atoms")
DESCRIPTOR_DEPENDENCIES(FractionRingCovalentRadii) { return {}; }
DescriptorResult FractionRingCovalentRadiiDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsRadiusCache.sumCovalentRadii > 0 ? tlsRadiusCache.sumRingCovalentRadii / tlsRadiusCache.sumCovalentRadii : 0.0;
}

// 42. Most common atomic radius
// 43. Most common covalent radius
// 44. Most common vdW radius
// 45. Ratio of max to min atomic radius
// 46. Ratio of max to min covalent radius
// 47. Ratio of max to min vdW radius

// 48. Bonding to atomic radius ratio
DECLARE_DESCRIPTOR(BondingToAtomicRadiusRatio, RadiusStatsDescriptor, "Ratio of mean bonding radius to mean atomic radius")
DESCRIPTOR_DEPENDENCIES(BondingToAtomicRadiusRatio) { return {}; }
DescriptorResult BondingToAtomicRadiusRatioDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsRadiusCache.bondingToAtomicRatio;
}

// 49. Bonding to vdW radius ratio
// 50. Mean heavy to mean all atom vdW ratio

void register_SumHeteroCovalentRadiiDescriptor() {
    auto descriptor = std::make_shared<SumHeteroCovalentRadiiDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_FractionHeteroVdwRadiiDescriptor() {
    auto descriptor = std::make_shared<FractionHeteroVdwRadiiDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_FractionHeteroCovalentRadiiDescriptor() {
    auto descriptor = std::make_shared<FractionHeteroCovalentRadiiDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_FractionAromaticVdwRadiiDescriptor() {
    auto descriptor = std::make_shared<FractionAromaticVdwRadiiDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_FractionAromaticCovalentRadiiDescriptor() {
    auto descriptor = std::make_shared<FractionAromaticCovalentRadiiDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_FractionRingVdwRadiiDescriptor() {
    auto descriptor = std::make_shared<FractionRingVdwRadiiDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_FractionRingCovalentRadiiDescriptor() {
    auto descriptor = std::make_shared<FractionRingCovalentRadiiDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_BondingToAtomicRadiusRatioDescriptor() {
    auto descriptor = std::make_shared<BondingToAtomicRadiusRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}
} // namespace desfact
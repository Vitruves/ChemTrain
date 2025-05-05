#include "../common.hpp"
#include <GraphMol/ROMol.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Atom.h>
#include <GraphMol/RingInfo.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <limits>  // For numeric_limits
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_scan.h>
#include <tbb/cache_aligned_allocator.h>
#include <tbb/task_arena.h>
#include <tbb/partitioner.h>

namespace desfact {

// Cache-aligned vectors for better memory access
template<typename T>
using aligned_vector = std::vector<T, tbb::cache_aligned_allocator<T>>;

// Helper: get bond polarity (difference in Pauling EN) - uses precomputed ENs
inline double bondPolarityPrecomputed(int an1, int an2) {
    double en1 = element::getElectronegativity(an1);
    double en2 = element::getElectronegativity(an2);
    return std::abs(en1 - en2);
}

// Thread-local counter structure to reduce atomic contention
struct alignas(64) BondCounters {  // Ensure 64-byte alignment for cache lines
    uint32_t singleBonds = 0;
    uint32_t doubleBonds = 0;
    uint32_t tripleBonds = 0;
    uint32_t aromaticBonds = 0;
    uint32_t rotatableBonds = 0;
    uint32_t bondsToHetero = 0;
    uint32_t bondsBetweenHetero = 0;
    uint32_t ccBonds = 0;
    uint32_t bondsToHydrogen = 0;
    uint32_t heavyAtomBonds = 0;
    uint32_t ringBonds = 0;
    uint32_t singleRingBonds = 0;
    uint32_t doubleRingBonds = 0;
    uint32_t rotatableRingBonds = 0;
    uint32_t bondsWithAromaticAtom = 0;
    uint32_t bondsBothAromaticAtoms = 0;
    uint32_t bondsWithSp3 = 0;
    uint32_t bondsBothSp3 = 0;
    uint32_t bondsWithSp2 = 0;
    uint32_t bondsBothSp2 = 0;
    uint32_t bondsWithSp = 0;
    uint32_t bondsBothSp = 0;
    uint32_t bondsWithHalogen = 0;
    uint32_t bondsBothHalogen = 0;
    uint32_t bondsWithMetal = 0;
    uint32_t bondsBothMetal = 0;
    uint32_t bondsWithRingAtom = 0;
    uint32_t bondsBothRingAtoms = 0;
    uint32_t bondsRingToNonRing = 0;

    void merge(const BondCounters& other) {
        singleBonds += other.singleBonds;
        doubleBonds += other.doubleBonds;
        tripleBonds += other.tripleBonds;
        aromaticBonds += other.aromaticBonds;
        rotatableBonds += other.rotatableBonds;
        bondsToHetero += other.bondsToHetero;
        bondsBetweenHetero += other.bondsBetweenHetero;
        ccBonds += other.ccBonds;
        bondsToHydrogen += other.bondsToHydrogen;
        heavyAtomBonds += other.heavyAtomBonds;
        ringBonds += other.ringBonds;
        singleRingBonds += other.singleRingBonds;
        doubleRingBonds += other.doubleRingBonds;
        rotatableRingBonds += other.rotatableRingBonds;
        bondsWithAromaticAtom += other.bondsWithAromaticAtom;
        bondsBothAromaticAtoms += other.bondsBothAromaticAtoms;
        bondsWithSp3 += other.bondsWithSp3;
        bondsBothSp3 += other.bondsBothSp3;
        bondsWithSp2 += other.bondsWithSp2;
        bondsBothSp2 += other.bondsBothSp2;
        bondsWithSp += other.bondsWithSp;
        bondsBothSp += other.bondsBothSp;
        bondsWithHalogen += other.bondsWithHalogen;
        bondsBothHalogen += other.bondsBothHalogen;
        bondsWithMetal += other.bondsWithMetal;
        bondsBothMetal += other.bondsBothMetal;
        bondsWithRingAtom += other.bondsWithRingAtom;
        bondsBothRingAtoms += other.bondsBothRingAtoms;
        bondsRingToNonRing += other.bondsRingToNonRing;
    }
};

// --- Thread-local cache for bond statistics ---
struct BondStatsCache {
    // --- Counts ---
    uint32_t numBonds = 0;
    uint32_t singleBonds = 0;
    uint32_t doubleBonds = 0;
    uint32_t tripleBonds = 0;
    uint32_t aromaticBonds = 0;
    uint32_t rotatableBonds = 0;
    uint32_t bondsToHetero = 0;
    uint32_t bondsBetweenHetero = 0;
    uint32_t ccBonds = 0;
    uint32_t bondsToHydrogen = 0;
    uint32_t heavyAtomBonds = 0;
    uint32_t ringBonds = 0;
    uint32_t singleRingBonds = 0;
    uint32_t doubleRingBonds = 0;
    uint32_t rotatableRingBonds = 0;
    uint32_t bondsWithAromaticAtom = 0;
    uint32_t bondsBothAromaticAtoms = 0;
    uint32_t bondsWithSp3 = 0;
    uint32_t bondsBothSp3 = 0;
    uint32_t bondsWithSp2 = 0;
    uint32_t bondsBothSp2 = 0;
    uint32_t bondsWithSp = 0;
    uint32_t bondsBothSp = 0;
    uint32_t bondsWithHalogen = 0;
    uint32_t bondsBothHalogen = 0;
    uint32_t bondsWithMetal = 0;
    uint32_t bondsBothMetal = 0;
    uint32_t bondsWithRingAtom = 0;
    uint32_t bondsBothRingAtoms = 0;
    uint32_t bondsRingToNonRing = 0;

    // --- Data for Stats ---
    aligned_vector<double> bondOrders;
    aligned_vector<double> bondPolarities;
    aligned_vector<int> bondAtomicNumDiffs;
    tbb::concurrent_unordered_map<double, int> bondOrderCounts;
    tbb::concurrent_unordered_map<int, int> bondPolarityCounts; // Key is polarity * 10, rounded

    // --- Calculated Stats ---
    double meanBondPolarity = 0.0;
    double maxBondPolarity = 0.0;
    double minBondPolarity = std::numeric_limits<double>::max(); // Initialize min properly
    double stddevBondPolarity = 0.0;
    double meanBondOrder = 0.0;
    double maxBondOrder = 0.0;
    double minBondOrder = std::numeric_limits<double>::max(); // Initialize min properly
    double stddevBondOrder = 0.0;
    double meanAtomicNumDiff = 0.0;
    double stddevAtomicNumDiff = 0.0;
    double mostCommonBondOrder = 0.0;
    double mostCommonBondPolarity = 0.0; // Value is rounded polarity / 10.0
    int uniqueBondOrders = 0;
    int uniqueBondPolarities = 0;

    // --- Precomputed Data (per molecule) ---
    aligned_vector<int> atomNums;
    aligned_vector<RDKit::Atom::HybridizationType> atomHybs;
    aligned_vector<bool> atomIsHetero;
    aligned_vector<bool> atomIsHeavy;
    aligned_vector<bool> atomInRing;
    aligned_vector<bool> atomIsAromatic;
    aligned_vector<bool> atomIsHalogen;
    aligned_vector<bool> atomIsMetal;
    aligned_vector<bool> bondIsRotatable;

    const RDKit::ROMol* populatedForMol = nullptr;

    void reset() {
        // Reset counts
        numBonds = 0; singleBonds = 0; doubleBonds = 0; tripleBonds = 0;
        aromaticBonds = 0; rotatableBonds = 0; bondsToHetero = 0;
        bondsBetweenHetero = 0; ccBonds = 0; bondsToHydrogen = 0;
        heavyAtomBonds = 0; ringBonds = 0; singleRingBonds = 0;
        doubleRingBonds = 0; rotatableRingBonds = 0; bondsWithAromaticAtom = 0;
        bondsBothAromaticAtoms = 0; bondsWithSp3 = 0; bondsBothSp3 = 0;
        bondsWithSp2 = 0; bondsBothSp2 = 0; bondsWithSp = 0; bondsBothSp = 0;
        bondsWithHalogen = 0; bondsBothHalogen = 0; bondsWithMetal = 0;
        bondsBothMetal = 0; bondsWithRingAtom = 0; bondsBothRingAtoms = 0;
        bondsRingToNonRing = 0;

        // Reset data vectors/maps
        bondOrders.clear();
        bondPolarities.clear();
        bondAtomicNumDiffs.clear();
        bondOrderCounts.clear();
        bondPolarityCounts.clear();

        // Reset calculated stats
        meanBondPolarity = 0.0; maxBondPolarity = 0.0;
        minBondPolarity = std::numeric_limits<double>::max(); // Re-initialize min
        stddevBondPolarity = 0.0; meanBondOrder = 0.0; maxBondOrder = 0.0;
        minBondOrder = std::numeric_limits<double>::max(); // Re-initialize min
        stddevBondOrder = 0.0; meanAtomicNumDiff = 0.0;
        stddevAtomicNumDiff = 0.0; mostCommonBondOrder = 0.0;
        mostCommonBondPolarity = 0.0; uniqueBondOrders = 0;
        uniqueBondPolarities = 0;

        // Clear precomputed data vectors
        atomNums.clear(); atomHybs.clear(); atomIsHetero.clear();
        atomIsHeavy.clear(); atomInRing.clear(); atomIsAromatic.clear();
        atomIsHalogen.clear(); atomIsMetal.clear(); bondIsRotatable.clear();

        populatedForMol = nullptr;
    }

    void computeFromMol(const RDKit::ROMol* mol) {
        if (!mol) return;
        reset();
        populatedForMol = mol;

        numBonds = mol->getNumBonds();
        uint32_t numAtoms = mol->getNumAtoms();
        if (numBonds == 0) return;

        // Pre-allocate with proper alignment and RESIZE (not just reserve)
        const size_t alignment = 64;  // Cache line size
        size_t alignedAtomSize = (numAtoms + alignment - 1) & ~(alignment - 1);
        size_t alignedBondSize = (numBonds + alignment - 1) & ~(alignment - 1);
        
        atomNums.resize(alignedAtomSize);
        atomHybs.resize(alignedAtomSize);
        atomIsHetero.resize(alignedAtomSize);
        atomIsHeavy.resize(alignedAtomSize);
        atomInRing.resize(alignedAtomSize);
        atomIsAromatic.resize(alignedAtomSize);
        atomIsHalogen.resize(alignedAtomSize);
        atomIsMetal.resize(alignedAtomSize);
        bondIsRotatable.resize(alignedBondSize);

        auto ringInfo = mol->getRingInfo();

        // Process atoms in cache-friendly chunks
        constexpr size_t ATOMS_PER_CHUNK = 4;  // Small chunk size for better load balancing
        
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, numAtoms, ATOMS_PER_CHUNK),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    const RDKit::Atom* atom = mol->getAtomWithIdx(i);
                    int an = atom->getAtomicNum();
                    
                    // Write to aligned locations
                    size_t aligned_idx = (i + alignment - 1) & ~(alignment - 1);
                    if (aligned_idx < alignedAtomSize) {  // Bounds check
                        atomNums[i] = an;  // Use original index for actual data
                        atomHybs[i] = atom->getHybridization();
                        atomIsHetero[i] = element::isHeteroatom(an);
                        atomIsHeavy[i] = an > 1;
                        atomInRing[i] = ringInfo->numAtomRings(i) > 0;
                        atomIsAromatic[i] = atom->getIsAromatic();
                        atomIsHalogen[i] = element::isHalogen(an);
                        atomIsMetal[i] = element::isMetal(an);
                    }
                }
            },
            tbb::static_partitioner()
        );

        // Process bonds in cache-friendly chunks
        constexpr size_t BONDS_PER_CHUNK = 4;
        
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, numBonds, BONDS_PER_CHUNK),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    const RDKit::Bond* bond = mol->getBondWithIdx(i);
                    size_t aligned_idx = (i + alignment - 1) & ~(alignment - 1);
                    if (aligned_idx < alignedBondSize) {  // Bounds check
                        bondIsRotatable[i] = mol::isRotatableBond(bond);  // Use original index
                    }
                }
            },
            tbb::static_partitioner()
        );

        // Use thread-local storage with proper alignment
        using TLS = tbb::enumerable_thread_specific<BondCounters, tbb::cache_aligned_allocator<BondCounters>>;
        TLS tls_counters;

        // Bond property structure with proper alignment
        struct alignas(64) BondProperty {
            double order;
            double polarity;
            int atomicNumDiff;
            int roundedPolarity;
            bool isSingle;
            bool isDouble;
            bool isTriple;
            bool isAromatic;
            bool toHetero;
            bool betweenHetero;
            bool isCC;
            bool toHydrogen;
            bool betweenHeavy;
            bool inRing;
            bool singleInRing;
            bool doubleInRing;
            bool rotatableInRing;
            bool withAromaticAtom;
            bool bothAromaticAtoms;
            bool withSp3;
            bool bothSp3;
            bool withSp2;
            bool bothSp2;
            bool withSp;
            bool bothSp;
            bool withHalogen;
            bool bothHalogen;
            bool withMetal;
            bool bothMetal;
            bool withRingAtom;
            bool bothRingAtoms;
            bool ringToNonRing;
        };

        // Allocate bond properties vector with cache alignment
        tbb::concurrent_vector<BondProperty, tbb::cache_aligned_allocator<BondProperty>> bondProps;
        bondProps.reserve(alignedBondSize);  // Pre-allocate with alignment

        // Process bonds in optimized chunks
        tbb::task_arena arena(tbb::task_arena::automatic);
        arena.execute([&] {
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, numBonds, BONDS_PER_CHUNK),
                [&](const tbb::blocked_range<size_t>& range) {
                    auto& local_counters = tls_counters.local();
                    
                    for (size_t i = range.begin(); i < range.end(); ++i) {
                        const RDKit::Bond* bond = mol->getBondWithIdx(i);
                        int idx1 = bond->getBeginAtomIdx();
                        int idx2 = bond->getEndAtomIdx();
                        
                        // Process bond properties in parallel
                        BondProperty prop;
                        
                        // Use precomputed atom data with bounds checking
                        if (idx1 < numAtoms && idx2 < numAtoms) {
                            int an1 = atomNums[idx1];
                            int an2 = atomNums[idx2];
                            auto hyb1 = atomHybs[idx1];
                            auto hyb2 = atomHybs[idx2];
                            bool a1IsHetero = atomIsHetero[idx1];
                            bool a2IsHetero = atomIsHetero[idx2];
                            bool a1IsHeavy = atomIsHeavy[idx1];
                            bool a2IsHeavy = atomIsHeavy[idx2];
                            bool a1InRing = atomInRing[idx1];
                            bool a2InRing = atomInRing[idx2];
                            bool a1IsAromatic = atomIsAromatic[idx1];
                            bool a2IsAromatic = atomIsAromatic[idx2];
                            bool a1IsHalogen = atomIsHalogen[idx1];
                            bool a2IsHalogen = atomIsHalogen[idx2];
                            bool a1IsMetal = atomIsMetal[idx1];
                            bool a2IsMetal = atomIsMetal[idx2];

                            // Use precomputed bond data
                            bool isRotBond = bondIsRotatable[i];
                            bool bondInRing = ringInfo->numBondRings(i) > 0;

                            // Calculate bond-specific values
                            RDKit::Bond::BondType bondType = bond->getBondType();
                            double currentOrder = bond->getBondTypeAsDouble();
                            double currentPolarity = bondPolarityPrecomputed(an1, an2);
                            int currentAtomicNumDiff = std::abs(an1 - an2);
                            int roundedPolarity = static_cast<int>(currentPolarity * 10.0 + 0.5);

                            // --- Accumulate Counts ---
                            if (bondType == RDKit::Bond::SINGLE) {
                                prop.isSingle = true;
                                local_counters.singleBonds++;
                                bondPolarityCounts[roundedPolarity]++;
                            } else if (bondType == RDKit::Bond::DOUBLE) {
                                prop.isDouble = true;
                                local_counters.doubleBonds++;
                                bondPolarityCounts[roundedPolarity]++;
                            } else if (bondType == RDKit::Bond::TRIPLE) {
                                prop.isTriple = true;
                                local_counters.tripleBonds++;
                                bondPolarityCounts[roundedPolarity]++;
                            } else if (bondType == RDKit::Bond::AROMATIC) {
                                prop.isAromatic = true;
                                local_counters.aromaticBonds++;
                                bondPolarityCounts[roundedPolarity]++;
                            }

                            // Store bond order count
                            bondOrderCounts[currentOrder]++;

                            // Heteroatom Counts
                            prop.toHetero = a1IsHetero || a2IsHetero;
                            prop.betweenHetero = a1IsHetero && a2IsHetero;
                            if (prop.toHetero) local_counters.bondsToHetero++;
                            if (prop.betweenHetero) local_counters.bondsBetweenHetero++;

                            // C/H/Heavy Counts
                            prop.isCC = an1 == 6 && an2 == 6;
                            if (prop.isCC) local_counters.ccBonds++;
                            prop.toHydrogen = an1 == 1 || an2 == 1;
                            if (prop.toHydrogen) local_counters.bondsToHydrogen++;
                            prop.betweenHeavy = a1IsHeavy && a2IsHeavy;
                            if (prop.betweenHeavy) local_counters.heavyAtomBonds++;

                            // Ring Counts
                            prop.inRing = bondInRing;
                            if (prop.inRing) {
                                local_counters.ringBonds++;
                                prop.singleInRing = bondType == RDKit::Bond::SINGLE;
                                prop.doubleInRing = bondType == RDKit::Bond::DOUBLE;
                                prop.rotatableInRing = isRotBond;
                                if (prop.singleInRing) local_counters.singleRingBonds++;
                                if (prop.doubleInRing) local_counters.doubleRingBonds++;
                                if (prop.rotatableInRing) local_counters.rotatableRingBonds++;
                            }

                            // Aromatic Atom Counts
                            prop.withAromaticAtom = a1IsAromatic || a2IsAromatic;
                            prop.bothAromaticAtoms = a1IsAromatic && a2IsAromatic;
                            if (prop.withAromaticAtom) local_counters.bondsWithAromaticAtom++;
                            if (prop.bothAromaticAtoms) local_counters.bondsBothAromaticAtoms++;

                            // Hybridization Counts
                            prop.withSp3 = hyb1 == RDKit::Atom::SP3 || hyb2 == RDKit::Atom::SP3;
                            prop.bothSp3 = hyb1 == RDKit::Atom::SP3 && hyb2 == RDKit::Atom::SP3;
                            prop.withSp2 = hyb1 == RDKit::Atom::SP2 || hyb2 == RDKit::Atom::SP2;
                            prop.bothSp2 = hyb1 == RDKit::Atom::SP2 && hyb2 == RDKit::Atom::SP2;
                            prop.withSp = hyb1 == RDKit::Atom::SP || hyb2 == RDKit::Atom::SP;
                            prop.bothSp = hyb1 == RDKit::Atom::SP && hyb2 == RDKit::Atom::SP;
                            if (prop.withSp3) local_counters.bondsWithSp3++;
                            if (prop.bothSp3) local_counters.bondsBothSp3++;
                            if (prop.withSp2) local_counters.bondsWithSp2++;
                            if (prop.bothSp2) local_counters.bondsBothSp2++;
                            if (prop.withSp) local_counters.bondsWithSp++;
                            if (prop.bothSp) local_counters.bondsBothSp++;

                            // Halogen Counts
                            prop.withHalogen = a1IsHalogen || a2IsHalogen;
                            prop.bothHalogen = a1IsHalogen && a2IsHalogen;
                            if (prop.withHalogen) local_counters.bondsWithHalogen++;
                            if (prop.bothHalogen) local_counters.bondsBothHalogen++;

                            // Metal Counts
                            prop.withMetal = a1IsMetal || a2IsMetal;
                            prop.bothMetal = a1IsMetal && a2IsMetal;
                            if (prop.withMetal) local_counters.bondsWithMetal++;
                            if (prop.bothMetal) local_counters.bondsBothMetal++;

                            // Atom Ring Membership Counts
                            prop.withRingAtom = a1InRing || a2InRing;
                            prop.bothRingAtoms = a1InRing && a2InRing;
                            if (prop.withRingAtom) local_counters.bondsWithRingAtom++;
                            if (prop.bothRingAtoms) local_counters.bondsBothRingAtoms++;
                            prop.ringToNonRing = (a1InRing && !a2InRing) || (!a1InRing && a2InRing);
                            if (prop.ringToNonRing) local_counters.bondsRingToNonRing++;

                            // Store for stats
                            prop.order = currentOrder;
                            prop.polarity = currentPolarity;
                            prop.atomicNumDiff = currentAtomicNumDiff;
                            prop.roundedPolarity = roundedPolarity;

                            // Update Min/Max
                            maxBondPolarity = std::max(maxBondPolarity, currentPolarity);
                            minBondPolarity = std::min(minBondPolarity, currentPolarity);
                            maxBondOrder = std::max(maxBondOrder, currentOrder);
                            minBondOrder = std::min(minBondOrder, currentOrder);

                            // Thread-safe updates using the local counter
                            if (prop.isSingle) local_counters.singleBonds++;
                            if (prop.isDouble) local_counters.doubleBonds++;
                            if (prop.isTriple) local_counters.tripleBonds++;
                            if (prop.isAromatic) local_counters.aromaticBonds++;
                            if (prop.toHetero) local_counters.bondsToHetero++;
                            if (prop.betweenHetero) local_counters.bondsBetweenHetero++;
                            if (prop.isCC) local_counters.ccBonds++;
                            if (prop.toHydrogen) local_counters.bondsToHydrogen++;
                            if (prop.betweenHeavy) local_counters.heavyAtomBonds++;
                            if (prop.inRing) local_counters.ringBonds++;
                            if (prop.singleInRing) local_counters.singleRingBonds++;
                            if (prop.doubleInRing) local_counters.doubleRingBonds++;
                            if (prop.rotatableInRing) local_counters.rotatableRingBonds++;
                            if (prop.withAromaticAtom) local_counters.bondsWithAromaticAtom++;
                            if (prop.bothAromaticAtoms) local_counters.bondsBothAromaticAtoms++;
                            if (prop.withSp3) local_counters.bondsWithSp3++;
                            if (prop.bothSp3) local_counters.bondsBothSp3++;
                            if (prop.withSp2) local_counters.bondsWithSp2++;
                            if (prop.bothSp2) local_counters.bondsBothSp2++;
                            if (prop.withSp) local_counters.bondsWithSp++;
                            if (prop.bothSp) local_counters.bondsBothSp++;
                            if (prop.withHalogen) local_counters.bondsWithHalogen++;
                            if (prop.bothHalogen) local_counters.bondsBothHalogen++;
                            if (prop.withMetal) local_counters.bondsWithMetal++;
                            if (prop.bothMetal) local_counters.bondsBothMetal++;
                            if (prop.withRingAtom) local_counters.bondsWithRingAtom++;
                            if (prop.bothRingAtoms) local_counters.bondsBothRingAtoms++;
                            if (prop.ringToNonRing) local_counters.bondsRingToNonRing++;

                            // Store bond property
                            bondProps.push_back(prop);
                        }
                    }
                },
                tbb::static_partitioner()
            );
        });

        // Merge counters using parallel reduction
        BondCounters final_counters = tbb::parallel_reduce(
            tbb::blocked_range<TLS::iterator>(tls_counters.begin(), tls_counters.end()),
            BondCounters(),
            [](const tbb::blocked_range<TLS::iterator>& range, BondCounters init) {
                for (auto it = range.begin(); it != range.end(); ++it) {
                    init.merge(*it);
                }
                return init;
            },
            [](BondCounters a, BondCounters b) {
                a.merge(b);
                return a;
            }
        );

        // Update member variables
        singleBonds = final_counters.singleBonds;
        doubleBonds = final_counters.doubleBonds;
        tripleBonds = final_counters.tripleBonds;
        aromaticBonds = final_counters.aromaticBonds;
        rotatableBonds = final_counters.rotatableBonds;
        bondsToHetero = final_counters.bondsToHetero;
        bondsBetweenHetero = final_counters.bondsBetweenHetero;
        ccBonds = final_counters.ccBonds;
        bondsToHydrogen = final_counters.bondsToHydrogen;
        heavyAtomBonds = final_counters.heavyAtomBonds;
        ringBonds = final_counters.ringBonds;
        singleRingBonds = final_counters.singleRingBonds;
        doubleRingBonds = final_counters.doubleRingBonds;
        rotatableRingBonds = final_counters.rotatableRingBonds;
        bondsWithAromaticAtom = final_counters.bondsWithAromaticAtom;
        bondsBothAromaticAtoms = final_counters.bondsBothAromaticAtoms;
        bondsWithSp3 = final_counters.bondsWithSp3;
        bondsBothSp3 = final_counters.bondsBothSp3;
        bondsWithSp2 = final_counters.bondsWithSp2;
        bondsBothSp2 = final_counters.bondsBothSp2;
        bondsWithSp = final_counters.bondsWithSp;
        bondsBothSp = final_counters.bondsBothSp;
        bondsWithHalogen = final_counters.bondsWithHalogen;
        bondsBothHalogen = final_counters.bondsBothHalogen;
        bondsWithMetal = final_counters.bondsWithMetal;
        bondsBothMetal = final_counters.bondsBothMetal;
        bondsWithRingAtom = final_counters.bondsWithRingAtom;
        bondsBothRingAtoms = final_counters.bondsBothRingAtoms;
        bondsRingToNonRing = final_counters.bondsRingToNonRing;

        // Calculate statistics using cache-friendly access
        auto stats = tbb::parallel_reduce(
            tbb::blocked_range<size_t>(0, bondProps.size()),
            std::make_tuple(0.0, 0.0, 0.0),
            [&](const tbb::blocked_range<size_t>& range, auto init) -> std::tuple<double, double, double> {
                double sum_p = std::get<0>(init);
                double sum_o = std::get<1>(init);
                double sum_d = std::get<2>(init);
                
                // Process in cache-line sized chunks
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    const auto& prop = bondProps[i];  // Access directly without alignment
                    sum_p += prop.polarity;
                    sum_o += prop.order;
                    sum_d += static_cast<double>(prop.atomicNumDiff);
                }
                return std::make_tuple(sum_p, sum_o, sum_d);
            },
            [](const std::tuple<double, double, double>& a, const std::tuple<double, double, double>& b) -> std::tuple<double, double, double> {
                return std::make_tuple(
                    std::get<0>(a) + std::get<0>(b),
                    std::get<1>(a) + std::get<1>(b),
                    std::get<2>(a) + std::get<2>(b)
                );
            }
        );

        // Calculate final statistics
        auto [sum_polarity, sum_order, sum_atomic_diff] = stats;
        if (numBonds > 0) {  // Guard against division by zero
            meanBondPolarity = sum_polarity / numBonds;
            meanBondOrder = sum_order / numBonds;
            meanAtomicNumDiff = sum_atomic_diff / numBonds;

            // Calculate variances in a separate parallel reduction to avoid memory issues
            auto variances = tbb::parallel_reduce(
                tbb::blocked_range<size_t>(0, bondProps.size()),
                std::make_tuple(0.0, 0.0, 0.0),
                [&](const tbb::blocked_range<size_t>& range, auto init) -> std::tuple<double, double, double> {
                    double var_p = std::get<0>(init);
                    double var_o = std::get<1>(init);
                    double var_d = std::get<2>(init);
                    
                    for (size_t i = range.begin(); i < range.end(); ++i) {
                        const auto& prop = bondProps[i];
                        double dp = prop.polarity - meanBondPolarity;
                        double do_ = prop.order - meanBondOrder;
                        double dd = static_cast<double>(prop.atomicNumDiff) - meanAtomicNumDiff;
                        var_p += dp * dp;
                        var_o += do_ * do_;
                        var_d += dd * dd;
                    }
                    return std::make_tuple(var_p, var_o, var_d);
                },
                [](const std::tuple<double, double, double>& a, const std::tuple<double, double, double>& b) -> std::tuple<double, double, double> {
                    return std::make_tuple(
                        std::get<0>(a) + std::get<0>(b),
                        std::get<1>(a) + std::get<1>(b),
                        std::get<2>(a) + std::get<2>(b)
                    );
                }
            );

            auto [varPolarity, varOrder, varAtomicNumDiff] = variances;
            stddevBondPolarity = std::sqrt(varPolarity / numBonds);
            stddevBondOrder = std::sqrt(varOrder / numBonds);
            stddevAtomicNumDiff = std::sqrt(varAtomicNumDiff / numBonds);
        }

        // Most Common
        int maxOrderCount = 0;
        for (const auto& pair : bondOrderCounts) {
            if (pair.second > maxOrderCount) {
                maxOrderCount = pair.second;
                mostCommonBondOrder = pair.first;
            }
        }

        int maxPolarityCount = 0;
        int maxPolarityValue = 0;
        for (const auto& pair : bondPolarityCounts) {
            if (pair.second > maxPolarityCount) {
                maxPolarityCount = pair.second;
                maxPolarityValue = pair.first;
            }
        }
        mostCommonBondPolarity = static_cast<double>(maxPolarityValue) / 10.0; // Cast to double

        // Unique Counts
        uniqueBondOrders = bondOrderCounts.size();
        uniqueBondPolarities = bondPolarityCounts.size();
    }
};

static thread_local BondStatsCache tlsBondCache;

// Helper to ensure cache is populated
inline void ensureCachePopulated(Context& context) {
    const RDKit::ROMol* mol = context.getMolecule();
    if (tlsBondCache.populatedForMol != mol || !mol) { // Also check if mol is valid
        tlsBondCache.computeFromMol(mol);
    }
}

// --- Base Class ---
class BondStatsDescriptor : public Descriptor {
public:
    using Descriptor::Descriptor; // Inherit constructor
    std::string getCategory() const override { return "BondStats"; }
};


// --- Descriptor Implementations ---

// 1. Total bond count
DECLARE_DESCRIPTOR(BondCount, BondStatsDescriptor, "Total number of bonds")
DESCRIPTOR_DEPENDENCIES(BondCount) { return {}; }
DescriptorResult BondCountDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return static_cast<int>(tlsBondCache.numBonds);
}

// 2. Single bond count
DECLARE_DESCRIPTOR(SingleBondCount, BondStatsDescriptor, "Number of single bonds")
DESCRIPTOR_DEPENDENCIES(SingleBondCount) { return {}; }
DescriptorResult SingleBondCountDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return static_cast<int>(tlsBondCache.singleBonds);
}

// 3. Double bond count
DECLARE_DESCRIPTOR(DoubleBondCount, BondStatsDescriptor, "Number of double bonds")
DESCRIPTOR_DEPENDENCIES(DoubleBondCount) { return {}; }
DescriptorResult DoubleBondCountDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return static_cast<int>(tlsBondCache.doubleBonds);
}

// 4. Triple bond count
DECLARE_DESCRIPTOR(TripleBondCount, BondStatsDescriptor, "Number of triple bonds")
DESCRIPTOR_DEPENDENCIES(TripleBondCount) { return {}; }
DescriptorResult TripleBondCountDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return static_cast<int>(tlsBondCache.tripleBonds);
}

// 5. Aromatic bond count
DECLARE_DESCRIPTOR(AromaticBondCount, BondStatsDescriptor, "Number of aromatic bonds")
DESCRIPTOR_DEPENDENCIES(AromaticBondCount) { return {}; }
DescriptorResult AromaticBondCountDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return static_cast<int>(tlsBondCache.aromaticBonds);
}

// 6. Rotatable bond count
DECLARE_DESCRIPTOR(RotatableBondCount, BondStatsDescriptor, "Number of rotatable bonds")
DESCRIPTOR_DEPENDENCIES(RotatableBondCount) { return {}; }
DescriptorResult RotatableBondCountDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return static_cast<int>(tlsBondCache.rotatableBonds);
}

// 7. Ratio double/single bonds
DECLARE_DESCRIPTOR(DoubleToSingleBondRatio, BondStatsDescriptor, "Ratio of double to single bonds")
DESCRIPTOR_DEPENDENCIES(DoubleToSingleBondRatio) { return {}; }
DescriptorResult DoubleToSingleBondRatioDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.singleBonds > 0 ? static_cast<double>(tlsBondCache.doubleBonds) / tlsBondCache.singleBonds : 0.0;
}

// 8. Ratio triple/single bonds
DECLARE_DESCRIPTOR(TripleToSingleBondRatio, BondStatsDescriptor, "Ratio of triple to single bonds")
DESCRIPTOR_DEPENDENCIES(TripleToSingleBondRatio) { return {}; }
DescriptorResult TripleToSingleBondRatioDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.singleBonds > 0 ? static_cast<double>(tlsBondCache.tripleBonds) / tlsBondCache.singleBonds : 0.0;
}

// 9. Ratio aromatic/total bonds
DECLARE_DESCRIPTOR(AromaticBondFraction, BondStatsDescriptor, "Fraction of aromatic bonds")
DESCRIPTOR_DEPENDENCIES(AromaticBondFraction) { return {}; }
DescriptorResult AromaticBondFractionDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.aromaticBonds) / tlsBondCache.numBonds : 0.0;
}

// 10. Max bond polarity
DECLARE_DESCRIPTOR(MaxBondPolarity, BondStatsDescriptor, "Maximum bond polarity (Pauling EN diff)")
DESCRIPTOR_DEPENDENCIES(MaxBondPolarity) { return {}; }
DescriptorResult MaxBondPolarityDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.maxBondPolarity;
}

// 11. Min bond polarity
DECLARE_DESCRIPTOR(MinBondPolarity, BondStatsDescriptor, "Minimum bond polarity (Pauling EN diff)")
DESCRIPTOR_DEPENDENCIES(MinBondPolarity) { return {}; }
DescriptorResult MinBondPolarityDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.minBondPolarity == std::numeric_limits<double>::max() ? 0.0 : tlsBondCache.minBondPolarity; // Handle case of no bonds
}

// 12. Stddev bond polarity - REMOVED (Zero Variance)

// 13. Mean bond order - REMOVED (Zero Variance)

// 14. Max bond order
DECLARE_DESCRIPTOR(MaxBondOrder, BondStatsDescriptor, "Maximum bond order")
DESCRIPTOR_DEPENDENCIES(MaxBondOrder) { return {}; }
DescriptorResult MaxBondOrderDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.maxBondOrder;
}

// 15. Min bond order
DECLARE_DESCRIPTOR(MinBondOrder, BondStatsDescriptor, "Minimum bond order")
DESCRIPTOR_DEPENDENCIES(MinBondOrder) { return {}; }
DescriptorResult MinBondOrderDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
     return tlsBondCache.minBondOrder == std::numeric_limits<double>::max() ? 0.0 : tlsBondCache.minBondOrder; // Handle case of no bonds
}

// 16. Stddev bond order - REMOVED (Zero Variance)

// 17. Fraction of bonds to heteroatoms
DECLARE_DESCRIPTOR(FractionBondsToHeteroatoms, BondStatsDescriptor, "Fraction of bonds to heteroatoms")
DESCRIPTOR_DEPENDENCIES(FractionBondsToHeteroatoms) { return {}; }
DescriptorResult FractionBondsToHeteroatomsDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.bondsToHetero) / tlsBondCache.numBonds : 0.0;
}

// 18. Fraction of bonds between heteroatoms
DECLARE_DESCRIPTOR(FractionBondsBetweenHeteroatoms, BondStatsDescriptor, "Fraction of bonds between heteroatoms")
DESCRIPTOR_DEPENDENCIES(FractionBondsBetweenHeteroatoms) { return {}; }
DescriptorResult FractionBondsBetweenHeteroatomsDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.bondsBetweenHetero) / tlsBondCache.numBonds : 0.0;
}

// 19. Fraction of bonds between carbons
DECLARE_DESCRIPTOR(FractionCCBonds, BondStatsDescriptor, "Fraction of C-C bonds")
DESCRIPTOR_DEPENDENCIES(FractionCCBonds) { return {}; }
DescriptorResult FractionCCBondsDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.ccBonds) / tlsBondCache.numBonds : 0.0;
}

// 20. Fraction of bonds in rings
DECLARE_DESCRIPTOR(FractionRingBonds, BondStatsDescriptor, "Fraction of bonds in rings")
DESCRIPTOR_DEPENDENCIES(FractionRingBonds) { return {}; }
DescriptorResult FractionRingBondsDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.ringBonds) / tlsBondCache.numBonds : 0.0;
}

// 21. Fraction of single bonds in rings
DECLARE_DESCRIPTOR(FractionSingleRingBonds, BondStatsDescriptor, "Fraction of single bonds in rings")
DESCRIPTOR_DEPENDENCIES(FractionSingleRingBonds) { return {}; }
DescriptorResult FractionSingleRingBondsDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.singleRingBonds) / tlsBondCache.numBonds : 0.0;
}

// 22. Fraction of double bonds in rings
DECLARE_DESCRIPTOR(FractionDoubleRingBonds, BondStatsDescriptor, "Fraction of double bonds in rings")
DESCRIPTOR_DEPENDENCIES(FractionDoubleRingBonds) { return {}; }
DescriptorResult FractionDoubleRingBondsDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.doubleRingBonds) / tlsBondCache.numBonds : 0.0;
}

// 23. Fraction of rotatable bonds in rings - REMOVED (Zero Variance)

// 24. Fraction of rotatable bonds
DECLARE_DESCRIPTOR(FractionRotatableBonds, BondStatsDescriptor, "Fraction of rotatable bonds")
DESCRIPTOR_DEPENDENCIES(FractionRotatableBonds) { return {}; }
DescriptorResult FractionRotatableBondsDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.rotatableBonds) / tlsBondCache.numBonds : 0.0;
}

// 25. Fraction of bonds with at least one aromatic atom
DECLARE_DESCRIPTOR(FractionBondsWithAromaticAtom, BondStatsDescriptor, "Fraction of bonds with at least one aromatic atom")
DESCRIPTOR_DEPENDENCIES(FractionBondsWithAromaticAtom) { return {}; }
DescriptorResult FractionBondsWithAromaticAtomDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.bondsWithAromaticAtom) / tlsBondCache.numBonds : 0.0;
}

// 26. Fraction of bonds with both aromatic atoms
DECLARE_DESCRIPTOR(FractionBondsBothAromaticAtoms, BondStatsDescriptor, "Fraction of bonds with both atoms aromatic")
DESCRIPTOR_DEPENDENCIES(FractionBondsBothAromaticAtoms) { return {}; }
DescriptorResult FractionBondsBothAromaticAtomsDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.bondsBothAromaticAtoms) / tlsBondCache.numBonds : 0.0;
}

// 27. Fraction of bonds with at least one sp3 atom
DECLARE_DESCRIPTOR(FractionBondsWithSp3Atom, BondStatsDescriptor, "Fraction of bonds with at least one sp3 atom")
DESCRIPTOR_DEPENDENCIES(FractionBondsWithSp3Atom) { return {}; }
DescriptorResult FractionBondsWithSp3AtomDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.bondsWithSp3) / tlsBondCache.numBonds : 0.0;
}

// 28. Fraction of bonds with both sp3 atoms
DECLARE_DESCRIPTOR(FractionBondsBothSp3Atoms, BondStatsDescriptor, "Fraction of bonds with both atoms sp3")
DESCRIPTOR_DEPENDENCIES(FractionBondsBothSp3Atoms) { return {}; }
DescriptorResult FractionBondsBothSp3AtomsDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.bondsBothSp3) / tlsBondCache.numBonds : 0.0;
}

// 29. Fraction of bonds with at least one sp2 atom
DECLARE_DESCRIPTOR(FractionBondsWithSp2Atom, BondStatsDescriptor, "Fraction of bonds with at least one sp2 atom")
DESCRIPTOR_DEPENDENCIES(FractionBondsWithSp2Atom) { return {}; }
DescriptorResult FractionBondsWithSp2AtomDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.bondsWithSp2) / tlsBondCache.numBonds : 0.0;
}

// 30. Fraction of bonds with both sp2 atoms
DECLARE_DESCRIPTOR(FractionBondsBothSp2Atoms, BondStatsDescriptor, "Fraction of bonds with both atoms sp2")
DESCRIPTOR_DEPENDENCIES(FractionBondsBothSp2Atoms) { return {}; }
DescriptorResult FractionBondsBothSp2AtomsDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.bondsBothSp2) / tlsBondCache.numBonds : 0.0;
}

// 31. Fraction of bonds with at least one sp atom
DECLARE_DESCRIPTOR(FractionBondsWithSpAtom, BondStatsDescriptor, "Fraction of bonds with at least one sp atom")
DESCRIPTOR_DEPENDENCIES(FractionBondsWithSpAtom) { return {}; }
DescriptorResult FractionBondsWithSpAtomDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.bondsWithSp) / tlsBondCache.numBonds : 0.0;
}

// 32. Fraction of bonds with both sp atoms
DECLARE_DESCRIPTOR(FractionBondsBothSpAtoms, BondStatsDescriptor, "Fraction of bonds with both atoms sp")
DESCRIPTOR_DEPENDENCIES(FractionBondsBothSpAtoms) { return {}; }
DescriptorResult FractionBondsBothSpAtomsDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.bondsBothSp) / tlsBondCache.numBonds : 0.0;
}

// 33. Fraction of bonds with at least one halogen
DECLARE_DESCRIPTOR(FractionBondsWithHalogen, BondStatsDescriptor, "Fraction of bonds with at least one halogen atom")
DESCRIPTOR_DEPENDENCIES(FractionBondsWithHalogen) { return {}; }
DescriptorResult FractionBondsWithHalogenDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.bondsWithHalogen) / tlsBondCache.numBonds : 0.0;
}

// 34. Fraction of bonds with at least one atom in a ring
DECLARE_DESCRIPTOR(FractionBondsWithRingAtom, BondStatsDescriptor, "Fraction of bonds with at least one atom in a ring")
DESCRIPTOR_DEPENDENCIES(FractionBondsWithRingAtom) { return {}; }
DescriptorResult FractionBondsWithRingAtomDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.bondsWithRingAtom) / tlsBondCache.numBonds : 0.0;
}

// 35. Fraction of bonds with both atoms in a ring
DECLARE_DESCRIPTOR(FractionBondsBothRingAtoms, BondStatsDescriptor, "Fraction of bonds with both atoms in a ring")
DESCRIPTOR_DEPENDENCIES(FractionBondsBothRingAtoms) { return {}; }
DescriptorResult FractionBondsBothRingAtomsDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.bondsBothRingAtoms) / tlsBondCache.numBonds : 0.0;
}

// 36. Fraction of bonds with at least one atom in a ring and one not
DECLARE_DESCRIPTOR(FractionBondsRingToNonRing, BondStatsDescriptor, "Fraction of bonds with one atom in a ring and one not")
DESCRIPTOR_DEPENDENCIES(FractionBondsRingToNonRing) { return {}; }
DescriptorResult FractionBondsRingToNonRingDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.numBonds > 0 ? static_cast<double>(tlsBondCache.bondsRingToNonRing) / tlsBondCache.numBonds : 0.0;
}

// 37. Most common bond order
DECLARE_DESCRIPTOR(MostCommonBondOrder, BondStatsDescriptor, "Most common bond order")
DESCRIPTOR_DEPENDENCIES(MostCommonBondOrder) { return {}; }
DescriptorResult MostCommonBondOrderDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.mostCommonBondOrder;
}

// 38. Most common bond polarity
DECLARE_DESCRIPTOR(MostCommonBondPolarity, BondStatsDescriptor, "Most common bond polarity (rounded 0.1)")
DESCRIPTOR_DEPENDENCIES(MostCommonBondPolarity) { return {}; }
DescriptorResult MostCommonBondPolarityDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return tlsBondCache.mostCommonBondPolarity;
}

// 39. Number of unique bond orders
DECLARE_DESCRIPTOR(UniqueBondOrderCount, BondStatsDescriptor, "Number of unique bond orders")
DESCRIPTOR_DEPENDENCIES(UniqueBondOrderCount) { return {}; }
DescriptorResult UniqueBondOrderCountDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return static_cast<int>(tlsBondCache.uniqueBondOrders);
}

// 40. Number of unique bond polarities (rounded 0.1)
DECLARE_DESCRIPTOR(UniqueBondPolarityCount, BondStatsDescriptor, "Number of unique bond polarities (rounded 0.1)")
DESCRIPTOR_DEPENDENCIES(UniqueBondPolarityCount) { return {}; }
DescriptorResult UniqueBondPolarityCountDescriptor::calculate(Context& context) const {
    ensureCachePopulated(context);
    return static_cast<int>(tlsBondCache.uniqueBondPolarities);
}

// 41. Mean atomic number difference across bonds - REMOVED (Zero Variance)

// 42. Stddev atomic number difference across bonds - REMOVED (Zero Variance)

// Consolidated registration function
void register_BondStatsDescriptors() {
    auto& registry = DescriptorRegistry::getInstance();

    registry.registerDescriptor(std::make_shared<BondCountDescriptor>());
    registry.registerDescriptor(std::make_shared<SingleBondCountDescriptor>());
    registry.registerDescriptor(std::make_shared<DoubleBondCountDescriptor>());
    registry.registerDescriptor(std::make_shared<TripleBondCountDescriptor>());
    registry.registerDescriptor(std::make_shared<AromaticBondCountDescriptor>());
    registry.registerDescriptor(std::make_shared<RotatableBondCountDescriptor>());
    registry.registerDescriptor(std::make_shared<DoubleToSingleBondRatioDescriptor>());
    registry.registerDescriptor(std::make_shared<TripleToSingleBondRatioDescriptor>());
    registry.registerDescriptor(std::make_shared<AromaticBondFractionDescriptor>());
    registry.registerDescriptor(std::make_shared<MaxBondPolarityDescriptor>());
    registry.registerDescriptor(std::make_shared<MinBondPolarityDescriptor>());
    registry.registerDescriptor(std::make_shared<MaxBondOrderDescriptor>());
    registry.registerDescriptor(std::make_shared<MinBondOrderDescriptor>());
    registry.registerDescriptor(std::make_shared<FractionBondsToHeteroatomsDescriptor>());
    registry.registerDescriptor(std::make_shared<FractionBondsBetweenHeteroatomsDescriptor>());
    registry.registerDescriptor(std::make_shared<FractionCCBondsDescriptor>());
    registry.registerDescriptor(std::make_shared<FractionRingBondsDescriptor>());
    registry.registerDescriptor(std::make_shared<FractionSingleRingBondsDescriptor>());
    registry.registerDescriptor(std::make_shared<FractionDoubleRingBondsDescriptor>());
    registry.registerDescriptor(std::make_shared<FractionRotatableBondsDescriptor>());
    registry.registerDescriptor(std::make_shared<FractionBondsWithAromaticAtomDescriptor>());
    registry.registerDescriptor(std::make_shared<FractionBondsBothAromaticAtomsDescriptor>());
    registry.registerDescriptor(std::make_shared<FractionBondsWithSp3AtomDescriptor>());
    registry.registerDescriptor(std::make_shared<FractionBondsBothSp3AtomsDescriptor>());
    registry.registerDescriptor(std::make_shared<FractionBondsWithSp2AtomDescriptor>());
    registry.registerDescriptor(std::make_shared<FractionBondsBothSp2AtomsDescriptor>());
    registry.registerDescriptor(std::make_shared<FractionBondsWithSpAtomDescriptor>());
    registry.registerDescriptor(std::make_shared<FractionBondsBothSpAtomsDescriptor>());
    registry.registerDescriptor(std::make_shared<FractionBondsWithHalogenDescriptor>());
    registry.registerDescriptor(std::make_shared<FractionBondsWithRingAtomDescriptor>());
    registry.registerDescriptor(std::make_shared<FractionBondsBothRingAtomsDescriptor>());
    registry.registerDescriptor(std::make_shared<FractionBondsRingToNonRingDescriptor>());
    registry.registerDescriptor(std::make_shared<MostCommonBondOrderDescriptor>());
    registry.registerDescriptor(std::make_shared<MostCommonBondPolarityDescriptor>());
    registry.registerDescriptor(std::make_shared<UniqueBondOrderCountDescriptor>());
    registry.registerDescriptor(std::make_shared<UniqueBondPolarityCountDescriptor>());
}


void register_BondCountDescriptor() {
    auto descriptor = std::make_shared<BondCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SingleBondCountDescriptor() {
    auto descriptor = std::make_shared<SingleBondCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_DoubleBondCountDescriptor() {
    auto descriptor = std::make_shared<DoubleBondCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_TripleBondCountDescriptor() {
    auto descriptor = std::make_shared<TripleBondCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_AromaticBondCountDescriptor() {
    auto descriptor = std::make_shared<AromaticBondCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_RotatableBondCountDescriptor() {
    auto descriptor = std::make_shared<RotatableBondCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_DoubleToSingleBondRatioDescriptor() {
    auto descriptor = std::make_shared<DoubleToSingleBondRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_TripleToSingleBondRatioDescriptor() {
    auto descriptor = std::make_shared<TripleToSingleBondRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_AromaticBondFractionDescriptor() {
    auto descriptor = std::make_shared<AromaticBondFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MaxBondPolarityDescriptor() {
    auto descriptor = std::make_shared<MaxBondPolarityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MinBondPolarityDescriptor() {
    auto descriptor = std::make_shared<MinBondPolarityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MaxBondOrderDescriptor() {
    auto descriptor = std::make_shared<MaxBondOrderDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MinBondOrderDescriptor() {
    auto descriptor = std::make_shared<MinBondOrderDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionBondsToHeteroatomsDescriptor() {
    auto descriptor = std::make_shared<FractionBondsToHeteroatomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionBondsBetweenHeteroatomsDescriptor() {
    auto descriptor = std::make_shared<FractionBondsBetweenHeteroatomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionCCBondsDescriptor() {
    auto descriptor = std::make_shared<FractionCCBondsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionRingBondsDescriptor() {
    auto descriptor = std::make_shared<FractionRingBondsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionSingleRingBondsDescriptor() {
    auto descriptor = std::make_shared<FractionSingleRingBondsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionDoubleRingBondsDescriptor() {
    auto descriptor = std::make_shared<FractionDoubleRingBondsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionRotatableBondsDescriptor() {
    auto descriptor = std::make_shared<FractionRotatableBondsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionBondsWithAromaticAtomDescriptor() {
    auto descriptor = std::make_shared<FractionBondsWithAromaticAtomDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionBondsBothAromaticAtomsDescriptor() {
    auto descriptor = std::make_shared<FractionBondsBothAromaticAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionBondsWithSp3AtomDescriptor() {
    auto descriptor = std::make_shared<FractionBondsWithSp3AtomDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionBondsBothSp3AtomsDescriptor() {
    auto descriptor = std::make_shared<FractionBondsBothSp3AtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionBondsWithSp2AtomDescriptor() {
    auto descriptor = std::make_shared<FractionBondsWithSp2AtomDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionBondsBothSp2AtomsDescriptor() {
    auto descriptor = std::make_shared<FractionBondsBothSp2AtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionBondsWithSpAtomDescriptor() {
    auto descriptor = std::make_shared<FractionBondsWithSpAtomDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionBondsBothSpAtomsDescriptor() {
    auto descriptor = std::make_shared<FractionBondsBothSpAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionBondsWithHalogenDescriptor() {
    auto descriptor = std::make_shared<FractionBondsWithHalogenDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionBondsWithRingAtomDescriptor() {
    auto descriptor = std::make_shared<FractionBondsWithRingAtomDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionBondsBothRingAtomsDescriptor() {
    auto descriptor = std::make_shared<FractionBondsBothRingAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FractionBondsRingToNonRingDescriptor() {
    auto descriptor = std::make_shared<FractionBondsRingToNonRingDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MostCommonBondOrderDescriptor() {
    auto descriptor = std::make_shared<MostCommonBondOrderDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MostCommonBondPolarityDescriptor() {
    auto descriptor = std::make_shared<MostCommonBondPolarityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_UniqueBondOrderCountDescriptor() {
    auto descriptor = std::make_shared<UniqueBondOrderCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_UniqueBondPolarityCountDescriptor() {
    auto descriptor = std::make_shared<UniqueBondPolarityCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

} // namespace desfact

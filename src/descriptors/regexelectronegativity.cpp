#include "../common.hpp"
#include <re2/re2.h>
#include <string>
#include <algorithm>
#include <memory>
#include <vector>
#include <numeric>
#include <map>
#include <cmath>
#include <set>

namespace desfact {

// --- Electronegativity Data Maps ---
// Using const maps for thread-safety
static const std::map<std::string, double> paulingElectronegativityMap = {
    {"H", 2.20}, {"Li", 0.98}, {"Be", 1.57}, {"B", 2.04}, {"C", 2.55}, {"N", 3.04}, {"O", 3.44}, {"F", 3.98},
    {"Na", 0.93}, {"Mg", 1.31}, {"Al", 1.61}, {"Si", 1.90}, {"P", 2.19}, {"S", 2.58}, {"Cl", 3.16},
    {"K", 0.82}, {"Ca", 1.00}, {"Sc", 1.36}, {"Ti", 1.54}, {"V", 1.63}, {"Cr", 1.66}, {"Mn", 1.55}, 
    {"Fe", 1.83}, {"Co", 1.88}, {"Ni", 1.91}, {"Cu", 1.90}, {"Zn", 1.65}, {"Ga", 1.81}, {"Ge", 2.01}, 
    {"As", 2.18}, {"Se", 2.55}, {"Br", 2.96}, {"Rb", 0.82}, {"Sr", 0.95}, {"Y", 1.22}, {"Zr", 1.33}, 
    {"Nb", 1.6}, {"Mo", 2.16}, {"Tc", 1.9}, {"Ru", 2.2}, {"Rh", 2.28}, {"Pd", 2.20}, {"Ag", 1.93}, 
    {"Cd", 1.69}, {"In", 1.78}, {"Sn", 1.96}, {"Sb", 2.05}, {"Te", 2.1}, {"I", 2.66},
    // Lowercase for aromatic atoms (RDKit convention)
    {"c", 2.55}, {"n", 3.04}, {"o", 3.44}, {"p", 2.19}, {"s", 2.58}, {"b", 2.04}, {"i", 2.66}
};

static const std::map<std::string, double> allredRochowENMap = {
    {"H", 2.20}, {"C", 2.50}, {"N", 3.07}, {"O", 3.50}, {"F", 4.10}, {"P", 2.06}, {"S", 2.44}, 
    {"Cl", 2.83}, {"Br", 2.74}, {"I", 2.21}, {"B", 2.01}, {"Si", 1.74}, {"As", 2.20}, {"Se", 2.48},
    {"c", 2.50}, {"n", 3.07}, {"o", 3.50}, {"p", 2.06}, {"s", 2.44}, {"b", 2.01}, {"i", 2.21}
};

static const std::map<std::string, double> mullikenENMap = {
    {"H", 2.25}, {"C", 2.63}, {"N", 3.06}, {"O", 3.61}, {"F", 4.43}, {"P", 2.32}, {"S", 2.69}, 
    {"Cl", 3.54}, {"Br", 3.24}, {"I", 2.37}, {"B", 2.05}, {"Si", 2.03}, {"As", 2.33}, {"Se", 2.75},
    {"c", 2.63}, {"n", 3.06}, {"o", 3.61}, {"p", 2.32}, {"s", 2.69}, {"b", 2.05}, {"i", 2.37}
};

static const std::map<std::string, double> sandersonENMap = {
    {"H", 2.59}, {"C", 2.75}, {"N", 3.19}, {"O", 3.65}, {"F", 4.00}, {"P", 2.55}, {"S", 2.96}, 
    {"Cl", 3.48}, {"Br", 3.22}, {"I", 2.14}, {"B", 2.28}, {"Si", 2.14}, {"As", 2.55}, {"Se", 2.84},
    {"c", 2.75}, {"n", 3.19}, {"o", 3.65}, {"p", 2.55}, {"s", 2.96}, {"b", 2.28}, {"i", 2.14}
};

static const std::map<std::string, double> atomicWeights = {
    {"H", 1.008}, {"Li", 6.941}, {"Be", 9.012}, {"B", 10.811}, {"C", 12.011}, {"N", 14.007}, 
    {"O", 15.999}, {"F", 18.998}, {"Na", 22.99}, {"Mg", 24.305}, {"Al", 26.982}, {"Si", 28.085}, 
    {"P", 30.974}, {"S", 32.065}, {"Cl", 35.453}, {"K", 39.098}, {"Ca", 40.078}, {"Sc", 44.956}, 
    {"Ti", 47.867}, {"V", 50.942}, {"Cr", 51.996}, {"Mn", 54.938}, {"Fe", 55.845}, {"Co", 58.933}, 
    {"Ni", 58.693}, {"Cu", 63.546}, {"Zn", 65.38}, {"Ga", 69.723}, {"Ge", 72.64}, {"As", 74.922}, 
    {"Se", 78.96}, {"Br", 79.904}, {"Rb", 85.468}, {"Sr", 87.62}, {"Y", 88.906}, {"Zr", 91.224}, 
    {"Nb", 92.906}, {"Mo", 95.94}, {"Tc", 98}, {"Ru", 101.07}, {"Rh", 102.906}, {"Pd", 106.42}, 
    {"Ag", 107.868}, {"Cd", 112.41}, {"In", 114.82}, {"Sn", 118.71}, {"Sb", 121.76}, {"Te", 127.6}, 
    {"I", 126.904},
    // Lowercase for aromatic atoms
    {"c", 12.011}, {"n", 14.007}, {"o", 15.999}, {"p", 30.974}, {"s", 32.065}, {"b", 10.811},
    {"i", 126.904}
};

// --- Optimized SMILES Element Extraction ---

// Thread-local cache: only create a new regex pattern once per thread
RE2& getElementRegex() {
    // Static per-thread instance, initialized once per thread 
    static thread_local RE2 pattern(
        "\\[(?:[0-9]+)?([A-Z][a-z]?|[cnospbih])(?:H[0-9]?)?(?:[+-][0-9]*)?(?::[0-9]+)?(?:@[@THALBOZRDGSP0-9]*)?\\]"
        "|(Cl|Br|Se|Si|As)"
        "|([A-Za-z])"
    );
    return pattern;
}

// Thread-local storage to avoid repeated allocations
struct ElementCache {
    std::vector<std::string> elements;
    std::vector<double> paulingValues;
    std::vector<double> arValues;
    std::vector<double> mullikenValues;
    std::vector<double> sandersonValues;
    
    // Reuse vectors to avoid reallocations
    void reset() {
        elements.clear();
        paulingValues.clear();
        arValues.clear();
        mullikenValues.clear();
        sandersonValues.clear();
    }
};

static thread_local ElementCache tlsElementCache;

// Single-pass extraction of all necessary data from SMILES
void extractElementData(const std::string& smiles, ElementCache& cache) {
    cache.reset();
    cache.elements.reserve(smiles.length() / 2);  // Reasonable estimate
    
    RE2& pattern = getElementRegex();
    re2::StringPiece input(smiles);
    std::string bracketed_elem, two_letter_elem, single_letter_elem;

    while (RE2::Consume(&input, pattern, &bracketed_elem, &two_letter_elem, &single_letter_elem)) {
        std::string element;
        
        if (!bracketed_elem.empty()) {
            element = bracketed_elem;
        } else if (!two_letter_elem.empty()) {
            element = two_letter_elem;
        } else if (!single_letter_elem.empty()) {
            // Only include known elements
            if (paulingElectronegativityMap.count(single_letter_elem) > 0) {
                element = single_letter_elem;
            }
        }
        
        if (!element.empty()) {
            cache.elements.push_back(element);
        }
        
        // Reset capture groups for next match
        bracketed_elem.clear();
        two_letter_elem.clear();
        single_letter_elem.clear();
    }
    
    // Pre-compute all values in one pass (avoids multiple loops later)
    cache.paulingValues.reserve(cache.elements.size());
    cache.arValues.reserve(cache.elements.size());
    cache.mullikenValues.reserve(cache.elements.size());
    cache.sandersonValues.reserve(cache.elements.size());
    
    for (const auto& elem : cache.elements) {
        auto it = paulingElectronegativityMap.find(elem);
        if (it != paulingElectronegativityMap.end()) {
            cache.paulingValues.push_back(it->second);
        }
        
        it = allredRochowENMap.find(elem);
        if (it != allredRochowENMap.end()) {
            cache.arValues.push_back(it->second);
        }
        
        it = mullikenENMap.find(elem);
        if (it != mullikenENMap.end()) {
            cache.mullikenValues.push_back(it->second);
        }
        
        it = sandersonENMap.find(elem);
        if (it != sandersonENMap.end()) {
            cache.sandersonValues.push_back(it->second);
        }
    }
}

// --- Calculation Helper Functions ---
// All functions below are simple and vectorizable without branches

inline double sumElectronegativity(const std::vector<double>& values) {
    return std::accumulate(values.begin(), values.end(), 0.0);
}

inline double meanElectronegativity(const std::vector<double>& values) {
    if (values.empty()) return 0.0;
    return sumElectronegativity(values) / static_cast<double>(values.size());
}

inline double maxElectronegativity(const std::vector<double>& values) {
    if (values.empty()) return 0.0;
    return *std::max_element(values.begin(), values.end());
}

inline double minElectronegativity(const std::vector<double>& values) {
    if (values.empty()) return 0.0;
    return *std::min_element(values.begin(), values.end());
}

inline double electronegativityRange(const std::vector<double>& values) {
    if (values.size() < 2) return 0.0;
    auto result = std::minmax_element(values.begin(), values.end());
    return *result.second - *result.first;
}

inline double electronegativityStdDev(const std::vector<double>& values) {
    size_t count = values.size();
    if (count < 2) return 0.0;

    double mean = meanElectronegativity(values);
    double variance = std::accumulate(values.begin(), values.end(), 0.0,
        [mean](double acc, double val) {
            double diff = val - mean;
            return acc + diff * diff;
        });

    return std::sqrt(variance / count);
}

inline double countElementsWithENAbove(const std::vector<double>& values, double threshold) {
    return static_cast<double>(std::count_if(values.begin(), values.end(),
        [threshold](double val) { return val > threshold; }
    ));
}

inline double ratioOfHighToLowEN(const std::vector<double>& values, double threshold) {
    if (values.empty()) return 0.0;
    double highCount = countElementsWithENAbove(values, threshold);
    double lowCount = static_cast<double>(values.size()) - highCount;
    return (lowCount > 0.0) ? highCount / lowCount : highCount;
}

// --- Element Filter Functions ---
// Avoid string comparisons by using type maps

// Fast element type check - bitsets would be faster but this works fine for small sets
struct ElementType {
    static bool isCarbon(const std::string& elem) { return elem == "C" || elem == "c"; }
    static bool isNitrogen(const std::string& elem) { return elem == "N" || elem == "n"; }
    static bool isOxygen(const std::string& elem) { return elem == "O" || elem == "o"; }
    static bool isSulfur(const std::string& elem) { return elem == "S" || elem == "s"; }
    static bool isAromatic(const std::string& elem) { return std::islower(elem[0]); }
    static bool isHalogen(const std::string& elem) { 
        return elem == "F" || elem == "Cl" || elem == "Br" || elem == "I"; 
    }
};

// Cache calculated results for element groups
struct ElementGroups {
    std::vector<size_t> carbon, nitrogen, oxygen, sulfur, aromatic, aliphatic, halogen;
    
    void analyze(const std::vector<std::string>& elements) {
        carbon.clear(); nitrogen.clear(); oxygen.clear(); sulfur.clear();
        aromatic.clear(); aliphatic.clear(); halogen.clear();
        
        for (size_t i = 0; i < elements.size(); ++i) {
            const auto& elem = elements[i];
            if (ElementType::isCarbon(elem)) carbon.push_back(i);
            if (ElementType::isNitrogen(elem)) nitrogen.push_back(i);
            if (ElementType::isOxygen(elem)) oxygen.push_back(i);
            if (ElementType::isSulfur(elem)) sulfur.push_back(i);
            if (ElementType::isHalogen(elem)) halogen.push_back(i);
            
            if (ElementType::isAromatic(elem)) {
                aromatic.push_back(i);
            } else {
                aliphatic.push_back(i);
            }
        }
    }
};

static thread_local ElementGroups tlsElementGroups;

// Calculate sum of EN values for a subset of elements
double sumENForElementIndices(const std::vector<double>& allValues, 
                              const std::vector<size_t>& indices) {
    double sum = 0.0;
    for (size_t idx : indices) {
        if (idx < allValues.size()) {
            sum += allValues[idx];
        }
    }
    return sum;
}

// --- Base Class ---
class ElectronegativityDescriptor : public Descriptor {
public:
    using Descriptor::Descriptor; // Inherit constructor
    std::string getCategory() const override { return "Electronegativity"; }
};

// --- Descriptor Implementations ---

// 1. Total Pauling Electronegativity Sum
DECLARE_DESCRIPTOR(SmilesPaulingENSum, ElectronegativityDescriptor, "Sum of Pauling electronegativity values of all atoms")
DESCRIPTOR_DEPENDENCIES(SmilesPaulingENSum) { return {}; }
DescriptorResult SmilesPaulingENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return sumElectronegativity(tlsElementCache.paulingValues);
}

// 2. Mean Pauling Electronegativity
DECLARE_DESCRIPTOR(SmilesPaulingENMean, ElectronegativityDescriptor, "Mean Pauling electronegativity of all atoms")
DESCRIPTOR_DEPENDENCIES(SmilesPaulingENMean) { return {}; }
DescriptorResult SmilesPaulingENMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return meanElectronegativity(tlsElementCache.paulingValues);
}

// 3. Maximum Pauling Electronegativity
DECLARE_DESCRIPTOR(SmilesPaulingENMax, ElectronegativityDescriptor, "Maximum Pauling electronegativity among all atoms")
DESCRIPTOR_DEPENDENCIES(SmilesPaulingENMax) { return {}; }
DescriptorResult SmilesPaulingENMaxDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return maxElectronegativity(tlsElementCache.paulingValues);
}

// 4. Minimum Pauling Electronegativity
DECLARE_DESCRIPTOR(SmilesPaulingENMin, ElectronegativityDescriptor, "Minimum Pauling electronegativity among all atoms")
DESCRIPTOR_DEPENDENCIES(SmilesPaulingENMin) { return {}; }
DescriptorResult SmilesPaulingENMinDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return minElectronegativity(tlsElementCache.paulingValues);
}

// 5. Pauling Electronegativity Range
DECLARE_DESCRIPTOR(SmilesPaulingENRange, ElectronegativityDescriptor, "Range of Pauling electronegativity values")
DESCRIPTOR_DEPENDENCIES(SmilesPaulingENRange) { return {}; }
DescriptorResult SmilesPaulingENRangeDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return electronegativityRange(tlsElementCache.paulingValues);
}

// 6. Pauling Electronegativity Standard Deviation
DECLARE_DESCRIPTOR(SmilesPaulingENStdDev, ElectronegativityDescriptor, "Standard deviation of Pauling electronegativity values")
DESCRIPTOR_DEPENDENCIES(SmilesPaulingENStdDev) { return {}; }
DescriptorResult SmilesPaulingENStdDevDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return electronegativityStdDev(tlsElementCache.paulingValues);
}

// 7. Count of High Electronegativity Atoms (Pauling > 3.0)
DECLARE_DESCRIPTOR(SmilesHighENCount, ElectronegativityDescriptor, "Count of atoms with Pauling EN > 3.0")
DESCRIPTOR_DEPENDENCIES(SmilesHighENCount) { return {}; }
DescriptorResult SmilesHighENCountDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return countElementsWithENAbove(tlsElementCache.paulingValues, 3.0);
}

// 8. Ratio of High to Low Electronegativity Atoms (threshold 3.0)
DECLARE_DESCRIPTOR(SmilesHighToLowENRatio, ElectronegativityDescriptor, "Ratio of atoms with EN > 3.0 to atoms with EN <= 3.0")
DESCRIPTOR_DEPENDENCIES(SmilesHighToLowENRatio) { return {}; }
DescriptorResult SmilesHighToLowENRatioDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return ratioOfHighToLowEN(tlsElementCache.paulingValues, 3.0);
}

// 9. Count of Very High Electronegativity Atoms (Pauling > 3.5)
DECLARE_DESCRIPTOR(SmilesVeryHighENCount, ElectronegativityDescriptor, "Count of atoms with Pauling EN > 3.5")
DESCRIPTOR_DEPENDENCIES(SmilesVeryHighENCount) { return {}; }
DescriptorResult SmilesVeryHighENCountDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return countElementsWithENAbove(tlsElementCache.paulingValues, 3.5);
}

// 10. Fraction of High Electronegativity Atoms
DECLARE_DESCRIPTOR(SmilesHighENFraction, ElectronegativityDescriptor, "Fraction of atoms with Pauling EN > 3.0")
DESCRIPTOR_DEPENDENCIES(SmilesHighENFraction) { return {}; }
DescriptorResult SmilesHighENFractionDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    if (tlsElementCache.elements.empty()) return 0.0;
    return countElementsWithENAbove(tlsElementCache.paulingValues, 3.0) / 
           static_cast<double>(tlsElementCache.elements.size());
}

// 11. Total Allred-Rochow Electronegativity Sum
DECLARE_DESCRIPTOR(SmilesAllredRochowENSum, ElectronegativityDescriptor, "Sum of Allred-Rochow electronegativity values")
DESCRIPTOR_DEPENDENCIES(SmilesAllredRochowENSum) { return {}; }
DescriptorResult SmilesAllredRochowENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return sumElectronegativity(tlsElementCache.arValues);
}

// 12. Mean Allred-Rochow Electronegativity
DECLARE_DESCRIPTOR(SmilesAllredRochowENMean, ElectronegativityDescriptor, "Mean Allred-Rochow electronegativity")
DESCRIPTOR_DEPENDENCIES(SmilesAllredRochowENMean) { return {}; }
DescriptorResult SmilesAllredRochowENMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return meanElectronegativity(tlsElementCache.arValues);
}

// 13. Allred-Rochow to Pauling Ratio
DECLARE_DESCRIPTOR(SmilesAllredToPaulingRatio, ElectronegativityDescriptor, "Ratio of Allred-Rochow to Pauling EN means")
DESCRIPTOR_DEPENDENCIES(SmilesAllredToPaulingRatio) { return {}; }
DescriptorResult SmilesAllredToPaulingRatioDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double arMean = meanElectronegativity(tlsElementCache.arValues);
    double paulingMean = meanElectronegativity(tlsElementCache.paulingValues);
    return (paulingMean != 0.0) ? arMean / paulingMean : 0.0;
}

// 14. Total Mulliken Electronegativity Sum
DECLARE_DESCRIPTOR(SmilesMullikenENSum, ElectronegativityDescriptor, "Sum of Mulliken electronegativity values")
DESCRIPTOR_DEPENDENCIES(SmilesMullikenENSum) { return {}; }
DescriptorResult SmilesMullikenENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return sumElectronegativity(tlsElementCache.mullikenValues);
}

// 15. Mean Mulliken Electronegativity
DECLARE_DESCRIPTOR(SmilesMullikenENMean, ElectronegativityDescriptor, "Mean Mulliken electronegativity")
DESCRIPTOR_DEPENDENCIES(SmilesMullikenENMean) { return {}; }
DescriptorResult SmilesMullikenENMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return meanElectronegativity(tlsElementCache.mullikenValues);
}

// 16. Total Sanderson Electronegativity Sum
DECLARE_DESCRIPTOR(SmilesSandersonENSum, ElectronegativityDescriptor, "Sum of Sanderson electronegativity values")
DESCRIPTOR_DEPENDENCIES(SmilesSandersonENSum) { return {}; }
DescriptorResult SmilesSandersonENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return sumElectronegativity(tlsElementCache.sandersonValues);
}

// 17. Mean Sanderson Electronegativity
DECLARE_DESCRIPTOR(SmilesSandersonENMean, ElectronegativityDescriptor, "Mean Sanderson electronegativity")
DESCRIPTOR_DEPENDENCIES(SmilesSandersonENMean) { return {}; }
DescriptorResult SmilesSandersonENMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return meanElectronegativity(tlsElementCache.sandersonValues);
}

// Helper for element fraction * mean EN calculations
double calculateElementFractionENProduct(const std::string& smiles, 
                                        std::function<bool(const std::string&)> elementFilter) {
    auto& cache = tlsElementCache;
    extractElementData(smiles, cache);
    if (cache.elements.empty()) return 0.0;

    // Count matching elements
    double matchCount = 0.0;
    for (const auto& elem : cache.elements) {
        if (elementFilter(elem)) matchCount += 1.0;
    }
    
    double fraction = matchCount / static_cast<double>(cache.elements.size());
    return fraction * meanElectronegativity(cache.paulingValues);
}

// 18. Carbon Fraction × Pauling EN
DECLARE_DESCRIPTOR(SmilesCarbonENProduct, ElectronegativityDescriptor, "Product of C fraction and mean Pauling EN")
DESCRIPTOR_DEPENDENCIES(SmilesCarbonENProduct) { return {}; }
DescriptorResult SmilesCarbonENProductDescriptor::calculate(Context& context) const {
    return calculateElementFractionENProduct(context.getSmiles(), ElementType::isCarbon);
}

// 19. Nitrogen Fraction × Pauling EN
DECLARE_DESCRIPTOR(SmilesNitrogenENProduct, ElectronegativityDescriptor, "Product of N fraction and mean Pauling EN")
DESCRIPTOR_DEPENDENCIES(SmilesNitrogenENProduct) { return {}; }
DescriptorResult SmilesNitrogenENProductDescriptor::calculate(Context& context) const {
    return calculateElementFractionENProduct(context.getSmiles(), ElementType::isNitrogen);
}

// 20. Oxygen Fraction × Pauling EN
DECLARE_DESCRIPTOR(SmilesOxygenENProduct, ElectronegativityDescriptor, "Product of O fraction and mean Pauling EN")
DESCRIPTOR_DEPENDENCIES(SmilesOxygenENProduct) { return {}; }
DescriptorResult SmilesOxygenENProductDescriptor::calculate(Context& context) const {
    return calculateElementFractionENProduct(context.getSmiles(), ElementType::isOxygen);
}

// 21. Sulfur Fraction × Pauling EN
DECLARE_DESCRIPTOR(SmilesSulfurENProduct, ElectronegativityDescriptor, "Product of S fraction and mean Pauling EN")
DESCRIPTOR_DEPENDENCIES(SmilesSulfurENProduct) { return {}; }
DescriptorResult SmilesSulfurENProductDescriptor::calculate(Context& context) const {
    return calculateElementFractionENProduct(context.getSmiles(), ElementType::isSulfur);
}

// 22. Halogen Fraction × Pauling EN
DECLARE_DESCRIPTOR(SmilesHalogenENProduct, ElectronegativityDescriptor, "Product of halogen fraction and mean Pauling EN")
DESCRIPTOR_DEPENDENCIES(SmilesHalogenENProduct) { return {}; }
DescriptorResult SmilesHalogenENProductDescriptor::calculate(Context& context) const {
    return calculateElementFractionENProduct(context.getSmiles(), ElementType::isHalogen);
}

// 23. EN Variance Across Scales
DECLARE_DESCRIPTOR(SmilesENScaleVariance, ElectronegativityDescriptor, "Variance between different EN scale means")
DESCRIPTOR_DEPENDENCIES(SmilesENScaleVariance) { return {}; }
DescriptorResult SmilesENScaleVarianceDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    
    // Store means of each scale (only include scales with values)
    std::vector<double> means;
    means.reserve(4);  // Pauling, Allred-Rochow, Mulliken, Sanderson
    
    if (!tlsElementCache.paulingValues.empty()) 
        means.push_back(meanElectronegativity(tlsElementCache.paulingValues));
    
    if (!tlsElementCache.arValues.empty()) 
        means.push_back(meanElectronegativity(tlsElementCache.arValues));
    
    if (!tlsElementCache.mullikenValues.empty()) 
        means.push_back(meanElectronegativity(tlsElementCache.mullikenValues));
    
    if (!tlsElementCache.sandersonValues.empty()) 
        means.push_back(meanElectronegativity(tlsElementCache.sandersonValues));
    
    if (means.size() < 2) return 0.0;  // Need at least 2 scales
    
    double overallMean = meanElectronegativity(means);
    
    // Calculate variance
    double sumSquares = std::accumulate(means.begin(), means.end(), 0.0,
        [overallMean](double acc, double mean) {
            double diff = mean - overallMean;
            return acc + diff * diff;
        });
    
    return sumSquares / static_cast<double>(means.size());
}

// 24. EN Consistency Index (max/min across scales)
DECLARE_DESCRIPTOR(SmilesENConsistencyIndex, ElectronegativityDescriptor, "Ratio of max to min mean EN across scales")
DESCRIPTOR_DEPENDENCIES(SmilesENConsistencyIndex) { return {}; }
DescriptorResult SmilesENConsistencyIndexDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    
    // Store means of each scale (only include scales with values)
    std::vector<double> means;
    means.reserve(4);  // Pauling, Allred-Rochow, Mulliken, Sanderson
    
    if (!tlsElementCache.paulingValues.empty()) 
        means.push_back(meanElectronegativity(tlsElementCache.paulingValues));
    
    if (!tlsElementCache.arValues.empty()) 
        means.push_back(meanElectronegativity(tlsElementCache.arValues));
    
    if (!tlsElementCache.mullikenValues.empty()) 
        means.push_back(meanElectronegativity(tlsElementCache.mullikenValues));
    
    if (!tlsElementCache.sandersonValues.empty()) 
        means.push_back(meanElectronegativity(tlsElementCache.sandersonValues));
    
    if (means.size() < 2) return 0.0;  // Need at least 2 scales
    
    auto result = std::minmax_element(means.begin(), means.end());
    double minMean = *result.first;
    double maxMean = *result.second;
    
    return (minMean != 0.0) ? maxMean / minMean : 0.0;
}

// Specific element type sums for electronegativity

// 25. Aromatic Carbon EN Sum
DECLARE_DESCRIPTOR(SmilesAromaticCarbonENSum, ElectronegativityDescriptor, "Sum of EN values for aromatic carbons")
DESCRIPTOR_DEPENDENCIES(SmilesAromaticCarbonENSum) { return {}; }
DescriptorResult SmilesAromaticCarbonENSumDescriptor::calculate(Context& context) const {
    auto& cache = tlsElementCache;
    extractElementData(context.getSmiles(), cache);
    
    double sum = 0.0;
    for (size_t i = 0; i < cache.elements.size(); ++i) {
        if (cache.elements[i] == "c" && i < cache.paulingValues.size()) {
            sum += cache.paulingValues[i];
        }
    }
    return sum;
}

// 26. Aliphatic Carbon EN Sum
DECLARE_DESCRIPTOR(SmilesAliphaticCarbonENSum, ElectronegativityDescriptor, "Sum of EN values for aliphatic carbons")
DESCRIPTOR_DEPENDENCIES(SmilesAliphaticCarbonENSum) { return {}; }
DescriptorResult SmilesAliphaticCarbonENSumDescriptor::calculate(Context& context) const {
    auto& cache = tlsElementCache;
    extractElementData(context.getSmiles(), cache);
    
    double sum = 0.0;
    for (size_t i = 0; i < cache.elements.size(); ++i) {
        if (cache.elements[i] == "C" && i < cache.paulingValues.size()) {
            sum += cache.paulingValues[i];
        }
    }
    return sum;
}

// 27. Aromatic to Aliphatic Carbon EN Ratio
DECLARE_DESCRIPTOR(SmilesAromaticToAliphaticCarbonENRatio, ElectronegativityDescriptor, "Ratio of aromatic to aliphatic carbon EN sums")
DESCRIPTOR_DEPENDENCIES(SmilesAromaticToAliphaticCarbonENRatio) { return {}; }
DescriptorResult SmilesAromaticToAliphaticCarbonENRatioDescriptor::calculate(Context& context) const {
    auto& cache = tlsElementCache;
    extractElementData(context.getSmiles(), cache);
    tlsElementGroups.analyze(cache.elements);
    
    // Find aromatic and aliphatic carbon indices
    std::vector<size_t> aromaticCarbons, aliphaticCarbons;
    
    for (size_t i = 0; i < cache.elements.size(); ++i) {
        if (cache.elements[i] == "c") aromaticCarbons.push_back(i);
        else if (cache.elements[i] == "C") aliphaticCarbons.push_back(i);
    }
    
    double aromaticSum = sumENForElementIndices(cache.paulingValues, aromaticCarbons);
    double aliphaticSum = sumENForElementIndices(cache.paulingValues, aliphaticCarbons);
    
    return (aliphaticSum != 0.0) ? aromaticSum / aliphaticSum : 0.0;
}

// Helper for element type EN ratio calculations
double calculateElementTypeENRatio(
    const std::string& smiles,
    std::function<bool(const std::string&)> numeratorFilter,
    std::function<bool(const std::string&)> denominatorFilter) {
    
    auto& cache = tlsElementCache;
    extractElementData(smiles, cache);
    
    double numeratorSum = 0.0;
    double denominatorSum = 0.0;
    
    for (size_t i = 0; i < cache.elements.size(); ++i) {
        if (i >= cache.paulingValues.size()) continue;
        
        if (numeratorFilter(cache.elements[i])) {
            numeratorSum += cache.paulingValues[i];
        }
        
        if (denominatorFilter(cache.elements[i])) {
            denominatorSum += cache.paulingValues[i];
        }
    }
    
    return (denominatorSum != 0.0) ? numeratorSum / denominatorSum : 0.0;
}

// 28. Nitrogen to Carbon EN Ratio
DECLARE_DESCRIPTOR(SmilesNitrogenToCarbonENRatio, ElectronegativityDescriptor, "Ratio of nitrogen to carbon EN sums")
DESCRIPTOR_DEPENDENCIES(SmilesNitrogenToCarbonENRatio) { return {}; }
DescriptorResult SmilesNitrogenToCarbonENRatioDescriptor::calculate(Context& context) const {
    return calculateElementTypeENRatio(
        context.getSmiles(),
        ElementType::isNitrogen,
        ElementType::isCarbon
    );
}

// 29. Oxygen to Carbon EN Ratio
DECLARE_DESCRIPTOR(SmilesOxygenToCarbonENRatio, ElectronegativityDescriptor, "Ratio of oxygen to carbon EN sums")
DESCRIPTOR_DEPENDENCIES(SmilesOxygenToCarbonENRatio) { return {}; }
DescriptorResult SmilesOxygenToCarbonENRatioDescriptor::calculate(Context& context) const {
    return calculateElementTypeENRatio(
        context.getSmiles(),
        ElementType::isOxygen,
        ElementType::isCarbon
    );
}

// 30. Heteroatom to Carbon EN Ratio
DECLARE_DESCRIPTOR(SmilesHeteroatomToCarbonENRatio, ElectronegativityDescriptor, "Ratio of heteroatom to carbon EN sums")
DESCRIPTOR_DEPENDENCIES(SmilesHeteroatomToCarbonENRatio) { return {}; }
DescriptorResult SmilesHeteroatomToCarbonENRatioDescriptor::calculate(Context& context) const {
    return calculateElementTypeENRatio(
        context.getSmiles(),
        [](const std::string& elem) { 
            return (elem != "C" && elem != "c" && elem != "H" && elem != "h"); 
        },
        ElementType::isCarbon
    );
}

// 31-34. Normalized EN Sums (identical to Mean values)
DECLARE_DESCRIPTOR(SmilesNormalizedPaulingENSum, ElectronegativityDescriptor, "Pauling EN sum normalized by atom count (same as Mean)")
DESCRIPTOR_DEPENDENCIES(SmilesNormalizedPaulingENSum) { return {}; }
DescriptorResult SmilesNormalizedPaulingENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return meanElectronegativity(tlsElementCache.paulingValues);
}

DECLARE_DESCRIPTOR(SmilesNormalizedAllredRochowENSum, ElectronegativityDescriptor, "Allred-Rochow EN sum normalized by atom count (same as Mean)")
DESCRIPTOR_DEPENDENCIES(SmilesNormalizedAllredRochowENSum) { return {}; }
DescriptorResult SmilesNormalizedAllredRochowENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return meanElectronegativity(tlsElementCache.arValues);
}

DECLARE_DESCRIPTOR(SmilesNormalizedMullikenENSum, ElectronegativityDescriptor, "Mulliken EN sum normalized by atom count (same as Mean)")
DESCRIPTOR_DEPENDENCIES(SmilesNormalizedMullikenENSum) { return {}; }
DescriptorResult SmilesNormalizedMullikenENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return meanElectronegativity(tlsElementCache.mullikenValues);
}

DECLARE_DESCRIPTOR(SmilesNormalizedSandersonENSum, ElectronegativityDescriptor, "Sanderson EN sum normalized by atom count (same as Mean)")
DESCRIPTOR_DEPENDENCIES(SmilesNormalizedSandersonENSum) { return {}; }
DescriptorResult SmilesNormalizedSandersonENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return meanElectronegativity(tlsElementCache.sandersonValues);
}

// 35. Electronegative Element Ratio (EN > 2.5)
DECLARE_DESCRIPTOR(SmilesElectronegativeElementRatio, ElectronegativityDescriptor, "Fraction of elements with Pauling EN > 2.5")
DESCRIPTOR_DEPENDENCIES(SmilesElectronegativeElementRatio) { return {}; }
DescriptorResult SmilesElectronegativeElementRatioDescriptor::calculate(Context& context) const {
    auto& cache = tlsElementCache;
    extractElementData(context.getSmiles(), cache);
    if (cache.elements.empty()) return 0.0;
    
    double electronegativeCount = countElementsWithENAbove(cache.paulingValues, 2.5);
    return electronegativeCount / static_cast<double>(cache.elements.size());
}

// 36. Electronegativity Diversity Index
DECLARE_DESCRIPTOR(SmilesENDiversityIndex, ElectronegativityDescriptor, "Fraction of unique Pauling EN values (rounded to 1 decimal)")
DESCRIPTOR_DEPENDENCIES(SmilesENDiversityIndex) { return {}; }
DescriptorResult SmilesENDiversityIndexDescriptor::calculate(Context& context) const {
    auto& cache = tlsElementCache;
    extractElementData(context.getSmiles(), cache);
    if (cache.elements.empty()) return 0.0;
    
    // Count unique EN values (rounded to 1 decimal place)
    std::set<int> uniqueValues;
    for (double val : cache.paulingValues) {
        uniqueValues.insert(static_cast<int>(std::round(val * 10.0)));
    }
    
    return static_cast<double>(uniqueValues.size()) / static_cast<double>(cache.elements.size());
}

// 37. Highly Electronegative Elements Weight (Pauling > 3.0 weighted double)
DECLARE_DESCRIPTOR(SmilesHighENElementsWeight, ElectronegativityDescriptor, "Mean EN with atoms EN>3.0 weighted double")
DESCRIPTOR_DEPENDENCIES(SmilesHighENElementsWeight) { return {}; }
DescriptorResult SmilesHighENElementsWeightDescriptor::calculate(Context& context) const {
    auto& cache = tlsElementCache;
    extractElementData(context.getSmiles(), cache);
    if (cache.paulingValues.empty()) return 0.0;
    
    double weightedSum = 0.0;
    double totalWeight = 0.0;
    for (double val : cache.paulingValues) {
        double weight = (val > 3.0) ? 2.0 : 1.0;
        weightedSum += val * weight;
        totalWeight += weight;
    }
    
    return (totalWeight > 0.0) ? weightedSum / totalWeight : 0.0;
}

// 38. EN Homogeneity Index (1 - Avg Absolute Deviation / Mean)
DECLARE_DESCRIPTOR(SmilesENHomogeneityIndex, ElectronegativityDescriptor, "Measure of Pauling electronegativity homogeneity")
DESCRIPTOR_DEPENDENCIES(SmilesENHomogeneityIndex) { return {}; }
DescriptorResult SmilesENHomogeneityIndexDescriptor::calculate(Context& context) const {
    auto& cache = tlsElementCache;
    extractElementData(context.getSmiles(), cache);
    if (cache.paulingValues.empty()) return 0.0;
    
    double mean = meanElectronegativity(cache.paulingValues);
    if (mean == 0.0) return 0.0;
    
    double sumAbsDifferences = std::accumulate(cache.paulingValues.begin(), cache.paulingValues.end(), 0.0,
        [mean](double acc, double val) { return acc + std::abs(val - mean); });
    
    double avgAbsDeviation = sumAbsDifferences / static_cast<double>(cache.paulingValues.size());
    return 1.0 - (avgAbsDeviation / mean);
}

// 39. Ratio of EN Sum to Molecular Size (using SMILES length as proxy)
DECLARE_DESCRIPTOR(SmilesENSumToSizeRatio, ElectronegativityDescriptor, "Ratio of Pauling EN sum to SMILES string length")
DESCRIPTOR_DEPENDENCIES(SmilesENSumToSizeRatio) { return {}; }
DescriptorResult SmilesENSumToSizeRatioDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    auto& cache = tlsElementCache;
    extractElementData(smiles, cache);
    double enSum = sumElectronegativity(cache.paulingValues);
    
    return enSum / static_cast<double>(smiles.length());
}

// 40. Squared Electronegativity Sum
DECLARE_DESCRIPTOR(SmilesSquaredENSum, ElectronegativityDescriptor, "Sum of squared Pauling electronegativity values")
DESCRIPTOR_DEPENDENCIES(SmilesSquaredENSum) { return {}; }
DescriptorResult SmilesSquaredENSumDescriptor::calculate(Context& context) const {
    auto& cache = tlsElementCache;
    extractElementData(context.getSmiles(), cache);
    
    return std::accumulate(cache.paulingValues.begin(), cache.paulingValues.end(), 0.0,
        [](double acc, double val) { return acc + val * val; });
}

// 41. Geometric Mean of Pauling EN
DECLARE_DESCRIPTOR(SmilesGeometricMeanEN, ElectronegativityDescriptor, "Geometric mean of Pauling electronegativity values")
DESCRIPTOR_DEPENDENCIES(SmilesGeometricMeanEN) { return {}; }
DescriptorResult SmilesGeometricMeanENDescriptor::calculate(Context& context) const {
    auto& cache = tlsElementCache;
    extractElementData(context.getSmiles(), cache);
    if (cache.paulingValues.empty()) return 0.0;
    
    // Filter out non-positive values before taking log
    double logSum = 0.0;
    size_t validCount = 0;
    for (double val : cache.paulingValues) {
        if (val > 0.0) {
            logSum += std::log(val);
            validCount++;
        }
    }
    
    if (validCount == 0) return 0.0;
    return std::exp(logSum / static_cast<double>(validCount));
}

// 42. Harmonic Mean of Pauling EN
DECLARE_DESCRIPTOR(SmilesHarmonicMeanEN, ElectronegativityDescriptor, "Harmonic mean of Pauling electronegativity values")
DESCRIPTOR_DEPENDENCIES(SmilesHarmonicMeanEN) { return {}; }
DescriptorResult SmilesHarmonicMeanENDescriptor::calculate(Context& context) const {
    auto& cache = tlsElementCache;
    extractElementData(context.getSmiles(), cache);
    if (cache.paulingValues.empty()) return 0.0;
    
    double reciprocalSum = 0.0;
    size_t validCount = 0;
    for (double val : cache.paulingValues) {
        if (val > 0.0) {
            reciprocalSum += 1.0 / val;
            validCount++;
        } else {
            // Harmonic mean is undefined if any value is zero or less
            return 0.0;
        }
    }
    
    if (validCount == 0 || reciprocalSum == 0.0) return 0.0;
    return static_cast<double>(validCount) / reciprocalSum;
}

// 43. Median Pauling EN
DECLARE_DESCRIPTOR(SmilesMedianEN, ElectronegativityDescriptor, "Median of Pauling electronegativity values")
DESCRIPTOR_DEPENDENCIES(SmilesMedianEN) { return {}; }
DescriptorResult SmilesMedianENDescriptor::calculate(Context& context) const {
    auto& cache = tlsElementCache;
    extractElementData(context.getSmiles(), cache);
    if (cache.paulingValues.empty()) return 0.0;
    
    // Make a copy for sorting (can't sort the original as it may be used by other threads)
    std::vector<double> sortedValues = cache.paulingValues;
    std::sort(sortedValues.begin(), sortedValues.end());
    
    size_t n = sortedValues.size();
    if (n % 2 == 0) {
        // Even number of elements: average of the two middle ones
        return (sortedValues[n/2 - 1] + sortedValues[n/2]) / 2.0;
    } else {
        // Odd number of elements: the middle one
        return sortedValues[n/2];
    }
}

// 44. Mode Pauling EN
DECLARE_DESCRIPTOR(SmilesModeEN, ElectronegativityDescriptor, "Mode of Pauling electronegativity values (rounded to 2 decimals)")
DESCRIPTOR_DEPENDENCIES(SmilesModeEN) { return {}; }
DescriptorResult SmilesModeENDescriptor::calculate(Context& context) const {
    auto& cache = tlsElementCache;
    extractElementData(context.getSmiles(), cache);
    if (cache.paulingValues.empty()) return 0.0;
    
    // Round to 2 decimal places and count frequencies
    std::unordered_map<int, int> valueCounts;
    for (double val : cache.paulingValues) {
        valueCounts[static_cast<int>(std::round(val * 100.0))]++;
    }
    
    int maxCount = 0;
    int modeValueRounded = 0;
    
    for (const auto& pair : valueCounts) {
        if (pair.second > maxCount) {
            maxCount = pair.second;
            modeValueRounded = pair.first;
        }
    }
    
    return static_cast<double>(modeValueRounded) / 100.0;
}

// 45. Carbon to Heteroatom EN Ratio
DECLARE_DESCRIPTOR(SmilesCarbonToHeteroatomENRatio, ElectronegativityDescriptor, "Ratio of carbon to heteroatom EN sums")
DESCRIPTOR_DEPENDENCIES(SmilesCarbonToHeteroatomENRatio) { return {}; }
DescriptorResult SmilesCarbonToHeteroatomENRatioDescriptor::calculate(Context& context) const {
    return calculateElementTypeENRatio(
        context.getSmiles(),
        ElementType::isCarbon,
        [](const std::string& elem) { 
            return (elem != "C" && elem != "c" && elem != "H" && elem != "h"); 
        }
    );
}

// 46. EN-Weighted Atom Count (Sum of Pauling EN values) - Equivalent to SmilesPaulingENSum
DECLARE_DESCRIPTOR(SmilesENWeightedAtomCount, ElectronegativityDescriptor, "Sum of atom counts weighted by EN (same as Sum)")
DESCRIPTOR_DEPENDENCIES(SmilesENWeightedAtomCount) { return {}; }
DescriptorResult SmilesENWeightedAtomCountDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    return sumElectronegativity(tlsElementCache.paulingValues);
}

// 47. EN Polarity Index ((MaxEN - MinEN) * MeanEN / AtomCount)
DECLARE_DESCRIPTOR(SmilesENPolarityIndex, ElectronegativityDescriptor, "Index of molecular polarity based on EN differences")
DESCRIPTOR_DEPENDENCIES(SmilesENPolarityIndex) { return {}; }
DescriptorResult SmilesENPolarityIndexDescriptor::calculate(Context& context) const {
    auto& cache = tlsElementCache;
    extractElementData(context.getSmiles(), cache);
    if (cache.paulingValues.size() < 2) return 0.0;
    
    double range = electronegativityRange(cache.paulingValues);
    double meanEN = meanElectronegativity(cache.paulingValues);
    
    return range * meanEN / static_cast<double>(cache.paulingValues.size());
}

// 48. EN-Weighted Mean Molecular Weight
DECLARE_DESCRIPTOR(SmilesENWeightedMolWeight, ElectronegativityDescriptor, "Mean molecular weight weighted by Pauling EN")
DESCRIPTOR_DEPENDENCIES(SmilesENWeightedMolWeight) { return {}; }
DescriptorResult SmilesENWeightedMolWeightDescriptor::calculate(Context& context) const {
    auto& cache = tlsElementCache;
    extractElementData(context.getSmiles(), cache);
    
    double totalWeightedWeight = 0.0;
    double totalEN = 0.0;
    
    for (size_t i = 0; i < cache.elements.size(); ++i) {
        if (i >= cache.paulingValues.size()) continue;
        
        auto weightIt = atomicWeights.find(cache.elements[i]);
        if (weightIt != atomicWeights.end()) {
            double en = cache.paulingValues[i];
            if (en > 0) {
                totalWeightedWeight += weightIt->second * en;
                totalEN += en;
            }
        }
    }
    
    return (totalEN > 0.0) ? totalWeightedWeight / totalEN : 0.0;
}

// 49. Multi-Scale EN Correlation (Pearson correlation between Pauling and Mulliken)
DECLARE_DESCRIPTOR(SmilesMultiScaleENCorrelation, ElectronegativityDescriptor, "Correlation between Pauling and Mulliken EN values")
DESCRIPTOR_DEPENDENCIES(SmilesMultiScaleENCorrelation) { return {}; }
DescriptorResult SmilesMultiScaleENCorrelationDescriptor::calculate(Context& context) const {
    auto& cache = tlsElementCache;
    extractElementData(context.getSmiles(), cache);

    std::vector<double> pauling, mulliken;
    for (size_t i = 0; i < cache.elements.size(); ++i) {
        const std::string& elem = cache.elements[i];
        auto pIt = paulingElectronegativityMap.find(elem);
        auto mIt = mullikenENMap.find(elem);
        if (pIt != paulingElectronegativityMap.end() && mIt != mullikenENMap.end()) {
            pauling.push_back(pIt->second);
            mulliken.push_back(mIt->second);
        }
    }
    if (pauling.size() < 2) return 0.0;

    double sumX = 0.0, sumY = 0.0, sumXY = 0.0, sumX2 = 0.0, sumY2 = 0.0;
    size_t n = pauling.size();
    for (size_t i = 0; i < n; ++i) {
        sumX += pauling[i];
        sumY += mulliken[i];
        sumXY += pauling[i] * mulliken[i];
        sumX2 += pauling[i] * pauling[i];
        sumY2 += mulliken[i] * mulliken[i];
    }
    double numerator = n * sumXY - sumX * sumY;
    double denominator = std::sqrt((n * sumX2 - sumX * sumX) * (n * sumY2 - sumY * sumY));
    if (denominator <= 0) return 1.0;
    return numerator / denominator;
}

// 50. EN Heterogeneity Index
DECLARE_DESCRIPTOR(SmilesENHeterogeneityIndex, ElectronegativityDescriptor, "Measure of Pauling electronegativity heterogeneity")
DESCRIPTOR_DEPENDENCIES(SmilesENHeterogeneityIndex) { return {}; }
DescriptorResult SmilesENHeterogeneityIndexDescriptor::calculate(Context& context) const {
    auto& cache = tlsElementCache;
    extractElementData(context.getSmiles(), cache);
    if (cache.elements.empty()) return 0.0;
    
    // Count distinct element types
    std::set<std::string> uniqueElements(cache.elements.begin(), cache.elements.end());
    
    // Get EN values for unique elements
    std::vector<double> uniqueENValues;
    for (const auto& element : uniqueElements) {
        auto it = paulingElectronegativityMap.find(element);
        if (it != paulingElectronegativityMap.end()) {
            uniqueENValues.push_back(it->second);
        }
    }
    
    if (uniqueENValues.empty()) return 0.0;
    
    // Calculate heterogeneity index
    double enRange = electronegativityRange(uniqueENValues);
    return enRange * static_cast<double>(uniqueElements.size()) / static_cast<double>(cache.elements.size());
}

// --- Additional Electronegativity Descriptors ---

// 51. Count of unique elements
DECLARE_DESCRIPTOR(SmilesUniqueElementCount, ElectronegativityDescriptor, "Count of unique elements in SMILES")
DESCRIPTOR_DEPENDENCIES(SmilesUniqueElementCount) { return {}; }
DescriptorResult SmilesUniqueElementCountDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::set<std::string> uniqueElems(tlsElementCache.elements.begin(), tlsElementCache.elements.end());
    return static_cast<double>(uniqueElems.size());
}

// 58. EN mean of unique elements
DECLARE_DESCRIPTOR(SmilesUniqueElementENMean, ElectronegativityDescriptor, "Mean Pauling EN for unique elements")
DESCRIPTOR_DEPENDENCIES(SmilesUniqueElementENMean) { return {}; }
DescriptorResult SmilesUniqueElementENMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::set<std::string> uniqueElems(tlsElementCache.elements.begin(), tlsElementCache.elements.end());
    double sum = 0.0;
    for (const auto& e : uniqueElems) {
        auto it = paulingElectronegativityMap.find(e);
        if (it != paulingElectronegativityMap.end()) sum += it->second;
    }
    return uniqueElems.empty() ? 0.0 : sum / uniqueElems.size();
}

// 59. EN range of unique elements
DECLARE_DESCRIPTOR(SmilesUniqueElementENRange, ElectronegativityDescriptor, "Range of Pauling EN for unique elements")
DESCRIPTOR_DEPENDENCIES(SmilesUniqueElementENRange) { return {}; }
DescriptorResult SmilesUniqueElementENRangeDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::set<double> uniqueEN;
    for (const auto& e : tlsElementCache.elements) {
        auto it = paulingElectronegativityMap.find(e);
        if (it != paulingElectronegativityMap.end()) uniqueEN.insert(it->second);
    }
    if (uniqueEN.size() < 2) return 0.0;
    return *uniqueEN.rbegin() - *uniqueEN.begin();
}

// 60. EN stddev of unique elements
DECLARE_DESCRIPTOR(SmilesUniqueElementENStdDev, ElectronegativityDescriptor, "Stddev of Pauling EN for unique elements")
DESCRIPTOR_DEPENDENCIES(SmilesUniqueElementENStdDev) { return {}; }
DescriptorResult SmilesUniqueElementENStdDevDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::set<double> uniqueEN;
    for (const auto& e : tlsElementCache.elements) {
        auto it = paulingElectronegativityMap.find(e);
        if (it != paulingElectronegativityMap.end()) uniqueEN.insert(it->second);
    }
    if (uniqueEN.size() < 2) return 0.0;
    double mean = std::accumulate(uniqueEN.begin(), uniqueEN.end(), 0.0) / uniqueEN.size();
    double var = 0.0;
    for (auto v : uniqueEN) var += (v - mean) * (v - mean);
    return std::sqrt(var / uniqueEN.size());
}

// 61. EN sum of aromatic atoms
DECLARE_DESCRIPTOR(SmilesAromaticENSum, ElectronegativityDescriptor, "Sum of Pauling EN for aromatic atoms")
DESCRIPTOR_DEPENDENCIES(SmilesAromaticENSum) { return {}; }
DescriptorResult SmilesAromaticENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (ElementType::isAromatic(tlsElementCache.elements[i]) && i < tlsElementCache.paulingValues.size())
            sum += tlsElementCache.paulingValues[i];
    return sum;
}

// 62. EN sum of aliphatic atoms
DECLARE_DESCRIPTOR(SmilesAliphaticENSum, ElectronegativityDescriptor, "Sum of Pauling EN for aliphatic atoms")
DESCRIPTOR_DEPENDENCIES(SmilesAliphaticENSum) { return {}; }
DescriptorResult SmilesAliphaticENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!ElementType::isAromatic(tlsElementCache.elements[i]) && i < tlsElementCache.paulingValues.size())
            sum += tlsElementCache.paulingValues[i];
    return sum;
}

// 63. Aromatic to aliphatic EN ratio
DECLARE_DESCRIPTOR(SmilesAromaticToAliphaticENRatio, ElectronegativityDescriptor, "Ratio of aromatic to aliphatic EN sums")
DESCRIPTOR_DEPENDENCIES(SmilesAromaticToAliphaticENRatio) { return {}; }
DescriptorResult SmilesAromaticToAliphaticENRatioDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double aromatic = 0.0, aliphatic = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i) {
        if (ElementType::isAromatic(tlsElementCache.elements[i]) && i < tlsElementCache.paulingValues.size())
            aromatic += tlsElementCache.paulingValues[i];
        else if (!ElementType::isAromatic(tlsElementCache.elements[i]) && i < tlsElementCache.paulingValues.size())
            aliphatic += tlsElementCache.paulingValues[i];
    }
    return aliphatic == 0.0 ? 0.0 : aromatic / aliphatic;
}

// 64. EN sum of halogens
DECLARE_DESCRIPTOR(SmilesHalogenENSum, ElectronegativityDescriptor, "Sum of Pauling EN for halogen atoms")
DESCRIPTOR_DEPENDENCIES(SmilesHalogenENSum) { return {}; }
DescriptorResult SmilesHalogenENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (ElementType::isHalogen(tlsElementCache.elements[i]) && i < tlsElementCache.paulingValues.size())
            sum += tlsElementCache.paulingValues[i];
    return sum;
}

// 65. EN mean of halogens
DECLARE_DESCRIPTOR(SmilesHalogenENMean, ElectronegativityDescriptor, "Mean Pauling EN for halogen atoms")
DESCRIPTOR_DEPENDENCIES(SmilesHalogenENMean) { return {}; }
DescriptorResult SmilesHalogenENMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0, count = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (ElementType::isHalogen(tlsElementCache.elements[i]) && i < tlsElementCache.paulingValues.size()) {
            sum += tlsElementCache.paulingValues[i];
            count += 1.0;
        }
    return count == 0.0 ? 0.0 : sum / count;
}

// 66. EN sum of non-halogens
DECLARE_DESCRIPTOR(SmilesNonHalogenENSum, ElectronegativityDescriptor, "Sum of Pauling EN for non-halogen atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonHalogenENSum) { return {}; }
DescriptorResult SmilesNonHalogenENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!ElementType::isHalogen(tlsElementCache.elements[i]) && i < tlsElementCache.paulingValues.size())
            sum += tlsElementCache.paulingValues[i];
    return sum;
}

// 67. EN mean of non-halogens
DECLARE_DESCRIPTOR(SmilesNonHalogenENMean, ElectronegativityDescriptor, "Mean Pauling EN for non-halogen atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonHalogenENMean) { return {}; }
DescriptorResult SmilesNonHalogenENMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0, count = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!ElementType::isHalogen(tlsElementCache.elements[i]) && i < tlsElementCache.paulingValues.size()) {
            sum += tlsElementCache.paulingValues[i];
            count += 1.0;
        }
    return count == 0.0 ? 0.0 : sum / count;
}

// 68. EN range of non-halogens
DECLARE_DESCRIPTOR(SmilesNonHalogenENRange, ElectronegativityDescriptor, "Range of Pauling EN for non-halogen atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonHalogenENRange) { return {}; }
DescriptorResult SmilesNonHalogenENRangeDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!ElementType::isHalogen(tlsElementCache.elements[i]) && i < tlsElementCache.paulingValues.size())
            vals.push_back(tlsElementCache.paulingValues[i]);
    if (vals.size() < 2) return 0.0;
    auto minmax = std::minmax_element(vals.begin(), vals.end());
    return *minmax.second - *minmax.first;
}

// 69. EN stddev of non-halogens
DECLARE_DESCRIPTOR(SmilesNonHalogenENStdDev, ElectronegativityDescriptor, "Stddev of Pauling EN for non-halogen atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonHalogenENStdDev) { return {}; }
DescriptorResult SmilesNonHalogenENStdDevDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!ElementType::isHalogen(tlsElementCache.elements[i]) && i < tlsElementCache.paulingValues.size())
            vals.push_back(tlsElementCache.paulingValues[i]);
    if (vals.size() < 2) return 0.0;
    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    double var = 0.0;
    for (auto v : vals) var += (v - mean) * (v - mean);
    return std::sqrt(var / vals.size());
}

// 70. EN sum of non-carbons
DECLARE_DESCRIPTOR(SmilesNonCarbonENSum, ElectronegativityDescriptor, "Sum of Pauling EN for non-carbon atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonCarbonENSum) { return {}; }
DescriptorResult SmilesNonCarbonENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!(tlsElementCache.elements[i] == "C" || tlsElementCache.elements[i] == "c") && i < tlsElementCache.paulingValues.size())
            sum += tlsElementCache.paulingValues[i];
    return sum;
}

// 71. EN mean of non-carbons
DECLARE_DESCRIPTOR(SmilesNonCarbonENMean, ElectronegativityDescriptor, "Mean Pauling EN for non-carbon atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonCarbonENMean) { return {}; }
DescriptorResult SmilesNonCarbonENMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0, count = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!(tlsElementCache.elements[i] == "C" || tlsElementCache.elements[i] == "c") && i < tlsElementCache.paulingValues.size()) {
            sum += tlsElementCache.paulingValues[i];
            count += 1.0;
        }
    return count == 0.0 ? 0.0 : sum / count;
}

// 72. EN range of non-carbons
DECLARE_DESCRIPTOR(SmilesNonCarbonENRange, ElectronegativityDescriptor, "Range of Pauling EN for non-carbon atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonCarbonENRange) { return {}; }
DescriptorResult SmilesNonCarbonENRangeDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!(tlsElementCache.elements[i] == "C" || tlsElementCache.elements[i] == "c") && i < tlsElementCache.paulingValues.size())
            vals.push_back(tlsElementCache.paulingValues[i]);
    if (vals.size() < 2) return 0.0;
    auto minmax = std::minmax_element(vals.begin(), vals.end());
    return *minmax.second - *minmax.first;
}

// 73. EN stddev of non-carbons
DECLARE_DESCRIPTOR(SmilesNonCarbonENStdDev, ElectronegativityDescriptor, "Stddev of Pauling EN for non-carbon atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonCarbonENStdDev) { return {}; }
DescriptorResult SmilesNonCarbonENStdDevDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!(tlsElementCache.elements[i] == "C" || tlsElementCache.elements[i] == "c") && i < tlsElementCache.paulingValues.size())
            vals.push_back(tlsElementCache.paulingValues[i]);
    if (vals.size() < 2) return 0.0;
    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    double var = 0.0;
    for (auto v : vals) var += (v - mean) * (v - mean);
    return std::sqrt(var / vals.size());
}

// 74. EN sum of non-nitrogens
DECLARE_DESCRIPTOR(SmilesNonNitrogenENSum, ElectronegativityDescriptor, "Sum of Pauling EN for non-nitrogen atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonNitrogenENSum) { return {}; }
DescriptorResult SmilesNonNitrogenENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!(tlsElementCache.elements[i] == "N" || tlsElementCache.elements[i] == "n") && i < tlsElementCache.paulingValues.size())
            sum += tlsElementCache.paulingValues[i];
    return sum;
}

// 75. EN mean of non-nitrogens
DECLARE_DESCRIPTOR(SmilesNonNitrogenENMean, ElectronegativityDescriptor, "Mean Pauling EN for non-nitrogen atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonNitrogenENMean) { return {}; }
DescriptorResult SmilesNonNitrogenENMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0, count = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!(tlsElementCache.elements[i] == "N" || tlsElementCache.elements[i] == "n") && i < tlsElementCache.paulingValues.size()) {
            sum += tlsElementCache.paulingValues[i];
            count += 1.0;
        }
    return count == 0.0 ? 0.0 : sum / count;
}

// 76. EN range of non-nitrogens
DECLARE_DESCRIPTOR(SmilesNonNitrogenENRange, ElectronegativityDescriptor, "Range of Pauling EN for non-nitrogen atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonNitrogenENRange) { return {}; }
DescriptorResult SmilesNonNitrogenENRangeDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!(tlsElementCache.elements[i] == "N" || tlsElementCache.elements[i] == "n") && i < tlsElementCache.paulingValues.size())
            vals.push_back(tlsElementCache.paulingValues[i]);
    if (vals.size() < 2) return 0.0;
    auto minmax = std::minmax_element(vals.begin(), vals.end());
    return *minmax.second - *minmax.first;
}

// 77. EN stddev of non-nitrogens
DECLARE_DESCRIPTOR(SmilesNonNitrogenENStdDev, ElectronegativityDescriptor, "Stddev of Pauling EN for non-nitrogen atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonNitrogenENStdDev) { return {}; }
DescriptorResult SmilesNonNitrogenENStdDevDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!(tlsElementCache.elements[i] == "N" || tlsElementCache.elements[i] == "n") && i < tlsElementCache.paulingValues.size())
            vals.push_back(tlsElementCache.paulingValues[i]);
    if (vals.size() < 2) return 0.0;
    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    double var = 0.0;
    for (auto v : vals) var += (v - mean) * (v - mean);
    return std::sqrt(var / vals.size());
}

// 78. EN sum of non-oxygens
DECLARE_DESCRIPTOR(SmilesNonOxygenENSum, ElectronegativityDescriptor, "Sum of Pauling EN for non-oxygen atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonOxygenENSum) { return {}; }
DescriptorResult SmilesNonOxygenENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!(tlsElementCache.elements[i] == "O" || tlsElementCache.elements[i] == "o") && i < tlsElementCache.paulingValues.size())
            sum += tlsElementCache.paulingValues[i];
    return sum;
}

// 79. EN mean of non-oxygens
DECLARE_DESCRIPTOR(SmilesNonOxygenENMean, ElectronegativityDescriptor, "Mean Pauling EN for non-oxygen atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonOxygenENMean) { return {}; }
DescriptorResult SmilesNonOxygenENMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0, count = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!(tlsElementCache.elements[i] == "O" || tlsElementCache.elements[i] == "o") && i < tlsElementCache.paulingValues.size()) {
            sum += tlsElementCache.paulingValues[i];
            count += 1.0;
        }
    return count == 0.0 ? 0.0 : sum / count;
}

// 80. EN range of non-oxygens
DECLARE_DESCRIPTOR(SmilesNonOxygenENRange, ElectronegativityDescriptor, "Range of Pauling EN for non-oxygen atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonOxygenENRange) { return {}; }
DescriptorResult SmilesNonOxygenENRangeDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!(tlsElementCache.elements[i] == "O" || tlsElementCache.elements[i] == "o") && i < tlsElementCache.paulingValues.size())
            vals.push_back(tlsElementCache.paulingValues[i]);
    if (vals.size() < 2) return 0.0;
    auto minmax = std::minmax_element(vals.begin(), vals.end());
    return *minmax.second - *minmax.first;
}

// 81. EN stddev of non-oxygens
DECLARE_DESCRIPTOR(SmilesNonOxygenENStdDev, ElectronegativityDescriptor, "Stddev of Pauling EN for non-oxygen atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonOxygenENStdDev) { return {}; }
DescriptorResult SmilesNonOxygenENStdDevDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!(tlsElementCache.elements[i] == "O" || tlsElementCache.elements[i] == "o") && i < tlsElementCache.paulingValues.size())
            vals.push_back(tlsElementCache.paulingValues[i]);
    if (vals.size() < 2) return 0.0;
    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    double var = 0.0;
    for (auto v : vals) var += (v - mean) * (v - mean);
    return std::sqrt(var / vals.size());
}

// 82. EN sum of non-sulfurs
DECLARE_DESCRIPTOR(SmilesNonSulfurENSum, ElectronegativityDescriptor, "Sum of Pauling EN for non-sulfur atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonSulfurENSum) { return {}; }
DescriptorResult SmilesNonSulfurENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!(tlsElementCache.elements[i] == "S" || tlsElementCache.elements[i] == "s") && i < tlsElementCache.paulingValues.size())
            sum += tlsElementCache.paulingValues[i];
    return sum;
}

// 83. EN mean of non-sulfurs
DECLARE_DESCRIPTOR(SmilesNonSulfurENMean, ElectronegativityDescriptor, "Mean Pauling EN for non-sulfur atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonSulfurENMean) { return {}; }
DescriptorResult SmilesNonSulfurENMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0, count = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!(tlsElementCache.elements[i] == "S" || tlsElementCache.elements[i] == "s") && i < tlsElementCache.paulingValues.size()) {
            sum += tlsElementCache.paulingValues[i];
            count += 1.0;
        }
    return count == 0.0 ? 0.0 : sum / count;
}

// 84. EN range of non-sulfurs
DECLARE_DESCRIPTOR(SmilesNonSulfurENRange, ElectronegativityDescriptor, "Range of Pauling EN for non-sulfur atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonSulfurENRange) { return {}; }
DescriptorResult SmilesNonSulfurENRangeDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!(tlsElementCache.elements[i] == "S" || tlsElementCache.elements[i] == "s") && i < tlsElementCache.paulingValues.size())
            vals.push_back(tlsElementCache.paulingValues[i]);
    if (vals.size() < 2) return 0.0;
    auto minmax = std::minmax_element(vals.begin(), vals.end());
    return *minmax.second - *minmax.first;
}

// 85. EN stddev of non-sulfurs
DECLARE_DESCRIPTOR(SmilesNonSulfurENStdDev, ElectronegativityDescriptor, "Stddev of Pauling EN for non-sulfur atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNonSulfurENStdDev) { return {}; }
DescriptorResult SmilesNonSulfurENStdDevDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (!(tlsElementCache.elements[i] == "S" || tlsElementCache.elements[i] == "s") && i < tlsElementCache.paulingValues.size())
            vals.push_back(tlsElementCache.paulingValues[i]);
    if (vals.size() < 2) return 0.0;
    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    double var = 0.0;
    for (auto v : vals) var += (v - mean) * (v - mean);
    return std::sqrt(var / vals.size());
}

// 86. EN sum of elements with EN > 2.5
DECLARE_DESCRIPTOR(SmilesHighENAtomSum, ElectronegativityDescriptor, "Sum of Pauling EN for atoms with EN > 2.5")
DESCRIPTOR_DEPENDENCIES(SmilesHighENAtomSum) { return {}; }
DescriptorResult SmilesHighENAtomSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v > 2.5) sum += v;
    return sum;
}

// 87. EN sum of elements with EN < 2.0
DECLARE_DESCRIPTOR(SmilesLowENAtomSum, ElectronegativityDescriptor, "Sum of Pauling EN for atoms with EN < 2.0")
DESCRIPTOR_DEPENDENCIES(SmilesLowENAtomSum) { return {}; }
DescriptorResult SmilesLowENAtomSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v < 2.0) sum += v;
    return sum;
}

// 88. EN mean of elements with EN > 2.5
DECLARE_DESCRIPTOR(SmilesHighENAtomMean, ElectronegativityDescriptor, "Mean Pauling EN for atoms with EN > 2.5")
DESCRIPTOR_DEPENDENCIES(SmilesHighENAtomMean) { return {}; }
DescriptorResult SmilesHighENAtomMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0, count = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v > 2.5) { sum += v; count += 1.0; }
    return count == 0.0 ? 0.0 : sum / count;
}

// 89. EN mean of elements with EN < 2.0
DECLARE_DESCRIPTOR(SmilesLowENAtomMean, ElectronegativityDescriptor, "Mean Pauling EN for atoms with EN < 2.0")
DESCRIPTOR_DEPENDENCIES(SmilesLowENAtomMean) { return {}; }
DescriptorResult SmilesLowENAtomMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0, count = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v < 2.0) { sum += v; count += 1.0; }
    return count == 0.0 ? 0.0 : sum / count;
}

// 90. EN range of elements with EN > 2.5
DECLARE_DESCRIPTOR(SmilesHighENAtomRange, ElectronegativityDescriptor, "Range of Pauling EN for atoms with EN > 2.5")
DESCRIPTOR_DEPENDENCIES(SmilesHighENAtomRange) { return {}; }
DescriptorResult SmilesHighENAtomRangeDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (double v : tlsElementCache.paulingValues) if (v > 2.5) vals.push_back(v);
    if (vals.size() < 2) return 0.0;
    auto minmax = std::minmax_element(vals.begin(), vals.end());
    return *minmax.second - *minmax.first;
}

// 91. EN stddev of elements with EN > 2.5
DECLARE_DESCRIPTOR(SmilesHighENAtomStdDev, ElectronegativityDescriptor, "Stddev of Pauling EN for atoms with EN > 2.5")
DESCRIPTOR_DEPENDENCIES(SmilesHighENAtomStdDev) { return {}; }
DescriptorResult SmilesHighENAtomStdDevDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (double v : tlsElementCache.paulingValues) if (v > 2.5) vals.push_back(v);
    if (vals.size() < 2) return 0.0;
    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    double var = 0.0;
    for (auto v : vals) var += (v - mean) * (v - mean);
    return std::sqrt(var / vals.size());
}

// 92. EN stddev of elements with EN < 2.0
DECLARE_DESCRIPTOR(SmilesLowENAtomStdDev, ElectronegativityDescriptor, "Stddev of Pauling EN for atoms with EN < 2.0")
DESCRIPTOR_DEPENDENCIES(SmilesLowENAtomStdDev) { return {}; }
DescriptorResult SmilesLowENAtomStdDevDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (double v : tlsElementCache.paulingValues) if (v < 2.0) vals.push_back(v);
    if (vals.size() < 2) return 0.0;
    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    double var = 0.0;
    for (auto v : vals) var += (v - mean) * (v - mean);
    return std::sqrt(var / vals.size());
}

// 93. EN sum of elements with EN between 2.0 and 3.0
DECLARE_DESCRIPTOR(SmilesMidENAtomSum, ElectronegativityDescriptor, "Sum of Pauling EN for atoms with 2.0 < EN < 3.0")
DESCRIPTOR_DEPENDENCIES(SmilesMidENAtomSum) { return {}; }
DescriptorResult SmilesMidENAtomSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v > 2.0 && v < 3.0) sum += v;
    return sum;
}

// 94. EN mean of elements with EN between 2.0 and 3.0
DECLARE_DESCRIPTOR(SmilesMidENAtomMean, ElectronegativityDescriptor, "Mean Pauling EN for atoms with 2.0 < EN < 3.0")
DESCRIPTOR_DEPENDENCIES(SmilesMidENAtomMean) { return {}; }
DescriptorResult SmilesMidENAtomMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0, count = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v > 2.0 && v < 3.0) { sum += v; count += 1.0; }
    return count == 0.0 ? 0.0 : sum / count;
}

// 95. EN range of elements with EN between 2.0 and 3.0
DECLARE_DESCRIPTOR(SmilesMidENAtomRange, ElectronegativityDescriptor, "Range of Pauling EN for atoms with 2.0 < EN < 3.0")
DESCRIPTOR_DEPENDENCIES(SmilesMidENAtomRange) { return {}; }
DescriptorResult SmilesMidENAtomRangeDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (double v : tlsElementCache.paulingValues) if (v > 2.0 && v < 3.0) vals.push_back(v);
    if (vals.size() < 2) return 0.0;
    auto minmax = std::minmax_element(vals.begin(), vals.end());
    return *minmax.second - *minmax.first;
}

// 96. EN stddev of elements with EN between 2.0 and 3.0
DECLARE_DESCRIPTOR(SmilesMidENAtomStdDev, ElectronegativityDescriptor, "Stddev of Pauling EN for atoms with 2.0 < EN < 3.0")
DESCRIPTOR_DEPENDENCIES(SmilesMidENAtomStdDev) { return {}; }
DescriptorResult SmilesMidENAtomStdDevDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (double v : tlsElementCache.paulingValues) if (v > 2.0 && v < 3.0) vals.push_back(v);
    if (vals.size() < 2) return 0.0;
    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    double var = 0.0;
    for (auto v : vals) var += (v - mean) * (v - mean);
    return std::sqrt(var / vals.size());
}

// 97. EN sum of elements with EN == 2.55 (carbon)
DECLARE_DESCRIPTOR(SmilesCarbonENSum, ElectronegativityDescriptor, "Sum of Pauling EN for carbon atoms")
DESCRIPTOR_DEPENDENCIES(SmilesCarbonENSum) { return {}; }
DescriptorResult SmilesCarbonENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if ((tlsElementCache.elements[i] == "C" || tlsElementCache.elements[i] == "c") && i < tlsElementCache.paulingValues.size())
            sum += tlsElementCache.paulingValues[i];
    return sum;
}

// 98. EN sum of elements with EN == 3.04 (nitrogen)
DECLARE_DESCRIPTOR(SmilesNitrogenENSum, ElectronegativityDescriptor, "Sum of Pauling EN for nitrogen atoms")
DESCRIPTOR_DEPENDENCIES(SmilesNitrogenENSum) { return {}; }
DescriptorResult SmilesNitrogenENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if ((tlsElementCache.elements[i] == "N" || tlsElementCache.elements[i] == "n") && i < tlsElementCache.paulingValues.size())
            sum += tlsElementCache.paulingValues[i];
    return sum;
}

// 99. EN sum of elements with EN == 3.44 (oxygen)
DECLARE_DESCRIPTOR(SmilesOxygenENSum, ElectronegativityDescriptor, "Sum of Pauling EN for oxygen atoms")
DESCRIPTOR_DEPENDENCIES(SmilesOxygenENSum) { return {}; }
DescriptorResult SmilesOxygenENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if ((tlsElementCache.elements[i] == "O" || tlsElementCache.elements[i] == "o") && i < tlsElementCache.paulingValues.size())
            sum += tlsElementCache.paulingValues[i];
    return sum;
}

// 100. EN sum of elements with EN == 2.58 (sulfur)
DECLARE_DESCRIPTOR(SmilesSulfurENSum, ElectronegativityDescriptor, "Sum of Pauling EN for sulfur atoms")
DESCRIPTOR_DEPENDENCIES(SmilesSulfurENSum) { return {}; }
DescriptorResult SmilesSulfurENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if ((tlsElementCache.elements[i] == "S" || tlsElementCache.elements[i] == "s") && i < tlsElementCache.paulingValues.size())
            sum += tlsElementCache.paulingValues[i];
    return sum;
}

// 101. EN sum of elements with EN == 3.98 (fluorine)
DECLARE_DESCRIPTOR(SmilesFluorineENSum, ElectronegativityDescriptor, "Sum of Pauling EN for fluorine atoms")
DESCRIPTOR_DEPENDENCIES(SmilesFluorineENSum) { return {}; }
DescriptorResult SmilesFluorineENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (tlsElementCache.elements[i] == "F" && i < tlsElementCache.paulingValues.size())
            sum += tlsElementCache.paulingValues[i];
    return sum;
}

// 102. EN sum of elements with EN == 3.16 (chlorine)
DECLARE_DESCRIPTOR(SmilesChlorineENSum, ElectronegativityDescriptor, "Sum of Pauling EN for chlorine atoms")
DESCRIPTOR_DEPENDENCIES(SmilesChlorineENSum) { return {}; }
DescriptorResult SmilesChlorineENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (tlsElementCache.elements[i] == "Cl" && i < tlsElementCache.paulingValues.size())
            sum += tlsElementCache.paulingValues[i];
    return sum;
}

// 103. EN sum of elements with EN == 2.96 (bromine)
DECLARE_DESCRIPTOR(SmilesBromineENSum, ElectronegativityDescriptor, "Sum of Pauling EN for bromine atoms")
DESCRIPTOR_DEPENDENCIES(SmilesBromineENSum) { return {}; }
DescriptorResult SmilesBromineENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (tlsElementCache.elements[i] == "Br" && i < tlsElementCache.paulingValues.size())
            sum += tlsElementCache.paulingValues[i];
    return sum;
}

// 104. EN sum of elements with EN == 2.66 (iodine)
DECLARE_DESCRIPTOR(SmilesIodineENSum, ElectronegativityDescriptor, "Sum of Pauling EN for iodine atoms")
DESCRIPTOR_DEPENDENCIES(SmilesIodineENSum) { return {}; }
DescriptorResult SmilesIodineENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if (tlsElementCache.elements[i] == "I" && i < tlsElementCache.paulingValues.size())
            sum += tlsElementCache.paulingValues[i];
    return sum;
}

// 105. EN sum of elements with EN == 3.98 or 3.16 or 2.96 or 2.66 (all halogens)
DECLARE_DESCRIPTOR(SmilesAllHalogenENSum, ElectronegativityDescriptor, "Sum of Pauling EN for all halogen atoms")
DESCRIPTOR_DEPENDENCIES(SmilesAllHalogenENSum) { return {}; }
DescriptorResult SmilesAllHalogenENSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if ((tlsElementCache.elements[i] == "F" || tlsElementCache.elements[i] == "Cl" ||
             tlsElementCache.elements[i] == "Br" || tlsElementCache.elements[i] == "I") && i < tlsElementCache.paulingValues.size())
            sum += tlsElementCache.paulingValues[i];
    return sum;
}

// 106. EN mean of all halogens
DECLARE_DESCRIPTOR(SmilesAllHalogenENMean, ElectronegativityDescriptor, "Mean Pauling EN for all halogen atoms")
DESCRIPTOR_DEPENDENCIES(SmilesAllHalogenENMean) { return {}; }
DescriptorResult SmilesAllHalogenENMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0, count = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i)
        if ((tlsElementCache.elements[i] == "F" || tlsElementCache.elements[i] == "Cl" ||
             tlsElementCache.elements[i] == "Br" || tlsElementCache.elements[i] == "I") && i < tlsElementCache.paulingValues.size()) {
            sum += tlsElementCache.paulingValues[i];
            count += 1.0;
        }
    return count == 0.0 ? 0.0 : sum / count;
}


// 108. EN sum of elements with EN < 1.5
DECLARE_DESCRIPTOR(SmilesVeryLowENAtomSum, ElectronegativityDescriptor, "Sum of Pauling EN for atoms with EN < 1.5")
DESCRIPTOR_DEPENDENCIES(SmilesVeryLowENAtomSum) { return {}; }
DescriptorResult SmilesVeryLowENAtomSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v < 1.5) sum += v;
    return sum;
}

// 109. EN mean of elements with EN < 1.5
DECLARE_DESCRIPTOR(SmilesVeryLowENAtomMean, ElectronegativityDescriptor, "Mean Pauling EN for atoms with EN < 1.5")
DESCRIPTOR_DEPENDENCIES(SmilesVeryLowENAtomMean) { return {}; }
DescriptorResult SmilesVeryLowENAtomMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0, count = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v < 1.5) { sum += v; count += 1.0; }
    return count == 0.0 ? 0.0 : sum / count;
}

// 110. EN stddev of elements with EN < 1.5
DECLARE_DESCRIPTOR(SmilesVeryLowENAtomStdDev, ElectronegativityDescriptor, "Stddev of Pauling EN for atoms with EN < 1.5")
DESCRIPTOR_DEPENDENCIES(SmilesVeryLowENAtomStdDev) { return {}; }
DescriptorResult SmilesVeryLowENAtomStdDevDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (double v : tlsElementCache.paulingValues) if (v < 1.5) vals.push_back(v);
    if (vals.size() < 2) return 0.0;
    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    double var = 0.0;
    for (auto v : vals) var += (v - mean) * (v - mean);
    return std::sqrt(var / vals.size());
}

// 111. EN sum of elements with EN > 3.0
DECLARE_DESCRIPTOR(SmilesVeryHighENAtomSum, ElectronegativityDescriptor, "Sum of Pauling EN for atoms with EN > 3.0")
DESCRIPTOR_DEPENDENCIES(SmilesVeryHighENAtomSum) { return {}; }
DescriptorResult SmilesVeryHighENAtomSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v > 3.0) sum += v;
    return sum;
}

// 112. EN mean of elements with EN > 3.0
DECLARE_DESCRIPTOR(SmilesVeryHighENAtomMean, ElectronegativityDescriptor, "Mean Pauling EN for atoms with EN > 3.0")
DESCRIPTOR_DEPENDENCIES(SmilesVeryHighENAtomMean) { return {}; }
DescriptorResult SmilesVeryHighENAtomMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0, count = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v > 3.0) { sum += v; count += 1.0; }
    return count == 0.0 ? 0.0 : sum / count;
}

// 113. EN stddev of elements with EN > 3.0
DECLARE_DESCRIPTOR(SmilesVeryHighENAtomStdDev, ElectronegativityDescriptor, "Stddev of Pauling EN for atoms with EN > 3.0")
DESCRIPTOR_DEPENDENCIES(SmilesVeryHighENAtomStdDev) { return {}; }
DescriptorResult SmilesVeryHighENAtomStdDevDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (double v : tlsElementCache.paulingValues) if (v > 3.0) vals.push_back(v);
    if (vals.size() < 2) return 0.0;
    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    double var = 0.0;
    for (auto v : vals) var += (v - mean) * (v - mean);
    return std::sqrt(var / vals.size());
}

// 114. EN sum of elements with EN between 1.5 and 2.5
DECLARE_DESCRIPTOR(SmilesLowMidENAtomSum, ElectronegativityDescriptor, "Sum of Pauling EN for atoms with 1.5 < EN < 2.5")
DESCRIPTOR_DEPENDENCIES(SmilesLowMidENAtomSum) { return {}; }
DescriptorResult SmilesLowMidENAtomSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v > 1.5 && v < 2.5) sum += v;
    return sum;
}

// 115. EN mean of elements with EN between 1.5 and 2.5
DECLARE_DESCRIPTOR(SmilesLowMidENAtomMean, ElectronegativityDescriptor, "Mean Pauling EN for atoms with 1.5 < EN < 2.5")
DESCRIPTOR_DEPENDENCIES(SmilesLowMidENAtomMean) { return {}; }
DescriptorResult SmilesLowMidENAtomMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0, count = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v > 1.5 && v < 2.5) { sum += v; count += 1.0; }
    return count == 0.0 ? 0.0 : sum / count;
}

// 116. EN stddev of elements with EN between 1.5 and 2.5
DECLARE_DESCRIPTOR(SmilesLowMidENAtomStdDev, ElectronegativityDescriptor, "Stddev of Pauling EN for atoms with 1.5 < EN < 2.5")
DESCRIPTOR_DEPENDENCIES(SmilesLowMidENAtomStdDev) { return {}; }
DescriptorResult SmilesLowMidENAtomStdDevDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (double v : tlsElementCache.paulingValues) if (v > 1.5 && v < 2.5) vals.push_back(v);
    if (vals.size() < 2) return 0.0;
    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    double var = 0.0;
    for (auto v : vals) var += (v - mean) * (v - mean);
    return std::sqrt(var / vals.size());
}

// 117. EN sum of elements with EN between 2.5 and 3.5
DECLARE_DESCRIPTOR(SmilesMidHighENAtomSum, ElectronegativityDescriptor, "Sum of Pauling EN for atoms with 2.5 < EN < 3.5")
DESCRIPTOR_DEPENDENCIES(SmilesMidHighENAtomSum) { return {}; }
DescriptorResult SmilesMidHighENAtomSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v > 2.5 && v < 3.5) sum += v;
    return sum;
}

// 118. EN mean of elements with EN between 2.5 and 3.5
DECLARE_DESCRIPTOR(SmilesMidHighENAtomMean, ElectronegativityDescriptor, "Mean Pauling EN for atoms with 2.5 < EN < 3.5")
DESCRIPTOR_DEPENDENCIES(SmilesMidHighENAtomMean) { return {}; }
DescriptorResult SmilesMidHighENAtomMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0, count = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v > 2.5 && v < 3.5) { sum += v; count += 1.0; }
    return count == 0.0 ? 0.0 : sum / count;
}

// 119. EN stddev of elements with EN between 2.5 and 3.5
DECLARE_DESCRIPTOR(SmilesMidHighENAtomStdDev, ElectronegativityDescriptor, "Stddev of Pauling EN for atoms with 2.5 < EN < 3.5")
DESCRIPTOR_DEPENDENCIES(SmilesMidHighENAtomStdDev) { return {}; }
DescriptorResult SmilesMidHighENAtomStdDevDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (double v : tlsElementCache.paulingValues) if (v > 2.5 && v < 3.5) vals.push_back(v);
    if (vals.size() < 2) return 0.0;
    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    double var = 0.0;
    for (auto v : vals) var += (v - mean) * (v - mean);
    return std::sqrt(var / vals.size());
}

// 120. EN sum of elements with EN between 3.0 and 4.0
DECLARE_DESCRIPTOR(SmilesHighVeryHighENAtomSum, ElectronegativityDescriptor, "Sum of Pauling EN for atoms with 3.0 < EN < 4.0")
DESCRIPTOR_DEPENDENCIES(SmilesHighVeryHighENAtomSum) { return {}; }
DescriptorResult SmilesHighVeryHighENAtomSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v > 3.0 && v < 4.0) sum += v;
    return sum;
}

// 121. EN mean of elements with EN between 3.0 and 4.0
DECLARE_DESCRIPTOR(SmilesHighVeryHighENAtomMean, ElectronegativityDescriptor, "Mean Pauling EN for atoms with 3.0 < EN < 4.0")
DESCRIPTOR_DEPENDENCIES(SmilesHighVeryHighENAtomMean) { return {}; }
DescriptorResult SmilesHighVeryHighENAtomMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0, count = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v > 3.0 && v < 4.0) { sum += v; count += 1.0; }
    return count == 0.0 ? 0.0 : sum / count;
}

// 122. EN stddev of elements with EN between 3.0 and 4.0
DECLARE_DESCRIPTOR(SmilesHighVeryHighENAtomStdDev, ElectronegativityDescriptor, "Stddev of Pauling EN for atoms with 3.0 < EN < 4.0")
DESCRIPTOR_DEPENDENCIES(SmilesHighVeryHighENAtomStdDev) { return {}; }
DescriptorResult SmilesHighVeryHighENAtomStdDevDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (double v : tlsElementCache.paulingValues) if (v > 3.0 && v < 4.0) vals.push_back(v);
    if (vals.size() < 2) return 0.0;
    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    double var = 0.0;
    for (auto v : vals) var += (v - mean) * (v - mean);
    return std::sqrt(var / vals.size());
}

// 123. EN sum of elements with EN between 1.0 and 2.0
DECLARE_DESCRIPTOR(SmilesVeryLowLowENAtomSum, ElectronegativityDescriptor, "Sum of Pauling EN for atoms with 1.0 < EN < 2.0")
DESCRIPTOR_DEPENDENCIES(SmilesVeryLowLowENAtomSum) { return {}; }
DescriptorResult SmilesVeryLowLowENAtomSumDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v > 1.0 && v < 2.0) sum += v;
    return sum;
}

// 124. EN mean of elements with EN between 1.0 and 2.0
DECLARE_DESCRIPTOR(SmilesVeryLowLowENAtomMean, ElectronegativityDescriptor, "Mean Pauling EN for atoms with 1.0 < EN < 2.0")
DESCRIPTOR_DEPENDENCIES(SmilesVeryLowLowENAtomMean) { return {}; }
DescriptorResult SmilesVeryLowLowENAtomMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0, count = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v > 1.0 && v < 2.0) { sum += v; count += 1.0; }
    return count == 0.0 ? 0.0 : sum / count;
}

// 125. EN stddev of elements with EN between 1.0 and 2.0
DECLARE_DESCRIPTOR(SmilesVeryLowLowENAtomStdDev, ElectronegativityDescriptor, "Stddev of Pauling EN for atoms with 1.0 < EN < 2.0")
DESCRIPTOR_DEPENDENCIES(SmilesVeryLowLowENAtomStdDev) { return {}; }
DescriptorResult SmilesVeryLowLowENAtomStdDevDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (double v : tlsElementCache.paulingValues) if (v > 1.0 && v < 2.0) vals.push_back(v);
    if (vals.size() < 2) return 0.0;
    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    double var = 0.0;
    for (auto v : vals) var += (v - mean) * (v - mean);
    return std::sqrt(var / vals.size());
}

// 126. Skewness of Pauling EN distribution
DECLARE_DESCRIPTOR(SmilesENSkewness, ElectronegativityDescriptor, "Skewness of Pauling EN values")
DESCRIPTOR_DEPENDENCIES(SmilesENSkewness) { return {}; }
DescriptorResult SmilesENSkewnessDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    const auto& v = tlsElementCache.paulingValues;
    if (v.size() < 3) return 0.0;
    double mean = meanElectronegativity(v);
    double m2 = 0.0, m3 = 0.0;
    for (double x : v) {
        double d = x - mean;
        m2 += d * d;
        m3 += d * d * d;
    }
    m2 /= v.size();
    m3 /= v.size();
    return m2 == 0.0 ? 0.0 : m3 / std::pow(m2, 1.5);
}

// 127. Kurtosis of Pauling EN distribution
DECLARE_DESCRIPTOR(SmilesENKurtosis, ElectronegativityDescriptor, "Kurtosis of Pauling EN values")
DESCRIPTOR_DEPENDENCIES(SmilesENKurtosis) { return {}; }
DescriptorResult SmilesENKurtosisDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    const auto& v = tlsElementCache.paulingValues;
    if (v.size() < 4) return 0.0;
    double mean = meanElectronegativity(v);
    double m2 = 0.0, m4 = 0.0;
    for (double x : v) {
        double d = x - mean;
        m2 += d * d;
        m4 += d * d * d * d;
    }
    m2 /= v.size();
    m4 /= v.size();
    return m2 == 0.0 ? 0.0 : m4 / (m2 * m2) - 3.0;
}

// 128. Shannon entropy of EN value histogram (rounded to 1 decimal)
DECLARE_DESCRIPTOR(SmilesENEntropy, ElectronegativityDescriptor, "Shannon entropy of Pauling EN value histogram")
DESCRIPTOR_DEPENDENCIES(SmilesENEntropy) { return {}; }
DescriptorResult SmilesENEntropyDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    if (tlsElementCache.paulingValues.empty()) return 0.0;
    std::map<int, int> hist;
    for (double v : tlsElementCache.paulingValues)
        hist[static_cast<int>(std::round(v * 10.0))]++;
    double entropy = 0.0, n = tlsElementCache.paulingValues.size();
    for (const auto& p : hist) {
        double f = p.second / n;
        entropy -= f * std::log2(f);
    }
    return entropy;
}

// 129. Product of mean EN and stddev
DECLARE_DESCRIPTOR(SmilesENMeanStdProduct, ElectronegativityDescriptor, "Product of mean and stddev of Pauling EN")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanStdProduct) { return {}; }
DescriptorResult SmilesENMeanStdProductDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double mean = meanElectronegativity(tlsElementCache.paulingValues);
    double stddev = electronegativityStdDev(tlsElementCache.paulingValues);
    return mean * stddev;
}

// 130. Ratio of max to min EN
DECLARE_DESCRIPTOR(SmilesENMaxMinRatio, ElectronegativityDescriptor, "Ratio of max to min Pauling EN")
DESCRIPTOR_DEPENDENCIES(SmilesENMaxMinRatio) { return {}; }
DescriptorResult SmilesENMaxMinRatioDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double minv = minElectronegativity(tlsElementCache.paulingValues);
    double maxv = maxElectronegativity(tlsElementCache.paulingValues);
    return minv == 0.0 ? 0.0 : maxv / minv;
}

// 131. Range divided by mean EN
DECLARE_DESCRIPTOR(SmilesENRangeOverMean, ElectronegativityDescriptor, "Range divided by mean Pauling EN")
DESCRIPTOR_DEPENDENCIES(SmilesENRangeOverMean) { return {}; }
DescriptorResult SmilesENRangeOverMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double mean = meanElectronegativity(tlsElementCache.paulingValues);
    double range = electronegativityRange(tlsElementCache.paulingValues);
    return mean == 0.0 ? 0.0 : range / mean;
}

// 132. Mean EN of first half of atoms
DECLARE_DESCRIPTOR(SmilesENMeanFirstHalf, ElectronegativityDescriptor, "Mean Pauling EN of first half of atoms")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanFirstHalf) { return {}; }
DescriptorResult SmilesENMeanFirstHalfDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    size_t n = tlsElementCache.paulingValues.size();
    if (n == 0) return 0.0;
    size_t half = n / 2;
    double sum = 0.0;
    for (size_t i = 0; i < half; ++i) sum += tlsElementCache.paulingValues[i];
    return half == 0 ? 0.0 : sum / half;
}

// 133. Mean EN of second half of atoms
DECLARE_DESCRIPTOR(SmilesENMeanSecondHalf, ElectronegativityDescriptor, "Mean Pauling EN of second half of atoms")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanSecondHalf) { return {}; }
DescriptorResult SmilesENMeanSecondHalfDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    size_t n = tlsElementCache.paulingValues.size();
    if (n == 0) return 0.0;
    size_t half = n / 2;
    double sum = 0.0;
    for (size_t i = half; i < n; ++i) sum += tlsElementCache.paulingValues[i];
    return (n - half) == 0 ? 0.0 : sum / (n - half);
}

// 134. Difference between mean EN of first and second half
DECLARE_DESCRIPTOR(SmilesENMeanHalfDiff, ElectronegativityDescriptor, "Difference between mean EN of first and second half")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanHalfDiff) { return {}; }
DescriptorResult SmilesENMeanHalfDiffDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    size_t n = tlsElementCache.paulingValues.size();
    if (n == 0) return 0.0;
    size_t half = n / 2;
    double sum1 = 0.0, sum2 = 0.0;
    for (size_t i = 0; i < half; ++i) sum1 += tlsElementCache.paulingValues[i];
    for (size_t i = half; i < n; ++i) sum2 += tlsElementCache.paulingValues[i];
    double mean1 = half == 0 ? 0.0 : sum1 / half;
    double mean2 = (n - half) == 0 ? 0.0 : sum2 / (n - half);
    return mean1 - mean2;
}

// 135. Mean absolute difference between consecutive EN values
DECLARE_DESCRIPTOR(SmilesENMeanAbsDiff, ElectronegativityDescriptor, "Mean absolute difference between consecutive Pauling EN values")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanAbsDiff) { return {}; }
DescriptorResult SmilesENMeanAbsDiffDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    const auto& v = tlsElementCache.paulingValues;
    if (v.size() < 2) return 0.0;
    double sum = 0.0;
    for (size_t i = 1; i < v.size(); ++i) sum += std::abs(v[i] - v[i-1]);
    return sum / (v.size() - 1);
}

// 136. Max absolute difference between consecutive EN values
DECLARE_DESCRIPTOR(SmilesENMaxAbsDiff, ElectronegativityDescriptor, "Max absolute difference between consecutive Pauling EN values")
DESCRIPTOR_DEPENDENCIES(SmilesENMaxAbsDiff) { return {}; }
DescriptorResult SmilesENMaxAbsDiffDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    const auto& v = tlsElementCache.paulingValues;
    if (v.size() < 2) return 0.0;
    double maxd = 0.0;
    for (size_t i = 1; i < v.size(); ++i) maxd = std::max(maxd, std::abs(v[i] - v[i-1]));
    return maxd;
}

// 137. Min absolute difference between consecutive EN values
DECLARE_DESCRIPTOR(SmilesENMinAbsDiff, ElectronegativityDescriptor, "Min absolute difference between consecutive Pauling EN values")
DESCRIPTOR_DEPENDENCIES(SmilesENMinAbsDiff) { return {}; }
DescriptorResult SmilesENMinAbsDiffDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    const auto& v = tlsElementCache.paulingValues;
    if (v.size() < 2) return 0.0;
    double mind = std::abs(v[1] - v[0]);
    for (size_t i = 2; i < v.size(); ++i) mind = std::min(mind, std::abs(v[i] - v[i-1]));
    return mind;
}

// 138. Mean EN of even-indexed atoms
DECLARE_DESCRIPTOR(SmilesENMeanEven, ElectronegativityDescriptor, "Mean Pauling EN of even-indexed atoms")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanEven) { return {}; }
DescriptorResult SmilesENMeanEvenDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0; size_t count = 0;
    for (size_t i = 0; i < tlsElementCache.paulingValues.size(); i += 2) {
        sum += tlsElementCache.paulingValues[i];
        count++;
    }
    return count == 0 ? 0.0 : sum / count;
}

// 139. Mean EN of odd-indexed atoms
DECLARE_DESCRIPTOR(SmilesENMeanOdd, ElectronegativityDescriptor, "Mean Pauling EN of odd-indexed atoms")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanOdd) { return {}; }
DescriptorResult SmilesENMeanOddDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0; size_t count = 0;
    for (size_t i = 1; i < tlsElementCache.paulingValues.size(); i += 2) {
        sum += tlsElementCache.paulingValues[i];
        count++;
    }
    return count == 0 ? 0.0 : sum / count;
}

// 140. Difference between mean EN of even and odd-indexed atoms
DECLARE_DESCRIPTOR(SmilesENMeanEvenOddDiff, ElectronegativityDescriptor, "Difference between mean EN of even and odd-indexed atoms")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanEvenOddDiff) { return {}; }
DescriptorResult SmilesENMeanEvenOddDiffDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sumEven = 0.0, sumOdd = 0.0; size_t countEven = 0, countOdd = 0;
    for (size_t i = 0; i < tlsElementCache.paulingValues.size(); ++i) {
        if (i % 2 == 0) { sumEven += tlsElementCache.paulingValues[i]; countEven++; }
        else { sumOdd += tlsElementCache.paulingValues[i]; countOdd++; }
    }
    double meanEven = countEven == 0 ? 0.0 : sumEven / countEven;
    double meanOdd = countOdd == 0 ? 0.0 : sumOdd / countOdd;
    return meanEven - meanOdd;
}

// 141. Product of max and min EN
DECLARE_DESCRIPTOR(SmilesENMaxMinProduct, ElectronegativityDescriptor, "Product of max and min Pauling EN")
DESCRIPTOR_DEPENDENCIES(SmilesENMaxMinProduct) { return {}; }
DescriptorResult SmilesENMaxMinProductDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double minv = minElectronegativity(tlsElementCache.paulingValues);
    double maxv = maxElectronegativity(tlsElementCache.paulingValues);
    return minv * maxv;
}

// 142. Ratio of mean EN to stddev
DECLARE_DESCRIPTOR(SmilesENMeanOverStd, ElectronegativityDescriptor, "Ratio of mean to stddev of Pauling EN")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanOverStd) { return {}; }
DescriptorResult SmilesENMeanOverStdDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double mean = meanElectronegativity(tlsElementCache.paulingValues);
    double stddev = electronegativityStdDev(tlsElementCache.paulingValues);
    return stddev == 0.0 ? 0.0 : mean / stddev;
}

// 143. Ratio of stddev to mean EN
DECLARE_DESCRIPTOR(SmilesENStdOverMean, ElectronegativityDescriptor, "Ratio of stddev to mean of Pauling EN")
DESCRIPTOR_DEPENDENCIES(SmilesENStdOverMean) { return {}; }
DescriptorResult SmilesENStdOverMeanDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double mean = meanElectronegativity(tlsElementCache.paulingValues);
    double stddev = electronegativityStdDev(tlsElementCache.paulingValues);
    return mean == 0.0 ? 0.0 : stddev / mean;
}

// 144. EN sum of atoms with atomic weight < 20
DECLARE_DESCRIPTOR(SmilesENSumLightAtoms, ElectronegativityDescriptor, "Sum of Pauling EN for atoms with atomic weight < 20")
DESCRIPTOR_DEPENDENCIES(SmilesENSumLightAtoms) { return {}; }
DescriptorResult SmilesENSumLightAtomsDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i) {
        auto it = atomicWeights.find(tlsElementCache.elements[i]);
        if (it != atomicWeights.end() && it->second < 20.0 && i < tlsElementCache.paulingValues.size())
            sum += tlsElementCache.paulingValues[i];
    }
    return sum;
}

// 145. EN sum of atoms with atomic weight > 50
DECLARE_DESCRIPTOR(SmilesENSumHeavyAtoms, ElectronegativityDescriptor, "Sum of Pauling EN for atoms with atomic weight > 50")
DESCRIPTOR_DEPENDENCIES(SmilesENSumHeavyAtoms) { return {}; }
DescriptorResult SmilesENSumHeavyAtomsDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i) {
        auto it = atomicWeights.find(tlsElementCache.elements[i]);
        if (it != atomicWeights.end() && it->second > 50.0 && i < tlsElementCache.paulingValues.size())
            sum += tlsElementCache.paulingValues[i];
    }
    return sum;
}

// 146. EN mean of atoms with atomic weight < 20
DECLARE_DESCRIPTOR(SmilesENMeanLightAtoms, ElectronegativityDescriptor, "Mean Pauling EN for atoms with atomic weight < 20")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanLightAtoms) { return {}; }
DescriptorResult SmilesENMeanLightAtomsDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0; size_t count = 0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i) {
        auto it = atomicWeights.find(tlsElementCache.elements[i]);
        if (it != atomicWeights.end() && it->second < 20.0 && i < tlsElementCache.paulingValues.size()) {
            sum += tlsElementCache.paulingValues[i];
            count++;
        }
    }
    return count == 0 ? 0.0 : sum / count;
}

// 147. EN mean of atoms with atomic weight > 50
DECLARE_DESCRIPTOR(SmilesENMeanHeavyAtoms, ElectronegativityDescriptor, "Mean Pauling EN for atoms with atomic weight > 50")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanHeavyAtoms) { return {}; }
DescriptorResult SmilesENMeanHeavyAtomsDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0; size_t count = 0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i) {
        auto it = atomicWeights.find(tlsElementCache.elements[i]);
        if (it != atomicWeights.end() && it->second > 50.0 && i < tlsElementCache.paulingValues.size()) {
            sum += tlsElementCache.paulingValues[i];
            count++;
        }
    }
    return count == 0 ? 0.0 : sum / count;
}

// 148. EN stddev of atoms with atomic weight < 20
DECLARE_DESCRIPTOR(SmilesENStdDevLightAtoms, ElectronegativityDescriptor, "Stddev of Pauling EN for atoms with atomic weight < 20")
DESCRIPTOR_DEPENDENCIES(SmilesENStdDevLightAtoms) { return {}; }
DescriptorResult SmilesENStdDevLightAtomsDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i) {
        auto it = atomicWeights.find(tlsElementCache.elements[i]);
        if (it != atomicWeights.end() && it->second < 20.0 && i < tlsElementCache.paulingValues.size())
            vals.push_back(tlsElementCache.paulingValues[i]);
    }
    if (vals.size() < 2) return 0.0;
    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    double var = 0.0;
    for (auto v : vals) var += (v - mean) * (v - mean);
    return std::sqrt(var / vals.size());
}

// 149. EN stddev of atoms with atomic weight > 50
DECLARE_DESCRIPTOR(SmilesENStdDevHeavyAtoms, ElectronegativityDescriptor, "Stddev of Pauling EN for atoms with atomic weight > 50")
DESCRIPTOR_DEPENDENCIES(SmilesENStdDevHeavyAtoms) { return {}; }
DescriptorResult SmilesENStdDevHeavyAtomsDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    std::vector<double> vals;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i) {
        auto it = atomicWeights.find(tlsElementCache.elements[i]);
        if (it != atomicWeights.end() && it->second > 50.0 && i < tlsElementCache.paulingValues.size())
            vals.push_back(tlsElementCache.paulingValues[i]);
    }
    if (vals.size() < 2) return 0.0;
    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    double var = 0.0;
    for (auto v : vals) var += (v - mean) * (v - mean);
    return std::sqrt(var / vals.size());
}

// 150. Ratio of mean EN of light to heavy atoms
DECLARE_DESCRIPTOR(SmilesENMeanLightHeavyRatio, ElectronegativityDescriptor, "Ratio of mean EN of light to heavy atoms")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanLightHeavyRatio) { return {}; }
DescriptorResult SmilesENMeanLightHeavyRatioDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sumLight = 0.0, countLight = 0.0, sumHeavy = 0.0, countHeavy = 0.0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i) {
        auto it = atomicWeights.find(tlsElementCache.elements[i]);
        if (it != atomicWeights.end() && i < tlsElementCache.paulingValues.size()) {
            if (it->second < 20.0) { sumLight += tlsElementCache.paulingValues[i]; countLight += 1.0; }
            if (it->second > 50.0) { sumHeavy += tlsElementCache.paulingValues[i]; countHeavy += 1.0; }
        }
    }
    double meanLight = countLight == 0.0 ? 0.0 : sumLight / countLight;
    double meanHeavy = countHeavy == 0.0 ? 0.0 : sumHeavy / countHeavy;
    return meanHeavy == 0.0 ? 0.0 : meanLight / meanHeavy;
}

// 51. Fraction of atoms with EN above the mean
DECLARE_DESCRIPTOR(SmilesENAboveMeanFraction, ElectronegativityDescriptor, "Fraction of atoms with Pauling EN above the mean")
DESCRIPTOR_DEPENDENCIES(SmilesENAboveMeanFraction) { return {}; }
DescriptorResult SmilesENAboveMeanFractionDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    if (tlsElementCache.paulingValues.empty()) return 0.0;
    double mean = meanElectronegativity(tlsElementCache.paulingValues);
    double count = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v > mean) count += 1.0;
    return count / tlsElementCache.paulingValues.size();
}

// 52. Fraction of atoms with EN below the mean
DECLARE_DESCRIPTOR(SmilesENBelowMeanFraction, ElectronegativityDescriptor, "Fraction of atoms with Pauling EN below the mean")
DESCRIPTOR_DEPENDENCIES(SmilesENBelowMeanFraction) { return {}; }
DescriptorResult SmilesENBelowMeanFractionDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    if (tlsElementCache.paulingValues.empty()) return 0.0;
    double mean = meanElectronegativity(tlsElementCache.paulingValues);
    double count = 0.0;
    for (double v : tlsElementCache.paulingValues) if (v < mean) count += 1.0;
    return count / tlsElementCache.paulingValues.size();
}

// 53. EN difference between most and least common element
DECLARE_DESCRIPTOR(SmilesENMostLeastCommonDiff, ElectronegativityDescriptor, "Difference in EN between most and least common element")
DESCRIPTOR_DEPENDENCIES(SmilesENMostLeastCommonDiff) { return {}; }
DescriptorResult SmilesENMostLeastCommonDiffDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    if (tlsElementCache.elements.empty()) return 0.0;
    std::map<std::string, int> freq;
    for (const auto& e : tlsElementCache.elements) freq[e]++;
    auto most = std::max_element(freq.begin(), freq.end(),
        [](const auto& a, const auto& b){ return a.second < b.second; });
    auto least = std::min_element(freq.begin(), freq.end(),
        [](const auto& a, const auto& b){ return a.second < b.second; });
    auto it1 = paulingElectronegativityMap.find(most->first);
    auto it2 = paulingElectronegativityMap.find(least->first);
    if (it1 == paulingElectronegativityMap.end() || it2 == paulingElectronegativityMap.end()) return 0.0;
    return std::abs(it1->second - it2->second);
}

// 54. EN mean of atoms with odd atomic number
DECLARE_DESCRIPTOR(SmilesENMeanOddZ, ElectronegativityDescriptor, "Mean Pauling EN for atoms with odd atomic number")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanOddZ) { return {}; }
DescriptorResult SmilesENMeanOddZDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0; int count = 0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i) {
        const auto& e = tlsElementCache.elements[i];
        int z = 0;
        if      (e == "H" || e == "h") z = 1;
        else if (e == "B" || e == "b") z = 5;
        else if (e == "N" || e == "n") z = 7;
        else if (e == "F") z = 9;
        else if (e == "P" || e == "p") z = 15;
        else if (e == "Cl") z = 17;
        else if (e == "K") z = 19;
        else if (e == "V") z = 23;
        else if (e == "Mn") z = 25;
        else if (e == "Co") z = 27;
        else if (e == "Cu") z = 29;
        else if (e == "Ga") z = 31;
        else if (e == "As") z = 33;
        else if (e == "Br") z = 35;
        else if (e == "I" || e == "i") z = 53;
        if (z % 2 == 1 && i < tlsElementCache.paulingValues.size()) {
            sum += tlsElementCache.paulingValues[i];
            count++;
        }
    }
    return count == 0 ? 0.0 : sum / count;
}

// 55. EN mean of atoms with even atomic number
DECLARE_DESCRIPTOR(SmilesENMeanEvenZ, ElectronegativityDescriptor, "Mean Pauling EN for atoms with even atomic number")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanEvenZ) { return {}; }
DescriptorResult SmilesENMeanEvenZDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0; int count = 0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i) {
        const auto& e = tlsElementCache.elements[i];
        int z = 0;
        if      (e == "He") z = 2;
        else if (e == "Be") z = 4;
        else if (e == "C" || e == "c") z = 6;
        else if (e == "O" || e == "o") z = 8;
        else if (e == "Ne") z = 10;
        else if (e == "Mg") z = 12;
        else if (e == "Si") z = 14;
        else if (e == "S" || e == "s") z = 16;
        else if (e == "Ca") z = 20;
        else if (e == "Ti") z = 22;
        else if (e == "Cr") z = 24;
        else if (e == "Fe") z = 26;
        else if (e == "Ni") z = 28;
        else if (e == "Zn") z = 30;
        else if (e == "Se") z = 34;
        if (z % 2 == 0 && z > 0 && i < tlsElementCache.paulingValues.size()) {
            sum += tlsElementCache.paulingValues[i];
            count++;
        }
    }
    return count == 0 ? 0.0 : sum / count;
}

// 56. EN mean of group 1/alkali metals
DECLARE_DESCRIPTOR(SmilesENMeanAlkali, ElectronegativityDescriptor, "Mean Pauling EN for alkali metals")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanAlkali) { return {}; }
DescriptorResult SmilesENMeanAlkaliDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0; int count = 0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i) {
        const auto& e = tlsElementCache.elements[i];
        if ((e == "Li" || e == "Na" || e == "K" || e == "Rb") && i < tlsElementCache.paulingValues.size()) {
            sum += tlsElementCache.paulingValues[i];
            count++;
        }
    }
    return count == 0 ? 0.0 : sum / count;
}

// 57. EN mean of group 2/alkaline earth metals
DECLARE_DESCRIPTOR(SmilesENMeanAlkalineEarth, ElectronegativityDescriptor, "Mean Pauling EN for alkaline earth metals")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanAlkalineEarth) { return {}; }
DescriptorResult SmilesENMeanAlkalineEarthDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0; int count = 0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i) {
        const auto& e = tlsElementCache.elements[i];
        if ((e == "Be" || e == "Mg" || e == "Ca" || e == "Sr") && i < tlsElementCache.paulingValues.size()) {
            sum += tlsElementCache.paulingValues[i];
            count++;
        }
    }
    return count == 0 ? 0.0 : sum / count;
}

// 58. EN mean of group 17/halogens
DECLARE_DESCRIPTOR(SmilesENMeanHalogen, ElectronegativityDescriptor, "Mean Pauling EN for halogens")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanHalogen) { return {}; }
DescriptorResult SmilesENMeanHalogenDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0; int count = 0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i) {
        const auto& e = tlsElementCache.elements[i];
        if ((e == "F" || e == "Cl" || e == "Br" || e == "I") && i < tlsElementCache.paulingValues.size()) {
            sum += tlsElementCache.paulingValues[i];
            count++;
        }
    }
    return count == 0 ? 0.0 : sum / count;
}

// 59. EN mean of group 16/chalcogens
DECLARE_DESCRIPTOR(SmilesENMeanChalcogen, ElectronegativityDescriptor, "Mean Pauling EN for chalcogens")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanChalcogen) { return {}; }
DescriptorResult SmilesENMeanChalcogenDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0; int count = 0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i) {
        const auto& e = tlsElementCache.elements[i];
        if ((e == "O" || e == "S" || e == "Se" || e == "Te") && i < tlsElementCache.paulingValues.size()) {
            sum += tlsElementCache.paulingValues[i];
            count++;
        }
    }
    return count == 0 ? 0.0 : sum / count;
}

// 60. EN mean of group 15/pnictogens
DECLARE_DESCRIPTOR(SmilesENMeanPnictogen, ElectronegativityDescriptor, "Mean Pauling EN for pnictogens")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanPnictogen) { return {}; }
DescriptorResult SmilesENMeanPnictogenDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0; int count = 0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i) {
        const auto& e = tlsElementCache.elements[i];
        if ((e == "N" || e == "P" || e == "As" || e == "Sb") && i < tlsElementCache.paulingValues.size()) {
            sum += tlsElementCache.paulingValues[i];
            count++;
        }
    }
    return count == 0 ? 0.0 : sum / count;
}

// 61. EN mean of group 14/tetrels
DECLARE_DESCRIPTOR(SmilesENMeanTetrel, ElectronegativityDescriptor, "Mean Pauling EN for tetrels")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanTetrel) { return {}; }
DescriptorResult SmilesENMeanTetrelDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0; int count = 0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i) {
        const auto& e = tlsElementCache.elements[i];
        if ((e == "C" || e == "Si" || e == "Ge" || e == "Sn") && i < tlsElementCache.paulingValues.size()) {
            sum += tlsElementCache.paulingValues[i];
            count++;
        }
    }
    return count == 0 ? 0.0 : sum / count;
}

// 62. EN mean of group 13/triels
DECLARE_DESCRIPTOR(SmilesENMeanTriel, ElectronegativityDescriptor, "Mean Pauling EN for triels")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanTriel) { return {}; }
DescriptorResult SmilesENMeanTrielDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0; int count = 0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i) {
        const auto& e = tlsElementCache.elements[i];
        if ((e == "B" || e == "Al" || e == "Ga" || e == "In") && i < tlsElementCache.paulingValues.size()) {
            sum += tlsElementCache.paulingValues[i];
            count++;
        }
    }
    return count == 0 ? 0.0 : sum / count;
}

// 63. EN mean of group 12/coinage metals
DECLARE_DESCRIPTOR(SmilesENMeanCoinage, ElectronegativityDescriptor, "Mean Pauling EN for coinage metals")
DESCRIPTOR_DEPENDENCIES(SmilesENMeanCoinage) { return {}; }
DescriptorResult SmilesENMeanCoinageDescriptor::calculate(Context& context) const {
    extractElementData(context.getSmiles(), tlsElementCache);
    double sum = 0.0; int count = 0;
    for (size_t i = 0; i < tlsElementCache.elements.size(); ++i) {
        const auto& e = tlsElementCache.elements[i];
        if ((e == "Cu" || e == "Ag" || e == "Au") && i < tlsElementCache.paulingValues.size()) {
            sum += tlsElementCache.paulingValues[i];
            count++;
        }
    }
    return count == 0 ? 0.0 : sum / count;
}


void register_SmilesPaulingENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesPaulingENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesPaulingENMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesPaulingENMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesPaulingENMaxDescriptor() {
    auto descriptor = std::make_shared<SmilesPaulingENMaxDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesPaulingENMinDescriptor() {
    auto descriptor = std::make_shared<SmilesPaulingENMinDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesPaulingENRangeDescriptor() {
    auto descriptor = std::make_shared<SmilesPaulingENRangeDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesPaulingENStdDevDescriptor() {
    auto descriptor = std::make_shared<SmilesPaulingENStdDevDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesHighENCountDescriptor() {
    auto descriptor = std::make_shared<SmilesHighENCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesHighToLowENRatioDescriptor() {
    auto descriptor = std::make_shared<SmilesHighToLowENRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesVeryHighENCountDescriptor() {
    auto descriptor = std::make_shared<SmilesVeryHighENCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesHighENFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesHighENFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAllredRochowENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesAllredRochowENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAllredRochowENMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesAllredRochowENMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAllredToPaulingRatioDescriptor() {
    auto descriptor = std::make_shared<SmilesAllredToPaulingRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesMullikenENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesMullikenENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesMullikenENMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesMullikenENMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesSandersonENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesSandersonENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesSandersonENMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesSandersonENMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesCarbonENProductDescriptor() {
    auto descriptor = std::make_shared<SmilesCarbonENProductDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNitrogenENProductDescriptor() {
    auto descriptor = std::make_shared<SmilesNitrogenENProductDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesOxygenENProductDescriptor() {
    auto descriptor = std::make_shared<SmilesOxygenENProductDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesSulfurENProductDescriptor() {
    auto descriptor = std::make_shared<SmilesSulfurENProductDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesHalogenENProductDescriptor() {
    auto descriptor = std::make_shared<SmilesHalogenENProductDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENScaleVarianceDescriptor() {
    auto descriptor = std::make_shared<SmilesENScaleVarianceDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENConsistencyIndexDescriptor() {
    auto descriptor = std::make_shared<SmilesENConsistencyIndexDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAromaticCarbonENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesAromaticCarbonENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAliphaticCarbonENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesAliphaticCarbonENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAromaticToAliphaticCarbonENRatioDescriptor() {
    auto descriptor = std::make_shared<SmilesAromaticToAliphaticCarbonENRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNitrogenToCarbonENRatioDescriptor() {
    auto descriptor = std::make_shared<SmilesNitrogenToCarbonENRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesOxygenToCarbonENRatioDescriptor() {
    auto descriptor = std::make_shared<SmilesOxygenToCarbonENRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesHeteroatomToCarbonENRatioDescriptor() {
    auto descriptor = std::make_shared<SmilesHeteroatomToCarbonENRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNormalizedPaulingENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesNormalizedPaulingENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNormalizedAllredRochowENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesNormalizedAllredRochowENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNormalizedMullikenENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesNormalizedMullikenENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNormalizedSandersonENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesNormalizedSandersonENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesElectronegativeElementRatioDescriptor() {
    auto descriptor = std::make_shared<SmilesElectronegativeElementRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENDiversityIndexDescriptor() {
    auto descriptor = std::make_shared<SmilesENDiversityIndexDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesHighENElementsWeightDescriptor() {
    auto descriptor = std::make_shared<SmilesHighENElementsWeightDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENHomogeneityIndexDescriptor() {
    auto descriptor = std::make_shared<SmilesENHomogeneityIndexDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENSumToSizeRatioDescriptor() {
    auto descriptor = std::make_shared<SmilesENSumToSizeRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesSquaredENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesSquaredENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesGeometricMeanENDescriptor() {
    auto descriptor = std::make_shared<SmilesGeometricMeanENDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesHarmonicMeanENDescriptor() {
    auto descriptor = std::make_shared<SmilesHarmonicMeanENDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesMedianENDescriptor() {
    auto descriptor = std::make_shared<SmilesMedianENDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesModeENDescriptor() {
    auto descriptor = std::make_shared<SmilesModeENDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesCarbonToHeteroatomENRatioDescriptor() {
    auto descriptor = std::make_shared<SmilesCarbonToHeteroatomENRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENWeightedAtomCountDescriptor() {
    auto descriptor = std::make_shared<SmilesENWeightedAtomCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENPolarityIndexDescriptor() {
    auto descriptor = std::make_shared<SmilesENPolarityIndexDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENWeightedMolWeightDescriptor() {
    auto descriptor = std::make_shared<SmilesENWeightedMolWeightDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesMultiScaleENCorrelationDescriptor() {
    auto descriptor = std::make_shared<SmilesMultiScaleENCorrelationDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENHeterogeneityIndexDescriptor() {
    auto descriptor = std::make_shared<SmilesENHeterogeneityIndexDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesUniqueElementCountDescriptor() {
    auto descriptor = std::make_shared<SmilesUniqueElementCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesUniqueElementENMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesUniqueElementENMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesUniqueElementENRangeDescriptor() {
    auto descriptor = std::make_shared<SmilesUniqueElementENRangeDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesUniqueElementENStdDevDescriptor() {
    auto descriptor = std::make_shared<SmilesUniqueElementENStdDevDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAromaticENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesAromaticENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAliphaticENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesAliphaticENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAromaticToAliphaticENRatioDescriptor() {
    auto descriptor = std::make_shared<SmilesAromaticToAliphaticENRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesHalogenENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesHalogenENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesHalogenENMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesHalogenENMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonHalogenENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesNonHalogenENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonHalogenENMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesNonHalogenENMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonHalogenENRangeDescriptor() {
    auto descriptor = std::make_shared<SmilesNonHalogenENRangeDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonHalogenENStdDevDescriptor() {
    auto descriptor = std::make_shared<SmilesNonHalogenENStdDevDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonCarbonENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesNonCarbonENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonCarbonENMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesNonCarbonENMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonCarbonENRangeDescriptor() {
    auto descriptor = std::make_shared<SmilesNonCarbonENRangeDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonCarbonENStdDevDescriptor() {
    auto descriptor = std::make_shared<SmilesNonCarbonENStdDevDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonNitrogenENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesNonNitrogenENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonNitrogenENMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesNonNitrogenENMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonNitrogenENRangeDescriptor() {
    auto descriptor = std::make_shared<SmilesNonNitrogenENRangeDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonNitrogenENStdDevDescriptor() {
    auto descriptor = std::make_shared<SmilesNonNitrogenENStdDevDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonOxygenENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesNonOxygenENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonOxygenENMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesNonOxygenENMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonOxygenENRangeDescriptor() {
    auto descriptor = std::make_shared<SmilesNonOxygenENRangeDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonOxygenENStdDevDescriptor() {
    auto descriptor = std::make_shared<SmilesNonOxygenENStdDevDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonSulfurENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesNonSulfurENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonSulfurENMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesNonSulfurENMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonSulfurENRangeDescriptor() {
    auto descriptor = std::make_shared<SmilesNonSulfurENRangeDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNonSulfurENStdDevDescriptor() {
    auto descriptor = std::make_shared<SmilesNonSulfurENStdDevDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesHighENAtomSumDescriptor() {
    auto descriptor = std::make_shared<SmilesHighENAtomSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesLowENAtomSumDescriptor() {
    auto descriptor = std::make_shared<SmilesLowENAtomSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesHighENAtomMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesHighENAtomMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesLowENAtomMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesLowENAtomMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesHighENAtomRangeDescriptor() {
    auto descriptor = std::make_shared<SmilesHighENAtomRangeDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesHighENAtomStdDevDescriptor() {
    auto descriptor = std::make_shared<SmilesHighENAtomStdDevDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesLowENAtomStdDevDescriptor() {
    auto descriptor = std::make_shared<SmilesLowENAtomStdDevDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesMidENAtomSumDescriptor() {
    auto descriptor = std::make_shared<SmilesMidENAtomSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesMidENAtomMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesMidENAtomMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesMidENAtomRangeDescriptor() {
    auto descriptor = std::make_shared<SmilesMidENAtomRangeDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesMidENAtomStdDevDescriptor() {
    auto descriptor = std::make_shared<SmilesMidENAtomStdDevDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesCarbonENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesCarbonENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNitrogenENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesNitrogenENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesOxygenENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesOxygenENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesSulfurENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesSulfurENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesFluorineENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesFluorineENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesChlorineENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesChlorineENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesBromineENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesBromineENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesIodineENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesIodineENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAllHalogenENSumDescriptor() {
    auto descriptor = std::make_shared<SmilesAllHalogenENSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAllHalogenENMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesAllHalogenENMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesVeryLowENAtomSumDescriptor() {
    auto descriptor = std::make_shared<SmilesVeryLowENAtomSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesVeryLowENAtomMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesVeryLowENAtomMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesVeryLowENAtomStdDevDescriptor() {
    auto descriptor = std::make_shared<SmilesVeryLowENAtomStdDevDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesVeryHighENAtomSumDescriptor() {
    auto descriptor = std::make_shared<SmilesVeryHighENAtomSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesVeryHighENAtomMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesVeryHighENAtomMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesVeryHighENAtomStdDevDescriptor() {
    auto descriptor = std::make_shared<SmilesVeryHighENAtomStdDevDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesLowMidENAtomSumDescriptor() {
    auto descriptor = std::make_shared<SmilesLowMidENAtomSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesLowMidENAtomMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesLowMidENAtomMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesLowMidENAtomStdDevDescriptor() {
    auto descriptor = std::make_shared<SmilesLowMidENAtomStdDevDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesMidHighENAtomSumDescriptor() {
    auto descriptor = std::make_shared<SmilesMidHighENAtomSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesMidHighENAtomMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesMidHighENAtomMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesMidHighENAtomStdDevDescriptor() {
    auto descriptor = std::make_shared<SmilesMidHighENAtomStdDevDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesHighVeryHighENAtomSumDescriptor() {
    auto descriptor = std::make_shared<SmilesHighVeryHighENAtomSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesHighVeryHighENAtomMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesHighVeryHighENAtomMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesHighVeryHighENAtomStdDevDescriptor() {
    auto descriptor = std::make_shared<SmilesHighVeryHighENAtomStdDevDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesVeryLowLowENAtomSumDescriptor() {
    auto descriptor = std::make_shared<SmilesVeryLowLowENAtomSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesVeryLowLowENAtomMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesVeryLowLowENAtomMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesVeryLowLowENAtomStdDevDescriptor() {
    auto descriptor = std::make_shared<SmilesVeryLowLowENAtomStdDevDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENSkewnessDescriptor() {
    auto descriptor = std::make_shared<SmilesENSkewnessDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENKurtosisDescriptor() {
    auto descriptor = std::make_shared<SmilesENKurtosisDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENEntropyDescriptor() {
    auto descriptor = std::make_shared<SmilesENEntropyDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanStdProductDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanStdProductDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMaxMinRatioDescriptor() {
    auto descriptor = std::make_shared<SmilesENMaxMinRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENRangeOverMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesENRangeOverMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanFirstHalfDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanFirstHalfDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanSecondHalfDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanSecondHalfDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanHalfDiffDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanHalfDiffDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanAbsDiffDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanAbsDiffDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMaxAbsDiffDescriptor() {
    auto descriptor = std::make_shared<SmilesENMaxAbsDiffDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMinAbsDiffDescriptor() {
    auto descriptor = std::make_shared<SmilesENMinAbsDiffDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanEvenDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanEvenDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanOddDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanOddDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanEvenOddDiffDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanEvenOddDiffDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMaxMinProductDescriptor() {
    auto descriptor = std::make_shared<SmilesENMaxMinProductDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanOverStdDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanOverStdDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENStdOverMeanDescriptor() {
    auto descriptor = std::make_shared<SmilesENStdOverMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENSumLightAtomsDescriptor() {
    auto descriptor = std::make_shared<SmilesENSumLightAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENSumHeavyAtomsDescriptor() {
    auto descriptor = std::make_shared<SmilesENSumHeavyAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanLightAtomsDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanLightAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanHeavyAtomsDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanHeavyAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENStdDevLightAtomsDescriptor() {
    auto descriptor = std::make_shared<SmilesENStdDevLightAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENStdDevHeavyAtomsDescriptor() {
    auto descriptor = std::make_shared<SmilesENStdDevHeavyAtomsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanLightHeavyRatioDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanLightHeavyRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENAboveMeanFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesENAboveMeanFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENBelowMeanFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesENBelowMeanFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMostLeastCommonDiffDescriptor() {
    auto descriptor = std::make_shared<SmilesENMostLeastCommonDiffDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanOddZDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanOddZDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanEvenZDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanEvenZDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanAlkaliDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanAlkaliDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanAlkalineEarthDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanAlkalineEarthDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanHalogenDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanHalogenDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanChalcogenDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanChalcogenDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanPnictogenDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanPnictogenDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanTetrelDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanTetrelDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanTrielDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanTrielDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesENMeanCoinageDescriptor() {
    auto descriptor = std::make_shared<SmilesENMeanCoinageDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}
}